"""Optional, physically grounded PDB2PQR/PROPKA -> APBS electrostatics.

No silent fallback. Canonical 11-map scores and Pareto ranks are unchanged:
electrostatics is a separate, explicitly requested, residue-matched layer.
"""
from __future__ import annotations

import math
import hashlib
import time
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

import numpy as np
from surface_potential import sample_residue_sas_potential


class ElectrostaticsError(ValueError):
    """An explicit electrostatic run could not be completed reliably."""


def availability() -> dict[str, Any]:
    return {
        "pdb2pqr": shutil.which("pdb2pqr") is not None,
        "apbs": shutil.which("apbs") is not None,
        "propka_via_pdb2pqr": True,
        "policy": "Requested electrostatics fail explicitly; no approximate substitute",
        "default_solver": "npbe", "solver_modes": ["npbe", "lpbe", "both"],
    }


def _run(args: list[str], cwd: Path, timeout: int = 240) -> str:
    try:
        result = subprocess.run(args, cwd=str(cwd), capture_output=True,
                                text=True, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ElectrostaticsError(f"External command failed: {Path(args[0]).name}: {exc}") from exc
    if result.returncode:
        snippet = (result.stderr or result.stdout)[-1500:].strip()
        raise ElectrostaticsError(
            f"{Path(args[0]).name} exited with code {result.returncode}: {snippet}"
        )
    return (result.stdout or "") + "\n" + (result.stderr or "")


_RESIDUE_ALIASES = {
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS",
    "CYM": "CYS", "CYX": "CYS", "TYM": "TYR",
}


def _pqr_atoms(path: Path) -> list[dict[str, Any]]:
    """Read whitespace PQR (including optional chain and insertion code)."""
    atoms = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM ", "HETATM ")):
            continue
        tokens = line.split()
        if len(tokens) < 10:
            raise ElectrostaticsError("Malformed PQR ATOM record")
        try:
            x, y, z, charge, radius = map(float, tokens[-5:])
        except ValueError as exc:
            raise ElectrostaticsError("PQR atomic coordinates/charges are invalid") from exc
        prefix = tokens[2:-5]
        if len(prefix) < 3:
            raise ElectrostaticsError("PQR record lacks residue identity")
        atom_name, res_name = prefix[:2]
        if len(prefix) >= 4:
            chain, resid = prefix[2], prefix[3]
        else:
            # PDB2PQR can concatenate a one-character chain ID and a four-
            # digit residue number, e.g. A1368, instead of separating them
            # with whitespace. Retain that chain for strict PDB/PQR matching.
            compact = re.fullmatch(r"([A-Za-z])(-?\d+)([A-Za-z]?)", prefix[2])
            if compact:
                chain = compact.group(1)
                resid = compact.group(2) + compact.group(3)
            else:
                chain, resid = "", prefix[2]
        m = re.fullmatch(r"(-?\d+)([A-Za-z]?)", resid)
        if not m:
            raise ElectrostaticsError(f"Cannot parse PQR residue index: {resid}")
        if not all(math.isfinite(v) for v in (x, y, z, charge, radius)) or radius < 0 or (radius == 0 and not atom_name.upper().startswith("H")):
            raise ElectrostaticsError(f"PQR has nonfinite values or invalid radius: atom={atom_name} residue={res_name} charge={charge} radius={radius}; record={line[:180]!r}")
        atoms.append({"name": atom_name, "res_name": res_name, "chain": chain,
                      "seq": int(m.group(1)), "icode": m.group(2),
                      "xyz": (x, y, z), "charge": charge, "radius": radius})
    if not atoms:
        raise ElectrostaticsError("PDB2PQR produced no atoms")
    return atoms


def _apbs_input(atoms: list[dict[str, Any]], ionic_mM: float,
                temp_K: float, solver: str = "lpbe") -> str:
    # Keep the helper default for archived LPBE audit scripts. Production compute passes its solver explicitly.
    if solver not in ("lpbe", "npbe"):
        raise ElectrostaticsError("Unsupported PB solver")
    coords = np.asarray([a["xyz"] for a in atoms], dtype=float)
    low, high = coords.min(axis=0), coords.max(axis=0)
    span = np.maximum(high - low, 1.0)
    fine = span + 20.0  # >=10 A padding on every side
    coarse = span + 40.0
    # APBS grid dimensions must be 32*n+1. Cap resolution/memory explicitly.
    dims = [min(257, max(65, 32 * math.ceil((v / 0.8 - 1) / 32) + 1))
            for v in fine]
    salt = float(ionic_mM) / 1000.0
    if salt < 0 or not math.isfinite(salt) or not (0 < temp_K < 1000):
        raise ElectrostaticsError("Invalid ionic strength or temperature")
    salt_lines = (
        f"    ion charge 1 conc {salt:.9g} radius 2.0\n"
        f"    ion charge -1 conc {salt:.9g} radius 2.0\n"
    ) if salt > 0 else ""
    return (
        "read\n    mol pqr structure.pqr\nend\n"
        "elec name protein_potential\n    mg-auto\n"
        f"    dime {' '.join(str(d) for d in dims)}\n"
        f"    cglen {' '.join(f'{v:.4f}' for v in coarse)}\n"
        f"    fglen {' '.join(f'{v:.4f}' for v in fine)}\n"
        "    cgcent mol 1\n    fgcent mol 1\n    mol 1\n"
        f"    {solver}\n    bcfl sdh\n    pdie 2.0\n    sdie 78.54\n"
        "    chgm spl2\n    srfm smol\n    srad 1.4\n    swin 0.3\n"
        "    sdens 10.0\n"
        + salt_lines
        + f"    temp {float(temp_K):.4f}\n"
        + "    calcenergy no\n    calcforce no\n"
        + "    write pot dx potential\nend\nquit\n"
    )


def _read_dx(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = path.read_text(errors="replace").splitlines()
    counts = origin = None
    deltas = []
    values = []
    inside = False
    expected = None
    for line in lines:
        parts = line.split()
        if "gridpositions counts" in line:
            counts = tuple(map(int, parts[-3:]))
        elif line.startswith("origin "):
            origin = np.asarray(list(map(float, parts[1:4])), dtype=float)
        elif line.startswith("delta "):
            deltas.append(list(map(float, parts[1:4])))
        elif "data follows" in line:
            m = re.search(r"items\s+(\d+)", line)
            if not m:
                raise ElectrostaticsError("DX array missing items count")
            expected = int(m.group(1))
            inside = True
        elif inside and len(values) < expected:
            try:
                values.extend(float(v) for v in parts[:expected - len(values)])
            except ValueError as exc:
                raise ElectrostaticsError("Invalid DX potential data") from exc
    if counts is None or origin is None or len(deltas) != 3:
        raise ElectrostaticsError("DX grid metadata missing")
    if expected != math.prod(counts) or len(values) != expected:
        raise ElectrostaticsError("Truncated or malformed DX grid")
    delta = np.asarray(deltas, dtype=float)
    if not np.allclose(delta, np.diag(np.diag(delta)), atol=1e-8):
        raise ElectrostaticsError("Non-axis-aligned DX grids are unsupported")
    spacing = np.diag(delta)
    if np.any(spacing <= 0):
        raise ElectrostaticsError("Invalid DX grid spacing")
    array = np.asarray(values, dtype=float).reshape(counts)
    if not np.isfinite(array).all():
        raise ElectrostaticsError("DX potential contains nonfinite values")
    return array, origin, spacing


def _interpolate(grid: np.ndarray, origin: np.ndarray, spacing: np.ndarray,
                 xyz: tuple[float, float, float]) -> float | None:
    pos = (np.asarray(xyz, dtype=float) - origin) / spacing
    lo = np.floor(pos).astype(int)
    hi = lo + 1
    if np.any(lo < 0) or np.any(hi >= grid.shape):
        return None
    t = pos - lo
    value = 0.0
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                w = (t[0] if a else 1-t[0]) * (t[1] if b else 1-t[1]) * (t[2] if c else 1-t[2])
                value += w * grid[lo[0]+a, lo[1]+b, lo[2]+c]
    return float(value)


# Some archival PDBs include isolated amino-acid records with only an N atom.
# Never fabricate missing coordinates or modify the original 11-map PDB.
_ORPHAN_MAX_COUNT = 4
_ORPHAN_MAX_FRACTION = 0.01
_CANONICAL_PDB_RESIDUES = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def _pdb_orphan_n_only_residues(
    prepared_pdb: str, valid: dict[str, str]
) -> tuple[str, list[dict[str, Any]]]:
    """Prepare an APBS-only candidate excluding strictly single-N fragments.

    Detection alone never changes input. Caller may retry this candidate only
    if PDB2PQR fails on the unchanged original. Any residue retaining a CA,
    C, O, side-chain atom, or more than one heavy atom is not auto-removed.
    """
    by_key: dict[str, dict[str, Any]] = {}
    lines = prepared_pdb.splitlines(keepends=True)
    for line in lines:
        if not line.startswith("ATOM  ") or len(line) < 27:
            continue
        try:
            chain = line[21].strip()
            seq = int(line[22:26].strip())
            icode = line[26].strip()
        except ValueError:
            continue
        key = f"{chain}:{seq}:{icode}"
        # A single-N fragment can be absent from the canonical CA-based
        # profile. Scan raw ATOM records, not only existing map residues.
        res_name = line[17:20].strip()
        if res_name not in _CANONICAL_PDB_RESIDUES:
            continue
        name = line[12:16].strip().upper()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element in ("H", "D") or (not element and name.startswith(("H", "D"))):
            continue
        entry = by_key.setdefault(key, {"names": set(), "res_names": set(), "lines": []})
        entry["names"].add(name)
        entry["res_names"].add(line[17:20].strip())
        entry["lines"].append(line)

    omitted = []
    for key, entry in by_key.items():
        if (entry["names"] == {"N"} and
                len(entry["res_names"]) == 1 and
                (key not in valid or entry["res_names"] == {valid[key]}) and
                len(entry["lines"]) == 1):
            omitted.append({
                "key": key, "res_name": next(iter(entry["res_names"])),
                "already_in_canonical_residue_profile": key in valid,
                "observed_heavy_atoms": ["N"],
                "missing_backbone_atoms": ["CA", "C", "O"],
                "reason": "orphan single N atom; C-alpha and carbonyl group absent",
            })
    omitted.sort(key=lambda x: x["key"])
    if not omitted:
        return prepared_pdb, []
    if (len(omitted) > _ORPHAN_MAX_COUNT or
            len(omitted) / max(len(valid), 1) > _ORPHAN_MAX_FRACTION):
        raise ElectrostaticsError(
            "PDB2PQR cannot safely recover this structure automatically: "
            f"{len(omitted)} N-only fragments among {len(valid)} residues "
            f"(maximum {_ORPHAN_MAX_COUNT} and {_ORPHAN_MAX_FRACTION:.0%}). "
            "Use a better experimental structure or independently validated remodeling."
        )
    excluded_keys = {v["key"] for v in omitted}
    remaining = []
    for line in lines:
        if line.startswith("ATOM  ") and len(line) >= 27:
            try:
                key = f"{line[21].strip()}:{int(line[22:26].strip())}:{line[26].strip()}"
            except ValueError:
                key = None
            if key in excluded_keys:
                continue
        remaining.append(line)
    if sum(line.startswith("ATOM  ") for line in remaining) < 4:
        raise ElectrostaticsError("APBS recovery would leave too few protein atoms")
    return "".join(remaining), omitted


def compute(prepared_pdb: str, residue_profiles: list[dict[str, Any]], *,
            pH: float, ionic_mM: float, temp_K: float, solver: str = "npbe") -> dict[str, Any]:
    """Run external programs and project APBS potential onto side-chain atoms.

    APBS writes potentials in kT/e. Surface-sidechain mean is descriptive, not
    adsorption free energy. Unmatched structures cause an explicit failure.
    """
    if solver not in ("npbe", "lpbe", "both"):
        raise ElectrostaticsError("Unsupported PB solver mode")
    tools = availability()
    if not tools["pdb2pqr"] or not tools["apbs"]:
        missing = ", ".join(n for n in ("pdb2pqr", "apbs") if not tools[n])
        raise ElectrostaticsError(
            f"Electrostatics requested but missing executable(s): {missing}. "
            "Install APBS and PDB2PQR with PROPKA support, then retry."
        )
    valid = {str(r["key"]): str(r["res_name"]) for r in residue_profiles}
    with tempfile.TemporaryDirectory(prefix="interfacescout_apbs_") as directory:
        wd = Path(directory)
        pdb_path = wd / "structure.pdb"
        pdb_path.write_text(prepared_pdb)
        command = [
            shutil.which("pdb2pqr"), "--ff=AMBER", "--keep-chain",
            "--titration-state-method=propka", "--whitespace", f"--with-ph={float(pH):.4f}",
            "structure.pdb", "structure.pqr",
        ]
        omitted = []
        original_pqr_failure = None
        try:
            _run(command, wd)
        except ElectrostaticsError as original_exc:
            original_pqr_failure = str(original_exc)
            cleaned, omitted = _pdb_orphan_n_only_residues(prepared_pdb, valid)
            if not omitted:
                raise ElectrostaticsError(
                    original_pqr_failure
                    + " No safely removable single-N fragment found; "
                      "automatic backbone reconstruction is not performed."
                ) from original_exc
            pdb_path.write_text(cleaned)
            (wd / "structure.pqr").unlink(missing_ok=True)
            try:
                _run(command, wd)
            except ElectrostaticsError as retry_exc:
                raise ElectrostaticsError(
                    "Original PDB2PQR failed; APBS-only omission of "
                    f"{[x['key'] for x in omitted]} also failed. "
                    "Original maps remain unchanged. "
                    f"Original: {original_pqr_failure[-500:]}; "
                    f"retry: {str(retry_exc)[-500:]}"
                ) from retry_exc
        pqr_path = wd / "structure.pqr"
        if not pqr_path.is_file():
            raise ElectrostaticsError("PDB2PQR did not produce structure.pqr")
        atoms = _pqr_atoms(pqr_path)
        pqr_sha = hashlib.sha256(pqr_path.read_bytes()).hexdigest()
        primary = "npbe" if solver == "both" else solver
        result = _solve_prepared(wd, atoms, residue_profiles, valid, omitted,
                                 original_pqr_failure, pH, ionic_mM, temp_K, primary, pqr_sha)
        result["requested_solver"] = solver
        if solver == "both":
            try:
                other = _solve_prepared(wd, atoms, residue_profiles, valid, omitted,
                                        original_pqr_failure, pH, ionic_mM, temp_K, "lpbe", pqr_sha)
                result["comparison"] = _paired_comparison(result, other)
            except (ElectrostaticsError, ValueError) as exc:
                result["comparison"] = {"status": "unavailable", "secondary_solver": "lpbe",
                                        "reason": str(exc), "fallback_used": False}
        return result


def _quality_assessment(log: str, *, solver: str, surface_meta: dict, records: list,
                        omitted: list, spacing: np.ndarray) -> dict[str, Any]:
    """Auditable numerical flags; not a claim of independent physical validation."""
    flags = []
    evidence = []
    # APBS can exit zero even when nonlinear overflow or focus warnings are emitted.
    # Do not flag generic asc_getToken EOF messages, common in successful APBS runs.
    for line in log.splitlines():
        if re.search(r"trapped exp overflows|exponential overflow|convergence not guaranteed|focusFillBound.*warning|failed to converge|did not converge|non.?finite|numerical divergence", line, re.I):
            evidence.append(line.strip()[:400])
    if evidence:
        flags.append("solver_diagnostic_review")
    total = int(surface_meta.get("n_exposed_residues_requested", 0))
    n = int(surface_meta.get("n_exposed_residues_with_surface_potential", 0))
    coverage = n / total if total else 0.
    if not total or coverage < .90:
        flags.append("low_surface_coverage")
    if omitted:
        flags.append("partial_input_recovery")
    if max(spacing) > .8001:
        flags.append("grid_spacing_above_target")
    values = [float(r["surface_potential_kT_per_e"]) for r in records
              if r.get("surface_potential_kT_per_e") is not None]
    if not values or not all(math.isfinite(v) for v in values):
        raise ElectrostaticsError("No finite solvent-accessible residue potentials")
    high = sum(abs(v) > 10 for v in values)
    if high:
        flags.append("high_absolute_potential")
        if solver == "lpbe":
            flags.append("review_linearization_applicability")
    return {"status": "review_recommended" if flags else "checks_passed",
            "flags": flags, "solver_log_evidence": evidence[:50],
            "surface_residue_coverage_fraction": coverage,
            "n_abs_potential_gt10": high, "abs_potential_review_threshold_kT_per_e": 10,
            "threshold_is_review_flag_not_validity_cutoff": True,
            "scope": "Output and numerical diagnostics only; independent grid convergence is not certified."}


def _paired_comparison(primary: dict, other: dict) -> dict[str, Any]:
    from scipy.stats import spearmanr
    if primary["pqr_sha256"] != other["pqr_sha256"] or primary["grid"] != other["grid"]:
        raise ElectrostaticsError("Paired PB inputs or grids differ")
    x = {r["key"]: r for r in primary["residues"]}
    y = {r["key"]: r for r in other["residues"]}
    if x.keys() != y.keys():
        raise ElectrostaticsError("Paired PB residue keys differ")
    a, b = [], []
    for key, r in x.items():
        q = y[key]
        if (r["surface_potential_samples"] != q["surface_potential_samples"] or
                r["surface_sampling_area_A2"] != q["surface_sampling_area_A2"]):
            raise ElectrostaticsError("Paired PB shell sampling differs")
        v, u = r["surface_potential_kT_per_e"], q["surface_potential_kT_per_e"]
        if (v is None) != (u is None):
            raise ElectrostaticsError("Paired PB potential availability differs")
        r["lpbe_surface_potential_kT_per_e"] = u
        r["lpbe_sidechain_mean_potential_kT_per_e"] = q["sidechain_mean_potential_kT_per_e"]
        r["npbe_minus_lpbe_kT_per_e"] = None if v is None else v - u
        if v is not None:
            a.append(v); b.append(u)
    if not a:
        raise ElectrostaticsError("No paired solvent-accessible potentials")
    a, b = np.asarray(a), np.asarray(b)
    rho = float(spearmanr(a, b).statistic) if len(a) > 1 and np.ptp(a) and np.ptp(b) else float("nan")
    return {"status": "computed", "primary_solver": "npbe", "secondary_solver": "lpbe",
            "same_pqr_grid_shell": True, "n_paired": len(a),
            "spearman": rho if math.isfinite(rho) else None,
            "median_abs_difference_kT_per_e": float(np.median(abs(a-b))),
            "p95_abs_difference_kT_per_e": float(np.quantile(abs(a-b), .95)),
            "sign_disagreement_fraction": float(np.mean(np.sign(a) != np.sign(b))),
            "lpbe_quality": other["quality"], "lpbe_elapsed_seconds": other["solver_elapsed_seconds"],
            "fallback_used": False}


def _solve_prepared(wd, atoms, residue_profiles, valid, omitted, original_pqr_failure,
                    pH, ionic_mM, temp_K, solver, pqr_sha):
    for previous in wd.glob("potential*.dx"):
        previous.unlink()
    (wd / "potential.in").write_text(_apbs_input(atoms, ionic_mM, temp_K, solver=solver))
    started = time.monotonic()
    apbs_log = _run([shutil.which("apbs"), "potential.in"], wd, timeout=2400)
    elapsed = time.monotonic() - started
    dx_path = wd / "potential.dx"
    if not dx_path.is_file():
        candidates = sorted(wd.glob("potential*.dx"))
        if len(candidates) == 1:
            dx_path = candidates[0]
        else:
            produced = ", ".join(sorted(p.name for p in wd.iterdir()))
            raise ElectrostaticsError(
                f"APBS did not write potential.dx (other grids={len(candidates)}). "
                f"Produced files: {produced}; APBS log tail: {apbs_log[-2500:]}"
            )
    grid, origin, spacing = _read_dx(dx_path)
    grouped: dict[str, list[dict]] = {}
    for atom in atoms:
        key = f'{atom["chain"]}:{atom["seq"]}:{atom["icode"]}'
        if key in valid and valid[key] == _RESIDUE_ALIASES.get(atom["res_name"], atom["res_name"]):
            grouped.setdefault(key, []).append(atom)
    omitted_keys = {x["key"] for x in omitted}
    expected_keys = set(valid) - omitted_keys
    if not grouped or set(grouped) != expected_keys:
        extra = sorted(set(grouped) - expected_keys)
        missing = sorted(expected_keys - set(grouped))
        raise ElectrostaticsError(
            f"PQR/PDB residue mapping incomplete: {len(grouped)}/{len(valid)}; "
            f"missing expected={missing[:12]}, extra={extra[:12]}. "
            "Only explicitly recorded N-only fragments may be absent."
        )
    # The principal potential descriptor is sampled *outside* the
    # molecule on an area-weighted probe-center accessible shell.
    # The legacy atom-coordinate estimator remains a separate field
    # for explicit comparison and never influences the 11 maps.
    exposed_keys = {
        str(row["key"]) for row in residue_profiles
        if row.get("surface_exposed", False) and str(row["key"]) in grouped
    }
    surface_by_residue, surface_meta = sample_residue_sas_potential(
        atoms, grid, origin, spacing, exposed_keys
    )
    backbone = {"N", "CA", "C", "O", "OXT"}
    records = []
    for row in residue_profiles:
        key = str(row["key"])
        group = grouped.get(key)
        if group is None:
            if key not in omitted_keys:
                raise ElectrostaticsError(f"Unaccounted PQR mapping gap: {key}")
            records.append({
                "key": key, "res_name": row["res_name"],
                "net_atomic_charge_e": None,
                "surface_potential_kT_per_e": None,
                "surface_sampling_area_A2": None,
                "surface_potential_samples": 0,
                "sidechain_mean_potential_kT_per_e": None,
                "sidechain_potential_samples": 0,
                "surface_exposed": bool(row.get("surface_exposed", False)),
                "pqr_mapped": False,
                "pqr_omission_reason": "missing backbone (orphan N-only residue)",
            })
            continue
        sidechain = [a for a in group if a["name"].upper() not in backbone
                     and not a["name"].upper().startswith("H")]
        samples = [_interpolate(grid, origin, spacing, a["xyz"]) for a in sidechain]
        samples = [v for v in samples if v is not None]
        if sidechain and len(samples) != len(sidechain):
            raise ElectrostaticsError(f"APBS grid excludes side-chain atom(s) at {key}")
        sas = surface_by_residue.get(key, {})
        records.append({
            "key": key, "res_name": row["res_name"],
            "net_atomic_charge_e": round(sum(a["charge"] for a in group), 6),
            "surface_potential_kT_per_e": sas.get("surface_potential_kT_per_e"),
            "surface_sampling_area_A2": sas.get("surface_sampling_area_A2"),
            "surface_potential_samples": sas.get("surface_potential_samples", 0),
            "sidechain_mean_potential_kT_per_e":
                round(sum(samples) / len(samples), 6) if samples else None,
            "sidechain_potential_samples": len(samples),
            "surface_exposed": bool(row.get("surface_exposed", False)),
            "pqr_mapped": True,
            "pqr_omission_reason": None,
        })
    quality = _quality_assessment(apbs_log, solver=solver, surface_meta=surface_meta,
                                  records=records, omitted=omitted, spacing=spacing)
    return {
        "status": "computed", "method": "PDB2PQR PROPKA protonation + APBS " + ("nonlinear" if solver == "npbe" else "linearized") + " Poisson-Boltzmann",
        "solver": solver, "pqr_sha256": pqr_sha, "solver_elapsed_seconds": elapsed,
        "grid": {"shape": list(grid.shape), "origin_A": origin.tolist(), "spacing_A": spacing.tolist()},
        "quality": quality, "fallback_used": False,
        "potential_unit": "kT/e", "pH": float(pH),
        "ionic_mM": float(ionic_mM), "temperature_K": float(temp_K),
        "n_pqr_atoms": len(atoms),
        "n_matched_residues": len(grouped),
        "n_input_residues": len(records),
        "n_unmatched_residues": len(set(valid) & omitted_keys),
        "n_orphan_records_removed_from_apbs_input": len(omitted),
        "quality_status": "partial_input_recovery" if omitted else "complete",
        "input_recovery": {
            "applied": bool(omitted),
            "policy": "APBS-only removal of strictly single-N orphan records after original PDB2PQR failure; no coordinates synthesized",
            "omitted_residues": omitted,
            "n_omitted": len(omitted),
            "original_PDB2PQR_failure": original_pqr_failure[-700:] if omitted else None,
            "canonical_maps_unchanged": True,
            "caveat": ("Real APBS is computed on a minimally altered PQR; "
                       "omitted residues have undefined potentials, and "
                       "nearby electrostatic values may be uncertain.")
                      if omitted else None,
        },
        "surface_sampling": surface_meta,
        "n_surface_potential_residues": surface_meta["n_exposed_residues_with_surface_potential"],
        "residues": records,
        "interpretation": "Primary surface_potential_kT_per_e is area-weighted APBS potential "
                          "on a discrete solvent-accessible shell (+0.25 A outside probe-center SAS); "
                          "sidechain_mean_potential_kT_per_e is a legacy atom-coordinate comparator. "
                          "Neither is affinity, binding energy, nor a recalibrated canonical map score.",
    }
