"""Optional, physically grounded PDB2PQR/PROPKA -> APBS electrostatics.

No silent fallback. Canonical 11-map scores and Pareto ranks are unchanged:
electrostatics is a separate, explicitly requested, residue-matched layer.
"""
from __future__ import annotations

import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

import numpy as np


class ElectrostaticsError(ValueError):
    """An explicit electrostatic run could not be completed reliably."""


def availability() -> dict[str, Any]:
    return {
        "pdb2pqr": shutil.which("pdb2pqr") is not None,
        "apbs": shutil.which("apbs") is not None,
        "propka_via_pdb2pqr": True,
        "policy": "Requested electrostatics fail explicitly; no approximate substitute",
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
                temp_K: float) -> str:
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
        "    lpbe\n    bcfl sdh\n    pdie 2.0\n    sdie 78.54\n"
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


def compute(prepared_pdb: str, residue_profiles: list[dict[str, Any]], *,
            pH: float, ionic_mM: float, temp_K: float) -> dict[str, Any]:
    """Run external programs and project APBS potential onto side-chain atoms.

    APBS writes potentials in kT/e. Surface-sidechain mean is descriptive, not
    adsorption free energy. Unmatched structures cause an explicit failure.
    """
    tools = availability()
    if not tools["pdb2pqr"] or not tools["apbs"]:
        missing = ", ".join(n for n in ("pdb2pqr", "apbs") if not tools[n])
        raise ElectrostaticsError(
            f"Electrostatics requested but missing executable(s): {missing}. "
            "Install APBS and PDB2PQR with PROPKA support, then retry."
        )
    with tempfile.TemporaryDirectory(prefix="interfacescout_apbs_") as directory:
        wd = Path(directory)
        (wd / "structure.pdb").write_text(prepared_pdb)
        command = [
            shutil.which("pdb2pqr"), "--ff=AMBER", "--keep-chain",
            "--titration-state-method=propka", "--whitespace", f"--with-ph={float(pH):.4f}",
            "structure.pdb", "structure.pqr",
        ]
        _run(command, wd)
        pqr_path = wd / "structure.pqr"
        if not pqr_path.is_file():
            raise ElectrostaticsError("PDB2PQR did not produce structure.pqr")
        atoms = _pqr_atoms(pqr_path)
        (wd / "potential.in").write_text(_apbs_input(atoms, ionic_mM, temp_K))
        apbs_log = _run([shutil.which("apbs"), "potential.in"], wd, timeout=600)
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
        valid = {str(r["key"]): str(r["res_name"]) for r in residue_profiles}
        grouped: dict[str, list[dict]] = {}
        for atom in atoms:
            key = f'{atom["chain"]}:{atom["seq"]}:{atom["icode"]}'
            if key in valid and valid[key] == _RESIDUE_ALIASES.get(atom["res_name"], atom["res_name"]):
                grouped.setdefault(key, []).append(atom)
        if not grouped or len(grouped) < len(valid):
            raise ElectrostaticsError(
                f"PQR/PDB residue mapping incomplete: {len(grouped)}/{len(valid)}. "
                "Check chain IDs, insertion codes and modified/missing residues."
            )
        backbone = {"N", "CA", "C", "O", "OXT"}
        records = []
        for row in residue_profiles:
            key = str(row["key"])
            group = grouped[key]
            sidechain = [a for a in group if a["name"].upper() not in backbone
                         and not a["name"].upper().startswith("H")]
            samples = [_interpolate(grid, origin, spacing, a["xyz"]) for a in sidechain]
            samples = [v for v in samples if v is not None]
            if sidechain and len(samples) != len(sidechain):
                raise ElectrostaticsError(f"APBS grid excludes side-chain atom(s) at {key}")
            records.append({
                "key": key, "res_name": row["res_name"],
                "net_atomic_charge_e": round(sum(a["charge"] for a in group), 6),
                "sidechain_mean_potential_kT_per_e":
                    round(sum(samples) / len(samples), 6) if samples else None,
                "sidechain_potential_samples": len(samples),
                "surface_exposed": bool(row.get("surface_exposed", False)),
            })
        return {
            "status": "computed", "method": "PDB2PQR PROPKA protonation + APBS linearized Poisson-Boltzmann",
            "potential_unit": "kT/e", "pH": float(pH),
            "ionic_mM": float(ionic_mM), "temperature_K": float(temp_K),
            "n_pqr_atoms": len(atoms), "n_matched_residues": len(records),
            "residues": records,
            "interpretation": "Side-chain sampled potential and PQR net atomic charge; "
                              "not affinity, binding energy, or a recalibrated 11-map score.",
        }
