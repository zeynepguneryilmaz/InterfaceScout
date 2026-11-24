"""Independent APBS solvent-accessible-shell potential audit for frozen TV48.

Experimental-science scope:
- One real PDB2PQR/PROPKA/APBS LPBE calculation per TV48 protein at pH
  7.4, 150 mM and 298 K.
- Sample a deterministic 1.4-A probe-center accessible shell (outward
  displacement +0.25 A) using PQR heavy-atom radii and 42 Fibonacci rays
  per heavy atom. A cKDTree excludes occluded rays. Surface-area-weighted
  residue means are independent of number of atoms and of mesh density.
- Eight prespecified size-stratified proteins have actual grid, salt and
  pH sensitivity checks; no grid-sensitivity claim if finer and coarse
  runs share effective spacing because of the explicit grid cap.
- Original side-chain-atom potential and all 11 maps remain unchanged.
- This is a technical continuum-model study, NOT adsorption validation.

Run: --panel-shard i --shards 8 | --sensitivity-shard i --shards 8
     --aggregate --shards 8
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
import traceback
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "backend"), str(ROOT / "validation")]
from electrostatics import _pqr_atoms, _apbs_input, _read_dx, _run
import strict_apo_compat_48 as prep
from interface_engine import analyze_all_maps
from model_settings import MODEL_VERSION
from export_bundle import MAP_NAMES

MANIFEST = list(csv.DictReader((ROOT / "validation" / "strict_apo_48_selected_candidates_v2_final.csv")
                                .open(encoding="utf-8", newline="")))
assert len(MANIFEST) == 48
OUT = Path("tv48_apbs_sas")
OUT.mkdir(exist_ok=True)

PROBE = 1.4
OFFSET = 0.25
N_DIRECTIONS = 42
SALT_REF = 150.0
TEMP = 298.0
PH_REF = 7.4


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def directions(n=N_DIRECTIONS):
    # Deterministic Fibonacci unit sphere sampling; no stochastic point loss.
    i = np.arange(n, dtype=float)
    z = 1 - 2 * (i + .5) / n
    theta = i * (math.pi * (3 - math.sqrt(5)))
    xy = np.sqrt(1 - z*z)
    return np.stack((np.cos(theta)*xy, np.sin(theta)*xy, z), axis=-1)


def build_sas_points(atoms, valid_keys, directions_count=N_DIRECTIONS):
    """Return (coordinates, residue keys, area-weights) on expanded vdW union.

    A sampled point is accepted only if outside every OTHER atom's
    expanded sphere. Hydrogens are excluded from the geometry.
    The outward 0.25-A shift places samples in a solvent-side shell.
    """
    heavy = [a for a in atoms if a["radius"] >= 0.7 and
             not re.match(r"^[0-9]*H", a["name"].upper())]
    if len(heavy) < 10:
        raise ValueError("PQR has too few heavy atoms")
    centers = np.asarray([a["xyz"] for a in heavy], dtype=np.float64)
    radii = np.asarray([a["radius"] + PROBE + OFFSET for a in heavy], dtype=np.float64)
    if not np.isfinite(centers).all() or not np.isfinite(radii).all():
        raise ValueError("Invalid PQR geometry")
    keys = [f'{a["chain"]}:{a["seq"]}:{a["icode"]}' for a in heavy]
    tree = cKDTree(centers)
    angular = directions(directions_count)
    points = []
    pointkeys = []
    pointweights = []
    maxrad = float(max(radii))
    for atom_index, (center, radius, key) in enumerate(zip(centers, radii, keys)):
        if key not in valid_keys:
            continue
        candidate = center + radius * angular
        neighbors = tree.query_ball_point(candidate, r=maxrad + 1e-5)
        weights = 4.0 * math.pi * radius * radius / directions_count
        for j, p in enumerate(candidate):
            idxs = [idx for idx in neighbors[j] if idx != atom_index]
            if idxs:
                delta = centers[idxs] - p
                squared = np.einsum("ij,ij->i", delta, delta)
                if np.any(squared < np.square(radii[idxs] - 1e-6)):
                    continue
            points.append(p)
            pointkeys.append(key)
            pointweights.append(weights)
    if not points:
        raise ValueError("No exposed solvent-accessible shell points")
    return (np.asarray(points, dtype=np.float64), pointkeys,
            np.asarray(pointweights, dtype=np.float64),
            {"n_heavy_atoms": len(heavy), "n_shell_points": len(points),
             "area_weighted_shell_A2": round(float(np.sum(pointweights)), 3),
             "directions_per_heavy_atom": directions_count,
             "probe_A": PROBE, "outward_offset_A": OFFSET})


def sample_grid(grid_meta, points):
    grid, origin, steps = grid_meta
    frac = (points - origin) / steps
    inside = np.all((frac >= 0) & (frac < (np.asarray(grid.shape) - 1)), axis=1)
    if not inside.all():
        raise ValueError(f"Grid does not cover shell: {int((~inside).sum())}/{len(inside)} points outside")
    sampled = map_coordinates(grid, frac.T, order=1, mode="constant",
                              cval=np.nan, prefilter=False)
    if not np.isfinite(sampled).all():
        raise ValueError("Nonfinite interpolated APBS shell potential")
    return sampled


def summarize_by_residue(keys, weights, values, permitted):
    sums = defaultdict(lambda: [0., 0., 0.])
    for key, area, phi in zip(keys, weights, values):
        if key not in permitted:
            continue
        record = sums[key]
        record[0] += float(area)
        record[1] += float(area) * float(phi)
        record[2] += 1
    return {key: {"surface_area_A2": round(v[0], 4),
                  "SAS_potential_kT_per_e": round(v[1] / v[0], 7),
                  "n_surface_samples": int(v[2])}
            for key, v in sums.items() if v[0] > 0}


def prepare_case(item):
    pid = item["pdb_id"].upper()
    prep.prepare(pid)
    cif = prep.OUT / f"{pid}-assembly1.cif"
    h = hashlib.sha256(cif.read_bytes()).hexdigest()
    if h != item["assembly_sha256"]:
        raise ValueError("Frozen source assembly SHA mismatch")
    pdb = (prep.OUT / f"{pid}_prepared.pdb").read_text(encoding="utf-8")
    a = analyze_all_maps(pdb_text=pdb, pH=PH_REF, ionic_mM=SALT_REF, temp_K=TEMP)
    if a["version"] != MODEL_VERSION or list(a["map_order"]) != list(MAP_NAMES):
        raise ValueError("Unexpected canonical model or maps")
    return pid, pdb, a, h


def make_pqr(workdir, pdb, pH):
    (workdir / "structure.pdb").write_text(pdb, encoding="utf-8")
    command = [shutil.which("pdb2pqr"), "--ff=AMBER", "--keep-chain",
               "--titration-state-method=propka", "--whitespace",
               f"--with-ph={pH:.4f}", "structure.pdb", "structure.pqr"]
    if not command[0]:
        raise RuntimeError("Real PDB2PQR unavailable")
    _run(command, workdir, timeout=400)
    atoms = _pqr_atoms(workdir / "structure.pqr")
    return atoms


def run_apbs(workdir, atoms, ionic, target_grid=0.8, label="baseline"):
    """Run a new real APBS job with explicit and recorded effective grid step."""
    apbs = shutil.which("apbs")
    if not apbs:
        raise RuntimeError("Real APBS unavailable")
    source = _apbs_input(atoms, ionic, TEMP)
    if abs(target_grid - .8) > 1e-8:
        match = re.search(r"(?m)^    dime (\d+) (\d+) (\d+)$", source)
        if not match:
            raise RuntimeError("Unable to locate APBS dime in production input")
        locations = np.asarray([a["xyz"] for a in atoms])
        span = np.maximum(np.ptp(locations, axis=0), 1.)
        fine = span + 20.
        dims = [min(257, max(65, 32 * math.ceil((v / target_grid - 1) / 32) + 1))
                for v in fine]
        source = source[:match.start()] + "    dime " + " ".join(map(str, dims)) + source[match.end():]
    match = re.search(r"(?m)^    dime (\d+) (\d+) (\d+)$", source)
    dims = [int(q) for q in match.groups()]
    fine_match = re.search(r"(?m)^    fglen ([0-9.]+) ([0-9.]+) ([0-9.]+)$", source)
    fine = list(map(float, fine_match.groups()))
    effective = [round(f / (d - 1), 6) for f, d in zip(fine, dims)]
    stem = f"potential_{label}"
    source = source.replace("write pot dx potential", f"write pot dx {stem}")
    input_path = workdir / f"run_{label}.in"
    input_path.write_text(source, encoding="utf-8")
    apbs_log = _run([apbs, input_path.name], workdir, timeout=780)
    # APBS can suffix DX files with the named ELEC calculation, e.g.
    # 'potential_baseline-protein_potential.dx'. Do not accept a stale
    # unrelated potential from another variant in the same work directory.
    expected = workdir / (stem + ".dx")
    matches = [expected] if expected.is_file() else sorted(
        candidate for candidate in workdir.glob(stem + "-*.dx")
        if candidate.is_file()
    )
    if len(matches) != 1:
        available = sorted(path.name for path in workdir.glob("*.dx"))
        raise RuntimeError(
            f"APBS DX selection failed for {stem}: found {len(matches)} valid "
            f"candidate(s). Available DX files: {available}; "
            f"APBS log tail: {apbs_log[-1200:]}"
        )
    dx = matches[0]
    result = _read_dx(dx)
    dx.unlink()
    return result, {"grid_dims": dims, "effective_spacing_A": effective,
                    "grid_target_A": target_grid, "fine_grid_lengths_A": fine,
                    "ionic_mM": ionic, "temperature_K": TEMP,
                    "actual_APBS_execution": True}


def surface_case(item):
    pid, pdb, analysis, h = prepare_case(item)
    all_rows = analysis["residue_profiles"]
    exposed = {x["key"] for x in all_rows if x["surface_exposed"]}
    valid = {x["key"] for x in all_rows}
    with tempfile.TemporaryDirectory(prefix="tv48_apbs_sas_") as td:
        wd = Path(td)
        atoms = make_pqr(wd, pdb, PH_REF)
        present = {f'{a["chain"]}:{a["seq"]}:{a["icode"]}' for a in atoms}
        if not valid.issubset(present):
            raise ValueError("Incomplete PQR-to-InterfaceScout residue identity")
        xyz, key, weight, geometry = build_sas_points(atoms, exposed)
        grid, gridmeta = run_apbs(wd, atoms, SALT_REF)
        sampled = sample_grid(grid, xyz)
        mapped = summarize_by_residue(key, weight, sampled, exposed)
        rows = [
            {"pdb_id": pid, "key": r["key"], "res_name": r["res_name"],
             "trsa_raw": r["trsa_raw"], "scrsa_raw": r["scrsa_raw"],
             "SAS_potential_kT_per_e": mapped.get(r["key"], {}).get("SAS_potential_kT_per_e"),
             "estimated_SAS_area_A2": mapped.get(r["key"], {}).get("surface_area_A2"),
             "n_shell_samples": mapped.get(r["key"], {}).get("n_surface_samples", 0)}
            for r in all_rows if r["surface_exposed"]
        ]
        # Compare to side-chain-atom sampling from the same PQR and same APBS grid.
        # This original estimator excludes glycine (no heavy side-chain).
        grid, origin, steps = grid
        atom_vals = defaultdict(list)
        backbone = {"N", "CA", "C", "O", "OXT"}
        for a in atoms:
            k = f'{a["chain"]}:{a["seq"]}:{a["icode"]}'
            if k not in exposed or a["name"].upper() in backbone or a["radius"] < 0.7:
                continue
            frac = (np.asarray(a["xyz"]) - origin) / steps
            if np.all((frac >= 0) & (frac < np.asarray(grid.shape) - 1)):
                atom_vals[k].append(float(map_coordinates(grid, frac.reshape(3, 1),
                                                          order=1, prefilter=False)[0]))
        comparable = [(mapped[k]["SAS_potential_kT_per_e"], float(np.mean(v)))
                      for k, v in atom_vals.items() if k in mapped and v]
        if len(comparable) < 20:
            raise ValueError("Too few paired SAS and atom-sampled potentials")
        pair = np.asarray(comparable)
        surface = np.asarray([x["SAS_potential_kT_per_e"] for x in rows
                              if x["SAS_potential_kT_per_e"] is not None])
        summary = {
            "pdb_id": pid, "status": "PASS", "source_sha256": h,
            "n_protein_residues": len(all_rows), "n_exposed_trsa_residues": len(exposed),
            "n_exposed_with_SAS_shell": len(mapped),
            "n_missing_SAS_shell": len(exposed) - len(mapped),
            "n_gly_exposed_with_SAS": sum(
                r["res_name"] == "GLY" and r["key"] in mapped for r in all_rows if r["surface_exposed"]),
            "n_compared_atom_sampling": len(comparable),
            "mean_abs_SAS_potential_kT_per_e": round(float(np.mean(np.abs(surface))), 7),
            "mean_abs_atom_sampling_kT_per_e": round(float(np.mean(np.abs(pair[:, 1]))), 7),
            "median_abs_difference_SAS_vs_atom_kT_per_e":
                round(float(np.median(np.abs(pair[:, 0] - pair[:, 1]))), 7),
            "SAS_vs_atom_spearman": round(float(spearmanr(pair[:, 0], pair[:, 1]).statistic), 7),
            "n_maps_unchanged": len(MAP_NAMES), "geometry": geometry, "grid": gridmeta,
            "method_note": "Discrete SAS-probe sphere approximate shell, not analytic SES; no APBS reranking."
        }
        write_csv(OUT / "per_protein" / f"{pid}_SAS_residues.csv", rows,
                  ["pdb_id", "key", "res_name", "trsa_raw", "scrsa_raw",
                   "SAS_potential_kT_per_e", "estimated_SAS_area_A2", "n_shell_samples"])
        write_json(OUT / "per_protein" / f"{pid}_SAS_summary.json", summary)
    return summary


# Frozen BEFORE running the first sensitivity job: eight evenly spaced size
# quantiles, rather than choosing pleasing/APBS-favorable protein results.
SENSITIVITY_SET = [
    MANIFEST[MANIFEST.index(x)] for x in
    [sorted(MANIFEST, key=lambda v: (int(v["total_residues"]), v["pdb_id"]))[i]
     for i in (0, 7, 14, 21, 28, 35, 42, 47)]
]


def sensitivity_case(item):
    pid, pdb, analysis, h = prepare_case(item)
    exposed = {r["key"] for r in analysis["residue_profiles"] if r["surface_exposed"]}
    variants = [
        ("baseline", PH_REF, SALT_REF, .8),
        ("grid_finer", PH_REF, SALT_REF, .6),
        ("grid_coarser", PH_REF, SALT_REF, 1.1),
        ("ionic_50", PH_REF, 50., .8),
        ("ionic_300", PH_REF, 300., .8),
        ("ph_65", 6.5, SALT_REF, .8),
        ("ph_85", 8.5, SALT_REF, .8),
    ]
    rows = []
    baseline_values = None
    baseline_grid = None
    with tempfile.TemporaryDirectory(prefix="tv48_apbs_sensitivity_") as td:
        wd = Path(td)
        for name, ph, ionic, target in variants:
            atoms = make_pqr(wd, pdb, ph)
            xyz, keys, weights, geometry = build_sas_points(atoms, exposed)
            grid, metadata = run_apbs(wd, atoms, ionic, target, label=name)
            values = summarize_by_residue(keys, weights, sample_grid(grid, xyz), exposed)
            if name == "baseline":
                baseline_values = values
                baseline_grid = metadata
            shared = sorted(set(baseline_values) & set(values))
            if len(shared) < 20:
                raise ValueError("Insufficient paired SAS residue coverage")
            ref = np.asarray([baseline_values[k]["SAS_potential_kT_per_e"] for k in shared])
            v = np.asarray([values[k]["SAS_potential_kT_per_e"] for k in shared])
            corr = float(spearmanr(ref, v).statistic) if np.ptp(ref) > 1e-9 and np.ptp(v) > 1e-9 else None
            row = {
                "pdb_id": pid, "variant": name, "pH": ph, "ionic_mM": ionic,
                "target_spacing_A": target, "effective_spacing_A": metadata["effective_spacing_A"],
                "grid_dims": metadata["grid_dims"],
                "n_shared_residues": len(shared),
                "fraction_shared_vs_baseline": round(len(shared) / len(baseline_values), 6),
                "median_abs_potential_change_kT_per_e": round(float(np.median(np.abs(v-ref))), 6),
                "mean_abs_potential_change_kT_per_e": round(float(np.mean(np.abs(v-ref))), 6),
                "spearman_vs_baseline": round(corr, 6) if corr is not None and np.isfinite(corr) else None,
                "identical_effective_grid":
                    metadata["grid_dims"] == baseline_grid["grid_dims"]
                    and metadata["fine_grid_lengths_A"] == baseline_grid["fine_grid_lengths_A"],
                "n_exposed_SAS_shell": len(values),
            }
            rows.append(row)
            print("APBS_SAS_SENSITIVITY " + json.dumps({
                "protein": pid, "variant": name, "median_abs_change": row["median_abs_potential_change_kT_per_e"],
                "n": len(shared)}), flush=True)
    write_csv(OUT / "sensitivity" / f"{pid}_conditions.csv", rows,
              list(rows[0]))
    return {"pdb_id": pid, "status": "PASS", "variants": len(rows), "assembly_sha256": h}


def run_shard(kind, shard, shards):
    items = MANIFEST if kind == "panel" else SENSITIVITY_SET
    report = []
    for i, item in enumerate(items):
        if i % shards != shard:
            continue
        pid = item["pdb_id"]
        try:
            result = surface_case(item) if kind == "panel" else sensitivity_case(item)
            state = "PASS"
        except Exception as exc:
            state = "FAIL"
            result = {"pdb_id": pid, "status": "FAIL", "error": str(exc),
                      "traceback_tail": traceback.format_exc()[-5000:]}
            write_json(OUT / "failures" / f"{kind}_{pid}.json", result)
        report.append(result)
        print("APBS_SAS_CASE " + json.dumps(
            {"mode": kind, "pdb_id": pid, "status": state,
             "error": result.get("error", "")[:250]}), flush=True)
    write_json(OUT / f"{kind}_shard_{shard}.json", report)
    # Fail the individual CI shard while retaining diagnostics and artifacts.
    # The aggregate will additionally report exact missing/failed proteins.
    nfail = sum(entry["status"] != "PASS" for entry in report)
    print("APBS_SAS_SHARD " + json.dumps(
        {"kind": kind, "shard": shard, "attempted": len(report), "failed": nfail}
    ), flush=True)
    if nfail:
        raise SystemExit(f"{nfail}/{len(report)} {kind} cases failed; see shard artifacts")


def aggregate(shards):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels, sensitivities, missing = [], [], []
    for group, records in (("panel", panels), ("sensitivity", sensitivities)):
        for i in range(shards):
            p = OUT / f"{group}_shard_{i}.json"
            if not p.is_file():
                missing.append(p.name)
            else:
                records.extend(json.loads(p.read_text(encoding="utf-8")))
    if missing:
        raise ValueError("Missing shard files: " + str(missing))
    failures = [x for x in panels + sensitivities if x["status"] != "PASS"]
    if failures or len(panels) != 48 or len(sensitivities) != len(SENSITIVITY_SET):
        write_json(OUT / "ERROR_SUMMARY.json", failures)
        raise ValueError("SAS scientific audit incomplete: " + json.dumps(failures[:2]))
    if {x["pdb_id"] for x in panels} != {x["pdb_id"] for x in MANIFEST}:
        raise ValueError("SAS protein set mismatch")
    rows = []
    for item in MANIFEST:
        p = OUT / "per_protein" / f'{item["pdb_id"]}_SAS_residues.csv'
        rows.append(pd.read_csv(p))
    residue_table = pd.concat(rows, ignore_index=True)
    all_sas = residue_table[residue_table.SAS_potential_kT_per_e.notna()]
    if len(all_sas) < 7000:
        raise ValueError("Unexpectedly low surface sampling coverage")
    p_summary = pd.DataFrame([{k: v for k, v in x.items() if not isinstance(v, dict)} for x in panels])
    p_summary.to_csv(OUT / "TV48_SAS_PER_PROTEIN.csv", index=False)
    residue_table.to_csv(OUT / "TV48_SAS_RESIDUE_POTENTIALS.csv", index=False)
    sensitivity_rows = pd.concat(
        [pd.read_csv(OUT / "sensitivity" / f'{x["pdb_id"]}_conditions.csv')
         for x in sensitivities], ignore_index=True)
    sensitivity_rows.to_csv(OUT / "TV48_SAS_PILOT_SENSITIVITY.csv", index=False)
    figroot = OUT / "figures"
    figroot.mkdir(exist_ok=True)
    plt.figure(figsize=(7.5, 5))
    plt.scatter(p_summary.mean_abs_atom_sampling_kT_per_e, p_summary.mean_abs_SAS_potential_kT_per_e)
    top = max(p_summary.mean_abs_atom_sampling_kT_per_e.max(),
              p_summary.mean_abs_SAS_potential_kT_per_e.max())
    plt.plot([0, top], [0, top], ls="--")
    plt.xlabel("Mean absolute atomic-sidechain APBS potential (kT/e)")
    plt.ylabel("Mean absolute SAS shell APBS potential (kT/e)")
    plt.title("TV48: change in physical sampling location")
    plt.tight_layout()
    plt.savefig(figroot / "01_surface_vs_atom_magnitude.png", dpi=240, bbox_inches="tight")
    plt.close()
    plt.figure(figsize=(9, 5))
    data = [sensitivity_rows.loc[sensitivity_rows.variant == x,
             "median_abs_potential_change_kT_per_e"].values
            for x in ["grid_finer", "grid_coarser", "ionic_50", "ionic_300", "ph_65", "ph_85"]]
    plt.boxplot(data, tick_labels=["Grid 0.6A", "Grid 1.1A", "50 mM", "300 mM", "pH 6.5", "pH 8.5"])
    plt.ylabel("Median absolute SAS potential change (kT/e)")
    plt.title("Eight prespecified proteins: grid / solution sensitivity")
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig(figroot / "02_grid_solution_sensitivity.png", dpi=240, bbox_inches="tight")
    plt.close()
    ratio = (p_summary.n_exposed_with_SAS_shell / p_summary.n_exposed_trsa_residues).values
    summary = {
        "status": "PASS", "n_proteins": 48, "n_sensitivity_proteins": len(sensitivities),
        "n_grid_solution_variants_per_pilot": 7, "n_exposed_tRSA_residues":
            int(p_summary.n_exposed_trsa_residues.sum()),
        "n_exposed_SAS_shell_residues": len(all_sas),
        "median_protein_SAS_coverage": round(float(np.median(ratio)), 6),
        "median_abs_surface_vs_atom_difference_kT_per_e":
            round(float(p_summary.median_abs_difference_SAS_vs_atom_kT_per_e.median()), 6),
        "median_SAS_vs_atom_spearman":
            round(float(p_summary.SAS_vs_atom_spearman.median()), 6),
        "SAS_sampling": {"probe_A": PROBE, "outward_shell_offset_A": OFFSET,
                         "directions_per_heavy_atom": N_DIRECTIONS,
                         "method": "fibonacci-expanded-PQR-vdW-spheres with occlusion and area weighting",
                         "surface_type": "discrete solvent-accessible shell (not SES)"},
        "sensitivity_selection": [x["pdb_id"] for x in SENSITIVITY_SET],
        "cautions": [
            "Neither atomic potential nor sampled SAS potential is adsorption free energy.",
            "SAS shell sampled at +0.25 A from probe-center accessible surface.",
            "A fixed-charge continuum model and finite angular sampling are approximations.",
            "Actual effective grid spacing must be inspected; grid cap can defeat attempted refinement.",
            "Only eight prespecified proteins have explicit pH/ionic/grid sensitivity.",
            "This technical panel is not independent experimental adsorption validation."
        ]
    }
    write_json(OUT / "SUMMARY.json", summary)
    print("APBS_SAS_AUDIT_SUMMARY " + json.dumps(summary), flush=True)


def main():
    pa = argparse.ArgumentParser()
    pa.add_argument("--panel-shard", type=int)
    pa.add_argument("--sensitivity-shard", type=int)
    pa.add_argument("--shards", type=int, default=8)
    pa.add_argument("--aggregate", action="store_true")
    a = pa.parse_args()
    if a.aggregate:
        aggregate(a.shards)
    elif a.panel_shard is not None:
        run_shard("panel", a.panel_shard, a.shards)
    elif a.sensitivity_shard is not None:
        run_shard("sensitivity", a.sensitivity_shard, a.shards)
    else:
        pa.error("Select panel shard, sensitivity shard or aggregate")


if __name__ == "__main__":
    main()
