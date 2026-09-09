"""Geometry-matched baseline evaluation on the external experimental manifest.

B1-B4 use the same exposed-residue coordinates, 8 A non-transitive same-face
coarse-patch geometry, redundancy suppression, and top-k prediction budget.
Only the center-ranking signal changes. B5 is the full InterfaceScout V2 rank.
Experimental labels are applied only after predictions are generated.
"""
from __future__ import annotations

import json
import urllib.request
from io import StringIO
from pathlib import Path
from typing import Dict, List

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

import main as v1
from v2.chemistry_freeze import apply_publication_chemistry
from v2.geometry import build_surface_geometry, ca_distance, same_face
from v2.gnm import solve_gnm
from v2.interface_engine import analyze_interface_v2
from v2.model_settings import COARSE_PATCH_RADIUS_A, MULTISCALE_RADII_A
from v2.prepare import prepare_pdb_text
from v2.surface_modes import get_surface_mode

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "benchmark_strict.json"
PATCH_A = COARSE_PATCH_RADIUS_A
SINGLE_RADIUS_A = MULTISCALE_RADII_A[1]


def fetch_pdb(pid: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=30) as r:
        return r.read().decode("utf-8", errors="replace")


def ca_map(pdb_text: str) -> Dict[str, np.ndarray]:
    structure = PDBParser(QUIET=True).get_structure("x", StringIO(pdb_text))
    model = next(structure.get_models())
    out = {}
    for chain in model:
        for res in chain:
            if is_aa(res, standard=True) and "CA" in res:
                out[f"{chain.id}:{int(res.id[1])}:{str(res.id[2]).strip()}"] = np.asarray(res["CA"].coord, float)
    return out


def gt_keys(coords: Dict[str, np.ndarray], case: dict) -> List[str]:
    out = []
    exact = {int(x) for x in case.get("gt_exact", [])}
    ranges = [(int(a), int(b)) for a, b in case.get("gt_ranges", [])]
    for key in coords:
        n = int(key.split(":")[1])
        if case["gt_kind"] == "exact" and n in exact:
            out.append(key)
        elif case["gt_kind"] == "ranges" and any(a <= n <= b for a, b in ranges):
            out.append(key)
    return sorted(out)


def suppress(scores: Dict[str, float], geometry: dict) -> List[str]:
    ordered = sorted(scores, key=lambda k: (-float(scores[k]), str(k)))
    kept = []
    for key in ordered:
        if float(scores[key]) <= 0:
            continue
        if any(ca_distance(key, prev, geometry) <= PATCH_A and same_face(key, prev, geometry) for prev in kept):
            continue
        kept.append(key)
    return kept


def patch_members(center: str, geometry: dict) -> List[str]:
    return sorted(k for k in geometry["coords"] if ca_distance(center, k, geometry) <= PATCH_A and same_face(center, k, geometry))


def topk_metric(rows: List[dict], gt: List[str], coords: Dict[str, np.ndarray], surface_keys: set[str], k: int, near_A: float) -> dict:
    selected = rows[:k]
    union_members = {m for p in selected for m in p.get("members", []) if m in coords}
    gset = set(gt)
    direct = union_members & gset
    overlap_recall = len(direct) / len(gset) if gset else 0.0
    overlap_precision = len(direct) / len(union_members) if union_members else 0.0

    surface_gt = gset & surface_keys
    background_fraction = len(surface_gt) / len(surface_keys) if surface_keys else 0.0
    enrichment = (overlap_precision / background_fraction) if background_fraction > 0 else None

    near_hits = 0
    min_distances = []
    for g in gt:
        if union_members:
            d = min(float(np.linalg.norm(coords[g] - coords[m])) for m in union_members)
            min_distances.append(d)
            near_hits += int(d <= near_A)
    near_recall = near_hits / len(gset) if gset else 0.0

    return {
        "n_union_members": len(union_members),
        "overlap_n": len(direct),
        "overlap_recall": overlap_recall,
        "overlap_precision": overlap_precision,
        "exposed_surface_background_fraction": background_fraction,
        "exposed_surface_enrichment": enrichment,
        f"near_{int(near_A)}A_recall": near_recall,
        "minimum_gt_to_topk_union_A": min(min_distances) if min_distances else None,
    }


def evaluate(case: dict) -> dict:
    raw = fetch_pdb(case["pdb_id"])
    prepared, _ = prepare_pdb_text(raw, chain=case.get("chain"))
    coords = ca_map(prepared)
    gt = gt_keys(coords, case)
    if not gt:
        raise RuntimeError(f"No mapped GT for {case['id']}")

    apply_publication_chemistry(v1)
    req = v1.AnalyzeRequest(pdb_text=prepared, chain=None, env=v1.EnvParams(pH=float(case["pH"]), ionic=150.0, temp=298.0))
    base = v1.analyze(req)
    mode = get_surface_mode(case["surface"])
    channel = base["chemistries"][mode.chemistry]

    geometry = build_surface_geometry(base, solve_gnm(prepared, cutoff_A=7.3))
    surface = {str(r["key"]): r for r in base.get("surface_residues", []) if str(r.get("key")) in geometry["coords"]}
    surface_keys = set(surface)
    chem_rows = {str(r["key"]): r for r in channel.get("residues", [])}
    center_rows = {str(r["center_key"]): r for r in channel.get("patch_centers", [])}
    outer_tag = f"{int(SINGLE_RADIUS_A)}A"

    signals = {
        "B1_scRSA": {key: float(r.get("scrsa", 0.0)) for key, r in surface.items()},
        "B2_membership": {key: (1.0 if key in chem_rows and float(chem_rows[key].get("local_score", 0.0)) > 0 else 0.0) for key in surface},
        "B3_local_compatibility": {key: float(chem_rows.get(key, {}).get("local_score", 0.0)) for key in surface},
        "B4_single_radius_9A": {key: float(center_rows.get(key, {}).get(f"density_{outer_tag}_norm", 0.0)) for key in surface},
    }

    variants = {}
    for name, scores in signals.items():
        centers = suppress(scores, geometry)
        variants[name] = [{"center_key": c, "members": patch_members(c, geometry), "center_score": scores[c]} for c in centers]

    full = analyze_interface_v2(surface=case["surface"], pdb_text=raw, chain=case.get("chain"), pH=float(case["pH"]))
    variants["B5_InterfaceScout"] = [{"center_key": p["center_key"], "members": p["members"], "center_score": p.get("patch_coherence", 0.0)} for p in full.get("patches", [])]

    out = {"id": case["id"], "protein": case["protein"], "tier": case["tier"], "gt_kind": case["gt_kind"], "n_gt": len(gt), "n_surface": len(surface_keys), "variants": {}}
    for name, rows in variants.items():
        vals = {"n_predictions": len(rows)}
        for k in (1, 3, 5):
            vals[f"top{k}_5A"] = topk_metric(rows, gt, coords, surface_keys, k, 5.0)
            vals[f"top{k}_8A"] = topk_metric(rows, gt, coords, surface_keys, k, 8.0)
        out["variants"][name] = vals
    return out


def summarize_group(rows: List[dict]) -> dict:
    variants = ["B1_scRSA", "B2_membership", "B3_local_compatibility", "B4_single_radius_9A", "B5_InterfaceScout"]
    summary = {}
    for name in variants:
        summary[name] = {}
        for k in (1, 3, 5):
            overlap = [float(r["variants"][name][f"top{k}_8A"]["overlap_recall"]) for r in rows]
            enrich = [r["variants"][name][f"top{k}_8A"]["exposed_surface_enrichment"] for r in rows]
            enrich = [float(v) for v in enrich if v is not None]
            summary[name][f"median_top{k}_overlap_recall"] = float(np.median(overlap)) if overlap else None
            summary[name][f"median_top{k}_enrichment"] = float(np.median(enrich)) if enrich else None
            for near in (5, 8):
                vals = [float(r["variants"][name][f"top{k}_{near}A"][f"near_{near}A_recall"]) for r in rows]
                summary[name][f"median_top{k}_near_{near}A_recall"] = float(np.median(vals)) if vals else None
                summary[name][f"n_top{k}_near_{near}A_nonzero"] = int(sum(v > 0 for v in vals))
    return summary


def main() -> None:
    manifest = json.loads(MANIFEST.read_text())
    rows = []
    for case in manifest["cases"]:
        print("RUN_MATCHED", case["id"], flush=True)
        rows.append(evaluate(case))
    payload = {
        "design": "geometry- and prediction-budget-matched center-ranking baselines; B4 uses the selected outer 9 A aggregation scale and top-k metrics use the union of the top-k patches",
        "enrichment_background": "all exposed surface residues",
        "results": rows,
        "summary": summarize_group(rows),
        "summary_exact": summarize_group([r for r in rows if r["gt_kind"] == "exact"]),
        "summary_regional": summarize_group([r for r in rows if r["gt_kind"] == "ranges"]),
    }
    Path("matched_baselines.json").write_text(json.dumps(payload, indent=2, sort_keys=True))
    print("MATCHED_SUMMARY " + json.dumps(payload["summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
