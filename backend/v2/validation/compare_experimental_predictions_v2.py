"""Compare saved blinded InterfaceScout predictions with audited experimental ground truth.

Predictions are never regenerated here. This stage reads the previously saved
InterfaceScout output and only then loads experimental localization labels.

Supported experimental ground-truth forms
-----------------------------------------
exact
    One residue-number set on the selected structural model.
ranges
    One or more residue-number ranges on the selected structural model.
equivalent_chain_ranges
    Homooligomeric residue-number regions that may occur on any symmetry-equivalent
    subunit. Each chain is scored as an alternative experimental realization; the best
    matching chain is reported for a prediction instead of incorrectly requiring all
    symmetry-equivalent copies to be contacted simultaneously.
"""
from __future__ import annotations

import csv
import json
import math
import urllib.request
from io import StringIO
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Sequence

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from v2.model_settings import COARSE_PATCH_RADIUS_A
from v2.prepare import prepare_pdb_text

HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = HERE / "benchmark_experimental_final.json"
PRED_DIR = Path("experimental_predictions")
OUTDIR = Path("experimental_comparison")
PROXIMITY_THRESHOLDS_A = (5.0, 8.0)
TOP_K = (1, 3, 5, 10)


def _fetch_pdb(pid: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")


def _ca_map(prepared_pdb: str) -> Dict[str, np.ndarray]:
    structure = PDBParser(QUIET=True).get_structure("comparison", StringIO(prepared_pdb))
    model = next(structure.get_models())
    coords: Dict[str, np.ndarray] = {}
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=True) or "CA" not in res:
                continue
            icode = str(res.id[2]).strip()
            coords[f"{chain.id}:{int(res.id[1])}:{icode}"] = np.asarray(res["CA"].coord, dtype=float)
    return coords


def _parts(key: str) -> tuple[str, int]:
    p = str(key).split(":")
    return p[0], int(p[1])


def _ranges_keys(coords: Dict[str, np.ndarray], ranges: Sequence[Sequence[int]], chain: str | None = None) -> List[str]:
    out = []
    for key in coords:
        ch, n = _parts(key)
        if chain is not None and ch != chain:
            continue
        if any(int(a) <= n <= int(b) for a, b in ranges):
            out.append(key)
    return sorted(out)


def _gt_realizations(coords: Dict[str, np.ndarray], case: dict) -> List[dict]:
    kind = case.get("gt_kind")
    if kind == "exact":
        wanted = {int(x) for x in case.get("gt_exact", [])}
        keys = sorted(k for k in coords if _parts(k)[1] in wanted)
        return [{"label": "single", "keys": keys}]
    if kind == "ranges":
        keys = _ranges_keys(coords, case.get("gt_ranges", []))
        return [{"label": "single", "keys": keys}]
    if kind == "equivalent_chain_ranges":
        ranges = case.get("gt_ranges", [])
        chains = sorted({_parts(k)[0] for k in coords})
        rows = []
        for chain in chains:
            keys = _ranges_keys(coords, ranges, chain=chain)
            if keys:
                rows.append({"label": f"chain_{chain}", "keys": keys})
        return rows
    raise ValueError(f"Unsupported gt_kind={kind!r} for {case.get('id')}")


def _members(patch: dict, coords: Dict[str, np.ndarray]) -> List[str]:
    return [k for k in patch.get("members", []) if k in coords]


def _patch_metrics_for_gt(patch: dict, gt: List[str], coords: Dict[str, np.ndarray]) -> dict:
    members = _members(patch, coords)
    mset, gset = set(members), set(gt)
    direct = sorted(mset & gset)
    distances = []
    near = {thr: [] for thr in PROXIMITY_THRESHOLDS_A}
    for g in gt:
        if not members:
            distances.append(math.inf)
            continue
        dmin = min(float(np.linalg.norm(coords[g] - coords[m])) for m in members)
        distances.append(dmin)
        for thr in PROXIMITY_THRESHOLDS_A:
            if dmin <= thr:
                near[thr].append(g)
    center = patch.get("center_key")
    center_min_gt = None
    if center in coords and gt:
        center_min_gt = min(float(np.linalg.norm(coords[center] - coords[g])) for g in gt)
    return {
        "display_rank": patch.get("display_rank"),
        "pareto_front": patch.get("pareto_front"),
        "classification": patch.get("classification"),
        "center_key": center,
        "n_members": len(members),
        "overlap_n": len(direct),
        "overlap_recall": len(direct) / len(gset) if gset else 0.0,
        "overlap_precision": len(direct) / len(mset) if mset else 0.0,
        "min_gt_to_patch_A": min(distances) if distances else None,
        "median_gt_to_patch_A": float(np.median(distances)) if distances else None,
        "center_min_gt_A": center_min_gt,
        "overlap_keys": direct,
        **{f"near_{int(thr)}A_recall": len(near[thr]) / len(gset) if gset else 0.0 for thr in PROXIMITY_THRESHOLDS_A},
        **{f"near_{int(thr)}A_keys": sorted(near[thr]) for thr in PROXIMITY_THRESHOLDS_A},
    }


def _choose_realization(patch: dict, realizations: List[dict], coords: Dict[str, np.ndarray]) -> dict:
    candidates = []
    for r in realizations:
        m = _patch_metrics_for_gt(patch, r["keys"], coords)
        m["gt_realization"] = r["label"]
        candidates.append(m)
    candidates.sort(key=lambda m: (
        -float(m.get("near_8A_recall", 0.0)),
        -float(m.get("overlap_recall", 0.0)),
        float(m.get("median_gt_to_patch_A") if m.get("median_gt_to_patch_A") is not None else math.inf),
        str(m.get("gt_realization")),
    ))
    return candidates[0] if candidates else {}


def _union_topk_metrics(patches: List[dict], realizations: List[dict], coords: Dict[str, np.ndarray], k: int) -> dict:
    selected = patches[:k]
    members = sorted({m for p in selected for m in _members(p, coords)})
    metrics = _choose_realization({"members": members, "center_key": None}, realizations, coords)
    metrics["k"] = k
    metrics["n_patches_used"] = len(selected)
    return metrics


def _surface_null(pred: dict, realizations: List[dict], coords: Dict[str, np.ndarray]) -> dict:
    surface = [k for k in pred.get("diagnostics", {}).get("surface_residue_keys", []) if k in coords]
    if not surface:
        return {"n_null_centers": 0}
    protein_centroid = np.mean(np.vstack(list(coords.values())), axis=0)

    def unit(v):
        n = float(np.linalg.norm(v))
        return v / n if n > 1e-12 else np.zeros(3)

    normals = {k: unit(coords[k] - protein_centroid) for k in surface}
    recalls = []
    for center in surface:
        members = [
            k for k in surface
            if float(np.linalg.norm(coords[center] - coords[k])) <= float(COARSE_PATCH_RADIUS_A)
            and float(np.dot(normals[center], normals[k])) > 0.0
        ]
        m = _choose_realization({"members": members, "center_key": center}, realizations, coords)
        recalls.append(float(m.get("near_8A_recall", 0.0)))
    return {
        "n_null_centers": len(recalls),
        "median_near_8A_recall": float(np.median(recalls)) if recalls else None,
        "q90_near_8A_recall": float(np.percentile(recalls, 90)) if recalls else None,
        "max_near_8A_recall": max(recalls) if recalls else None,
    }


def compare_case(case: dict) -> dict:
    pred_path = PRED_DIR / "raw" / f"{case['id']}.json"
    if not pred_path.exists():
        raise FileNotFoundError(f"Prediction file missing: {pred_path}")
    envelope = json.loads(pred_path.read_text(encoding="utf-8"))
    pred = envelope["interfacescout_output"]

    raw = _fetch_pdb(case["pdb_id"])
    prepared, prep = prepare_pdb_text(raw, chain=case.get("chain"))
    coords = _ca_map(prepared)
    realizations = _gt_realizations(coords, case)
    if not realizations or not any(r["keys"] for r in realizations):
        raise RuntimeError(f"No experimental GT residues map to prepared PDB for {case['id']}")

    patches = pred.get("patches", [])
    patch_metrics = [_choose_realization(p, realizations, coords) for p in patches]
    topk = {str(k): _union_topk_metrics(patches, realizations, coords, k) for k in TOP_K}
    primary = [m for m in patch_metrics if int(m.get("pareto_front") or 999) == 1]

    def best(rows: Iterable[dict], field: str):
        vals = [float(r[field]) for r in rows if r.get(field) is not None]
        return max(vals) if vals else None

    return {
        "id": case["id"],
        "analysis_set": case.get("analysis_set", "primary"),
        "protein": case.get("protein"),
        "pdb_id": case.get("pdb_id"),
        "chain": prep.get("selected_chain"),
        "surface": case.get("surface"),
        "pH": case.get("pH"),
        "tier": case.get("tier"),
        "gt_kind": case.get("gt_kind"),
        "experimental_evidence": case.get("evidence"),
        "doi": case.get("doi"),
        "gt_realizations": [{"label": r["label"], "n_residues": len(r["keys"]), "keys": r["keys"]} for r in realizations],
        "prediction_file": str(pred_path),
        "prediction_generated_without_gt": bool(envelope.get("ground_truth_visible_to_predictor") is False),
        "n_predicted_patches": len(patches),
        "n_primary_patches": len(primary),
        "topk": topk,
        "best_primary_overlap_recall": best(primary, "overlap_recall"),
        "best_primary_near_5A_recall": best(primary, "near_5A_recall"),
        "best_primary_near_8A_recall": best(primary, "near_8A_recall"),
        "best_any_overlap_recall": best(patch_metrics, "overlap_recall"),
        "best_any_near_5A_recall": best(patch_metrics, "near_5A_recall"),
        "best_any_near_8A_recall": best(patch_metrics, "near_8A_recall"),
        "matched_surface_null": _surface_null(pred, realizations, coords),
        "patch_metrics": patch_metrics,
        "notes": case.get("notes"),
    }


def _summary(rows: List[dict]) -> dict:
    def med(path):
        vals = []
        for r in rows:
            x = r
            for key in path:
                x = x[key]
            if x is not None:
                vals.append(float(x))
        return float(median(vals)) if vals else None

    return {
        "n_conditions": len(rows),
        "n_unique_proteins": len({r["protein"] for r in rows}),
        "median_top1_overlap_recall": med(["topk", "1", "overlap_recall"]),
        "median_top3_overlap_recall": med(["topk", "3", "overlap_recall"]),
        "median_top5_overlap_recall": med(["topk", "5", "overlap_recall"]),
        "median_top1_near_5A_recall": med(["topk", "1", "near_5A_recall"]),
        "median_top3_near_5A_recall": med(["topk", "3", "near_5A_recall"]),
        "median_top5_near_5A_recall": med(["topk", "5", "near_5A_recall"]),
        "median_top1_near_8A_recall": med(["topk", "1", "near_8A_recall"]),
        "median_top3_near_8A_recall": med(["topk", "3", "near_8A_recall"]),
        "median_top5_near_8A_recall": med(["topk", "5", "near_8A_recall"]),
    }


def main(manifest_path: Path = DEFAULT_MANIFEST) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    OUTDIR.mkdir(parents=True, exist_ok=True)
    per_case_dir = OUTDIR / "per_case"
    per_case_dir.mkdir(exist_ok=True)

    results = []
    patch_rows = []
    summary_rows = []
    for case in manifest.get("cases", []):
        print(f"COMPARE {case['id']}", flush=True)
        r = compare_case(case)
        results.append(r)
        (per_case_dir / f"{case['id']}_comparison.json").write_text(json.dumps(r, indent=2, sort_keys=True), encoding="utf-8")
        for p in r["patch_metrics"]:
            patch_rows.append({
                "case_id": r["id"], "analysis_set": r["analysis_set"], "protein": r["protein"], "tier": r["tier"],
                "display_rank": p.get("display_rank"), "pareto_front": p.get("pareto_front"),
                "center_key": p.get("center_key"), "gt_realization": p.get("gt_realization"),
                "n_members": p.get("n_members"), "overlap_n": p.get("overlap_n"),
                "overlap_recall": p.get("overlap_recall"), "overlap_precision": p.get("overlap_precision"),
                "near_5A_recall": p.get("near_5A_recall"), "near_8A_recall": p.get("near_8A_recall"),
                "min_gt_to_patch_A": p.get("min_gt_to_patch_A"), "center_min_gt_A": p.get("center_min_gt_A"),
            })
        summary_rows.append({
            "case_id": r["id"], "analysis_set": r["analysis_set"], "protein": r["protein"], "pdb_id": r["pdb_id"],
            "surface": r["surface"], "pH": r["pH"], "tier": r["tier"], "gt_kind": r["gt_kind"],
            "n_patches": r["n_predicted_patches"],
            "top1_overlap_recall": r["topk"]["1"]["overlap_recall"],
            "top3_overlap_recall": r["topk"]["3"]["overlap_recall"],
            "top5_overlap_recall": r["topk"]["5"]["overlap_recall"],
            "top1_near_5A_recall": r["topk"]["1"]["near_5A_recall"],
            "top3_near_5A_recall": r["topk"]["3"]["near_5A_recall"],
            "top5_near_5A_recall": r["topk"]["5"]["near_5A_recall"],
            "top1_near_8A_recall": r["topk"]["1"]["near_8A_recall"],
            "top3_near_8A_recall": r["topk"]["3"]["near_8A_recall"],
            "top5_near_8A_recall": r["topk"]["5"]["near_8A_recall"],
            "top1_gt_realization": r["topk"]["1"].get("gt_realization"),
            "best_primary_overlap_recall": r["best_primary_overlap_recall"],
            "best_primary_near_8A_recall": r["best_primary_near_8A_recall"],
            "null_median_near_8A_recall": r["matched_surface_null"].get("median_near_8A_recall"),
            "prediction_generated_without_gt": r["prediction_generated_without_gt"],
            "doi": r["doi"],
        })

    primary = [r for r in results if r["analysis_set"] == "primary"]
    secondary = [r for r in results if r["analysis_set"] == "secondary"]
    payload = {
        "stage": "comparison_only_after_saved_blinded_predictions",
        "prediction_directory": str(PRED_DIR),
        "comparison_tolerances_A": list(PROXIMITY_THRESHOLDS_A),
        "comparison_tolerances_are_not_model_radii": True,
        "top_k": list(TOP_K),
        "summary_all": _summary(results),
        "summary_primary": _summary(primary),
        "summary_secondary": _summary(secondary),
        "results": results,
    }
    (OUTDIR / "experimental_comparison.json").write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    if summary_rows:
        with (OUTDIR / "comparison_summary.csv").open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary_rows[0]))
            w.writeheader(); w.writerows(summary_rows)
    if patch_rows:
        with (OUTDIR / "all_patch_comparisons.csv").open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(patch_rows[0]))
            w.writeheader(); w.writerows(patch_rows)

    print("COMPARISON_SUMMARY " + json.dumps({"primary": payload["summary_primary"], "secondary": payload["summary_secondary"]}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
