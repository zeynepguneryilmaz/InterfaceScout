"""Compare frozen InterfaceScout predictions with experimental ground truth.

Stage 2 never reruns or modifies the predictor. It reads raw prediction JSON
files generated without access to experimental labels and only then exposes the
locked benchmark ground truth.

Ground-truth resolution is respected:
- exact / anchor evidence: exact overlap plus distance-based proximity;
- regional evidence: overlap recall/precision plus distance-based proximity;
- equivalent_chain_ranges: a sequence-defined region on a homooligomer is
  scored against each symmetry-equivalent chain separately and the best matching
  chain is retained; the predictor is not required to contact every copy;
- orientation/domain evidence requires a dedicated scorer and is not forced
  into residue-overlap metrics.
"""
from __future__ import annotations

import csv
import json
import math
import urllib.request
from io import StringIO
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List

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


def _resseq(key: str) -> int:
    return int(str(key).split(":")[1])


def _chain(key: str) -> str:
    return str(key).split(":", 1)[0]


def _keys_for_ranges(keys: Iterable[str], ranges: List[List[int]]) -> List[str]:
    out = []
    for key in keys:
        n = _resseq(key)
        if any(int(a) <= n <= int(b) for a, b in ranges):
            out.append(key)
    return sorted(out)


def _gt_groups(coords: Dict[str, np.ndarray], case: dict) -> List[dict]:
    kind = case.get("gt_kind")
    if kind == "exact":
        wanted = {int(x) for x in case.get("gt_exact", [])}
        keys = sorted(k for k in coords if _resseq(k) in wanted)
        return [{"group": "all", "keys": keys}]
    if kind == "ranges":
        keys = _keys_for_ranges(coords.keys(), case.get("gt_ranges", []))
        return [{"group": "all", "keys": keys}]
    if kind == "equivalent_chain_ranges":
        ranges = case.get("gt_ranges", [])
        groups = []
        for chain_id in sorted({_chain(k) for k in coords}):
            keys = _keys_for_ranges((k for k in coords if _chain(k) == chain_id), ranges)
            if keys:
                groups.append({"group": chain_id, "keys": keys})
        return groups
    raise ValueError(f"Unsupported residue-level gt_kind={kind!r} for {case.get('id')}")


def _members(patch: dict, coords: Dict[str, np.ndarray]) -> List[str]:
    return [k for k in patch.get("members", []) if k in coords]


def _patch_metrics(patch: dict, gt: List[str], coords: Dict[str, np.ndarray]) -> dict:
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


def _metric_choice_key(m: dict) -> tuple:
    meddist = m.get("median_gt_to_patch_A")
    meddist = float(meddist) if meddist is not None and math.isfinite(float(meddist)) else 1e12
    return (
        float(m.get("near_8A_recall", 0.0)),
        float(m.get("near_5A_recall", 0.0)),
        float(m.get("overlap_recall", 0.0)),
        float(m.get("overlap_precision", 0.0)),
        -meddist,
    )


def _best_group_metrics(patch: dict, gt_groups: List[dict], coords: Dict[str, np.ndarray]) -> dict:
    candidates = []
    for group in gt_groups:
        m = _patch_metrics(patch, group["keys"], coords)
        m["matched_gt_group"] = group["group"]
        m["n_gt_in_matched_group"] = len(group["keys"])
        candidates.append(m)
    if not candidates:
        raise RuntimeError("No experimental ground-truth groups mapped to the prepared structure")
    return max(candidates, key=_metric_choice_key)


def _union_topk_metrics(patches: List[dict], gt_groups: List[dict], coords: Dict[str, np.ndarray], k: int) -> dict:
    selected = patches[:k]
    members = sorted({m for p in selected for m in _members(p, coords)})
    metrics = _best_group_metrics({"members": members, "center_key": None}, gt_groups, coords)
    metrics["k"] = k
    metrics["n_patches_used"] = len(selected)
    return metrics


def _surface_null(pred: dict, gt_groups: List[dict], coords: Dict[str, np.ndarray], patch_radius_A: float = COARSE_PATCH_RADIUS_A) -> dict:
    surface = [k for k in pred.get("diagnostics", {}).get("surface_residue_keys", []) if k in coords]
    if not surface:
        return {"n_null_centers": 0}
    protein_centroid = np.mean(np.vstack(list(coords.values())), axis=0)
    def unit(v):
        n = float(np.linalg.norm(v)); return v / n if n > 1e-12 else np.zeros(3)
    normals = {k: unit(coords[k] - protein_centroid) for k in surface}
    recalls = []
    for center in surface:
        members = [k for k in surface if float(np.linalg.norm(coords[center]-coords[k])) <= patch_radius_A and float(np.dot(normals[center], normals[k])) > 0.0]
        m = _best_group_metrics({"members": members}, gt_groups, coords)
        recalls.append(float(m["near_8A_recall"]))
    return {
        "patch_radius_A": float(patch_radius_A),
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

    prepared_path = envelope.get("archived_prepared_pdb")
    if prepared_path and Path(prepared_path).exists():
        prepared = Path(prepared_path).read_text(encoding="utf-8")
        prep = envelope.get("structure_preparation", {})
    else:
        raw = _fetch_pdb(case["pdb_id"])
        prepared, prep = prepare_pdb_text(raw, chain=case.get("chain"))
    coords = _ca_map(prepared)
    gt_groups = [g for g in _gt_groups(coords, case) if g["keys"]]
    if not gt_groups:
        raise RuntimeError(f"No experimental GT residues map to prepared PDB for {case['id']}")

    patches = pred.get("patches", [])
    patch_metrics = [_best_group_metrics(p, gt_groups, coords) for p in patches]
    topk = {str(k): _union_topk_metrics(patches, gt_groups, coords, k) for k in TOP_K}
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
        "source_url": case.get("source_url"),
        "n_gt_groups": len(gt_groups),
        "gt_groups": gt_groups,
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
        "matched_surface_null": _surface_null(pred, gt_groups, coords),
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
            vals.append(float(x))
        return float(median(vals)) if vals else None
    return {
        "n_conditions": len(rows),
        "n_unique_proteins": len({r['protein'] for r in rows}),
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
        if case.get("gt_kind") not in {"exact", "ranges", "equivalent_chain_ranges"}:
            print(f"SKIP_NONRESIDUE {case.get('id')} gt_kind={case.get('gt_kind')}", flush=True)
            continue
        print(f"COMPARE {case['id']}", flush=True)
        r = compare_case(case)
        results.append(r)
        (per_case_dir / f"{case['id']}_comparison.json").write_text(json.dumps(r, indent=2, sort_keys=True), encoding="utf-8")
        for p in r["patch_metrics"]:
            patch_rows.append({
                "case_id": r["id"], "analysis_set": r["analysis_set"], "protein": r["protein"], "tier": r["tier"],
                "display_rank": p.get("display_rank"), "pareto_front": p.get("pareto_front"),
                "center_key": p.get("center_key"), "matched_gt_group": p.get("matched_gt_group"),
                "n_members": p.get("n_members"), "overlap_n": p.get("overlap_n"),
                "overlap_recall": p.get("overlap_recall"), "overlap_precision": p.get("overlap_precision"),
                "near_5A_recall": p.get("near_5A_recall"), "near_8A_recall": p.get("near_8A_recall"),
                "min_gt_to_patch_A": p.get("min_gt_to_patch_A"), "center_min_gt_A": p.get("center_min_gt_A"),
            })
        summary_rows.append({
            "case_id": r["id"], "analysis_set": r["analysis_set"], "protein": r["protein"], "pdb_id": r["pdb_id"], "surface": r["surface"],
            "pH": r["pH"], "tier": r["tier"], "gt_kind": r["gt_kind"], "n_gt_groups": r["n_gt_groups"],
            "n_patches": r["n_predicted_patches"], "top1_overlap_recall": r["topk"]["1"]["overlap_recall"],
            "top3_overlap_recall": r["topk"]["3"]["overlap_recall"], "top5_overlap_recall": r["topk"]["5"]["overlap_recall"],
            "top1_near_5A_recall": r["topk"]["1"]["near_5A_recall"], "top3_near_5A_recall": r["topk"]["3"]["near_5A_recall"],
            "top5_near_5A_recall": r["topk"]["5"]["near_5A_recall"], "top1_near_8A_recall": r["topk"]["1"]["near_8A_recall"],
            "top3_near_8A_recall": r["topk"]["3"]["near_8A_recall"], "top5_near_8A_recall": r["topk"]["5"]["near_8A_recall"],
            "best_primary_overlap_recall": r["best_primary_overlap_recall"], "best_primary_near_8A_recall": r["best_primary_near_8A_recall"],
            "null_median_near_8A_recall": r["matched_surface_null"].get("median_near_8A_recall"),
            "prediction_generated_without_gt": r["prediction_generated_without_gt"], "doi": r["doi"], "source_url": r["source_url"],
        })

    primary = [r for r in results if r.get("analysis_set") == "primary"]
    secondary = [r for r in results if r.get("analysis_set") == "secondary"]
    payload = {
        "stage": "comparison_only_after_frozen_predictions",
        "prediction_directory": str(PRED_DIR),
        "proximity_thresholds_A": list(PROXIMITY_THRESHOLDS_A),
        "top_k": list(TOP_K),
        "summary_primary": _summary(primary),
        "summary_secondary": _summary(secondary),
        "summary_all": _summary(results),
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

    print("COMPARISON_SUMMARY " + json.dumps(payload["summary_primary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
