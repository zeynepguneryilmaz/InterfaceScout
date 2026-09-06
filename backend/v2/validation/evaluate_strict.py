"""Evaluate InterfaceScout V2 against locked experimental coarse-interface GT.

Ground-truth labels are never passed to the predictor. They are used only after
predictions are generated. Validation reports continuous metrics rather than a
post-hoc fitted success threshold.
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

from v2.interface_engine import analyze_interface_v2
from v2.prepare import prepare_pdb_text

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "benchmark_strict.json"
COARSE_NEAR_A = 8.0  # same frozen spatial scale as the V2 coarse patch definition


def _fetch_pdb(pdb_id: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=30) as r:
        return r.read().decode("utf-8", errors="replace")


def _ca_map(pdb_text: str) -> Dict[str, np.ndarray]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("validation", StringIO(pdb_text))
    model = next(structure.get_models())
    out: Dict[str, np.ndarray] = {}
    for chain in model:
        for residue in chain:
            if not is_aa(residue, standard=True) or "CA" not in residue:
                continue
            icode = str(residue.id[2]).strip()
            key = f"{chain.id}:{int(residue.id[1])}:{icode}"
            out[key] = np.asarray(residue["CA"].coord, dtype=float)
    return out


def _is_gt_resseq(resseq: int, case: dict) -> bool:
    if case["gt_kind"] == "exact":
        return int(resseq) in {int(x) for x in case.get("gt_exact", [])}
    if case["gt_kind"] == "ranges":
        return any(int(a) <= int(resseq) <= int(b) for a, b in case.get("gt_ranges", []))
    raise ValueError(f"Unsupported GT kind: {case['gt_kind']}")


def _gt_keys(coords: Dict[str, np.ndarray], case: dict) -> List[str]:
    keys = []
    for key in coords:
        try:
            resseq = int(key.split(":")[1])
        except Exception:
            continue
        if _is_gt_resseq(resseq, case):
            keys.append(key)
    return sorted(keys)


def _patch_metrics(patch: dict, gt: List[str], coords: Dict[str, np.ndarray]) -> dict:
    members = [k for k in patch.get("members", []) if k in coords]
    mset = set(members)
    gset = set(gt)
    overlap = sorted(mset & gset)
    overlap_recall = len(overlap) / len(gset) if gset else 0.0
    overlap_precision = len(overlap) / len(mset) if mset else 0.0

    near_hits = []
    for g in gt:
        if not members:
            continue
        dmin = min(float(np.linalg.norm(coords[g] - coords[m])) for m in members)
        if dmin <= COARSE_NEAR_A:
            near_hits.append(g)
    near_recall = len(near_hits) / len(gset) if gset else 0.0

    center = patch.get("center_key")
    center_min_gt_A = None
    if center in coords and gt:
        center_min_gt_A = min(float(np.linalg.norm(coords[center] - coords[g])) for g in gt)

    return {
        "display_rank": patch.get("display_rank"),
        "pareto_front": patch.get("pareto_front"),
        "center_key": center,
        "center_min_gt_A": center_min_gt_A,
        "n_members": len(members),
        "overlap_n": len(overlap),
        "overlap_recall": overlap_recall,
        "overlap_precision": overlap_precision,
        "near_8A_recall": near_recall,
        "overlap_keys": overlap,
        "near_8A_keys": sorted(near_hits),
        "chemistry_support": patch.get("chemistry_support"),
        "mean_accessibility": patch.get("mean_accessibility"),
        "patch_coherence": patch.get("patch_coherence"),
        "dynamic_coupling_abs": patch.get("dynamic_coupling_abs"),
        "orientation_coherence": patch.get("orientation_coherence"),
    }


def _surface_null_patches(surface_keys: List[str], coords: Dict[str, np.ndarray]) -> List[dict]:
    """Every exposed residue is used once as a geometry-matched null centre.

    The null uses the same non-transitive 8 A neighbourhood and same-hemisphere
    rule as V2, but ignores material chemistry, V1 patch coherence and GNM.
    Therefore it asks whether the V2 ranking localizes experiment better than an
    arbitrary solvent-exposed surface location of the same geometric scale.
    """
    valid = [k for k in surface_keys if k in coords]
    if not valid:
        return []
    protein_centroid = np.mean(np.vstack(list(coords.values())), axis=0)

    def unit(v: np.ndarray) -> np.ndarray:
        n = float(np.linalg.norm(v))
        return v / n if n > 1e-12 else np.zeros(3)

    normals = {k: unit(coords[k] - protein_centroid) for k in valid}
    out = []
    for center in valid:
        members = []
        for key in valid:
            d = float(np.linalg.norm(coords[center] - coords[key]))
            if d <= COARSE_NEAR_A and float(np.dot(normals[center], normals[key])) > 0.0:
                members.append(key)
        out.append({"center_key": center, "members": members})
    return out


def _best_top_k(metrics: List[dict], k: int, field: str) -> float:
    return max((float(m[field]) for m in metrics[:k]), default=0.0)


def evaluate_case(case: dict) -> dict:
    raw = _fetch_pdb(case["pdb_id"])
    prepared, prep = prepare_pdb_text(raw, chain=case.get("chain"))
    coords = _ca_map(prepared)
    gt = _gt_keys(coords, case)
    if not gt:
        raise RuntimeError(f"No GT residues mapped into prepared PDB for {case['id']}")

    pred = analyze_interface_v2(
        surface=case["surface"],
        pdb_text=raw,
        chain=case.get("chain"),
        pH=float(case["pH"]),
    )
    metrics = [_patch_metrics(p, gt, coords) for p in pred.get("patches", [])]
    primary = [m for m in metrics if int(m.get("pareto_front") or 999) == 1]

    def best(rows: List[dict], field: str) -> float:
        return max((float(r[field]) for r in rows), default=0.0)

    primary_union = set()
    for p in pred.get("primary_patches", []):
        primary_union.update(p.get("members", []))

    surface_keys = list(pred.get("diagnostics", {}).get("surface_residue_keys", []))
    null_patches = _surface_null_patches(surface_keys, coords)
    null_metrics = [_patch_metrics(p, gt, coords) for p in null_patches]
    top1_near = _best_top_k(metrics, 1, "near_8A_recall")
    null_near = [float(m["near_8A_recall"]) for m in null_metrics]
    null_ge_top1_fraction = (
        float(np.mean([v >= top1_near for v in null_near])) if null_near else None
    )
    null_median = float(np.median(null_near)) if null_near else None

    patch_sizes = [int(m["n_members"]) for m in metrics]

    return {
        "id": case["id"],
        "protein": case["protein"],
        "pdb_id": case["pdb_id"],
        "chain": prep.get("selected_chain"),
        "surface": case["surface"],
        "pH": case["pH"],
        "tier": case["tier"],
        "gt_kind": case["gt_kind"],
        "n_gt_structure_residues": len(gt),
        "n_surface_residues": len(surface_keys),
        "n_predicted_patches": len(metrics),
        "n_primary_pareto_patches": len(primary),
        "median_patch_size": float(np.median(patch_sizes)) if patch_sizes else None,
        "max_patch_size": max(patch_sizes) if patch_sizes else None,
        "top1_near_8A_recall": _best_top_k(metrics, 1, "near_8A_recall"),
        "top3_near_8A_recall": _best_top_k(metrics, 3, "near_8A_recall"),
        "top5_near_8A_recall": _best_top_k(metrics, 5, "near_8A_recall"),
        "top10_near_8A_recall": _best_top_k(metrics, 10, "near_8A_recall"),
        "top1_overlap_recall": _best_top_k(metrics, 1, "overlap_recall"),
        "top3_overlap_recall": _best_top_k(metrics, 3, "overlap_recall"),
        "top5_overlap_recall": _best_top_k(metrics, 5, "overlap_recall"),
        "best_primary_overlap_recall": best(primary, "overlap_recall"),
        "best_primary_near_8A_recall": best(primary, "near_8A_recall"),
        "best_any_overlap_recall": best(metrics, "overlap_recall"),
        "best_any_near_8A_recall": best(metrics, "near_8A_recall"),
        "primary_union_n_members": len(primary_union),
        "matched_surface_null": {
            "n_null_centres": len(null_metrics),
            "median_near_8A_recall": null_median,
            "fraction_null_at_least_as_good_as_v2_top1": null_ge_top1_fraction,
        },
        "patch_metrics": metrics,
        "gt_keys": gt,
        "notes": case.get("notes"),
    }


def summarize(results: List[dict]) -> dict:
    unique = sorted({r["protein"] for r in results})
    exact = [r for r in results if r["gt_kind"] == "exact"]
    regional = [r for r in results if r["gt_kind"] == "ranges"]

    def med(rows: List[dict], field: str) -> float | None:
        if not rows:
            return None
        return float(np.median([float(r[field]) for r in rows]))

    return {
        "n_conditions": len(results),
        "n_unique_proteins": len(unique),
        "unique_proteins": unique,
        "median_top1_near_8A_recall": med(results, "top1_near_8A_recall"),
        "median_top3_near_8A_recall": med(results, "top3_near_8A_recall"),
        "median_top5_near_8A_recall": med(results, "top5_near_8A_recall"),
        "median_best_primary_near_8A_recall_all": med(results, "best_primary_near_8A_recall"),
        "median_best_primary_near_8A_recall_exact": med(exact, "best_primary_near_8A_recall"),
        "median_best_primary_near_8A_recall_regional": med(regional, "best_primary_near_8A_recall"),
        "n_conditions_top3_near_recall_nonzero": sum(float(r["top3_near_8A_recall"]) > 0.0 for r in results),
        "n_conditions_top5_near_recall_nonzero": sum(float(r["top5_near_8A_recall"]) > 0.0 for r in results),
        "n_conditions_primary_near_recall_nonzero": sum(float(r["best_primary_near_8A_recall"]) > 0.0 for r in results),
        "interpretation": "Continuous coarse-patch metrics plus an all-surface matched geometric null; no adsorption-label-fitted threshold."
    }


def main() -> None:
    manifest = json.loads(MANIFEST.read_text())
    results = []
    for case in manifest["cases"]:
        print(f"RUN {case['id']}", flush=True)
        result = evaluate_case(case)
        results.append(result)
        print(
            f"RESULT {case['id']} top1={result['top1_near_8A_recall']:.3f} "
            f"top3={result['top3_near_8A_recall']:.3f} top5={result['top5_near_8A_recall']:.3f} "
            f"null_ge_top1={result['matched_surface_null']['fraction_null_at_least_as_good_as_v2_top1']}",
            flush=True,
        )

    out = {
        "policy": manifest["policy"],
        "summary": summarize(results),
        "results": results,
        "pending_high_quality_cases": manifest.get("pending_high_quality_cases", []),
    }
    target = Path("v2_strict_validation.json")
    target.write_text(json.dumps(out, indent=2, sort_keys=True))
    print("SUMMARY " + json.dumps(out["summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
