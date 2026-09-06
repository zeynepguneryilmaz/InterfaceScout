"""Evaluate InterfaceScout V2 against locked experimental coarse-interface GT.

This script is deliberately separate from the predictor.  Ground-truth labels
are never passed to ``analyze_interface_v2`` and are used only after predictions
have been generated.
"""

from __future__ import annotations

import json
import urllib.request
from io import StringIO
from pathlib import Path
from typing import Dict, Iterable, List

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

    return {
        "display_rank": patch.get("display_rank"),
        "pareto_front": patch.get("pareto_front"),
        "n_members": len(members),
        "overlap_n": len(overlap),
        "overlap_recall": overlap_recall,
        "overlap_precision": overlap_precision,
        "near_8A_recall": near_recall,
        "overlap_keys": overlap,
        "near_8A_keys": sorted(near_hits),
        "chemistry_support": patch.get("chemistry_support"),
        "mean_accessibility": patch.get("mean_accessibility"),
        "dynamic_coupling_abs": patch.get("dynamic_coupling_abs"),
        "orientation_coherence": patch.get("orientation_coherence"),
    }


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

    n_surface = len({r["key"] for p in pred.get("patches", []) for r in p.get("member_residues", [])})
    # n_surface above is only the union represented by generated patches; report
    # it as predicted-region universe rather than pretending it is total SASA surface.

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
        "n_predicted_patches": len(metrics),
        "n_primary_pareto_patches": len(primary),
        "best_primary_overlap_recall": best(primary, "overlap_recall"),
        "best_primary_near_8A_recall": best(primary, "near_8A_recall"),
        "best_any_overlap_recall": best(metrics, "overlap_recall"),
        "best_any_near_8A_recall": best(metrics, "near_8A_recall"),
        "primary_union_n_members": len(primary_union),
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
        "median_best_primary_near_8A_recall_all": med(results, "best_primary_near_8A_recall"),
        "median_best_primary_near_8A_recall_exact": med(exact, "best_primary_near_8A_recall"),
        "median_best_primary_near_8A_recall_regional": med(regional, "best_primary_near_8A_recall"),
        "n_conditions_primary_near_recall_nonzero": sum(float(r["best_primary_near_8A_recall"]) > 0.0 for r in results),
        "n_conditions_any_near_recall_nonzero": sum(float(r["best_any_near_8A_recall"]) > 0.0 for r in results),
        "interpretation": "Coarse near-recall is the primary pilot metric; no success threshold is fitted from these labels."
    }


def main() -> None:
    manifest = json.loads(MANIFEST.read_text())
    results = []
    for case in manifest["cases"]:
        print(f"RUN {case['id']}", flush=True)
        result = evaluate_case(case)
        results.append(result)
        print(
            f"RESULT {case['id']} primary_near={result['best_primary_near_8A_recall']:.3f} "
            f"primary_overlap={result['best_primary_overlap_recall']:.3f} "
            f"any_near={result['best_any_near_8A_recall']:.3f}",
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
