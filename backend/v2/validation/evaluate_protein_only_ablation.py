"""Evaluate material-blind InterfaceScout ablations on the locked development set.

The benchmark manifest may contain material labels because they describe the
published experiments, but those labels are never passed to the predictor.
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

from v2.prepare import prepare_pdb_text
from v2.protein_only import analyze_protein_only

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "benchmark_strict.json"
COARSE_NEAR_A = 8.0
VARIANTS = (
    "M0_surface",
    "M1_surface_rin",
    "M2_surface_gnm",
    "M3_surface_rin_gnm",
)


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
            key = f"{chain.id}:{int(residue.id[1])}:{str(residue.id[2]).strip()}"
            out[key] = np.asarray(residue["CA"].coord, dtype=float)
    return out


def _is_gt_resseq(resseq: int, case: dict) -> bool:
    if case["gt_kind"] == "exact":
        return int(resseq) in {int(x) for x in case.get("gt_exact", [])}
    if case["gt_kind"] == "ranges":
        return any(int(a) <= int(resseq) <= int(b) for a, b in case.get("gt_ranges", []))
    raise ValueError(f"Unsupported GT kind: {case['gt_kind']}")


def _gt_keys(coords: Dict[str, np.ndarray], case: dict) -> List[str]:
    return sorted(
        key for key in coords
        if _is_gt_resseq(int(key.split(":")[1]), case)
    )


def _metrics(patch: dict, gt: List[str], coords: Dict[str, np.ndarray]) -> dict:
    members = [k for k in patch.get("members", []) if k in coords]
    gset = set(gt)
    mset = set(members)
    overlap = mset & gset
    near = 0
    for g in gt:
        if members and min(float(np.linalg.norm(coords[g] - coords[m])) for m in members) <= COARSE_NEAR_A:
            near += 1
    return {
        "rank": patch.get("display_rank"),
        "center_key": patch.get("center_key"),
        "n_members": len(members),
        "overlap_recall": len(overlap) / len(gset) if gset else 0.0,
        "near_8A_recall": near / len(gset) if gset else 0.0,
    }


def _best_top(rows: List[dict], k: int, field: str) -> float:
    return max((float(r[field]) for r in rows[:k]), default=0.0)


def evaluate_case(case: dict) -> dict:
    raw = _fetch_pdb(case["pdb_id"])
    prepared, _ = prepare_pdb_text(raw, chain=case.get("chain"))
    coords = _ca_map(prepared)
    gt = _gt_keys(coords, case)
    if not gt:
        raise RuntimeError(f"No GT mapped for {case['id']}")

    pred = analyze_protein_only(
        pdb_text=raw,
        chain=case.get("chain"),
        pH=float(case["pH"]),
    )

    out = {
        "id": case["id"],
        "protein": case["protein"],
        "pH": case["pH"],
        "experimental_surface_for_reference_only": case.get("surface"),
        "material_was_predictor_input": False,
        "tier": case["tier"],
        "gt_kind": case["gt_kind"],
        "n_gt_structure_residues": len(gt),
        "variants": {},
    }
    for variant in VARIANTS:
        rows = [_metrics(p, gt, coords) for p in pred["variants"][variant]]
        out["variants"][variant] = {
            "n_patches": len(rows),
            "top1_near_8A_recall": _best_top(rows, 1, "near_8A_recall"),
            "top3_near_8A_recall": _best_top(rows, 3, "near_8A_recall"),
            "top5_near_8A_recall": _best_top(rows, 5, "near_8A_recall"),
            "top1_overlap_recall": _best_top(rows, 1, "overlap_recall"),
            "top3_overlap_recall": _best_top(rows, 3, "overlap_recall"),
            "top5_overlap_recall": _best_top(rows, 5, "overlap_recall"),
            "patch_metrics": rows,
        }
    return out


def _summary(results: List[dict]) -> dict:
    out = {
        "n_conditions": len(results),
        "n_unique_proteins": len({r["protein"] for r in results}),
        "material_blind": True,
        "variants": {},
    }
    for variant in VARIANTS:
        vals = [r["variants"][variant] for r in results]
        out["variants"][variant] = {
            "median_top1_near_8A_recall": float(np.median([v["top1_near_8A_recall"] for v in vals])),
            "median_top3_near_8A_recall": float(np.median([v["top3_near_8A_recall"] for v in vals])),
            "median_top5_near_8A_recall": float(np.median([v["top5_near_8A_recall"] for v in vals])),
            "median_top1_overlap_recall": float(np.median([v["top1_overlap_recall"] for v in vals])),
            "n_top3_near_nonzero": int(sum(v["top3_near_8A_recall"] > 0 for v in vals)),
            "n_top5_near_nonzero": int(sum(v["top5_near_8A_recall"] > 0 for v in vals)),
        }
    return out


def main() -> None:
    manifest = json.loads(MANIFEST.read_text())
    results = []
    for case in manifest["cases"]:
        print(f"RUN {case['id']}", flush=True)
        r = evaluate_case(case)
        results.append(r)
        msg = " ".join(
            f"{v}:T1={r['variants'][v]['top1_near_8A_recall']:.3f},T3={r['variants'][v]['top3_near_8A_recall']:.3f},T5={r['variants'][v]['top5_near_8A_recall']:.3f}"
            for v in VARIANTS
        )
        print(f"RESULT {case['id']} {msg}", flush=True)

    payload = {
        "experiment": "protein-only material-blind ablation",
        "development_set_warning": "These systems were previously inspected; use results for ablation/development only, not final blind external validation.",
        "summary": _summary(results),
        "results": results,
    }
    Path("protein_only_ablation.json").write_text(json.dumps(payload, indent=2, sort_keys=True))
    print("SUMMARY " + json.dumps(payload["summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
