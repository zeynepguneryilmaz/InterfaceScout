"""Create evidence-resolution-matched summaries from already saved experimental comparisons.

This script does not alter predictions, model parameters, surface modes, ground truth,
or patch ranking. It only avoids interpreting a broad experimental domain/region as if
InterfaceScout were expected to recover every residue in that region with one local patch.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import median

COMP_PATH = Path("experimental_comparison/experimental_comparison.json")
OUTDIR = Path("experimental_comparison")


def first_rank(rows: list[dict], field: str) -> int | None:
    ranks = [int(r["display_rank"]) for r in rows if r.get("display_rank") is not None and float(r.get(field) or 0.0) > 0.0]
    return min(ranks) if ranks else None


def summarize_case(r: dict) -> dict:
    tier = str(r.get("tier") or "")
    resolution = "residue_or_anchor" if tier.startswith("A") else "region_or_peptide"
    f_overlap = first_rank(r["patch_metrics"], "overlap_recall")
    f_near5 = first_rank(r["patch_metrics"], "near_5A_recall")
    f_near8 = first_rank(r["patch_metrics"], "near_8A_recall")
    return {
        "case_id": r["id"],
        "analysis_set": r["analysis_set"],
        "protein": r["protein"],
        "tier": tier,
        "evidence_resolution": resolution,
        "first_direct_overlap_rank": f_overlap,
        "first_near_5A_rank": f_near5,
        "first_near_8A_rank": f_near8,
        "top3_near_8A_hit": bool(f_near8 is not None and f_near8 <= 3),
        "top5_near_8A_hit": bool(f_near8 is not None and f_near8 <= 5),
        "top5_direct_overlap_hit": bool(f_overlap is not None and f_overlap <= 5),
        "top3_near_8A_recall": r["topk"]["3"]["near_8A_recall"],
        "top5_near_8A_recall": r["topk"]["5"]["near_8A_recall"],
        "top5_overlap_recall": r["topk"]["5"]["overlap_recall"],
        "best_primary_near_8A_recall": r.get("best_primary_near_8A_recall"),
        "null_median_near_8A_recall": r["matched_surface_null"].get("median_near_8A_recall"),
        "null_q90_near_8A_recall": r["matched_surface_null"].get("q90_near_8A_recall"),
        "doi": r.get("doi"),
    }


def set_summary(rows: list[dict]) -> dict:
    if not rows:
        return {"n_conditions": 0, "n_unique_proteins": 0}
    f8 = [x["first_near_8A_rank"] for x in rows if x["first_near_8A_rank"] is not None]
    return {
        "n_conditions": len(rows),
        "n_unique_proteins": len({x["protein"] for x in rows}),
        "n_top3_near_8A_hits": sum(bool(x["top3_near_8A_hit"]) for x in rows),
        "n_top5_near_8A_hits": sum(bool(x["top5_near_8A_hit"]) for x in rows),
        "n_top5_direct_overlap_hits": sum(bool(x["top5_direct_overlap_hit"]) for x in rows),
        "median_first_near_8A_rank_among_hits": float(median(f8)) if f8 else None,
    }


def main() -> None:
    payload = json.loads(COMP_PATH.read_text(encoding="utf-8"))
    rows = [summarize_case(r) for r in payload["results"]]
    primary = [r for r in rows if r["analysis_set"] == "primary"]
    secondary = [r for r in rows if r["analysis_set"] == "secondary"]
    primary_A = [r for r in primary if r["evidence_resolution"] == "residue_or_anchor"]
    primary_B = [r for r in primary if r["evidence_resolution"] == "region_or_peptide"]

    if rows:
        with (OUTDIR / "resolution_matched_case_summary.csv").open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)

    summaries = {
        "primary_all": set_summary(primary),
        "primary_residue_or_anchor": set_summary(primary_A),
        "primary_region_or_peptide": set_summary(primary_B),
        "secondary": set_summary(secondary),
    }
    (OUTDIR / "resolution_matched_summary.json").write_text(json.dumps({"summaries": summaries, "cases": rows}, indent=2), encoding="utf-8")

    summary_rows = []
    for name, s in summaries.items():
        summary_rows.append({"group": name, **s})
    with (OUTDIR / "resolution_matched_group_summary.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0]))
        w.writeheader(); w.writerows(summary_rows)

    print("RESOLUTION_MATCHED_SUMMARY " + json.dumps(summaries, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
