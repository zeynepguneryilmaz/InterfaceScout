"""Create resolution-aware summaries from saved experimental comparisons."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import median

COMP_PATH = Path("experimental_comparison/experimental_comparison.json")
OUTDIR = Path("experimental_comparison")


def first_position(rows: list[dict], field: str) -> int | None:
    positions = [
        int(r["display_rank"])
        for r in rows
        if r.get("display_rank") is not None and float(r.get(field) or 0.0) > 0.0
    ]
    return min(positions) if positions else None


def summarize_case(r: dict) -> dict:
    f_overlap = first_position(r["patch_metrics"], "overlap_recall")
    f_near5 = first_position(r["patch_metrics"], "near_5A_recall")
    f_near8 = first_position(r["patch_metrics"], "near_8A_recall")
    return {
        "case_id": r["id"],
        "analysis_set": r["analysis_set"],
        "protein": r["protein"],
        "evidence_resolution": r.get("evidence_resolution"),
        "first_direct_overlap_position": f_overlap,
        "first_near_5A_position": f_near5,
        "first_near_8A_position": f_near8,
        "top3_near_8A_hit": bool(f_near8 is not None and f_near8 <= 3),
        "top5_near_8A_hit": bool(f_near8 is not None and f_near8 <= 5),
        "top5_direct_overlap_hit": bool(f_overlap is not None and f_overlap <= 5),
        "top3_near_8A_recall": r["topk"]["3"]["near_8A_recall"],
        "top5_near_8A_recall": r["topk"]["5"]["near_8A_recall"],
        "top5_overlap_recall": r["topk"]["5"]["overlap_recall"],
        "best_front1_near_8A_recall": r.get("best_primary_near_8A_recall"),
        "null_median_near_8A_recall": r["matched_surface_null"].get("median_near_8A_recall"),
        "null_q90_near_8A_recall": r["matched_surface_null"].get("q90_near_8A_recall"),
        "doi": r.get("doi"),
    }


def set_summary(rows: list[dict]) -> dict:
    if not rows:
        return {"n_conditions": 0, "n_unique_proteins": 0}
    f8 = [x["first_near_8A_position"] for x in rows if x["first_near_8A_position"] is not None]
    top5_recall = [float(x["top5_near_8A_recall"]) for x in rows]
    return {
        "n_conditions": len(rows),
        "n_unique_proteins": len({x["protein"] for x in rows}),
        "n_top3_near_8A_hits": sum(bool(x["top3_near_8A_hit"]) for x in rows),
        "n_top5_near_8A_hits": sum(bool(x["top5_near_8A_hit"]) for x in rows),
        "n_top5_direct_overlap_hits": sum(bool(x["top5_direct_overlap_hit"]) for x in rows),
        "median_first_near_8A_position_among_hits": float(median(f8)) if f8 else None,
        "median_top5_near_8A_recall": float(median(top5_recall)) if top5_recall else None,
    }


def main() -> None:
    payload = json.loads(COMP_PATH.read_text(encoding="utf-8"))
    rows = [summarize_case(r) for r in payload["results"]]
    primary = [r for r in rows if r["analysis_set"] == "primary"]
    secondary = [r for r in rows if r["analysis_set"] == "secondary"]

    if rows:
        with (OUTDIR / "case_summary.csv").open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    summaries = {
        "primary": set_summary(primary),
        "secondary": set_summary(secondary),
    }
    (OUTDIR / "summary.json").write_text(
        json.dumps({"summaries": summaries, "cases": rows}, indent=2), encoding="utf-8"
    )

    summary_rows = [{"group": name, **values} for name, values in summaries.items()]
    with (OUTDIR / "group_summary.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0]))
        w.writeheader()
        w.writerows(summary_rows)

    print("BENCHMARK_SUMMARY " + json.dumps(summaries, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
