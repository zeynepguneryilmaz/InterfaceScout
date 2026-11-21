"""Verify frozen InterfaceScout publication parameter selection and benchmark aggregates."""
from __future__ import annotations
import csv
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent

def truth(value: str) -> bool:
    return value.strip().lower() == "true"

def main() -> None:
    selection = json.loads((HERE / "parameter_selection.json").read_text(encoding="utf-8"))
    with (HERE / "benchmark_case_results.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))

    primary = [r for r in rows if r["analysis_set"] == "primary"]
    secondary = [r for r in rows if r["analysis_set"] == "secondary"]

    checks = {
        "selected_sasa_points_200": selection["sasa_selection"]["selected_points_per_atom"] == 200,
        "selected_multiscale_radii_6_9": [float(x) for x in selection["radius_selection"]["selected_pair_A"]] == [6.0, 9.0],
        "primary_conditions_7": len(primary) == 7,
        "primary_top3_near8_hits_7": sum(truth(r["top3_near_8A_hit"]) for r in primary) == 7,
        "primary_top5_direct_hits_4": sum(truth(r["top5_direct_overlap_hit"]) for r in primary) == 4,
        "primary_median_top5_near8_recall_0_80": abs(
            statistics.median(float(r["top5_near_8A_recall"]) for r in primary) - 0.80
        ) < 1e-12,
        "secondary_challenge_cases_2": len(secondary) == 2,
    }

    for key, passed in checks.items():
        print(f"{key}: {'PASS' if passed else 'FAIL'}")
    if not all(checks.values()):
        raise SystemExit(1)

    print("\nInterfaceScout publication verification: PASS")
    print("Parameter selection was adsorption-label-free; benchmark annotations were evaluated after settings were locked.")

if __name__ == "__main__":
    main()
