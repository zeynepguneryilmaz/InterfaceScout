"""Run the strict V2 evaluator on the frozen core plus independently verified additions."""

from __future__ import annotations

import json
from pathlib import Path

from v2.validation.evaluate_strict import evaluate_case, summarize

HERE = Path(__file__).resolve().parent
CORE = HERE / "benchmark_strict.json"
ADDITIONAL = HERE / "benchmark_additional.json"


def main() -> None:
    core = json.loads(CORE.read_text())
    additional = json.loads(ADDITIONAL.read_text()) if ADDITIONAL.exists() else {"cases": [], "pending_high_quality_cases": []}
    cases = list(core.get("cases", [])) + list(additional.get("cases", []))
    results = []
    for case in cases:
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
        "policy": core["policy"],
        "summary": summarize(results),
        "results": results,
        "pending_high_quality_cases": list(core.get("pending_high_quality_cases", [])) + list(additional.get("pending_high_quality_cases", [])),
    }
    target = Path("v2_extended_validation.json")
    target.write_text(json.dumps(out, indent=2, sort_keys=True))
    print("SUMMARY " + json.dumps(out["summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
