"""Reproduce the frozen InterfaceScout publication development and benchmark package."""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

from v2.validation import select_core_parameters_final
from v2.validation.compare_experimental_predictions_v2 import main as compare_benchmark
from v2.validation.generate_experimental_predictions import main as generate_predictions
from v2.validation.summarize_experimental_validation import main as summarize_benchmark

HERE = Path(__file__).resolve().parent
BACKEND = HERE.parents[1]
REPO = HERE.parents[2]
PUBLICATION_DATA = REPO / "publication_data"
MANIFEST = HERE / "benchmark_experimental_final.json"


def _remove(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def _copytree(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def _git_sha() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return None


def _figure4_source(comparison: dict, output_path: Path) -> None:
    rows = []
    for result in comparison.get("results", []):
        if result.get("analysis_set") != "primary":
            continue
        case_id = result["id"]
        top5 = result["topk"]["5"]
        chosen_realization = top5.get("gt_realization")
        gt = next((x for x in result.get("gt_realizations", []) if x["label"] == chosen_realization), None)
        for key in (gt or {}).get("keys", []):
            rows.append({"case_id": case_id, "residue_key": key, "role": "experimental_region"})

        pred = json.loads((BACKEND / "experimental_predictions" / "raw" / f"{case_id}.json").read_text(encoding="utf-8"))
        patches = pred["interfacescout_output"].get("patches", [])
        top1 = set(patches[0].get("members", [])) if patches else set()
        top5_members = {m for patch in patches[:5] for m in patch.get("members", [])}
        for key in sorted(top1):
            rows.append({"case_id": case_id, "residue_key": key, "role": "top1_patch"})
        for key in sorted(top5_members - top1):
            rows.append({"case_id": case_id, "residue_key": key, "role": "top2_to_top5_added"})

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with output_path.open("w", newline="", encoding="utf-8-sig") as fh:
            writer = csv.DictWriter(fh, fieldnames=["case_id", "residue_key", "role"])
            writer.writeheader()
            writer.writerows(rows)


def _verification(development: dict, benchmark_summary: dict) -> dict:
    primary = benchmark_summary["summaries"]["primary"]
    checks = {
        "development_sasa_points_200": development["sasa_selection"]["selected_points_per_atom"] == 200,
        "development_multiscale_radii_6_9": [float(x) for x in development["radius_selection"]["selected_pair_A"]] == [6.0, 9.0],
        "primary_conditions_7": primary["n_conditions"] == 7,
        "primary_unique_proteins_6": primary["n_unique_proteins"] == 6,
        "primary_top3_near8_hits_7": primary["n_top3_near_8A_hits"] == 7,
        "primary_top5_direct_hits_4": primary["n_top5_direct_overlap_hits"] == 4,
        "primary_median_top5_near8_recall_0_80": abs(float(primary["median_top5_near_8A_recall"]) - 0.80) < 1e-12,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "n_mismatches": sum(not value for value in checks.values()),
        "checks": checks,
    }


def main() -> None:
    old_cwd = Path.cwd()
    os.chdir(BACKEND)
    try:
        for path in [
            BACKEND / "core_parameter_selection_consensus.json",
            BACKEND / "selection_pdbs",
            BACKEND / "experimental_predictions",
            BACKEND / "experimental_comparison",
            PUBLICATION_DATA,
        ]:
            _remove(path)

        print("[1/4] Development-panel parameter robustness", flush=True)
        select_core_parameters_final.main()

        print("[2/4] Prediction-first literature benchmark", flush=True)
        generate_predictions(MANIFEST)

        print("[3/4] Experimental-region comparison", flush=True)
        compare_benchmark(MANIFEST)
        summarize_benchmark()

        print("[4/4] Publication package and verification", flush=True)
        dev_dir = PUBLICATION_DATA / "development"
        bench_dir = PUBLICATION_DATA / "benchmark"
        figure_dir = PUBLICATION_DATA / "figure_source_data"
        dev_dir.mkdir(parents=True, exist_ok=True)
        bench_dir.mkdir(parents=True, exist_ok=True)
        figure_dir.mkdir(parents=True, exist_ok=True)

        shutil.copy2(BACKEND / "core_parameter_selection_consensus.json", dev_dir / "core_parameter_selection_consensus.json")
        _copytree(BACKEND / "selection_pdbs", dev_dir / "pdbs")
        shutil.copy2(MANIFEST, bench_dir / "manifest.json")
        _copytree(BACKEND / "experimental_predictions", bench_dir / "predictions")
        _copytree(BACKEND / "experimental_comparison", bench_dir / "comparison")

        comparison = json.loads((BACKEND / "experimental_comparison" / "experimental_comparison.json").read_text(encoding="utf-8"))
        benchmark_summary = json.loads((BACKEND / "experimental_comparison" / "summary.json").read_text(encoding="utf-8"))
        development = json.loads((BACKEND / "core_parameter_selection_consensus.json").read_text(encoding="utf-8"))

        shutil.copy2(BACKEND / "experimental_comparison" / "case_summary.csv", figure_dir / "figure3_case_metrics.csv")
        shutil.copy2(BACKEND / "experimental_comparison" / "all_patch_comparisons.csv", figure_dir / "patch_comparison_metrics.csv")
        _figure4_source(comparison, figure_dir / "figure4_structural_mapping.csv")

        verification = _verification(development, benchmark_summary)
        verification.update({
            "repository_commit": _git_sha(),
            "model_version": "2.4.0-selected-parameters",
            "public_application_version": "1.0-publication",
        })
        (PUBLICATION_DATA / "verification.json").write_text(
            json.dumps(verification, indent=2, sort_keys=True), encoding="utf-8"
        )
        (PUBLICATION_DATA / "README.md").write_text(
            "# InterfaceScout publication data\n\n"
            "This directory is generated by `python -m v2.validation.reproduce_publication` from the frozen publication workflow.\n\n"
            "- `development/`: numerical convergence, radius-pair selection, and exact development PDB snapshots.\n"
            "- `benchmark/`: final benchmark manifest, all canonical map outputs, residue values, patch tables, prepared PDBs, propensity tracks, and comparison results.\n"
            "- `figure_source_data/`: compact source tables for manuscript figures.\n"
            "- `verification.json`: manuscript-level reproducibility checks.\n",
            encoding="utf-8",
        )

        if verification["status"] != "PASS":
            raise RuntimeError("Publication verification failed: " + json.dumps(verification["checks"], sort_keys=True))
        print("PUBLICATION_REPRODUCTION PASS: 0 mismatches", flush=True)
    finally:
        os.chdir(old_cwd)


if __name__ == "__main__":
    main()
