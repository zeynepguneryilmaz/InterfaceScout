"""Generate blinded InterfaceScout predictions for the final experimental benchmark.

Only the following manifest fields are permitted to enter prediction:
    id, protein, pdb_id, chain, surface, pH

Experimental contact labels, evidence tiers, DOI metadata and comparison fields
are never passed to the predictor. Raw predictions and exact PDB inputs are
written before any experimental comparison is performed.
"""
from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path

from v2.interface_engine import analyze_interface_v2
from v2.model_settings import (
    MODEL_VERSION,
    SASA_POINTS,
    SASA_PROBE_A,
    SC_RSA_THRESHOLD,
    MULTISCALE_RADII_A,
    COARSE_PATCH_RADIUS_A,
)
from v2.prepare import prepare_pdb_text

HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = HERE / "benchmark_experimental_final.json"
OUTDIR = Path("experimental_predictions")

ALLOWED_INPUT_FIELDS = ("id", "protein", "pdb_id", "chain", "surface", "pH")


def prediction_input(case: dict) -> dict:
    return {k: case.get(k) for k in ALLOWED_INPUT_FIELDS}


def _fetch_pdb(pdb_id: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


def flatten_patch(case_id: str, patch: dict) -> dict:
    center = patch.get("center_residue") or {}
    return {
        "case_id": case_id,
        "display_rank": patch.get("display_rank"),
        "pareto_front": patch.get("pareto_front"),
        "classification": patch.get("classification"),
        "center_key": patch.get("center_key"),
        "center_res_name": center.get("res_name"),
        "center_res_seq": center.get("res_seq"),
        "center_chain": center.get("chain"),
        "n_members": patch.get("n_members"),
        "n_chemistry_seeds": patch.get("n_chemistry_seeds"),
        "diameter_A": patch.get("diameter_A"),
        "chemistry_support": patch.get("chemistry_support"),
        "mean_accessibility": patch.get("mean_accessibility"),
        "patch_coherence": patch.get("patch_coherence"),
        "orientation_coherence": patch.get("orientation_coherence"),
        "members": ";".join(str(x) for x in patch.get("members", [])),
        "seed_members": ";".join(str(x) for x in patch.get("seed_members", [])),
    }


def main(manifest_path: Path = DEFAULT_MANIFEST) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    OUTDIR.mkdir(parents=True, exist_ok=True)
    raw_dir = OUTDIR / "raw"
    csv_dir = OUTDIR / "patch_csv"
    pdb_full_dir = OUTDIR / "pdb" / "full_rcsb"
    pdb_prepared_dir = OUTDIR / "pdb" / "prepared_for_analysis"
    for directory in (raw_dir, csv_dir, pdb_full_dir, pdb_prepared_dir):
        directory.mkdir(parents=True, exist_ok=True)

    generated = []
    all_patch_rows = []
    pdb_rows = []

    for full_case in manifest.get("cases", []):
        case = prediction_input(full_case)
        missing = [k for k in ("id", "pdb_id", "surface", "pH") if case.get(k) is None]
        if missing:
            raise ValueError(f"Missing prediction inputs for {case.get('id')}: {missing}")

        print(f"PREDICT {case['id']}", flush=True)

        # Archive the exact coordinate input independently of ground truth.
        pdb_text = _fetch_pdb(str(case["pdb_id"]))
        prepared, prep_report = prepare_pdb_text(pdb_text, chain=case.get("chain"))
        full_path = pdb_full_dir / f"{case['id']}_{case['pdb_id']}.pdb"
        prepared_path = pdb_prepared_dir / f"{case['id']}_{case['pdb_id']}_prepared.pdb"
        full_path.write_text(pdb_text, encoding="utf-8")
        prepared_path.write_text(prepared, encoding="utf-8")

        pred = analyze_interface_v2(
            surface=str(case["surface"]),
            pdb_text=pdb_text,
            chain=case.get("chain"),
            pH=float(case["pH"]),
            ionic_mM=150.0,
            temp_K=298.0,
        )

        envelope = {
            "prediction_input": case,
            "ground_truth_visible_to_predictor": False,
            "model_version": MODEL_VERSION,
            "fixed_model_settings": {
                "sasa_points_per_atom": SASA_POINTS,
                "sasa_probe_A": SASA_PROBE_A,
                "scrsa_threshold": SC_RSA_THRESHOLD,
                "multiscale_radii_A": list(MULTISCALE_RADII_A),
                "coarse_patch_radius_A": COARSE_PATCH_RADIUS_A,
            },
            "structure_preparation": prep_report,
            "archived_full_pdb": str(full_path.as_posix()),
            "archived_prepared_pdb": str(prepared_path.as_posix()),
            "interfacescout_output": pred,
        }
        (raw_dir / f"{case['id']}.json").write_text(
            json.dumps(envelope, indent=2, sort_keys=True), encoding="utf-8"
        )

        patch_rows = [flatten_patch(case["id"], p) for p in pred.get("patches", [])]
        all_patch_rows.extend(patch_rows)
        if patch_rows:
            with (csv_dir / f"{case['id']}_patches.csv").open("w", newline="", encoding="utf-8-sig") as fh:
                w = csv.DictWriter(fh, fieldnames=list(patch_rows[0]))
                w.writeheader(); w.writerows(patch_rows)

        pdb_rows.append({
            "case_id": case["id"],
            "protein": case.get("protein"),
            "pdb_id": case["pdb_id"],
            "requested_chain": case.get("chain"),
            "selected_chain": prep_report.get("selected_chain"),
            "full_pdb": str(full_path.as_posix()),
            "prepared_pdb": str(prepared_path.as_posix()),
        })
        generated.append({
            **case,
            "model_version": pred.get("version"),
            "n_surface_residues": pred.get("diagnostics", {}).get("n_surface_residues"),
            "n_patches": pred.get("n_patches"),
            "n_primary_patches": pred.get("n_pareto_primary_patches"),
            "raw_output": str((raw_dir / f"{case['id']}.json").as_posix()),
        })

    if generated:
        with (OUTDIR / "prediction_index.csv").open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(generated[0]))
            w.writeheader(); w.writerows(generated)
    if all_patch_rows:
        with (OUTDIR / "all_predicted_patches.csv").open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(all_patch_rows[0]))
            w.writeheader(); w.writerows(all_patch_rows)
    if pdb_rows:
        with (OUTDIR / "pdb_manifest.csv").open("w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(pdb_rows[0]))
            w.writeheader(); w.writerows(pdb_rows)

    metadata = {
        "stage": "prediction_only",
        "manifest_file": str(manifest_path),
        "allowed_prediction_fields": list(ALLOWED_INPUT_FIELDS),
        "ground_truth_visible_to_predictor": False,
        "n_cases": len(generated),
        "model_version": MODEL_VERSION,
        "fixed_model_settings": {
            "sasa_points_per_atom": SASA_POINTS,
            "sasa_probe_A": SASA_PROBE_A,
            "scrsa_threshold": SC_RSA_THRESHOLD,
            "multiscale_radii_A": list(MULTISCALE_RADII_A),
            "coarse_patch_radius_A": COARSE_PATCH_RADIUS_A,
        },
    }
    (OUTDIR / "prediction_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("PREDICTION_SUMMARY " + json.dumps(metadata, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
