"""Generate prediction-first InterfaceScout outputs for the final benchmark.

Only id, protein, pdb_id, chain, surface, and pH are permitted to enter the
benchmark predictor. Experimental localization is loaded only by the later
comparison stage.
"""
from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path

from v2.interface_engine import analyze_all_maps, analyze_interface_v2
from v2.model_settings import (
    MODEL_VERSION,
    SASA_POINTS,
    SASA_PROBE_A,
    SC_RSA_THRESHOLD,
    MULTISCALE_RADII_A,
    COARSE_PATCH_RADIUS_A,
)
from v2.prepare import prepare_pdb_text
from v2.surface_modes import get_surface_mode

HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = HERE / "benchmark_experimental_final.json"
OUTDIR = Path("experimental_predictions")
ALLOWED_INPUT_FIELDS = ("id", "protein", "pdb_id", "chain", "surface", "pH")


def prediction_input(case: dict) -> dict:
    return {k: case.get(k) for k in ALLOWED_INPUT_FIELDS}


def _fetch_pdb(pdb_id: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60) as response:
        return response.read().decode("utf-8", errors="replace")


def flatten_patch(case_id: str, map_key: str, patch: dict) -> dict:
    center = patch.get("center_residue") or {}
    return {
        "case_id": case_id,
        "map_key": map_key,
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


def _residue_rows(case_id: str, map_key: str, map_payload: dict) -> list[dict]:
    rows = []
    for r in map_payload.get("residues", []):
        rows.append({
            "case_id": case_id,
            "map_key": map_key,
            "internal_chemistry_key": map_payload.get("internal_chemistry_key"),
            "chain": r.get("chain"),
            "res_seq": r.get("res_seq"),
            "icode": r.get("icode", ""),
            "res_name": r.get("res_name"),
            "scrsa": r.get("scrsa"),
            "local_score": r.get("local_score"),
            "propensity": r.get("propensity"),
            "multiscale_persistence": r.get("multiscale_persistence"),
        })
    return rows


def _selected_public_map(all_maps: dict, surface: str) -> str:
    internal = get_surface_mode(surface).chemistry
    for key, payload in all_maps.get("maps", {}).items():
        if payload.get("internal_chemistry_key") == internal:
            return key
    raise KeyError(f"No canonical public map found for validation surface {surface!r} -> {internal!r}")


def _write_propensity_pdb(prepared: str, map_payload: dict, output_path: Path, case_id: str) -> None:
    values = {}
    for row in map_payload.get("residues", []):
        key = f"{row.get('chain')}:{int(row.get('res_seq'))}:{str(row.get('icode') or '').strip()}"
        values[key] = float(row.get("propensity") or 0.0)

    out = [
        "REMARK 900 INTERFACESCOUT PUBLICATION PROPERTY TRACK",
        f"REMARK 900 CASE {case_id}",
        f"REMARK 900 MAP {map_payload.get('key')}",
        "REMARK 900 B-FACTOR FIELD CONTAINS NORMALIZED INTERFACESCOUT RESIDUE PROPENSITY 0-100",
        "REMARK 900 VALUES ARE NOT EXPERIMENTAL CRYSTALLOGRAPHIC B-FACTORS",
    ]
    for raw in prepared.splitlines():
        line = raw
        if raw.startswith(("ATOM  ", "HETATM")) and len(raw) >= 27:
            chain = raw[21].strip()
            try:
                seq = int(raw[22:26])
            except ValueError:
                out.append(raw)
                continue
            icode = raw[26].strip()
            key = f"{chain}:{seq}:{icode}"
            value = values.get(key, 0.0)
            line = raw.ljust(66)
            line = line[:60] + f"{value:6.2f}" + line[66:]
        out.append(line)
    output_path.write_text("\n".join(out) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(manifest_path: Path = DEFAULT_MANIFEST) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    OUTDIR.mkdir(parents=True, exist_ok=True)
    raw_dir = OUTDIR / "raw"
    all_maps_dir = OUTDIR / "all_maps"
    residue_dir = OUTDIR / "residues"
    patch_dir = OUTDIR / "patches"
    pdb_full_dir = OUTDIR / "pdb" / "full_rcsb"
    pdb_prepared_dir = OUTDIR / "pdb" / "prepared_for_analysis"
    pdb_track_dir = OUTDIR / "pdb" / "propensity_tracks"
    for directory in (raw_dir, all_maps_dir, residue_dir, patch_dir, pdb_full_dir, pdb_prepared_dir, pdb_track_dir):
        directory.mkdir(parents=True, exist_ok=True)

    generated = []
    selected_patch_rows = []
    pdb_rows = []

    for full_case in manifest.get("cases", []):
        case = prediction_input(full_case)
        missing = [k for k in ("id", "pdb_id", "surface", "pH") if case.get(k) is None]
        if missing:
            raise ValueError(f"Missing prediction inputs for {case.get('id')}: {missing}")

        print(f"PREDICT {case['id']}", flush=True)
        pdb_text = _fetch_pdb(str(case["pdb_id"]))
        prepared, prep_report = prepare_pdb_text(pdb_text, chain=case.get("chain"))
        full_path = pdb_full_dir / f"{case['id']}_{case['pdb_id']}.pdb"
        prepared_path = pdb_prepared_dir / f"{case['id']}_{case['pdb_id']}_prepared.pdb"
        full_path.write_text(pdb_text, encoding="utf-8")
        prepared_path.write_text(prepared, encoding="utf-8")

        all_maps = analyze_all_maps(
            pdb_text=pdb_text,
            chain=case.get("chain"),
            pH=float(case["pH"]),
            ionic_mM=150.0,
            temp_K=298.0,
        )
        selected_map_key = _selected_public_map(all_maps, str(case["surface"]))
        selected_map = all_maps["maps"][selected_map_key]

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
            "selected_public_map": selected_map_key,
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
        (all_maps_dir / f"{case['id']}_all_maps.json").write_text(
            json.dumps(all_maps, indent=2, sort_keys=True), encoding="utf-8"
        )

        residue_rows = []
        all_patch_rows = []
        for map_key in all_maps.get("map_order", []):
            payload = all_maps["maps"][map_key]
            residue_rows.extend(_residue_rows(case["id"], map_key, payload))
            all_patch_rows.extend(flatten_patch(case["id"], map_key, p) for p in payload.get("patches", []))
        _write_csv(residue_dir / f"{case['id']}_residue_values.csv", residue_rows)
        _write_csv(patch_dir / f"{case['id']}_all_maps_patches.csv", all_patch_rows)

        selected_rows = [flatten_patch(case["id"], selected_map_key, p) for p in pred.get("patches", [])]
        selected_patch_rows.extend(selected_rows)
        _write_csv(patch_dir / f"{case['id']}_selected_map_patches.csv", selected_rows)

        track_path = pdb_track_dir / f"{case['id']}_{selected_map_key}_propensity.pdb"
        _write_propensity_pdb(prepared, selected_map, track_path, case["id"])

        pdb_rows.append({
            "case_id": case["id"],
            "protein": case.get("protein"),
            "pdb_id": case["pdb_id"],
            "requested_chain": case.get("chain"),
            "selected_chain": prep_report.get("selected_chain"),
            "full_pdb": str(full_path.as_posix()),
            "prepared_pdb": str(prepared_path.as_posix()),
            "propensity_track": str(track_path.as_posix()),
        })
        generated.append({
            **case,
            "selected_public_map": selected_map_key,
            "model_version": pred.get("version"),
            "n_surface_residues": pred.get("diagnostics", {}).get("n_surface_residues"),
            "n_patches": pred.get("n_patches"),
            "n_front1_patches": pred.get("n_pareto_primary_patches"),
            "raw_output": str((raw_dir / f"{case['id']}.json").as_posix()),
            "all_maps_output": str((all_maps_dir / f"{case['id']}_all_maps.json").as_posix()),
        })

    _write_csv(OUTDIR / "prediction_index.csv", generated)
    _write_csv(OUTDIR / "all_selected_map_patches.csv", selected_patch_rows)
    _write_csv(OUTDIR / "pdb_manifest.csv", pdb_rows)

    metadata = {
        "stage": "prediction_first",
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
