"""TV48: run the existing InterfaceScout production engine + real APBS on all 48.

Inputs are the pre-frozen quality-screened biological assembly 1 files, with
assembly SHA-256 verified before any calculation. No fitted electrostatic
corrections or modifications to the canonical 11-map/Pareto algorithm.

Usage:
 python validation/tv48_apbs_full.py --shard 0 --shards 8
 python validation/tv48_apbs_full.py --aggregate tv48_apbs_results
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "backend"), str(ROOT / "validation")]
from interface_engine import analyze_all_maps
from export_bundle import (
    _bfactor_pdb, _csv_bytes, _map_scores, _topk_scores,
    _xlsx_bytes, structure_tables,
)
from model_settings import MODEL_VERSION, PUBLIC_VERSION
import strict_apo_compat_48 as preparation

MANIFEST_PATH = ROOT / "validation" / "strict_apo_48_selected_candidates_v2_final.csv"
MANIFEST = list(csv.DictReader(MANIFEST_PATH.open(encoding="utf-8")))
EXPECTED = [item["pdb_id"].upper() for item in MANIFEST]
assert len(EXPECTED) == 48 and len(set(EXPECTED)) == 48, "TV48 manifest must contain 48 unique entries"
OUT = Path("tv48_apbs_results")

SETTINGS = dict(pH=7.4, ionic_mM=150.0, temp_K=298.0)
STATUS_FIELDS = [
    "pdb_id", "status", "assembly_sha256", "prepared_sha256", "model_version",
    "software_version", "n_maps", "n_residues", "n_surface_residues",
    "n_pqr_atoms", "n_matched_residues", "n_patches", "canonical_maps_identical",
    "residue_profiles_identical", "elapsed_s", "error_type", "error",
]
MAP_FIELDS = [
    "pdb_id", "map", "n_residues", "n_patches",
    "n_primary_patches", "top1_center_key", "top5_center_keys",
]


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_one(item: dict) -> tuple[dict, list[dict]]:
    pid = item["pdb_id"].upper()
    start = time.monotonic()
    details = {"pdb_id": pid, "status": "FAIL"}
    permap = []
    try:
        preparation.prepare(pid)
        cif = preparation.OUT / f"{pid}-assembly1.cif"
        actual = hashlib.sha256(cif.read_bytes()).hexdigest()
        if actual != item["assembly_sha256"]:
            raise ValueError(f"Frozen assembly SHA mismatch: {actual}")
        pdb_path = preparation.OUT / f"{pid}_prepared.pdb"
        pdb = pdb_path.read_text(encoding="utf-8")
        prepared_sha = hashlib.sha256(pdb.encode("utf-8")).hexdigest()
        input_args = dict(pdb_text=pdb, chain=None, **SETTINGS)

        # Two real independent production runs: verify the optional layer has
        # no side effects on any map, patch, Pareto rank, or residue descriptor.
        baseline = analyze_all_maps(**input_args, electrostatics=False)
        result = analyze_all_maps(**input_args, electrostatics=True)
        if baseline["maps"] != result["maps"]:
            raise AssertionError("APBS changed a canonical 11-map result")
        if baseline["residue_profiles"] != result["residue_profiles"]:
            raise AssertionError("APBS changed canonical residue descriptors")
        if result["n_maps"] != 11 or len(result["map_order"]) != 11:
            raise AssertionError("Missing one or more canonical chemistry maps")
        if result["version"] != MODEL_VERSION:
            raise AssertionError("Unexpected production model version")

        profiles = result["residue_profiles"]
        electro = result["electrostatics"]
        if electro["status"] != "computed":
            raise AssertionError("APBS electrostatics not computed")
        if electro["n_matched_residues"] != len(profiles):
            raise AssertionError("PQR residue coverage is incomplete")
        if {v["key"] for v in electro["residues"]} != {v["key"] for v in profiles}:
            raise AssertionError("PQR residue identities do not match the original profile")
        potentials = [
            row["sidechain_mean_potential_kT_per_e"]
            for row in electro["residues"]
            if row["sidechain_mean_potential_kT_per_e"] is not None
        ]
        if not potentials:
            raise AssertionError("No interpolated APBS side-chain potential")

        work = OUT / "proteins" / pid
        work.mkdir(parents=True, exist_ok=True)
        (work / f"{pid}_prepared.pdb").write_text(pdb, encoding="utf-8")
        with gzip.open(work / "analysis_complete.json.gz", "wt", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, separators=(",", ":"))
        (work / "electrostatics.csv").write_bytes(_csv_bytes(electro["residues"]))
        tables = structure_tables(pid, result)
        tables["APBS_Electrostatics"] = electro["residues"]
        (work / "complete_numeric_results.xlsx").write_bytes(_xlsx_bytes(tables))

        # All 11 0-100 B-factor visualizations plus Top1/3/5 membership PDBs,
        # following the production exporter. No synthetic browser screenshots.
        bfactor = work / "pdb_bfactor"
        bfactor.mkdir(exist_ok=True)
        n_patches = 0
        for k in result["map_order"]:
            m = result["maps"][k]
            n_patches += len(m["patches"])
            permap.append({
                "pdb_id": pid, "map": k, "n_residues": len(m["residues"]),
                "n_patches": len(m["patches"]),
                "n_primary_patches": m["n_pareto_primary_patches"],
                "top1_center_key": m["patches"][0]["center_key"] if m["patches"] else "",
                "top5_center_keys": ";".join(x["center_key"] for x in m["patches"][:5]),
            })
            (bfactor / f"{k}_Bfactor.pdb").write_bytes(
                _bfactor_pdb(pdb, _map_scores(result, k), f"{k}: propensity")
            )
            for top in (1, 3, 5):
                (bfactor / f"{k}_Top{top}_Bfactor.pdb").write_bytes(
                    _bfactor_pdb(pdb, _topk_scores(m, top), f"{k}: Top-{top} patch membership")
                )
        details.update({
            "status": "PASS", "assembly_sha256": actual,
            "prepared_sha256": prepared_sha, "model_version": MODEL_VERSION,
            "software_version": PUBLIC_VERSION, "n_maps": 11,
            "n_residues": len(profiles),
            "n_surface_residues": sum(bool(v["surface_exposed"]) for v in profiles),
            "n_pqr_atoms": electro["n_pqr_atoms"],
            "n_matched_residues": electro["n_matched_residues"],
            "n_patches": n_patches, "canonical_maps_identical": True,
            "residue_profiles_identical": True,
        })
        (work / "provenance.json").write_text(json.dumps(details, indent=2) + "\n", encoding="utf-8")
    except Exception as error:
        details["error_type"] = type(error).__name__
        details["error"] = str(error)[:2000]
        failure_dir = OUT / "proteins" / pid
        failure_dir.mkdir(parents=True, exist_ok=True)
        (failure_dir / "failure.json").write_text(
            json.dumps({"status": details, "traceback": traceback.format_exc()}, indent=2)
            + "\n", encoding="utf-8"
        )
    details["elapsed_s"] = round(time.monotonic() - start, 2)
    return details, permap


def run_shard(index: int, shards: int) -> None:
    if shards <= 0 or not 0 <= index < shards:
        raise ValueError("Invalid shard selection")
    OUT.mkdir(parents=True, exist_ok=True)
    statuses = []
    maps = []
    for i, item in enumerate(MANIFEST):
        if i % shards != index:
            continue
        status, permap = run_one(item)
        statuses.append(status)
        maps.extend(permap)
        print("TV48_APBS_" + status["status"], json.dumps(status), flush=True)
    write_csv(OUT / f"status_{index}.csv", statuses, STATUS_FIELDS)
    write_csv(OUT / f"maps_{index}.csv", maps, MAP_FIELDS)
    print("TV48_APBS_SHARD", json.dumps({
        "shard": index, "total": len(statuses),
        "passes": sum(s["status"] == "PASS" for s in statuses),
        "failures": [s["pdb_id"] for s in statuses if s["status"] != "PASS"],
    }), flush=True)
    # Persist every individual failure and aggregate it in the final CI gate.
    # Do not hide failures by dropping a member or substituting a PDB.


def aggregate(path: Path, shards: int) -> None:
    rows = []
    map_rows = []
    for i in range(shards):
        status_file = path / f"status_{i}.csv"
        map_file = path / f"maps_{i}.csv"
        if not status_file.exists() or not map_file.exists():
            raise RuntimeError(f"Missing shard evidence for index {i}")
        rows.extend(csv.DictReader(status_file.open(encoding="utf-8")))
        map_rows.extend(csv.DictReader(map_file.open(encoding="utf-8")))
    if len(rows) != 48 or set(r["pdb_id"] for r in rows) != set(EXPECTED):
        raise RuntimeError("TV48 dataset coverage or ID uniqueness check failed")
    statuses_ok = [r for r in rows if r["status"] == "PASS"]
    failures = [r for r in rows if r["status"] != "PASS"]
    complete = len(statuses_ok) == 48 and len(map_rows) == 528
    summary = {
        "status": "PASS" if complete else "INCOMPLETE",
        "dataset": "TV48", "software": PUBLIC_VERSION, "model": MODEL_VERSION,
        "conditions": SETTINGS, "n_structures_expected": 48,
        "n_structures_pass": len(statuses_ok),
        "n_structures_fail": len(failures), "n_maps_computed": len(map_rows),
        "n_residues": sum(int(r["n_residues"]) for r in statuses_ok),
        "n_pqr_atoms": sum(int(r["n_pqr_atoms"]) for r in statuses_ok),
        "n_patches": sum(int(r["n_patches"]) for r in statuses_ok),
        "failed": [{"pdb_id": r["pdb_id"], "reason": r["error"]} for r in failures],
        "interpretation": "Technical reproducibility and PQR/APBS execution, not experimental adsorption-site accuracy",
    }
    write_csv(path / "all_status_48.csv", rows, STATUS_FIELDS)
    write_csv(path / "all_maps_528.csv", map_rows, MAP_FIELDS)
    (path / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("TV48_APBS_SUMMARY", json.dumps(summary), flush=True)
    if not complete:
        raise SystemExit("TV48 APBS did not pass 48/48 or 528/528; see SUMMARY.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int)
    parser.add_argument("--shards", type=int, default=8)
    parser.add_argument("--aggregate", type=Path)
    args = parser.parse_args()
    if args.aggregate:
        aggregate(args.aggregate, args.shards)
    elif args.shard is not None:
        run_shard(args.shard, args.shards)
    else:
        parser.error("Supply --shard or --aggregate")
