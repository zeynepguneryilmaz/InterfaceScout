"""Run the frozen TV48 panel through the public InterfaceScout API with real APBS.

Uses the same SHA-checked RCSB biological-assembly input as final_v110_48.py.
No synthetic, cached or substituted electrostatics. Each case is recorded as
PASS or FAIL; all failures survive into the uploaded aggregate report.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "backend"), str(ROOT / "validation")]
import strict_apo_compat_48 as prep
from interface_engine import analyze_all_maps
from model_settings import MODEL_VERSION, PUBLIC_VERSION

MANIFEST = list(csv.DictReader(
    (ROOT / "validation" / "strict_apo_48_selected_candidates_v2_final.csv").open(
        encoding="utf-8", newline=""
    )
))
assert len(MANIFEST) == 48 and len({r["pdb_id"] for r in MANIFEST}) == 48
OUT = Path("tv48_apbs")
OUT.mkdir(exist_ok=True)


def _write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _number(value):
    return round(float(value), 6) if value is not None else None


def analyze_one(item: dict) -> dict:
    pdb_id = item["pdb_id"].upper()
    assembly = prep.prepare(pdb_id)
    cif_path = prep.OUT / f"{pdb_id}-assembly1.cif"
    source_hash = hashlib.sha256(cif_path.read_bytes()).hexdigest()
    if source_hash != item["assembly_sha256"]:
        raise ValueError(f"Frozen assembly SHA-256 mismatch for {pdb_id}: {source_hash}")

    pdb = (prep.OUT / f"{pdb_id}_prepared.pdb").read_text(encoding="utf-8")
    settings = dict(pdb_text=pdb, chain=None, pH=7.4, ionic_mM=150.0, temp_K=298.0)
    baseline = analyze_all_maps(**settings)
    enriched = analyze_all_maps(**settings, electrostatics=True)

    unchanged = (
        baseline["maps"] == enriched["maps"]
        and baseline["residue_profiles"] == enriched["residue_profiles"]
        and baseline["map_order"] == enriched["map_order"]
    )
    if not unchanged:
        raise AssertionError(f"Canonical InterfaceScout results changed with APBS for {pdb_id}")
    if enriched["n_maps"] != 11 or len(enriched["map_order"]) != 11:
        raise AssertionError(f"Not all eleven maps were produced for {pdb_id}")
    profiles = enriched["residue_profiles"]
    electro = enriched["electrostatics"]
    if electro["status"] != "computed" or electro["n_matched_residues"] != len(profiles):
        raise AssertionError(f"Incomplete PDB/PQR correspondence for {pdb_id}")
    if len(electro["residues"]) != len(profiles):
        raise AssertionError(f"Incomplete electrostatics residue table for {pdb_id}")

    sidechain_values = [
        row["sidechain_mean_potential_kT_per_e"]
        for row in electro["residues"]
        if row["sidechain_mean_potential_kT_per_e"] is not None
    ]
    exposed_values = [
        row["sidechain_mean_potential_kT_per_e"]
        for row in electro["residues"]
        if row["surface_exposed"] and row["sidechain_mean_potential_kT_per_e"] is not None
    ]
    if not sidechain_values or not exposed_values:
        raise AssertionError(f"No sampled side-chain potentials on surface for {pdb_id}")

    meta = {
        "pdb_id": pdb_id,
        "status": "PASS",
        "software_version": PUBLIC_VERSION,
        "model_version": MODEL_VERSION,
        "assembly_sha256": source_hash,
        "prepared_pdb_sha256": hashlib.sha256(pdb.encode("utf-8")).hexdigest(),
        "n_chains": assembly["n_chains"],
        "n_maps": 11,
        "n_residues": len(profiles),
        "n_surface_residues": sum(bool(x["surface_exposed"]) for x in profiles),
        "n_pqr_atoms": electro["n_pqr_atoms"],
        "n_potential_residues": len(sidechain_values),
        "n_surface_potential_residues": len(exposed_values),
        "mean_surface_potential_kT_per_e": _number(statistics.mean(exposed_values)),
        "mean_abs_surface_potential_kT_per_e": _number(statistics.mean(abs(v) for v in exposed_values)),
        "min_surface_potential_kT_per_e": _number(min(exposed_values)),
        "max_surface_potential_kT_per_e": _number(max(exposed_values)),
        "sum_pqr_residue_charge_e": _number(sum(x["net_atomic_charge_e"] for x in electro["residues"])),
        "canonical_unchanged": True,
        "conditions": {"pH": 7.4, "ionic_mM": 150.0, "temperature_K": 298.0},
    }
    _write_json(OUT / "proteins" / f"{pdb_id}.json", {
        "summary": meta,
        "electrostatics": electro,
        "note": "APBS potential is a protein-side descriptive annotation (kT/e), NOT adsorption affinity."
    })
    return meta


def run_shard(shard: int, shards: int) -> None:
    rows = []
    for i, item in enumerate(MANIFEST):
        if i % shards != shard:
            continue
        pid = item["pdb_id"].upper()
        try:
            record = analyze_one(item)
        except Exception as exc:
            record = {
                "pdb_id": pid, "status": "FAIL", "error_type": type(exc).__name__,
                "error": str(exc), "traceback_tail": traceback.format_exc()[-4500:],
            }
            _write_json(OUT / "failures" / f"{pid}.json", record)
        rows.append(record)
        print("TV48_APBS_CASE " + json.dumps({
            "pdb_id": pid, "status": record["status"],
            "n_residues": record.get("n_residues"),
            "n_pqr_atoms": record.get("n_pqr_atoms"),
            "error": record.get("error", "")[:400],
        }), flush=True)
    _write_json(OUT / f"shard_{shard}.json", rows)
    print("TV48_APBS_SHARD " + json.dumps({
        "shard": shard, "attempted": len(rows),
        "passed": sum(x["status"] == "PASS" for x in rows),
        "failed": sum(x["status"] != "PASS" for x in rows),
    }), flush=True)
    # Do not block the aggregate/report: failures must be preserved as evidence.


def aggregate(shards: int) -> None:
    files = [OUT / f"shard_{i}.json" for i in range(shards)]
    missing = [p.name for p in files if not p.is_file()]
    records = sum((json.loads(p.read_text(encoding="utf-8")) for p in files if p.is_file()), [])
    by_id = {r["pdb_id"]: r for r in records}
    duplicate_ids = sorted({r["pdb_id"] for r in records if sum(x["pdb_id"] == r["pdb_id"] for x in records) > 1})
    expected = {r["pdb_id"].upper() for r in MANIFEST}
    absent = sorted(expected - by_id.keys())
    unexpected = sorted(by_id.keys() - expected)
    passed = [x for x in records if x["status"] == "PASS"]
    failed = [x for x in records if x["status"] != "PASS"]
    summary = {
        "purpose": "Real external-solver APBS/PROPKA/PDB2PQR on frozen TV48",
        "software_version": PUBLIC_VERSION, "model_version": MODEL_VERSION,
        "conditions": {"pH": 7.4, "ionic_mM": 150.0, "temperature_K": 298.0},
        "n_expected": 48, "n_attempted": len(records),
        "n_pass": len(passed), "n_fail": len(failed),
        "all_48_pass": len(passed) == 48 and not missing and not absent and not duplicate_ids and not unexpected,
        "n_successful_map_instances": 11 * len(passed),
        "n_matched_residues_successful": sum(x["n_residues"] for x in passed),
        "n_surface_potential_residues_successful": sum(x["n_surface_potential_residues"] for x in passed),
        "missing_shard_files": missing, "missing_pdb_ids": absent,
        "unexpected_pdb_ids": unexpected, "duplicate_pdb_ids": duplicate_ids,
        "failed_structures": [{"pdb_id": r["pdb_id"], "error_type": r.get("error_type"),
                               "error": r.get("error")} for r in failed],
        "interpretation": "Technical APBS processing only; no affinity, adsorption-site, or experimental validation.",
    }
    _write_json(OUT / "SUMMARY.json", summary)
    cols = ["pdb_id", "status", "n_chains", "n_residues", "n_surface_residues",
            "n_pqr_atoms", "n_potential_residues", "n_surface_potential_residues",
            "mean_surface_potential_kT_per_e",
            "mean_abs_surface_potential_kT_per_e",
            "min_surface_potential_kT_per_e", "max_surface_potential_kT_per_e",
            "sum_pqr_residue_charge_e", "canonical_unchanged", "error_type", "error"]
    with (OUT / "ALL_48_RESULTS.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for pid in (x["pdb_id"].upper() for x in MANIFEST):
            if pid in by_id:
                w.writerow(by_id[pid])
            else:
                w.writerow({"pdb_id": pid, "status": "NOT_EXECUTED", "error": "Missing result"})
    print("TV48_APBS_SUMMARY " + json.dumps(summary), flush=True)
    if not summary["all_48_pass"]:
        raise SystemExit("TV48 APBS run incomplete or has real solver failures; see artifacts.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int)
    parser.add_argument("--shards", type=int, default=8)
    parser.add_argument("--aggregate", action="store_true")
    args = parser.parse_args()
    if args.aggregate:
        aggregate(args.shards)
    elif args.shard is not None and 0 <= args.shard < args.shards:
        run_shard(args.shard, args.shards)
    else:
        parser.error("Provide --shard n or --aggregate")


if __name__ == "__main__":
    main()
