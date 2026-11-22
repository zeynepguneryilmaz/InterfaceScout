"""Current 33 structure selections + four explicitly non-canonical comparison structures.
Technical InterfaceScout compatibility only: NO experimental adsorption validation.
Uses exact installed InterfaceScout backend and established assembly preparation.
"""
from __future__ import annotations
import argparse
import csv
import json
import subprocess
import sys
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit_literature_33_engine_20261122 as engine

CANONICAL = """1A5R 1AKI 1AML 1AO6 1BEB 1F6R 1FNF 1H3X 1HRC 1HZH
1IGT 1JNJ 1MBN 1TGU 1UBQ 2LZM 2CBA 2GB1 2HAV 2HHB
2OED 2PTN 2VH7 4EUL 3ECA 3GHG 3STN 4F5S 5J17 6O5E
7KWI 7RSA 8ONG""".split()
COMPARISON_ONLY = ["1F6S", "3V83", "1DEQ", "1JY3"]
IDS = CANONICAL + COMPARISON_ONLY
assert len(CANONICAL) == 33 and len(IDS) == len(set(IDS)) == 37
OUT = Path("selected33_engine")
OUT.mkdir(exist_ok=True)
engine.OUT = OUT

def run_worker(pid: str):
    result = engine.work(pid)
    if result["n_maps"] != 11:
        raise RuntimeError("Not all 11 chemistry maps returned")
    if len(result["n_patches_per_map"]) != 11:
        raise RuntimeError("11-map patch evidence missing")
    return result

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", help="Compute prepared PDB by ID")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--shards", type=int, default=6)
    opts = ap.parse_args()
    if opts.worker:
        try:
            print("RESULT_JSON " + json.dumps(run_worker(opts.worker)), flush=True)
        except BaseException:
            traceback.print_exc()
            sys.exit(1)
        return
    assert opts.shards >= 1 and 0 <= opts.shard < opts.shards
    subset = IDS[opts.shard::opts.shards]
    rows = []
    for pid in subset:
        entry = {"pdb": pid, "selection": "CANONICAL" if pid in CANONICAL else "COMPARISON_ONLY",
                 "environment": "pH 7.4; ionic strength 150 mM; 298 K; fixed baseline, NOT assay matched",
                 "interpretation": "Technical maps only; first model, altloc collapsed, protein-only cofactors stripped"}
        try:
            entry.update(engine.prepare(pid))
            proc = subprocess.run(
                [sys.executable, __file__, "--worker", pid],
                text=True, capture_output=True, timeout=900
            )
            if proc.returncode or "RESULT_JSON " not in proc.stdout:
                raise RuntimeError((proc.stderr or proc.stdout or "worker_failed")[-2500:])
            entry.update(json.loads(proc.stdout.split("RESULT_JSON ", 1)[1].splitlines()[0]))
            entry["status"] = "PASS_11MAP_TECHNICAL_ONLY"
        except Exception as e:
            entry["status"] = "FAIL_PREP_OR_ENGINE"
            entry["error"] = f"{type(e).__name__}: {e}"
        rows.append(entry)
        print("SELECTED33", pid, entry["status"], entry.get("n_residues"), str(entry.get("error", ""))[:100], flush=True)
    target = OUT / f"shard_{opts.shard}"
    target.mkdir(exist_ok=True)
    (target / "selected_structure_results.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    cols = ["pdb","selection","status","n_maps","n_residues","n_surface","n_chains",
            "n_patches_per_map","combined_checks","original_sha256","protein_only_sha256",
            "environment","interpretation","preparation_impact","error"]
    with (target / "selected_structure_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, cols)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(row.get(k), ensure_ascii=False)
                             if isinstance(row.get(k), (dict, list)) else row.get(k, "")
                             for k in cols})
    if len(rows) != len(subset):
        raise RuntimeError("Partial shard output")

if __name__ == "__main__":
    main()
