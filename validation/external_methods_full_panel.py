"""Fixed ISV76 external-method full-panel audit.

The real P2Rank/PEP-Patch commands and resource-monitoring function are shared
with preflight; no new competitor scoring, no adsorption-site accuracy claims.

One GitHub Actions shard handles 4-5 structures and retains each outcome,
including failures. Maximal external runner time is an operational limitation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

import external_methods_preflight as compare

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--tool",choices=("p2rank","pep_patch"),required=True)
    ap.add_argument("--shard",type=int,required=True)
    ap.add_argument("--shards",type=int,default=16)
    args=ap.parse_args()
    assert 0<=args.shard<args.shards
    ids=list(compare.ACTIVE)
    assert len(ids)==76 and len(set(ids))==76
    compare.PIDS=tuple(pid for i,pid in enumerate(ids) if i%args.shards==args.shard)
    compare.OUT=Path("external_methods_full_panel").resolve()
    original_run=compare.run
    def bounded_run(cmd,cwd,timeout):
        # A 240s method limit is explicit, not a fake zero or an omitted case.
        return original_run(cmd,cwd,min(timeout,240))
    compare.run=bounded_run
    sys.argv=[str(Path(__file__).resolve()),"--tool",args.tool]
    try:
        compare.main()
    finally:
        d=compare.OUT/args.tool
        assert d.exists()
        # Preserve the metric records, exact input hashes and area/residue
        # reports. Giant APBS DX/mesh files are original solver intermediate
        # artifacts and their sizes are already in per-case inventory.
        # We do not archive them (large storage for 76 proteins).
        for path in d.rglob("*"):
            if path.is_file() and path.suffix.lower() in (
                ".dx",".ply",".npz",".pdb",".gz"
            ):
                path.unlink()
        records=[]
        for path in d.glob("*_result.json"):
            records.append(json.loads(path.read_text()))
        summary={
          "tool":args.tool,"shard":args.shard,
          "n_attempted":len(compare.PIDS),
          "n_pass":sum(q.get("overall_status")=="PASS" for q in records),
          "n_partial_or_fail":sum(q.get("overall_status")!="PASS" for q in records),
          "pdb_ids":list(compare.PIDS),
          "all_failed_cases_retained":True,
          "science_scope":"P2Rank predicts small-molecule pockets; PEP-Patch measures electrostatic/hydrophobic patches; neither predicts material adsorption affinity",
          "resources":"75-ms sampled lower-bound RSS; wall-clock for respective tasks",
          "per_tool_time_limit_s":240,
          "source_comparator_script":"validation/external_methods_preflight.py",
        }
        (d/f"FULL_SHARD_{args.shard}.json").write_text(
            json.dumps(summary,indent=2),encoding="utf-8")
        print("EXTERNAL_FULL_SHARD "+json.dumps(summary),flush=True)

if __name__=="__main__":
    main()
