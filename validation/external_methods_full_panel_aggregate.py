"""Aggregate real 76x2 independent competitor runs, including failures.

No synthetic metrics, no adsorption site accuracy. Outputs descriptive task
availability, wall time, sampled RSS and tool-native patch/pocket quantities.
"""
from pathlib import Path
import csv
import json
import math
import numpy as np
import pandas as pd

ROOT=Path("external_methods_full_panel")
OUT=ROOT/"final"
OUT.mkdir(parents=True,exist_ok=True)
manifest=list(csv.DictReader(Path("validation/ISV_active_candidates_76_v3.csv").open(encoding="utf-8")))
expected={v["proposed_pdb_id"] for v in manifest}
assert len(expected)==76
records=[]
summary={"expected_PDBs":76,"tools":{},"interpretation":{
  "no_common_accuracy":"P2Rank pockets are small-ligand sites; PEP-Patch chemistry is different from protein-material adsorption",
  "workload_difference":"Do not calculate direct speed advantage or disadvantage across these distinct tasks",
  "measured_memory":"75-ms sampled process-tree RSS lower bound, not OS-wide peak RAM",
  "failures":"Real nonzero exit, exception or timeout counted as failed, never as zero or omitted result"}}
for tool in ["p2rank","pep_patch"]:
    t=ROOT/tool
    manifest_ids=[]
    shardfiles=list(t.glob("FULL_SHARD_*.json"))
    per=[]
    for file in shardfiles:
        manifest_ids+=json.loads(file.read_text())["pdb_ids"]
    for p in expected:
        f=t/(p+"_result.json")
        if not f.exists():
            per.append({"pdb_id":p,"tool":tool,"overall_status":"MISSING_JOB_OUTPUT",
                        "note":"No result JSON, not zero activity"})
            continue
        item=json.loads(f.read_text())
        taskresults=item.get("attempts",[])
        d={"pdb_id":p,"tool":tool,"overall_status":item.get("overall_status","UNKNOWN"),
           "input_sha256":item.get("prepared_PDB_sha256"),
           "n_output_files":item.get("output_inventory",{}).get("n_files"),
           "total_output_bytes":item.get("output_inventory",{}).get("total_bytes")}
        for attempt in taskresults:
            key=attempt["task"]
            run=attempt["runtime"]
            d[key+"_status"]=run.get("status")
            d[key+"_wall_s"]=run.get("wall_s")
            d[key+"_sampled_peak_RSS_MiB"]=run.get("peak_RSS_sampled_MiB")
            d[key+"_error_tail"]=(run.get("stderr_tail","") or "")[-500:] if run.get("status")!="PASS" else ""
        per.append(d)
    if len(manifest_ids)!=76 or set(manifest_ids)!=expected or len(set(manifest_ids))!=76:
        summary["tools"][tool]={"status":"INCOMPLETE_SHARD_MEMBERSHIP",
                                "n_shards":len(shardfiles),"n_members":len(manifest_ids)}
    df=pd.DataFrame(per).sort_values("pdb_id")
    df.to_csv(OUT/(tool+"_76_CASE_METRICS.csv"),index=False)
    stats={"n_expected":76,"n_explicit_cases":int((df.overall_status!="MISSING_JOB_OUTPUT").sum()),
           "n_pass_complete":int((df.overall_status=="PASS").sum()),
           "n_partial_or_fail":int((df.overall_status=="PARTIAL_OR_FAIL").sum()),
           "n_missing":int((df.overall_status=="MISSING_JOB_OUTPUT").sum())}
    for c in df.columns:
        if c.endswith("_status") and c!="overall_status":
            stats[c+"_counts"]=df[c].value_counts(dropna=False).to_dict()
        if c.endswith("_wall_s") or c.endswith("_sampled_peak_RSS_MiB"):
            series=pd.to_numeric(df.loc[df[c.replace("_wall_s","_status").replace("_sampled_peak_RSS_MiB","_status")]=="PASS",c],errors="coerce")
            if len(series.dropna()):
                stats[c+"_median_PASS"]=round(float(series.median()),4)
                stats[c+"_p90_PASS"]=round(float(series.quantile(.9)),4)
    summary["tools"][tool]={**summary["tools"].get(tool,{}),**stats}
    print("EXTERNAL_FULL_TOOL "+json.dumps(summary["tools"][tool],allow_nan=False),flush=True)
(OUT/"SUMMARY.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf-8")
if any(summary["tools"][x].get("n_explicit_cases",0)!=76 for x in ["p2rank","pep_patch"]):
    raise SystemExit("Not all independent method outcomes available; see final/SUMMARY.json")
