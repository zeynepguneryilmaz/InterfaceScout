"""External FreeSASA vs InterfaceScout: independent SASA and runtime audit.

This compares matching protein coordinates and residue IDs. Not adsorption
accuracy or a head-to-head speed competition: FreeSASA is SASA-only;
InterfaceScout runs eleven maps and Pareto patches (optional APBS disabled).
Different default atomic radii are retained/documented as confounders.

Per-case evidence: successful/failed command, subprocess peak RSS, elapsed
wall time, actual FreeSASA JSON parameters, residue-level SASA deviations,
exposed-residue agreement. Failure never becomes an invented score.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
import numpy as np
import psutil
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
from isv_full_80_qa import prepare
from core import TOTAL_RESIDUE_REF_ASA

OUT=Path("isv76_external_freesasa_audit")
OUT.mkdir(exist_ok=True)
ACTIVE=list(csv.DictReader(
    (ROOT/"validation"/"ISV_active_candidates_76_v3.csv").open(
        encoding="utf-8",newline="")))
assert len(ACTIVE)==76 and len({r["proposed_pdb_id"] for r in ACTIVE})==76


def write_json(p,value):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                 encoding="utf-8")
def write_csv(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    if not rows:raise ValueError("Empty residue comparison rows")
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def monitored(cmd,cwd=None,timeout=900):
    """Independent subprocess wall clock and peak resident memory.

    The process tree is polled, so peak RSS is a sampled lower bound, not
    a hardware-independent guarantee. No GPU tasks involved.
    """
    started=time.perf_counter()
    # Output is written to files: a large FreeSASA JSON must not deadlock
    # in a subprocess PIPE while the parent is sampling process memory.
    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as so, \
         tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as se:
        proc=psutil.Popen(cmd,cwd=str(cwd) if cwd else None,
                          stdout=so,stderr=se,text=True)
        peak=0
        while proc.poll() is None:
            try:
                children=proc.children(recursive=True)
                mem=proc.memory_info().rss+sum(
                    p.memory_info().rss for p in children if p.is_running())
                peak=max(peak,mem)
            except (psutil.Error,OSError):pass
            if time.perf_counter()-started>timeout:
                proc.kill()
                proc.wait()
                raise TimeoutError("Tool timed out after "+str(timeout)+" seconds")
            time.sleep(.075)
        proc.wait()
        so.seek(0);se.seek(0)
        stdout=so.read();stderr=se.read()
    return {
        "returncode":proc.returncode,
        "wall_seconds":round(time.perf_counter()-started,5),
        "sampled_peak_rss_MiB":round(peak/1048576,3),
        "stdout":stdout,
        "stderr":stderr[-2400:],
        "max_rss_is_sampled_lower_bound":True,
    }


def worker(pdb_path, out_path):
    from interface_engine import analyze_all_maps
    prepared=pdb_path.read_text(encoding="utf-8")
    result=analyze_all_maps(pdb_text=prepared,pH=7.4,ionic_mM=150.,
                            temp_K=298.,electrostatics=False)
    profiles=result["residue_profiles"]
    if result["n_maps"]!=11:raise ValueError("Missing map in production engine")
    write_json(out_path, {
        "n_maps":result["n_maps"],"model_version":result["version"],
        "profiles":[{
            "key":r["key"],"res_name":r["res_name"],
            "sasa_A2":r["residue_sasa"],"sidechain_sasa_A2":r["sidechain_sasa"],
            "exposed":r["surface_exposed"],
        } for r in profiles],
        "n_all_patches":sum(len(m["patches"]) for m in result["maps"].values()),
    })


def freesasa_parse(payload):
    if not isinstance(payload.get("results"),list) or not payload["results"]:
        raise ValueError("FreeSASA JSON has no results")
    data={}
    for structure in payload["results"][0].get("structures",[]):
        for chain in structure.get("chains",[]):
            for rec in chain.get("residues",[]):
                m=re.fullmatch(r"(-?\d+)([A-Za-z]?)",str(rec["number"]).strip())
                if m is None:raise ValueError("Unexpected FreeSASA residue number")
                key=f'{chain["label"]}:{int(m.group(1))}:{m.group(2)}'
                if key in data:raise ValueError("Duplicate FreeSASA residue "+key)
                area=rec["area"]
                data[key]={
                    "res_name":rec["name"],
                    "sasa_A2":float(area["total"]),
                    "sidechain_sasa_A2":float(area["side-chain"]),
                }
    if not data:
        structures=payload.get("results",[{}])[0].get("structures",[])
        raise ValueError(
            "FreeSASA JSON has no parsed residue records; "
            "result_keys="+str(sorted(payload["results"][0].keys()))+
            " structure_keys="+str([sorted(q.keys()) for q in structures[:1]])+
            " chain_keys="+str([[sorted(ch.keys()) for ch in q.get("chains",[])[:1]]
                                for q in structures[:1]])+
            " preview="+json.dumps(payload["results"][0],ensure_ascii=False)[:1800]
        )
    return data


def one(row):
    pid=row["proposed_pdb_id"]
    meta,pdbpath=prepare(row)
    from prepare import prepare_pdb_text
    src=pdbpath.read_text(encoding="utf-8")
    normalized,_=prepare_pdb_text(src)
    root=OUT/"proteins"/pid
    root.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        inp=Path(td)/"structure.pdb"
        inp.write_text(normalized,encoding="utf-8")
        exported=Path(td)/"interface.json"
        own=monitored([sys.executable,str(Path(__file__).resolve()),"--worker",
                       str(inp),str(exported)],timeout=1200)
        if own["returncode"]:
            raise RuntimeError("InterfaceScout failed: "+own["stderr"][-1000:])
        data=json.loads(exported.read_text())
        # Ubuntu apt FreeSASA 2.1.2 implements --depth=residue. Older
        # upstream online docs mention --output-depth; do not use that
        # unsupported alias with this pinned system-package environment.
        cmd=["freesasa","--shrake-rupley","-n","200","--probe-radius","1.4",
             "--n-threads","1","--format=json","--depth=residue",str(inp)]
        external=monitored(cmd,timeout=1200)
        if external["returncode"]:
            raise RuntimeError("FreeSASA failed: "+external["stderr"][-1200:])
        payload=json.loads(external["stdout"])
    fs=freesasa_parse(payload)
    params=payload["results"][0].get("parameters",{})
    if "Shrake" not in str(params.get("algorithm","")):
        raise ValueError("External solver algorithm mismatch "+str(params))
    if abs(float(params.get("probe-radius",-1))-1.4)>.001:
        raise ValueError("External solver probe mismatch")
    rows=[]
    miss=[]
    for p in data["profiles"]:
        key=p["key"]
        b=fs.get(key)
        if b is None:
            miss.append(key);continue
        if p["res_name"] != b["res_name"]:
            raise ValueError(pid+" mismatched residue identity "+key)
        ref=TOTAL_RESIDUE_REF_ASA[p["res_name"]]
        classification=b["sasa_A2"]/ref >= .05
        rows.append({
            "pdb_id":pid,"key":key,"res_name":p["res_name"],
            "InterfaceScout_total_SASA_A2":p["sasa_A2"],
            "FreeSASA_total_SASA_A2":b["sasa_A2"],
            "InterfaceScout_sidechain_SASA_A2":p["sidechain_sasa_A2"],
            "FreeSASA_sidechain_SASA_A2":b["sidechain_sasa_A2"],
            "InterfaceScout_tRSA_exposed":bool(p["exposed"]),
            "FreeSASA_tRSA_exposed_same_Tien_reference":bool(classification),
        })
    if len(rows)<.95*len(data["profiles"]):
        raise ValueError(pid+": matched "+str(len(rows))+" of "+
                         str(len(data["profiles"]))+" residues; missing "+str(miss[:20]))
    a=np.array([r["InterfaceScout_total_SASA_A2"] for r in rows])
    b=np.array([r["FreeSASA_total_SASA_A2"] for r in rows])
    aa=np.array([r["InterfaceScout_sidechain_SASA_A2"] for r in rows])
    bb=np.array([r["FreeSASA_sidechain_SASA_A2"] for r in rows])
    agree=np.mean([r["InterfaceScout_tRSA_exposed"]==
                   r["FreeSASA_tRSA_exposed_same_Tien_reference"] for r in rows])
    result={
        "pdb_id":pid,"category":row["category"],"status":"PASS",
        "input_sha256":hashlib.sha256(normalized.encode()).hexdigest(),
        "n_original_residues":len(data["profiles"]),
        "n_matched_residues":len(rows),
        "n_missing_from_FreeSASA":len(miss),
        "n_maps":data["n_maps"],"model_version":data["model_version"],
        "n_patches":data["n_all_patches"],
        "FreeSASA_source":payload.get("source"),
        "FreeSASA_classifier":payload["results"][0].get("classifier"),
        "FreeSASA_parameters":params,
        "residue_SASA_spearman":float(spearmanr(a,b).statistic),
        "residue_SASA_median_absolute_difference_A2":float(np.median(abs(a-b))),
        "residue_SASA_95p_absolute_difference_A2":float(np.quantile(abs(a-b),.95)),
        "sidechain_SASA_spearman":float(spearmanr(aa,bb).statistic)
                                  if np.std(aa)>1e-12 and np.std(bb)>1e-12 else None,
        "same_tRSA_exposed_classification_fraction":float(agree),
        "InterfaceScout_11map_no_APBS_wall_s":own["wall_seconds"],
        "FreeSASA_SASA_only_wall_s":external["wall_seconds"],
        "InterfaceScout_sampled_peak_RSS_MiB":own["sampled_peak_rss_MiB"],
        "FreeSASA_sampled_peak_RSS_MiB":external["sampled_peak_rss_MiB"],
        "timing_note":"Different workloads: do not calculate head-to-head speed ratio",
        "memory_note":"75 ms polling of process RSS; lower bound on peaks",
    }
    write_csv(root/"MATCHED_RESIDUE_SASA.csv",rows)
    write_json(root/"SUMMARY.json",result)
    return result


def shard(i,n):
    records=[]
    for idx,row in enumerate(ACTIVE):
        if idx%n!=i:continue
        try:r=one(row)
        except Exception as exc:
            r={"pdb_id":row["proposed_pdb_id"],"status":"FAIL",
               "error":str(exc),"traceback":traceback.format_exc()[-4000:]}
            write_json(OUT/"failures"/(r["pdb_id"]+".json"),r)
        records.append(r)
        print("EXTERNAL_SASA_CASE "+json.dumps({
            "pdb_id":r["pdb_id"],"status":r["status"],
            "n_residues":r.get("n_matched_residues",0),
            "rho":r.get("residue_SASA_spearman"),
            "seconds_11maps":r.get("InterfaceScout_11map_no_APBS_wall_s"),
            "seconds_external_SASA":r.get("FreeSASA_SASA_only_wall_s"),
            "error":r.get("error","")[:150]}),flush=True)
    write_json(OUT/f"shard_{i}.json",records)
    if any(r["status"]!="PASS" for r in records):
        raise SystemExit("Comparison shard incomplete")


def aggregate(shards):
    import pandas as pd
    rows=[]
    for idx in range(shards):
        p=OUT/f"shard_{idx}.json"
        if not p.exists():raise ValueError("Missing shard "+str(idx))
        rows+=json.loads(p.read_text())
    good=[r for r in rows if r["status"]=="PASS"]
    bad=[r for r in rows if r["status"]!="PASS"]
    write_json(OUT/"COMPLETENESS.json",{
        "n_expected":76,"n_run":len(rows),"n_pass":len(good),
        "n_fail":len(bad),"failures":bad})
    if len(rows)!=76 or len({r["pdb_id"] for r in rows})!=76:
        raise ValueError("Not all 76 PDBs accounted for")
    df=pd.DataFrame(good)
    if not df.empty:
        df.to_csv(OUT/"ALL_PASS_PROTEIN_COMPARISONS.csv",index=False)
        def med(k):
            return float(df[k].median())
        summary={
            "status":"PASS" if not bad else "PARTIAL",
            "n_structures":len(rows),"n_pass":len(good),
            "n_fail":len(bad),
            "median_protein_residue_SASA_spearman":med("residue_SASA_spearman"),
            "median_protein_total_SASA_median_absolute_delta_A2":
                med("residue_SASA_median_absolute_difference_A2"),
            "median_protein_exposure_agreement":med("same_tRSA_exposed_classification_fraction"),
            "median_11map_wall_s_no_APBS":med("InterfaceScout_11map_no_APBS_wall_s"),
            "median_freesasa_SASA_only_wall_s":med("FreeSASA_SASA_only_wall_s"),
            "median_11map_sampled_peak_RSS_MiB":med("InterfaceScout_sampled_peak_RSS_MiB"),
            "median_freesasa_sampled_peak_RSS_MiB":med("FreeSASA_sampled_peak_RSS_MiB"),
            "runtime_caveat":"Different workloads! 11 maps and Pareto versus independent SASA only",
            "radii_caveat":"Biopython element radii versus FreeSASA ProtOr default",
            "surface_metrics_only":"Not affinity or contact-site prediction benchmark",
        }
        write_json(OUT/"SUMMARY.json",summary)
        print("EXTERNAL_SASA_FINAL "+json.dumps(summary),flush=True)
    if bad:raise SystemExit("Some external tool runs incomplete; preserve failures")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--worker",nargs=2)
    p.add_argument("--shard",type=int)
    p.add_argument("--shards",type=int,default=16)
    p.add_argument("--aggregate",action="store_true")
    a=p.parse_args()
    if a.worker:worker(Path(a.worker[0]),Path(a.worker[1]))
    elif a.aggregate:aggregate(a.shards)
    elif a.shard is not None and 0<=a.shard<a.shards:shard(a.shard,a.shards)
    else:p.error("provide worker, shard, or aggregate")


if __name__=="__main__":
    main()
