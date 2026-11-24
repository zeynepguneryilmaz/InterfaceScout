"""Run the original SURFMAP native CLI on the fixed 76-case ISV panel.

Real SURFMAP matrices/maps are required for PASS. A command exit code alone
is insufficient. Failed cases retain stdout/stderr and exact command.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, subprocess, sys, time, traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"validation"),str(ROOT/"backend")]
from isv_full_80_qa import prepare
from prepare import prepare_pdb_text

MANIFEST=list(csv.DictReader((ROOT/"validation"/"ISV_active_candidates_76_v3.csv").open(encoding="utf-8",newline="")))
assert len(MANIFEST)==len({r["proposed_pdb_id"] for r in MANIFEST})==76
OUT=Path("surfmap_ISV76_native")
FEATURES=("kyte_doolittle","wimley_white","stickiness","circular_variance")
def save(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")

def run(pid, row, modes, per_task_timeout):
    home=OUT/"proteins"/pid;home.mkdir(parents=True,exist_ok=True)
    record={"pdb_id":pid,"status":"FAIL","modes":{},"source":"Original SURFMAP CLI, native (no Docker)"}
    try:
        meta, raw=prepare(row)
        prepared, _=prepare_pdb_text(raw.read_text(encoding="utf-8"))
        inp=home/"input_prepared.pdb";inp.write_text(prepared,encoding="utf-8")
        record["input_sha256"]=hashlib.sha256(prepared.encode()).hexdigest()
        record["original_assembly_sha256"]=meta.get("assembly_cif_sha256")
        for mode in modes:
            destination=home/mode
            cmd=["surfmap","-pdb",str(inp.resolve()),"-tomap",mode,
                 "-d",str(destination.resolve()),"--png"]
            t0=time.perf_counter()
            info={"command":cmd,"status":"FAIL"}
            try:
                result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                    text=True,errors="replace",timeout=per_task_timeout,check=False)
                info.update(returncode=result.returncode,
                    wall_s=round(time.perf_counter()-t0,3))
                (home/(mode+".stdout.log")).write_text(result.stdout,encoding="utf-8")
                (home/(mode+".stderr.log")).write_text(result.stderr,encoding="utf-8")
                files=[{"path":str(p.relative_to(home)),"bytes":p.stat().st_size}
                    for p in destination.rglob("*") if p.is_file()]
                matrices=[x for x in files if "matrix" in x["path"].lower() and
                    x["path"].endswith(".txt") and x["bytes"]>100]
                maps=[x for x in files if x["path"].lower().endswith((".png",".pdf"))
                    and x["bytes"]>100]
                info.update(output_inventory=files,
                    n_matrices=len(matrices),n_maps=len(maps))
                info["status"]="PASS" if result.returncode==0 and matrices and maps else "INVALID_OR_MISSING_SCIENTIFIC_OUTPUT"
            except subprocess.TimeoutExpired as exc:
                info.update(status="TIMEOUT",wall_s=round(time.perf_counter()-t0,3),
                    error=str(exc))
            except Exception as exc:
                info.update(status="EXCEPTION",error=str(exc),
                    wall_s=round(time.perf_counter()-t0,3))
            record["modes"][mode]=info
        record["status"]="PASS" if all(
            record["modes"].get(x,{}).get("status")=="PASS" for x in modes) else "PARTIAL_OR_FAIL"
    except Exception as exc:
        record.update(status="PREPARATION_FAILURE",error=str(exc),
            traceback=traceback.format_exc()[-3000:])
    save(home/"SUMMARY.json",record)
    print("SURFMAP_CASE "+json.dumps({"pdb":pid,"status":record["status"],
        "modes":{k:v["status"] for k,v in record["modes"].items()}},ensure_ascii=False),flush=True)
    return record

def aggregate(shards, modes):
    expected={r["proposed_pdb_id"] for r in MANIFEST}
    records=[]
    for shard in range(shards):
        file=OUT/("shard_"+str(shard)+".json")
        if not file.exists():
            raise SystemExit("Missing results from shard "+str(shard))
        records.extend(json.loads(file.read_text(encoding="utf-8")))
    observed=[r["pdb_id"] for r in records]
    counts={m:{s:sum(r.get("modes",{}).get(m,{}).get("status","MISSING")==s
                   for r in records) for s in ("PASS","INVALID_OR_MISSING_SCIENTIFIC_OUTPUT","TIMEOUT","EXCEPTION","MISSING")} for m in modes}
    summary={"expected":76,"observed":len(records),
             "unique":len(set(observed)),
             "same_frozen_manifest":set(observed)==expected,
             "n_fully_passed":sum(r["status"]=="PASS" for r in records),
             "n_failed_or_partial":sum(r["status"]!="PASS" for r in records),
             "mode_counts":counts,
             "interpretation":"Surface feature projection, not material adsorption accuracy"}
    save(OUT/"final"/"SUMMARY.json",summary)
    out=OUT/"final"/"ALL_76_RESULTS.csv";out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["pdb_id","status"]+list(modes)+["error"])
        w.writeheader()
        for r in sorted(records,key=lambda v:v["pdb_id"]):
            w.writerow({"pdb_id":r["pdb_id"],"status":r["status"],
                **{m:r.get("modes",{}).get(m,{}).get("status","MISSING") for m in modes},
                "error":r.get("error","")})
    print("SURFMAP_ISV76_FINAL "+json.dumps(summary),flush=True)
    if not summary["same_frozen_manifest"] or summary["unique"]!=76 or summary["n_fully_passed"]!=76:
        raise SystemExit("SURFMAP 76-case verification incomplete; see retained evidence")

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--shard",type=int)
    a.add_argument("--shards",type=int,default=16)
    a.add_argument("--aggregate",action="store_true")
    a.add_argument("--electrostatics",action="store_true",
        help="Run electrostatics separately from native four-map all mode")
    a.add_argument("--timeout",type=int,default=1200)
    opt=a.parse_args()
    modes=("all","electrostatics") if opt.electrostatics else ("all",)
    if opt.aggregate:return aggregate(opt.shards,modes)
    if opt.shard is None or not 0<=opt.shard<opt.shards:a.error("valid --shard required")
    selected=[r for i,r in enumerate(MANIFEST) if i%opt.shards==opt.shard]
    results=[run(r["proposed_pdb_id"],r,modes,opt.timeout) for r in selected]
    save(OUT/("shard_"+str(opt.shard)+".json"),results)
    if any(r["status"]!="PASS" for r in results):
        raise SystemExit("SURFMAP shard contains real failures; inspect artifact")

if __name__=="__main__":main()
