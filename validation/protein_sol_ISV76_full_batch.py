"""Protein-sol Patches: complete frozen ISV76 batch with real result downloads.

Submits public PDBs via the validated Protein-sol POST and extracts the
server-generated JavaScript timestamp, then fetches both genuine
potlfield.pdb and ratio.pdb as specified by the provider NGL viewer.
Runs sequentially inside each shard, with deliberately low concurrency.
Every structure produces a manifest row even after HTTP, parsing, timeout,
or scientific service failures. This is NOT local PB or simulated output.

Protein-sol results expire after seven days: archive immediately.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re,sys,time,traceback
from pathlib import Path
from urllib.parse import urlparse
from requests import Session
from bs4 import BeautifulSoup
sys.path[:0]=[str(Path(__file__).resolve().parents[1]/"validation"),
              str(Path(__file__).resolve().parents[1]/"backend")]
from isv_full_80_qa import prepare
from prepare import prepare_pdb_text

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=list(csv.DictReader((ROOT/"validation"/"ISV_active_candidates_76_v3.csv").open(encoding="utf-8",newline="")))
assert len(MANIFEST)==76 and len({x["proposed_pdb_id"] for x in MANIFEST})==76
BASE="https://protein-sol.manchester.ac.uk"
POST=BASE+"/cgi-bin/patches/patches.php"
OUT=Path("protein_sol_ISV76").resolve()
AGENT="InterfaceScout-research-comparison/1.0; public-structure benchmarking"
MAX_RETRIES=6
POLL_SECONDS=(2,4,6,9,13,18)

def atomic_json(path,obj):
    path.parent.mkdir(exist_ok=True,parents=True)
    tmp=path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    tmp.replace(path)

def safe_request(s,method,url,**kwargs):
    p=urlparse(url)
    if p.scheme!="https" or p.hostname!="protein-sol.manchester.ac.uk":
        raise RuntimeError("Refuse non-provider download or upload: "+url)
    resp=s.request(method,url,timeout=(15,60),allow_redirects=True,**kwargs)
    p2=urlparse(resp.url)
    if p2.hostname!="protein-sol.manchester.ac.uk":
        raise RuntimeError("Unexpected final host "+resp.url)
    return resp

def valid_pdb_bytes(data):
    if len(data)<200:return False
    lines=data.splitlines()
    return sum(x.startswith((b"ATOM  ",b"HETATM")) for x in lines)>=4

def parse_pdb(path):
    """Residue-aligned B factor distribution from provider PDB; no units assumed."""
    group={}
    atoms=0
    for line in path.read_text(encoding="utf-8",errors="replace").splitlines():
        if not line.startswith(("ATOM  ","HETATM")) or len(line)<66:continue
        try:
            rid=f"{line[21].strip()}:{int(line[22:26].strip())}:{line[26].strip()}"
            val=float(line[60:66].strip())
        except ValueError:continue
        entry=group.setdefault(rid,{"aa":line[17:20].strip(),"bfactor":[]})
        entry["bfactor"].append(val);atoms+=1
    if atoms<4:return {"valid":False,"n_atoms":atoms}
    import statistics
    values=[v for row in group.values() for v in row["bfactor"]]
    return {
        "valid":True,"n_atoms":atoms,"n_residues":len(group),
        "bfactor_min":min(values),"bfactor_max":max(values),
        "bfactor_median":statistics.median(values),
        "n_unique_bfactors":len({round(v,2) for v in values}),
        "bfactor_interpretation":"Raw exported PDB B factor field; physical units not assumed pending provenance",
    }

def process(row):
    pid=row["proposed_pdb_id"]
    d=OUT/"proteins"/pid;d.mkdir(parents=True,exist_ok=True)
    t0=time.monotonic()
    result={"pdb_id":pid,"status":"NOT_RUN",
        "service":"Protein-sol Patches",
        "parameters":"Service-controlled settings (not identical to InterfaceScout pH assumptions)",
        "service_server_RAM_CPU":"NOT_OBSERVABLE_EXTERNALLY",
        "scientific_accuracy":"NOT_AN_ADSORPTION_PREDICTION_METRIC"}
    try:
        metadata,raw=prepare(row)
        prepared,_=prepare_pdb_text(raw.read_text(encoding="utf-8"))
        b=prepared.encode()
        sourcehash=hashlib.sha256(b).hexdigest()
        result.update({"input_sha256":sourcehash,"input_bytes":len(b),
                       "original_assembly_sha256":metadata["assembly_cif_sha256"]})
        (d/"input_prepared.pdb").write_bytes(b)
        with Session() as s:
            s.headers.update({"User-Agent":AGENT})
            posted=safe_request(s,"POST",POST,
                files={"uploadedFile":(pid+".pdb",b,"chemical/x-pdb")},
                data={"multiprediction":"Upload"})
            result["post_status"]=posted.status_code
            result["post_final_url"]=posted.url
            result["post_elapsed_s"]=round(time.monotonic()-t0,3)
            result["post_response_bytes"]=len(posted.content)
            if posted.status_code!=200:
                raise RuntimeError("Upload HTTP "+str(posted.status_code))
            soup=BeautifulSoup(posted.content,"html.parser")
            (d/"upload_response.html").write_bytes(posted.content)
            rawhtml=posted.text
            m=re.search(r'var\s+timestamp\s*=\s*["\']([0-9a-zA-Z_-]+)["\']',rawhtml)
            if not m:
                raise RuntimeError("No server job timestamp in upload HTML; no inferred URL")
            jobid=m.group(1)
            if not re.fullmatch(r"[a-zA-Z0-9_-]{8,50}",jobid):
                raise ValueError("Invalid job token")
            result["job_id"]=jobid
            result_url=f"{BASE}/results/patches/run-{jobid}/results.html"
            result["result_url"]=result_url
            resolved=False
            result["poll_attempts"]=[]
            for attempt,delay in enumerate((0,)+POLL_SECONDS):
                if delay:time.sleep(delay)
                try:
                    r=safe_request(s,"GET",result_url)
                    valid=(r.status_code==200 and
                        "Potential field PDB" in r.text and
                        "NPP ratio PDB" in r.text and
                        jobid in r.text)
                    result["poll_attempts"].append({"attempt":attempt,
                        "status_code":r.status_code,"bytes":len(r.content),
                        "both_output_labels":valid})
                    if valid:
                        (d/"results.html").write_bytes(r.content)
                        resolved=True;break
                except Exception as e:
                    result["poll_attempts"].append({"attempt":attempt,"error":str(e)[:200]})
            if not resolved:
                result["status"]="SERVER_JOB_NOT_READY_OR_FAILED"
            else:
                results={}
                for mode in ("potlfield","ratio"):
                    target=f"{BASE}/results/patches/run-{jobid}/{mode}.pdb"
                    a=safe_request(s,"GET",target)
                    info={"http_status":a.status_code,"bytes":len(a.content),
                          "content_type":a.headers.get("Content-Type"),
                          "source_url":target}
                    if a.status_code==200 and valid_pdb_bytes(a.content):
                        p=d/(mode+".pdb");p.write_bytes(a.content)
                        info["sha256"]=hashlib.sha256(a.content).hexdigest()
                        info.update(parse_pdb(p))
                        info["status"]="VALID_PDB" if info.get("valid") else "INVALID_PDB"
                    else:
                        info["status"]="INVALID_OR_MISSING_PDB"
                        info["preview"]=a.text[:250]
                    results[mode]=info
                result["provider_outputs"]=results
                result["status"]=(
                    "SUCCESS_BOTH_PDB" if all(v["status"]=="VALID_PDB" for v in results.values())
                    else "PARTIAL_OR_BAD_PROVIDER_OUTPUT")
    except Exception as exc:
        result["error"]=str(exc)[:850]
        result["exception_type"]=type(exc).__name__
        if result["status"]=="NOT_RUN":
            result["status"]="FAILED"
        result["traceback_tail"]=traceback.format_exc()[-1400:]
    result["elapsed_s"]=round(time.monotonic()-t0,3)
    atomic_json(d/"SUMMARY.json",result)
    print("PROTEIN_SOL_76_CASE "+json.dumps({
        "pdb":pid,"status":result["status"],"elapsed_s":result["elapsed_s"],
        "job_id":result.get("job_id"),"error":result.get("error")},ensure_ascii=False),flush=True)
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--shard",type=int,required=True)
    ap.add_argument("--shards",type=int,default=4)
    a=ap.parse_args()
    if not 0<=a.shard<a.shards:ap.error("Invalid shard")
    selected=[x for i,x in enumerate(MANIFEST) if i%a.shards==a.shard]
    summary={"shard":a.shard,"n_selected":len(selected),"case_outcomes":[]}
    for i,row in enumerate(selected):
        if i:
            # Explicit politeness delay; jobs also use max-parallel=1.
            time.sleep(4)
        x=process(row)
        summary["case_outcomes"].append(x)
        atomic_json(OUT/("shard_"+str(a.shard)+".json"),summary)
    summary.update({"n_success":sum(x["status"]=="SUCCESS_BOTH_PDB" for x in summary["case_outcomes"]),
                    "n_partial_or_failed":sum(x["status"]!="SUCCESS_BOTH_PDB" for x in summary["case_outcomes"])})
    atomic_json(OUT/("shard_"+str(a.shard)+".json"),summary)
    print("PROTEIN_SOL_76_SHARD_SUMMARY "+json.dumps({
        "shard":a.shard,"completed":len(summary["case_outcomes"]),
        "success":summary["n_success"],"problems":summary["n_partial_or_failed"]}),flush=True)
    # Preserve all error evidence, but do not mislabel submitted jobs success.
    if summary["n_partial_or_failed"]:
        raise SystemExit("Some proteins failed or returned incomplete scientific output")

if __name__=="__main__":
    main()
