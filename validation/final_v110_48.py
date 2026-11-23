"""Frozen 48-structure verification of the final public InterfaceScout model.

Recomputes the *production* analyze_all_maps API from SHA-checked biological
assemblies; writes per-protein, per-map, per-residue outputs and a deterministic
Pareto-criterion ablation. No experimental labels or old-model cached maps.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
import sensitivity_48 as original
import strict_apo_compat_48 as prep
from interface_engine import analyze_all_maps
from model_settings import MODEL_VERSION,PUBLIC_VERSION
OUT=Path("final_v110_48");OUT.mkdir(exist_ok=True)
OBJECTIVES=["chemistry_support","mean_accessibility","patch_coherence","orientation_coherence"]
def putcsv(path,rows):
    if not rows:return
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def nondominated(rows,keys):
    result=set()
    for i,p in enumerate(rows):
        if not any(j!=i and all(float(q[k])>=float(p[k]) for k in keys)
          and any(float(q[k])>float(p[k]) for k in keys) for j,q in enumerate(rows)):
            result.add(p["center_key"])
    return result
def test(item):
    pid=item["pdb_id"].upper()
    prep.prepare(pid)
    cif=prep.OUT/f"{pid}-assembly1.cif"
    source_sha=hashlib.sha256(cif.read_bytes()).hexdigest()
    if source_sha!=item["assembly_sha256"]:raise ValueError("Frozen biological-assembly SHA changed")
    pdb=(prep.OUT/f"{pid}_prepared.pdb").read_text(encoding="utf-8")
    data=analyze_all_maps(pdb_text=pdb,pH=7.4,chain=None)
    if data["version"]!=MODEL_VERSION or data["n_maps"]!=11:raise ValueError("Incorrect production model or map count")
    if data["method"]["surface_gate"]!="total_residue_RSA":raise ValueError("Not dual accessibility")
    profiles=data["residue_profiles"]
    names={r["key"] for r in profiles}
    surface={r["key"] for r in profiles if r["surface_exposed"]}
    if len(names)!=len(profiles):raise ValueError("Duplicate residue key")
    if surface!=set(data["diagnostics"]["surface_residue_keys"]):raise ValueError("Surface mismatch")
    for r in profiles:
        if r["surface_exposed"]!=(r["trsa_raw"]>=.05):raise ValueError("tRSA gate mismatch")
        if r["res_name"]=="GLY" and (r["scrsa"]!=0 or r["sidechain_sasa"]!=0):
            raise ValueError("GLY chemical contribution nonzero")
    permap=[];ablation=[]
    for ch in data["map_order"]:
        entry=data["maps"][ch];patches=entry["patches"]
        front={p["center_key"] for p in patches if p["pareto_front"]==1}
        if front!=nondominated(patches,OBJECTIVES):raise ValueError("Pareto front mismatch "+ch)
        for p in patches:
            if not set(p["members"]).issubset(surface):raise ValueError("Patch extends beyond exposed geometry")
            vals=[next(r["trsa"] for r in profiles if r["key"]==k) for k in p["members"]]
            if abs(float(p["mean_accessibility"])-float(np.mean(vals)))>1e-7:
                raise ValueError("Patch uses incorrect accessibility "+ch)
        permap.append({"protein":pid,"channel":ch,"n_patches":len(patches),
            "n_primary":len(front),"n_chemistry_residues":len(entry["residues"]),
            "top1_center":patches[0]["center_key"] if patches else "",
            "top5_centers":";".join(p["center_key"] for p in patches[:5]),
            "top5_members":";".join(sorted(set(k for p in patches[:5] for k in p["members"]))),
            "pareto_front_centers":";".join(sorted(front))})
        for omitted in OBJECTIVES:
            selected=nondominated(patches,[k for k in OBJECTIVES if k!=omitted])
            ablation.append({"protein":pid,"channel":ch,"omitted":omitted,
                "baseline_front1_n":len(front),"ablated_front1_n":len(selected),
                "front1_jaccard":original.jacc(front,selected)})
    perres=[{"protein":pid,"key":r["key"],"res_name":r["res_name"],
      "tRSA":r["trsa_raw"],"scRSA":r["scrsa_raw"],
      "surface_exposed":int(r["surface_exposed"])} for r in profiles]
    provenance={"protein":pid,"source_sha256":source_sha,
      "prepared_sha256":hashlib.sha256(pdb.encode()).hexdigest(),
      "model_version":MODEL_VERSION,"software_version":PUBLIC_VERSION,
      "n_residues":len(profiles),"n_surface":len(surface),
      "n_exposed_gly":sum(r["res_name"]=="GLY" and r["surface_exposed"] for r in profiles),
      "n_maps":11,"status":"PASS"}
    return permap,perres,ablation,provenance
def shard(i,n):
    maps=[];res=[];ab=[];prov=[];errors=[]
    for j,item in enumerate(original.MANIFEST):
        if j%n!=i:continue
        try:
            a,b,c,d=test(item)
            maps+=a;res+=b;ab+=c;prov.append(d)
            print("FINAL48_PASS",json.dumps(d),flush=True)
        except Exception as e:
            err={"protein":item["pdb_id"],"error":str(e),"trace":traceback.format_exc()[-1500:]}
            errors.append(err);print("FINAL48_FAIL",json.dumps(err),flush=True)
    putcsv(OUT/f"maps_{i}.csv",maps);putcsv(OUT/f"residues_{i}.csv",res)
    putcsv(OUT/f"ablation_{i}.csv",ab)
    (OUT/f"provenance_{i}.json").write_text(json.dumps(prov,indent=2))
    (OUT/f"failures_{i}.json").write_text(json.dumps(errors,indent=2))
    if errors:raise SystemExit(2)
def aggregate(path):
    import pandas as pd
    def load(prefix):
        files=sorted(path.glob(prefix+"_[0-9].csv"))
        if len(files)!=4:raise ValueError(f"Missing {prefix} shards")
        return pd.concat([pd.read_csv(p) for p in files],ignore_index=True)
    m=load("maps");r=load("residues");a=load("ablation")
    if len(m)!=528 or m.protein.nunique()!=48 or m.channel.nunique()!=11 or len(a)!=2112:
        raise ValueError("Incomplete final 48x11 verification")
    m.to_csv(path/"all_528_maps.csv",index=False)
    r.to_csv(path/"all_residue_tRSA_scRSA.csv",index=False)
    a.to_csv(path/"all_pareto_ablations.csv",index=False)
    files=sorted(path.glob("provenance_[0-9].json"))
    prov=sum([json.loads(p.read_text()) for p in files],[])
    if len(prov)!=48 or any(p["status"]!="PASS" for p in prov):
        raise ValueError("Missing SHA-provenance passes")
    (path/"provenance.json").write_text(json.dumps(prov,indent=2))
    summary={"status":"PASS","model":MODEL_VERSION,"software":PUBLIC_VERSION,
      "n_proteins":48,"n_maps":528,"n_residues":len(r),"n_ablations":len(a),
      "n_exposed_gly":int(((r.res_name=="GLY") & (r.surface_exposed==1)).sum()),
      "per_channel_patch_counts":m.groupby("channel").n_patches.agg(["mean","min","max"]).to_dict("index"),
      "pareto_ablation_mean_jaccard":a.groupby("omitted").front1_jaccard.mean().to_dict(),
      "interpretation":"Computational correctness and descriptor changes only; no experimental predictive validation"}
    (path/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("FINAL48_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--shard",type=int)
    p.add_argument("--shards",type=int,default=4);p.add_argument("--aggregate",type=Path)
    args=p.parse_args()
    if args.aggregate:aggregate(args.aggregate)
    elif args.shard is not None:shard(args.shard,args.shards)
    else:p.error("Specify --shard or --aggregate")
