"""Frozen 48 × 11 glycine policy sensitivity, baseline vs GLY-excluded.
No production changes, no adsorption labels, SHA-verified assembly inputs.
"""
import sys,csv,json,hashlib,argparse,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
import core,geometry,strict_apo_compat_48 as prep,sensitivity_48 as sens
MANIFEST=sens.MANIFEST
OUT=Path("glycine_sensitivity_48")
OUT.mkdir(exist_ok=True)
def compute(row):
    pid=row["pdb_id"].upper()
    prep.prepare(pid)
    cif=prep.OUT/f"{pid}-assembly1.cif"
    sha=hashlib.sha256(cif.read_bytes()).hexdigest()
    assert sha==row["assembly_sha256"],"Source assembly hash mismatch"
    pdb=(prep.OUT/f"{pid}_prepared.pdb").read_text()
    nodes=geometry.extract_ca_nodes(pdb)
    original=core._surface_residues(pdb,7.4)
    gly=[r for r in original if r["res_name"]=="GLY"]
    non_gly=[r for r in original if r["res_name"]!="GLY"]
    assert len(original)==len(non_gly)+len(gly)
    before,bmask=sens.run_cached(pdb,nodes,original,.05,(6.,9.),8.,7.4)
    after,amask=sens.run_cached(pdb,nodes,non_gly,.05,(6.,9.),8.,7.4)
    assert all(not x.startswith("GLY:") for x in amask)
    out=[]
    for channel in sens.CHANNELS:
        a,b=before[channel],after[channel]
        out.append({"protein":pid,"channel":channel,"n_gly_total":len(gly),
          "n_gly_surface":sum(bool(x["surface_exposed"]) for x in gly),
          "n_surface_baseline":len(bmask),"n_surface_no_gly":len(amask),
          "surface_jaccard":sens.jacc(bmask,amask),
          "top5_members_jaccard":sens.jacc(a["members"],b["members"]),
          "top5_centers_jaccard":sens.jacc(a["centers"],b["centers"]),
          "top1_members_jaccard":sens.jacc(a["top1_members"],b["top1_members"]),
          "front1_centers_jaccard":sens.jacc(a["front1_centers"],b["front1_centers"]),
          "n_patches_baseline":a["n_patches"],"n_patches_no_gly":b["n_patches"],
          "n_top5_members_baseline":len(a["members"]),"n_top5_members_no_gly":len(b["members"])})
    return out,{"protein":pid,"assembly_sha256":sha,"n_gly":len(gly),
                "n_gly_surface":sum(bool(x["surface_exposed"]) for x in gly),"status":"PASS"}
def writecsv(p,rows):
    if not rows:return
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def shard(i,n):
    rows=[];logs=[];failed=[]
    for j,item in enumerate(MANIFEST):
        if j%n!=i:continue
        try:
            r,meta=compute(item);rows+=r;logs.append(meta);print("GLY_OK",json.dumps(meta),flush=True)
        except Exception as exc:
            fail={"protein":item["pdb_id"],"error":str(exc),"trace":traceback.format_exc()[-2000:]}
            failed.append(fail);print("GLY_FAIL",json.dumps(fail),flush=True)
    writecsv(OUT/f"channels_{i}.csv",rows)
    (OUT/f"provenance_{i}.json").write_text(json.dumps(logs,indent=2))
    (OUT/f"failures_{i}.json").write_text(json.dumps(failed,indent=2))
    print("GLY_SHARD",i,len(logs),len(failed),flush=True)
    if failed:raise SystemExit(2)
def aggregate(path):
    import pandas as pd
    files=sorted(path.glob("channels_[0-9].csv"))
    assert len(files)==4,files
    d=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
    assert len(d)==48*11 and d.protein.nunique()==48 and d.channel.nunique()==11
    d.to_csv(path/"all_528_comparisons.csv",index=False)
    metrics=["surface_jaccard","top5_members_jaccard","top5_centers_jaccard","top1_members_jaccard","front1_centers_jaccard"]
    by_channel=d.groupby("channel")[metrics].agg(["mean","median","min"])
    by_channel.to_csv(path/"channel_summary.csv")
    by_protein=d.groupby("protein")[metrics+["n_gly_total","n_gly_surface"]].mean().reset_index()
    by_protein.to_csv(path/"protein_summary.csv",index=False)
    rng=np.random.default_rng(20261123)
    boot={}
    for metric in metrics:
        arr=by_protein[metric].to_numpy()
        samples=np.mean(rng.choice(arr,size=(5000,len(arr)),replace=True),axis=1)
        boot[metric]={"mean":float(np.mean(arr)),"median":float(np.median(arr)),
                      "ci95":[float(np.quantile(samples,.025)),float(np.quantile(samples,.975))],
                      "min":float(np.min(arr))}
    summary={"n_proteins":48,"n_channels":11,"n_comparisons":len(d),
       "proteins_with_gly":int((by_protein.n_gly_total>0).sum()),
       "proteins_with_exposed_gly":int((by_protein.n_gly_surface>0).sum()),
       "map_comparisons_changed_top5":int((d.top5_members_jaccard<1-1e-10).sum()),
       "map_comparisons_changed_front1":int((d.front1_centers_jaccard<1-1e-10).sum()),
       "bootstrap_protein_cluster":boot,
       "note":"Deterministic leave-glycine-out sensitivity; NOT an independent SASA implementation or validation of glycine physical exposure."}
    (path/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("GLY_FINAL_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--shard",type=int);p.add_argument("--shards",type=int,default=4);p.add_argument("--aggregate",type=Path)
    a=p.parse_args()
    if a.aggregate:aggregate(a.aggregate)
    elif a.shard is not None:shard(a.shard,a.shards)
    else:p.error("Provide --shard or --aggregate")
