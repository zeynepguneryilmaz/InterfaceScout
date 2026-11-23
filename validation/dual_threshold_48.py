"""48 structure total-RSA threshold sensitivity. Prespecified grid, no labels.
Uses dual geometry/side-chain model, preserving non-GLY chemistry and GLY geometry.
"""
import sys,csv,json,argparse,hashlib,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
import core,geometry,strict_apo_compat_48 as prep,sensitivity_48 as sens,dual_accessibility_48 as dual
OUT=Path("dual_threshold_48");OUT.mkdir(exist_ok=True)
GRID=(.01,.03,.05,.075,.10)
def dump(p,items):
    if items:
        with p.open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
def one(item):
    pid=item["pdb_id"].upper()
    prep.prepare(pid)
    cif=prep.OUT/f"{pid}-assembly1.cif"
    if hashlib.sha256(cif.read_bytes()).hexdigest()!=item["assembly_sha256"]:
        raise ValueError("Frozen assembly sha mismatch")
    pdb=(prep.OUT/f"{pid}_prepared.pdb").read_text()
    nodes=geometry.extract_ca_nodes(pdb)
    raw=core._surface_residues(pdb,7.4)
    dual_rows,residues=dual.dual_rows(pdb,raw)
    base_maps,base_surf=dual.run_dual(pdb,dual_rows,nodes)
    comparisons=[]
    for t in GRID:
        alt=[dict(r,surface_exposed=bool(r["trsa_raw"]>=t)) for r in dual_rows]
        amap,asurf=dual.run_dual(pdb,alt,nodes)
        for ch in sens.CHANNELS:
            a,b=base_maps[ch],amap[ch]
            comparisons.append({"protein":pid,"channel":ch,"threshold":t,
                "n_surface":len(asurf),"n_gly_surface":sum(r["res_name"]=="GLY" and r["key"] in asurf for r in alt),
                "surface_jaccard":sens.jacc(base_surf,asurf),
                "top5_members_jaccard":sens.jacc(a["members"],b["members"]),
                "top5_centers_jaccard":sens.jacc(a["centers"],b["centers"]),
                "top1_members_jaccard":sens.jacc(a["top1_members"],b["top1_members"]),
                "front1_centers_jaccard":sens.jacc(a["front1_centers"],b["front1_centers"]),
                "n_patches":b["n_patches"],"n_baseline_patches":a["n_patches"]})
    return comparisons
def shard(i,n):
    allrows=[];failed=[];n_ok=0
    for j,item in enumerate(sens.MANIFEST):
        if j%n!=i:continue
        try:
            rows=one(item)
            allrows+=rows;n_ok+=1
            print("THRESHOLD_OK",item["pdb_id"],len(rows),flush=True)
        except Exception as e:
            err={"protein":item["pdb_id"],"error":str(e),"trace":traceback.format_exc()[-1600:]}
            failed.append(err);print("THRESHOLD_FAIL",json.dumps(err),flush=True)
    dump(OUT/f"rows_{i}.csv",allrows)
    (OUT/f"errors_{i}.json").write_text(json.dumps(failed,indent=2))
    if failed:raise SystemExit(2)
def aggregate(path):
    import pandas as pd
    files=sorted(path.glob("rows_[0-9].csv"))
    if len(files)!=4:raise RuntimeError("Missing shards")
    d=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
    assert len(d)==48*11*len(GRID) and d.protein.nunique()==48 and d.channel.nunique()==11
    d.to_csv(path/"all_2640_map_comparisons.csv",index=False)
    cols=["surface_jaccard","top5_members_jaccard","top5_centers_jaccard","top1_members_jaccard","front1_centers_jaccard"]
    d.groupby("threshold")[cols+["n_surface","n_gly_surface"]].mean().to_csv(path/"threshold_summary.csv")
    d.groupby(["threshold","channel"])[cols].mean().to_csv(path/"per_channel.csv")
    byprotein=d.groupby(["threshold","protein"])[cols].mean().reset_index()
    rng=np.random.default_rng(20261123)
    boots=[]
    for (t),sub in byprotein.groupby("threshold"):
        for col in cols:
            a=sub[col].to_numpy()
            b=np.mean(rng.choice(a,(2500,len(a)),replace=True),axis=1)
            boots.append({"threshold":float(t),"metric":col,"mean":float(a.mean()),
                "ci95_low":float(np.quantile(b,.025)),"ci95_high":float(np.quantile(b,.975))})
    pd.DataFrame(boots).to_csv(path/"protein_cluster_bootstrap.csv",index=False)
    result={"n_proteins":48,"n_channels":11,"n_map_comparisons":len(d),
      "threshold_grid":GRID,"reference_threshold":.05,
      "metrics_mean":d.groupby("threshold")[cols].mean().reset_index().to_dict("records"),
      "interpretation":"Sensitivity relative to tRSA=0.05, not experimental threshold optimization"}
    (path/"SUMMARY.json").write_text(json.dumps(result,indent=2))
    print("THRESHOLD_FINAL",json.dumps(result),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--shard",type=int)
    p.add_argument("--shards",type=int,default=4);p.add_argument("--aggregate",type=Path)
    a=p.parse_args()
    if a.aggregate:aggregate(a.aggregate)
    elif a.shard is not None:shard(a.shard,a.shards)
    else:p.error("shard or aggregate")
