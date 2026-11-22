"""48-development-panel direct SASA convergence and all-pairs radius consensus.
Runs original published InterfaceScout kernel; does not change code or the 48 IDs.
Coordinates are biological assembly 1 and SHA-256 checked against frozen manifest.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, sys, time, traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
sys.path.insert(0,str(ROOT/"validation"))
import core, geometry, strict_apo_compat_48 as prep, sensitivity_48 as s
from model_settings import RADIUS_SELECTION_CANDIDATES_A
OUT=Path("parameter_decision_48")
OUT.mkdir(exist_ok=True)
MANIFEST=s.MANIFEST
POINTS=[100,200,500,1000]
RADII=[tuple(map(float,x)) for x in RADIUS_SELECTION_CANDIDATES_A]
assert len(RADII)==10 and (6.,9.) in RADII
CHANNELS=s.CHANNELS

def dump_rows(name,rows):
    p=OUT/name
    if rows:
        with p.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def overlap(a,b,kind):
    return s.jacc(a[kind],b[kind])

def run_one(row):
    pid=row["pdb_id"].upper()
    prep.prepare(pid)
    cif=prep.OUT/f"{pid}-assembly1.cif"
    sha=hashlib.sha256(cif.read_bytes()).hexdigest()
    if sha!=row["assembly_sha256"]:raise ValueError("SOURCE SHA256 mismatch: "+pid)
    pdb=(prep.OUT/f"{pid}_prepared.pdb").read_text(encoding="utf-8")
    nodes=geometry.extract_ca_nodes(pdb)
    profiles={}; maps={}; masks={}
    for pts in POINTS:
        core.SASA_POINTS=pts
        r=core._surface_residues(pdb,7.4)
        d={v["key"]:float(v["scrsa_raw"]) for v in r}
        if profiles and set(d)!=set(next(iter(profiles.values()))):
            raise ValueError("Residue keys changed across SASA sampling")
        profiles[pts]=d
        v,mask=s.run_cached(pdb,nodes,r,.05,(6.,9.),8.,7.4)
        maps[pts]=v; masks[pts]=mask
    sasa_rows=[]; threshold_rows=[]
    for ia,pa in enumerate(POINTS):
        for pb in POINTS[ia+1:]:
            x,y=profiles[pa],profiles[pb]
            diffs=np.array([abs(x[k]-y[k]) for k in x],dtype=float)
            maskj=s.jacc(masks[pa],masks[pb])
            crossing=sum((x[k]>=.05)!=(y[k]>=.05) for k in x)
            threshold_rows.append({"protein":pid,"points_A":pa,"points_B":pb,
                "n_residues":len(x),"surface_count_A":len(masks[pa]),"surface_count_B":len(masks[pb]),
                "surface_jaccard":maskj,"threshold_crossings":crossing,
                "mae_scrsa_raw":float(diffs.mean()),"p95_abs_scrsa_raw":float(np.quantile(diffs,.95)),
                "max_abs_scrsa_raw":float(diffs.max())})
            for ch in CHANNELS:
                a,b=maps[pa][ch],maps[pb][ch]
                sasa_rows.append({"protein":pid,"channel":ch,"points_A":pa,"points_B":pb,
                  "top5_members_jaccard":overlap(a,b,"members"),
                  "top5_centers_jaccard":overlap(a,b,"centers"),
                  "top1_members_jaccard":overlap(a,b,"top1_members"),
                  "front1_centers_jaccard":overlap(a,b,"front1_centers"),
                  "n_patches_A":a["n_patches"],"n_patches_B":b["n_patches"],
                  "nonempty_either":int(a["n_patches"]>0 or b["n_patches"]>0),
                  "threshold_crossings":crossing})
    core.SASA_POINTS=200
    reference=core._surface_residues(pdb,7.4)
    candidate_maps={}
    for aa,bb in RADII:
        v,_=s.run_cached(pdb,nodes,reference,.05,(aa,bb),8.,7.4)
        candidate_maps[f"{int(aa)}/{int(bb)}"]=v
    radius_rows=[]; candidate_rows=[]
    names=list(candidate_maps)
    for ia,pa in enumerate(names):
        for pb in names[ia+1:]:
            for ch in CHANNELS:
                a,b=candidate_maps[pa][ch],candidate_maps[pb][ch]
                radius_rows.append({"protein":pid,"channel":ch,"radius_A":pa,"radius_B":pb,
                  "top5_members_jaccard":overlap(a,b,"members"),
                  "top5_centers_jaccard":overlap(a,b,"centers"),
                  "top1_members_jaccard":overlap(a,b,"top1_members"),
                  "front1_centers_jaccard":overlap(a,b,"front1_centers"),
                  "n_patches_A":a["n_patches"],"n_patches_B":b["n_patches"],
                  "nonempty_either":int(a["n_patches"]>0 or b["n_patches"]>0)})
    for pa in names:
        for ch in CHANNELS:
            a=candidate_maps[pa][ch]
            others=[candidate_maps[pb][ch] for pb in names if pb!=pa]
            candidate_rows.append({"protein":pid,"channel":ch,"candidate_radius":pa,
               "consensus_top5_members":float(np.mean([overlap(a,b,"members") for b in others])),
               "consensus_top5_centers":float(np.mean([overlap(a,b,"centers") for b in others])),
               "consensus_top1_members":float(np.mean([overlap(a,b,"top1_members") for b in others])),
               "mean_patches":a["n_patches"],
               "pairs_involving_nonempty":sum(a["n_patches"]>0 or b["n_patches"]>0 for b in others)})
    return sasa_rows,threshold_rows,radius_rows,candidate_rows,{"protein":pid,"source_sha256":sha,"status":"PASS","n_residues":len(profiles[200])}

def run_shard(shard,shards):
    out=[[],[],[],[],[]];fail=[]
    for i,row in enumerate(MANIFEST):
        if i%shards!=shard:continue
        pid=row["pdb_id"]
        started=time.time()
        try:
            data=run_one(row)
            for group,part in zip(out,data):
                if isinstance(part,list):group.extend(part)
                else:group.append(part)
            print("PROTEIN_OK "+json.dumps({"pdb":pid,"seconds":round(time.time()-started,1)}),flush=True)
        except Exception as exc:
            fail.append({"protein":pid,"error":str(exc),"trace":traceback.format_exc()[-2600:]})
            print("PROTEIN_FAIL "+json.dumps(fail[-1]),flush=True)
    for filename,rows in zip(["sasa_direct","sasa_residue_errors","radius_pairwise","radius_consensus","provenance"],out):
        dump_rows(f"{filename}_{shard}.csv",rows)
    (OUT/f"failures_{shard}.json").write_text(json.dumps(fail,indent=2))
    print("SHARD_SUMMARY "+json.dumps({"shard":shard,"success":len(out[-1]),"fail":len(fail),"sasa_pairs":len(out[0]),"radius_pairs":len(out[2])}),flush=True)
    if fail:raise SystemExit(2)

def agg(root):
    import pandas as pd
    def read(prefix):
        paths=sorted(root.glob(prefix+"_*.csv"))
        if len(paths)!=4:raise RuntimeError(f"missing shards for {prefix}: {paths}")
        return pd.concat([pd.read_csv(p) for p in paths],ignore_index=True)
    sas=read("sasa_direct");err=read("sasa_residue_errors");rpair=read("radius_pairwise");rc=read("radius_consensus");prov=read("provenance")
    if len(prov)!=48 or prov.protein.nunique()!=48:raise RuntimeError(f"Incomplete sample {len(prov)}")
    if len(sas)!=48*11*6 or len(rpair)!=48*11*45 or len(rc)!=48*11*10:raise RuntimeError("Incomplete pairwise grid")
    def interval_by_protein(df,metric,keys):
        protein=df.groupby(keys+["protein"],as_index=False)[metric].mean()
        rows=[];rng=np.random.default_rng(20261122)
        for v,g in protein.groupby(keys,sort=True):
            arr=g[metric].to_numpy()
            boot=np.mean(rng.choice(arr,(4000,len(arr)),replace=True),axis=1)
            rec=dict(zip(keys, v if isinstance(v,tuple) else (v,)))
            rec.update({"n_proteins":len(arr),"mean":float(arr.mean()),"median":float(np.median(arr)),
                "ci95_low":float(np.quantile(boot,.025)),"ci95_high":float(np.quantile(boot,.975))})
            rows.append(rec)
        return pd.DataFrame(rows)
    sasa_summary=interval_by_protein(sas,"top5_members_jaccard",["points_A","points_B"])
    sasa_centers=interval_by_protein(sas,"top5_centers_jaccard",["points_A","points_B"])
    sasa_summary=sasa_summary.merge(sasa_centers[["points_A","points_B","mean","ci95_low","ci95_high"]],on=["points_A","points_B"],suffixes=("","_centers"))
    scalar=err.groupby(["points_A","points_B"]).agg(mae_scrsa=("mae_scrsa_raw","mean"),
      avg_surface_jaccard=("surface_jaccard","mean"),median_crossing=("threshold_crossings","median"),
      sum_crossings=("threshold_crossings","sum")).reset_index()
    sasa_summary=sasa_summary.merge(scalar,on=["points_A","points_B"])
    radius_summary=interval_by_protein(rc,"consensus_top5_members",["candidate_radius"])
    radius_centers=interval_by_protein(rc,"consensus_top5_centers",["candidate_radius"])
    radius_summary=radius_summary.merge(radius_centers[["candidate_radius","mean","ci95_low","ci95_high"]],on="candidate_radius",suffixes=("","_centers"))
    radius_summary=radius_summary.sort_values("mean",ascending=False)
    pair_summary=interval_by_protein(rpair,"top5_members_jaccard",["radius_A","radius_B"])
    # Data for conditional metric ignoring comparisons where both maps are empty.
    relevant=rpair[rpair.nonempty_either==1]
    nonempty_summary=interval_by_protein(relevant,"top5_members_jaccard",["radius_A","radius_B"])
    radius_summary_nonempty=interval_by_protein(rc[rc.pairs_involving_nonempty>0],"consensus_top5_members",["candidate_radius"])
    for label,df in [("sasa_direct_all",sas),("sasa_residue_errors_all",err),
                     ("radius_pairwise_all",rpair),("radius_consensus_all",rc),
                     ("sasa_pair_summary",sasa_summary),("radius_consensus_summary",radius_summary),
                     ("radius_pair_summary",pair_summary),("radius_pair_nonempty_summary",nonempty_summary),
                     ("radius_consensus_nonempty_summary",radius_summary_nonempty),
                     ("source_provenance",prov)]:
        df.to_csv(root/f"{label}.csv",index=False)
    top=radius_summary.iloc[0]
    current=radius_summary.query("candidate_radius == '6/9'").iloc[0]
    result={"n_proteins":48,"n_channels":11,"n_sasa_pair_map_comparisons":len(sas),
       "n_radius_pair_map_comparisons":len(rpair),"n_radius_candidates":len(RADII),
       "radius_best_consensus":str(top["candidate_radius"]),"radius_best_mean":float(top["mean"]),
       "radius_6_9_rank":int(radius_summary.reset_index(drop=True).index[radius_summary.reset_index(drop=True).candidate_radius.eq("6/9")][0])+1,
       "radius_6_9_consensus":float(current["mean"])}
    (root/"SUMMARY.json").write_text(json.dumps(result,indent=2))
    print("SUMMARY "+json.dumps(result),flush=True)
    print("SASA_SUMMARY\n"+sasa_summary.round(4).to_string(index=False),flush=True)
    print("RADIUS_SUMMARY\n"+radius_summary.round(4).to_string(index=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--shard",type=int);p.add_argument("--shards",type=int,default=4);p.add_argument("--aggregate",type=Path)
    a=p.parse_args()
    if a.aggregate:agg(a.aggregate)
    elif a.shard is not None:run_shard(a.shard,a.shards)
    else:p.error("select shard or aggregate")
