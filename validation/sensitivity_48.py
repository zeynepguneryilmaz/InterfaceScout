"""InterfaceScout 48-structure parameter sensitivity. Recomputes using the frozen public kernel.
No adsorption labels; source biological-assembly SHA256 is checked against the frozen manifest.
Usage: python validation/sensitivity_48.py --shard 0 --shards 4
       python validation/sensitivity_48.py --aggregate output/sensitivity48
"""
from __future__ import annotations
import argparse, copy, csv, hashlib, json, math, os, sys, time, traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
sys.path.insert(0,str(ROOT/"validation"))
import core, coarse_patch, geometry, strict_apo_compat_48 as prep
from model_settings import RADIUS_SELECTION_CANDIDATES_A

MANIFEST=list(csv.DictReader(open(ROOT/"validation"/"strict_apo_48_selected_candidates_v2_final.csv",encoding="utf-8")))
assert len(MANIFEST)==48 and len(set(x["pdb_id"] for x in MANIFEST))==48
BASE_THRESHOLD=0.05
BASE_RADII=(6.0,9.0)
BASE_PATCH=8.0
BASE_PH=7.4
BASE_POINTS=200
CHANNELS=list(core.CHEMISTRIES)
assert len(CHANNELS)==11
FIELDS=["protein","channel","family","setting","n_surface","n_seeds","n_patches","n_front1","top5_center_jaccard","top5_members_jaccard","top1_members_jaccard","front1_center_jaccard","surface_jaccard","sasa_mae_scrsa","elapsed_s","error"]
OUT=Path("sensitivity48")
OUT.mkdir(exist_ok=True)

def jacc(a,b):
    a=set(a); b=set(b)
    return len(a&b)/len(a|b) if (a|b) else 1.0

def summarize_map(patches, seeds):
    top=patches[:5]
    return {"n_seeds":len(seeds),"n_patches":len(patches),
        "n_front1":sum(int(p.get("pareto_front",0))==1 for p in patches),
        "centers":set(str(p["center_key"]) for p in top),
        "members":set(m for p in top for m in p.get("members",[])),
        "top1_members":set(patches[0].get("members",[])) if patches else set(),
        "front1_centers":set(str(p["center_key"]) for p in patches if int(p.get("pareto_front",0))==1),
        "patches":patches}

def run_cached(pdb_text,ca_nodes,rows,threshold,radii,patch_radius,pH=BASE_PH):
    cached=[]
    for r in rows:
        r=r.copy()
        r["surface_exposed"]=bool(float(r["scrsa_raw"])>=threshold)
        r["charge_fraction"]=round(core.charged_fraction(r["res_name"],pH),5)
        r["charge_descriptor"]=round(core.residue_charge_descriptor(r["res_name"],pH),5)
        r["ionization_sensitive"]=bool(r["res_name"] in core.PKA and abs(pH-core.PKA[r["res_name"]])<=1.)
        cached.append(r)
    prior=core._surface_residues
    prior_radii=core.PATCH_RADII_A
    prior_patch=coarse_patch.PATCH_SCALE_A
    try:
        core._surface_residues=lambda text,pH: [r.copy() for r in cached]
        core.PATCH_RADII_A=tuple(radii)
        coarse_patch.PATCH_SCALE_A=float(patch_radius)
        request=core.AnalyzeRequest(pdb_text=pdb_text,env=core.EnvParams(pH=pH,ionic=150.,temp=298.))
        result=core.analyze(request)
        data={}
        for key in CHANNELS:
            patches=coarse_patch.build_coarse_patches(v1_result=result,chemistry=key,ca_nodes=ca_nodes)
            seeds=[r for r in result["chemistries"][key]["residues"] if float(r.get("local_score",0))>0]
            data[key]=summarize_map(patches,seeds)
        surf=set(str(r["key"]) for r in result["surface_residues"])
        return data,surf
    finally:
        core._surface_residues=prior
        core.PATCH_RADII_A=prior_radii
        coarse_patch.PATCH_SCALE_A=prior_patch

def calc_rows(pid,family,setting,data,surface,baseline,bsurface,elapsed,mae=None):
    result=[]
    for key in CHANNELS:
        v=data[key]; b=baseline[key]
        result.append({"protein":pid,"channel":key,"family":family,"setting":setting,
            "n_surface":len(surface),"n_seeds":v["n_seeds"],"n_patches":v["n_patches"],
            "n_front1":v["n_front1"],
            "top5_center_jaccard":jacc(v["centers"],b["centers"]),
            "top5_members_jaccard":jacc(v["members"],b["members"]),
            "top1_members_jaccard":jacc(v["top1_members"],b["top1_members"]),
            "front1_center_jaccard":jacc(v["front1_centers"],b["front1_centers"]),
            "surface_jaccard":jacc(surface,bsurface),
            "sasa_mae_scrsa":mae if mae is not None else "",
            "elapsed_s":round(elapsed/len(CHANNELS),4),"error":""})
    return result

def pareto_ablation(pid,baseline):
    obj=["chemistry_support","mean_accessibility","patch_coherence","orientation_coherence"]
    rows=[]
    for channel,b in baseline.items():
        patches=b["patches"]
        for omitted in obj:
            keep=[v for v in obj if v!=omitted]
            front=set()
            for i,p in enumerate(patches):
                pv=[float(p[k]) for k in keep]
                if not any(
                   all(float(q[k])>=pv[j] for j,k in enumerate(keep)) and
                   any(float(q[k])>pv[j] for j,k in enumerate(keep))
                   for z,q in enumerate(patches) if z!=i
                ):
                    front.add(str(p["center_key"]))
            rows.append({"protein":pid,"channel":channel,"omitted":omitted,
                "baseline_front1":b["n_front1"],"ablation_front1":len(front),
                "front1_jaccard":jacc(front,b["front1_centers"])})
    return rows

def run_shard(shard,shards):
    csvpath=OUT/f"shard_{shard}.csv"
    errors=[]; ablations=[]; timings=[]; per_protein=[]
    indices=[i for i in range(len(MANIFEST)) if i%shards==shard]
    with csvpath.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS);writer.writeheader()
        for i in indices:
            row=MANIFEST[i]
            pid=row["pdb_id"].upper()
            started=time.time()
            try:
                p=prep.prepare(pid)
                cif=prep.OUT/f"{pid}-assembly1.cif"
                actual=hashlib.sha256(cif.read_bytes()).hexdigest()
                expected=row["assembly_sha256"]
                if expected!=actual:raise ValueError(f"ASSEMBLY_SHA_MISMATCH {expected} {actual}")
                pdb=(prep.OUT/f"{pid}_prepared.pdb").read_text(encoding="utf-8")
                ca_nodes=geometry.extract_ca_nodes(pdb)
                core.SASA_POINTS=BASE_POINTS
                core.SC_RSA_THRESHOLD=BASE_THRESHOLD
                baseline_rows=core._surface_residues(pdb,BASE_PH)
                baseline,bsurf=run_cached(pdb,ca_nodes,baseline_rows,BASE_THRESHOLD,BASE_RADII,BASE_PATCH)
                writer.writerows(calc_rows(pid,"baseline","reference",baseline,bsurf,baseline,bsurf,0))
                ablations.extend(pareto_ablation(pid,baseline))
                cases=[]
                cases += [("scrsa_threshold",str(x),x,BASE_RADII,BASE_PATCH,BASE_PH,None) for x in [0.01,0.03,0.075,0.10]]
                cases += [("radii",f"{a:g}_{b:g}",BASE_THRESHOLD,(a,b),BASE_PATCH,BASE_PH,None) for a,b in RADIUS_SELECTION_CANDIDATES_A if (a,b)!=BASE_RADII]
                cases += [("patch_radius",str(x),BASE_THRESHOLD,BASE_RADII,float(x),BASE_PH,None) for x in [6,10,12]]
                cases += [("pH",str(x),BASE_THRESHOLD,BASE_RADII,BASE_PATCH,float(x),None) for x in [5.5,9.0]]
                for family,setting,thr,radii,pradius,ph,points in cases:
                    t0=time.time()
                    v,s=run_cached(pdb,ca_nodes,baseline_rows,thr,radii,pradius,ph)
                    writer.writerows(calc_rows(pid,family,setting,v,s,baseline,bsurf,time.time()-t0))
                # Independent Shrake-Rupley reruns; these are actual per-atom SASA resamplings.
                for points in [100,500,1000]:
                    t0=time.time()
                    core.SASA_POINTS=int(points)
                    raw=core._surface_residues(pdb,BASE_PH)
                    mae=float(np.mean([abs(float(a["scrsa_raw"])-float(b["scrsa_raw"])) for a,b in zip(raw,baseline_rows)]))
                    v,s=run_cached(pdb,ca_nodes,raw,BASE_THRESHOLD,BASE_RADII,BASE_PATCH)
                    writer.writerows(calc_rows(pid,"sasa_points",str(points),v,s,baseline,bsurf,time.time()-t0,mae))
                core.SASA_POINTS=BASE_POINTS
                f.flush()
                t=time.time()-started
                timings.append(t)
                per_protein.append({"pdb":pid,"seconds":round(t,2),"structure_sha256":actual,"prepared_sha256":p["pdb_sha256"],"status":"PASS"})
                print("PROTEIN_OK "+json.dumps({"pdb":pid,"seconds":round(t,1),"cases":1+len(cases)+3}),flush=True)
            except Exception as e:
                err={"pdb":pid,"error":str(e),"trace":traceback.format_exc()[-4000:]}
                errors.append(err)
                per_protein.append({"pdb":pid,"status":"FAIL","error":str(e)})
                print("PROTEIN_FAIL "+json.dumps(err),flush=True)
    with (OUT/f"ablation_{shard}.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["protein","channel","omitted","baseline_front1","ablation_front1","front1_jaccard"]);w.writeheader();w.writerows(ablations)
    (OUT/f"provenance_{shard}.json").write_text(json.dumps(per_protein,indent=2))
    (OUT/f"errors_{shard}.json").write_text(json.dumps(errors,indent=2))
    print("SHARD_SUMMARY "+json.dumps({"shard":shard,"n":len(indices),"passes":len(timings),"failures":len(errors),"runtime_sum_s":round(sum(timings),1)}),flush=True)
    if errors:sys.exit(2)

def aggregate(directory):
    import pandas as pd
    paths=sorted(directory.rglob("shard_*.csv"))
    if len(paths)!=4:raise RuntimeError(f"Expected four shard results, got {paths}")
    df=pd.concat([pd.read_csv(p) for p in paths],ignore_index=True)
    if df.protein.nunique()!=48 or df.channel.nunique()!=11:raise RuntimeError("Missing protein or channel")
    df.to_csv(directory/"sensitivity_48_all.csv",index=False)
    a=pd.concat([pd.read_csv(p) for p in directory.rglob("ablation_*.csv")],ignore_index=True)
    a.to_csv(directory/"pareto_ablation_48.csv",index=False)
    grp=df.query("family != 'baseline'").groupby(["family","setting"],sort=True,dropna=False)
    tab=grp.agg(n_entries=("protein","count"),n_proteins=("protein","nunique"),
        median_top5_members_jaccard=("top5_members_jaccard","median"),
        mean_top5_members_jaccard=("top5_members_jaccard","mean"),
        mean_top5_centers_jaccard=("top5_center_jaccard","mean"),
        mean_top1_members_jaccard=("top1_members_jaccard","mean"),
        mean_surface_jaccard=("surface_jaccard","mean"),
        median_n_patches=("n_patches","median"),
        mean_sasa_mae=("sasa_mae_scrsa","mean")).reset_index()
    tab.to_csv(directory/"sensitivity_summary.csv",index=False)
    channel=df.query("family != 'baseline'").groupby(["channel","family","setting"]).agg(
        n_proteins=("protein","nunique"),mean_top5_members_jaccard=("top5_members_jaccard","mean"),
        mean_top5_centers_jaccard=("top5_center_jaccard","mean"),
        mean_n_patches=("n_patches","mean")).reset_index()
    channel.to_csv(directory/"channel_sensitivity.csv",index=False)
    ab=a.groupby("omitted").agg(n_proteins=("protein","nunique"),
         mean_front1_jaccard=("front1_jaccard","mean"),
         median_front1_jaccard=("front1_jaccard","median"),
         mean_new_front1=("ablation_front1","mean")).reset_index()
    ab.to_csv(directory/"ablation_summary.csv",index=False)
    # Protein-level bootstrap: channels are repeated measures, not independent biological replicates.
    rng=np.random.default_rng(20261122); boot=[]
    p=df.query("family != 'baseline'").groupby(["family","setting","protein"],as_index=False)["top5_members_jaccard"].mean()
    for (family,setting),v in p.groupby(["family","setting"],sort=True):
        values=v["top5_members_jaccard"].to_numpy()
        means=[np.mean(rng.choice(values,len(values),replace=True)) for j in range(1000)]
        boot.append({"family":family,"setting":setting,"protein_n":len(values),
                     "mean_protein_avg_top5_jaccard":round(float(values.mean()),5),
                     "bootstrap95_low":round(float(np.quantile(means,.025)),5),
                     "bootstrap95_high":round(float(np.quantile(means,.975)),5)})
    pd.DataFrame(boot).to_csv(directory/"protein_level_bootstrap.csv",index=False)
    summary={"proteins":int(df.protein.nunique()),"channels":int(df.channel.nunique()),
       "conditions":int(df[["family","setting"]].drop_duplicates().shape[0]),
       "map_instance_count":len(df),"ablation_count":len(a),
       "parameter_sensitivity":tab.to_dict("records"),
       "pareto_ablation":ab.to_dict("records"),
       "limitations":["No holdout validation","No measured adsorption accuracy",
            "No independent local-normal/geodesic metric or perturbation of temperature/ionic strength"]}
    (directory/"FINAL_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("FINAL_SUMMARY "+json.dumps({"proteins":summary["proteins"],"channels":summary["channels"],
        "conditions":summary["conditions"],"map_instances":len(df),"ablation_records":len(a),
        "families":tab.groupby("family")["mean_top5_members_jaccard"].agg(["min","max"]).round(4).to_dict("index"),
        "ablation":ab.to_dict("records")},separators=(",",":")),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--shard",type=int)
    ap.add_argument("--shards",type=int,default=4)
    ap.add_argument("--aggregate",type=Path)
    args=ap.parse_args()
    if args.aggregate:aggregate(args.aggregate)
    elif args.shard is not None:run_shard(args.shard,args.shards)
    else:ap.error("Specify --shard or --aggregate")
