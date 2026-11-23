"""Production v1.1.0: criterion dependence, rank stability, independent baselines,
and prespecified descriptive examples. No adsorption labels enter ranking methods.

Usage: --panel-shard I --shards 4 | --cases | --aggregate DIR
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys,traceback
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
from interface_engine import analyze_all_maps
from research_workflows import combine_maps
from geometry import extract_ca_nodes
from model_settings import MODEL_VERSION
import sensitivity_48 as sens
import strict_apo_compat_48 as prepare48
import dual_literature_benchmark as literature
OUT=Path("extended_v110");OUT.mkdir(exist_ok=True)
CRITERIA=("chemistry_support","mean_accessibility","patch_coherence","orientation_coherence")
PAIRS=[(CRITERIA[i],CRITERIA[j]) for i in range(4) for j in range(i+1,4)]
CASE_PAIRS=[("hydrophobic","pi_carbon"),("cationic","hbond_donor"),("anionic","hbond_acceptor")]
CASE_IDS={"ubiquitin_citrate_au","hcaii_silica_ph63"}
def csvwrite(path,records):
    if records:
        with path.open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
def jacc(a,b):
    a,b=set(a),set(b)
    return len(a&b)/len(a|b) if a|b else 1.
def fronts(values):
    """Pareto layers computed independently of InterfaceScout production sorter."""
    n=len(values)
    if not n:return []
    x=np.asarray(values,dtype=float)
    remaining=np.ones(n,dtype=bool);level=np.zeros(n,dtype=int);rank=1
    while remaining.any():
        active=np.flatnonzero(remaining)
        v=x[active]
        dominates=(v[:,None,:]>=v[None,:,:]).all(axis=2)&(v[:,None,:]>v[None,:,:]).any(axis=2)
        top=active[~dominates.any(axis=0)]
        if not len(top):raise RuntimeError("Invalid Pareto layer")
        level[top]=rank;remaining[top]=False;rank+=1
    return level.tolist()
def rerank(patches,values):
    levels=fronts(values)
    order=sorted(range(len(patches)),key=lambda i:(levels[i],-values[i][2],-values[i][0],str(patches[i]["center_key"])))
    return order,levels
def panel_case(item):
    pid=item["pdb_id"].upper()
    prepare48.prepare(pid)
    cif=prepare48.OUT/f"{pid}-assembly1.cif"
    actual=hashlib.sha256(cif.read_bytes()).hexdigest()
    if actual!=item["assembly_sha256"]:raise ValueError("Frozen source hash mismatch")
    pdb=(prepare48.OUT/f"{pid}_prepared.pdb").read_text()
    a=analyze_all_maps(pdb_text=pdb,pH=7.4,chain=None)
    if a["version"]!=MODEL_VERSION or a["n_maps"]!=11:raise ValueError("Not final production")
    corr=[];stability=[];detail=[]
    for ch in a["map_order"]:
        patches=a["maps"][ch]["patches"]
        if not patches:
            detail.append({"protein":pid,"channel":ch,"n_patches":0,
                "front1_count":0,"top5_centers":""})
            continue
        vals=np.array([[float(p[k]) for k in CRITERIA] for p in patches],dtype=float)
        baseline,levels=rerank(patches,vals)
        if baseline != list(range(len(patches))):
            raise ValueError("Production rank mismatch "+pid+" "+ch)
        for x,y in PAIRS:
            v1=vals[:,CRITERIA.index(x)];v2=vals[:,CRITERIA.index(y)]
            constant=(len(patches)<3 or np.ptp(v1)<1e-12 or np.ptp(v2)<1e-12)
            rho=float(spearmanr(v1,v2).statistic) if not constant else None
            corr.append({"protein":pid,"channel":ch,"n_patches":len(patches),
                "criterion_a":x,"criterion_b":y,"spearman_rho":rho if rho is not None and np.isfinite(rho) else "",
                "evaluable":int(rho is not None and np.isfinite(rho))})
        basefront={patches[i]["center_key"] for i in range(len(patches)) if levels[i]==1}
        basetop={patches[i]["center_key"] for i in baseline[:5]}
        # Multiplicative bounded deterministic perturbations probe numerical ranking
        # sensitivity. They are *not* a measurement-error confidence model.
        for pct in (0.01,0.05):
            for trial in range(20):
                rng=np.random.default_rng(20261123+int(pid.encode().hex(),16)%1000000+sum(ord(c) for c in ch)*29+int(pct*100)*1000+trial)
                pert=vals*np.maximum(0.,1.+rng.normal(0,pct,size=vals.shape))
                pert=np.clip(pert,0,1)
                order,newlevels=rerank(patches,pert)
                nf={patches[i]["center_key"] for i in range(len(patches)) if newlevels[i]==1}
                nt={patches[i]["center_key"] for i in order[:5]}
                stability.append({"protein":pid,"channel":ch,"perturbation_sd_fraction":pct,"trial":trial,
                    "n_patches":len(patches),"front1_jaccard":jacc(basefront,nf),
                    "top5_center_jaccard":jacc(basetop,nt),
                    "top1_center_retained":int(bool(baseline and order and baseline[0]==order[0]))})
        detail.append({"protein":pid,"channel":ch,"n_patches":len(patches),
            "front1_count":len(basefront),"top5_centers":";".join(patches[i]["center_key"] for i in baseline[:5])})
    return corr,stability,detail,{"protein":pid,"source_sha256":actual,"model":MODEL_VERSION,"status":"PASS"}
def do_panel(i,n):
    corr=[];stab=[];details=[];provs=[];errors=[]
    for j,item in enumerate(sens.MANIFEST):
        if j%n!=i:continue
        try:
            x,y,z,p=panel_case(item)
            corr+=x;stab+=y;details+=z;provs.append(p)
            print("EXT_PANEL_PASS",item["pdb_id"],len(z),flush=True)
        except Exception as e:
            er={"protein":item["pdb_id"],"error":str(e),"trace":traceback.format_exc()[-1200:]}
            errors.append(er);print("EXT_PANEL_FAIL",json.dumps(er),flush=True)
    csvwrite(OUT/f"correlations_{i}.csv",corr)
    csvwrite(OUT/f"stability_{i}.csv",stab)
    csvwrite(OUT/f"patches_{i}.csv",details)
    (OUT/f"provenance_{i}.json").write_text(json.dumps(provs,indent=2))
    (OUT/f"errors_{i}.json").write_text(json.dumps(errors,indent=2))
    if errors:raise SystemExit(2)
def exposed_geometry(pdb,a):
    nodes={n["key"]:np.asarray(n["coord"],dtype=float) for n in extract_ca_nodes(pdb)}
    rows={r["key"]:r for r in a["residue_profiles"] if r["surface_exposed"] and r["key"] in nodes}
    keys=sorted(rows)
    centroid=np.mean(np.array(list(nodes.values())),axis=0)
    pts=np.array([nodes[k] for k in keys])
    norms=pts-centroid
    norms=norms/np.maximum(np.linalg.norm(norms,axis=1)[:,None],1e-12)
    distance=np.linalg.norm(pts[:,None,:]-pts[None,:,:],axis=-1)
    face=(norms@norms.T)>0
    return keys,rows,pts,distance,face
def independent_baselines(pdb,a,k=5):
    """Build patches *without* reading chemical maps or production patches.

    SASA peaks: descending tRSA local maxima, 8A same-face NMS.
    Geometry coverage: deterministic farthest-point selection over exposed CA
    with a geometry-only starting point (farthest from CA centroid).
    Both patch sets use the same 8A radius + hemisphere rule solely to make
    coarse regions/coverage interpretable. They are NOT 11 chemistry maps.
    """
    keys,rows,pts,distance,face=exposed_geometry(pdb,a)
    if not keys:return {"sasa_only":[],"geometry_fps":[]},{}
    dloc=distance<=8
    nms=dloc & face
    scores=np.asarray([float(rows[key]["trsa_raw"]) for key in keys])
    ranked=sorted(range(len(keys)),key=lambda i:(-scores[i],keys[i]))
    sasa=[]
    for i in ranked:
        if not any(nms[i,j] for j in sasa):sasa.append(i)
        if len(sasa)==k:break
    if len(sasa)<k:
        for i in ranked:
            if i not in sasa:sasa.append(i)
            if len(sasa)==k:break
    centroid=pts.mean(axis=0)
    first=max(range(len(keys)),key=lambda i:(float(np.linalg.norm(pts[i]-centroid)),keys[i]))
    geom=[first]
    while len(geom)<min(k,len(keys)):
        choose=max((i for i in range(len(keys)) if i not in geom),
                   key=lambda i:(float(np.min(distance[i,geom])),keys[i]))
        geom.append(choose)
    def make(indices):
        return [{"center_key":keys[i],"members":[keys[j] for j in range(len(keys)) if nms[i,j]]}
                for i in indices]
    baselines={"sasa_only":make(sasa),"geometry_fps":make(geom)}
    return baselines,{"n_surface":len(keys),"n_sasa_centers":len(sasa),
      "n_geometry_centers":len(geom)}
def do_cases():
    comparisons=[];case_details=[];provs=[];errors=[];case_examples=[];examples_details=[]
    # These two case studies and three map pairs are predeclared in source.
    for case in literature.MANIFEST["cases"]:
        try:
            pdb,prov=literature.prepare(case)
            a=analyze_all_maps(pdb_text=pdb,pH=float(case["pH"]),chain=None)
            if a["version"]!=MODEL_VERSION or a["n_maps"]!=11:raise ValueError("Not production")
            labels=literature.labels_for_case(case,a["residue_profiles"])
            coords={n["key"]:np.asarray(n["coord"],float) for n in extract_ca_nodes(pdb)}
            if labels-set(coords):raise ValueError("Missing evidence coordinates")
            from final_v110_literature import measure
            baselines,metadata=independent_baselines(pdb,a)
            for meth,patches in baselines.items():
                score=measure(patches,labels,coords)
                comparisons.append({"case":case["id"],"protein":case["pdb_id"],
                    "analysis_set":case["analysis_set"],"method":meth,"n_labels":len(labels),**score})
            # Independent baselines are evaluated ONCE per condition.
            # Production map scores are retained separately, not replicated.
            for ch in a["map_order"]:
                patches=a["maps"][ch]["patches"]
                score=measure(patches,labels,coords)
                case_details.append({"case":case["id"],"protein":case["pdb_id"],"channel":ch,
                    "method":"interfacescout_pareto","n_labels":len(labels),**score})
            if case["id"] in CASE_IDS:
                summary={"case":case["id"],"pdb_id":case["pdb_id"],"pH":float(case["pH"]),
                    "n_surface":metadata["n_surface"],"n_maps":a["n_maps"]}
                for x,y in CASE_PAIRS:
                    c=combine_maps(a,[x,y],max_regions=None)
                    example={"case":case["id"],"pdb_id":case["pdb_id"],"pH":float(case["pH"]),
                        "map_a":x,"map_b":y,"n_exposed":c["n_exposed_residues"],
                        "exact_overlap_residues":c["n_exact_overlap"],
                        "exact_overlap_keys":";".join(c["exact_same_residue_overlap"]),
                        "colocated_regions":sum(q["all_maps_colocated"] for q in c["regions"]),
                        "distinct_residue_only_regions":sum(q["distinct_residue_colocation_only"] for q in c["regions"]),
                        "n_total_regions":c["regions_total"]}
                    case_examples.append(example)
                    for index,region in enumerate(c["regions"][:3],1):
                        p=region["regional_profile"]
                        examples_details.append({"case":case["id"],"pair":x+"|"+y,
                            "rank":index,"source_center":region["source_center_key"],
                            "source_map":region["source_map"],"member_count":region["n_members"],
                            "exact_overlap_count":len(region["exact_all_map_seed_overlap"]),
                            "n_support_a":p["supporting_residue_count_by_map"][x],
                            "n_support_b":p["supporting_residue_count_by_map"][y],
                            "mean_support_a":p["mean_local_support_by_map"][x],
                            "mean_support_b":p["mean_local_support_by_map"][y],
                            "min_distinct_CA_A":p["min_distinct_intermap_CA_distance_A"],
                            "members":";".join(region["member_keys"])})
                summary["map_patch_count_by_channel"]={k:len(a["maps"][k]["patches"]) for k in a["map_order"]}
                (OUT/f"case_{case['id']}_summary.json").write_text(json.dumps(summary,indent=2))
            prov.update({"case":case["id"],"model":MODEL_VERSION,"status":"PASS"})
            provs.append(prov);print("EXT_CASE_PASS",case["id"],flush=True)
        except Exception as e:
            er={"case":case["id"],"error":str(e),"trace":traceback.format_exc()[-1200:]}
            errors.append(er);print("EXT_CASE_FAIL",json.dumps(er),flush=True)
    csvwrite(OUT/"independent_baselines.csv",comparisons)
    csvwrite(OUT/"production_per_channel.csv",case_details)
    csvwrite(OUT/"case_study_map_combinations.csv",case_examples)
    csvwrite(OUT/"case_study_region_details.csv",examples_details)
    (OUT/"case_provenance.json").write_text(json.dumps(provs,indent=2))
    (OUT/"case_errors.json").write_text(json.dumps(errors,indent=2))
    if errors:raise SystemExit(2)
def aggregate(path):
    import pandas as pd
    corr=pd.concat([pd.read_csv(x) for x in sorted(path.glob("correlations_[0-9].csv"))],ignore_index=True)
    stab=pd.concat([pd.read_csv(x) for x in sorted(path.glob("stability_[0-9].csv"))],ignore_index=True)
    patches=pd.concat([pd.read_csv(x) for x in sorted(path.glob("patches_[0-9].csv"))],ignore_index=True)
    ctrl=pd.read_csv(path/"independent_baselines.csv")
    prod=pd.read_csv(path/"production_per_channel.csv")
    case=pd.read_csv(path/"case_study_map_combinations.csv")
    if patches.protein.nunique()!=48 or len(patches)!=528 or len(corr)!=6*int((patches.n_patches>0).sum()) or len(stab)!=40*int((patches.n_patches>0).sum()) or len(ctrl)!=18 or len(prod)!=99 or len(case)!=6:
        raise ValueError("Incomplete extended production validation")
    provenance=[]
    for p in sorted(path.glob("provenance_[0-9].json")):provenance.extend(json.loads(p.read_text()))
    if len(provenance)!=48:raise ValueError("Incomplete provenance")
    corr.to_csv(path/"all_criterion_correlations.csv",index=False)
    stab.to_csv(path/"all_rank_stability_trials.csv",index=False)
    patches.to_csv(path/"all_patch_counts.csv",index=False)
    summary={"status":"PASS","model":MODEL_VERSION,"n_proteins":48,
        "n_maps":len(patches),"criterion_pairs":len(corr),
        "criterion_pairs_evaluable":int(corr.evaluable.sum()),
        "rank_perturbation_trials":len(stab),
        "mean_spearman_by_pair":[],
        "rank_stability_by_noise":[],
        "baseline_primary":{},
        "case_studies":case.to_dict("records"),
        "limitations":"Descriptive structural robustness. The 9 historical experimental conditions are not an independent holdout. Baselines share geometric exposure/patch-size choices, but independently generate surface centers; no material-specific interaction is modeled."}
    for (x,y),part in corr.groupby(["criterion_a","criterion_b"]):
        valid=part.spearman_rho.dropna()
        protein=part.groupby("protein").spearman_rho.mean().dropna()
        summary["mean_spearman_by_pair"].append({"a":x,"b":y,
          "mean_per_protein_rho":None if len(protein)==0 else float(protein.mean()),
          "median_per_map_rho":None if len(valid)==0 else float(valid.median()),
          "n_evaluable_maps":len(valid)})
    for n,part in stab.groupby("perturbation_sd_fraction"):
        summary["rank_stability_by_noise"].append({"noise_fraction":float(n),
          "mean_front1_jaccard":float(part.front1_jaccard.mean()),
          "median_front1_jaccard":float(part.front1_jaccard.median()),
          "mean_top5_jaccard":float(part.top5_center_jaccard.mean()),
          "mean_top1_retention":float(part.top1_center_retained.mean())})
    metrics=("top5_hit_8A","top5_recall_8A","top5_direct_hit","top5_direct_recall")
    # Report a condition-level average of the 11 original maps, but do not use
    # 11 channels as independent biological repetitions.
    # Attach frozen evidence-set labels from the per-condition controls; the\n    # per-channel CSV intentionally stores case/protein/channel/method only.\n    prod_labeled=prod.merge(ctrl[["case","analysis_set"]].drop_duplicates(),on="case",validate="many_to_one")\n    original=prod_labeled.groupby(["case","analysis_set","protein"],as_index=False)[list(metrics)].mean()
    methods=pd.concat([original.assign(method="interfacescout_11map_average"),ctrl],ignore_index=True)
    methods.to_csv(path/"per_condition_independent_baseline_comparison.csv",index=False)
    rng=np.random.default_rng(20261123)
    bootstrap=[]
    for method,part in methods[methods.analysis_set=="primary"].groupby("method"):
        summary["baseline_primary"][method]={k:float(part[k].mean()) for k in metrics}
    for baseline in ("sasa_only","geometry_fps"):
        paired=original.merge(ctrl[ctrl.method==baseline],on=["case","protein","analysis_set"],suffixes=("_is","_baseline"))
        for metric in metrics:
            primary=paired[paired.analysis_set=="primary"]
            deltas=primary.assign(diff=primary[metric+"_is"]-primary[metric+"_baseline"])
            byprotein=deltas.groupby("protein")["diff"].mean().to_numpy()
            draws=np.mean(rng.choice(byprotein,(4000,len(byprotein)),replace=True),axis=1)
            bootstrap.append({"baseline":baseline,"metric":metric,"n_proteins":len(byprotein),
                "mean_protein_cluster_delta":float(byprotein.mean()),
                "ci95_low":float(np.quantile(draws,.025)),
                "ci95_high":float(np.quantile(draws,.975))})
    csvwrite(path/"independent_baseline_protein_bootstrap.csv",bootstrap)
    (path/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("EXTENDED_VALIDATION_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--panel-shard",type=int)
    parser.add_argument("--shards",type=int,default=4)
    parser.add_argument("--cases",action="store_true")
    parser.add_argument("--aggregate",type=Path)
    args=parser.parse_args()
    if args.aggregate:aggregate(args.aggregate)
    elif args.cases:do_cases()
    elif args.panel_shard is not None:do_panel(args.panel_shard,args.shards)
    else:parser.error("Choose --panel-shard, --cases or --aggregate")
