"""Reproduce current production-model evidence concordance and fixed-candidate baselines.

The nine previously curated literature conditions are NOT an untouched holdout.
Same exact production patches are ranked by (1) public Pareto/display rank,
(2) chemistry support alone, or (3) mean tRSA alone. Comparisons isolate the
value of multi-criterion *prioritization*, NOT of discovering patches from a
random sequence location or a new material-specific adsorption prediction.
"""
from __future__ import annotations
import sys,csv,json,hashlib,traceback,argparse
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
import dual_literature_benchmark as bench
from interface_engine import analyze_all_maps
from model_settings import MODEL_VERSION
MANIFEST=bench.MANIFEST
OUT=Path("final_v110_literature");OUT.mkdir(exist_ok=True)
METRICS=("top5_hit_8A","top5_recall_8A","top5_direct_hit","top5_direct_recall")
METHODS=("pareto","chemistry_only","accessibility_only")
def put(path,items):
    if not items:return
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
def measure(patches,labels,coords):
    chosen=patches[:5]
    centers=[p["center_key"] for p in chosen if p["center_key"] in coords]
    members=set(k for p in chosen for k in p["members"])
    near={k for k in labels if k in coords and any(
        np.linalg.norm(coords[k]-coords[c])<=8.0 for c in centers)}
    return {"top5_hit_8A":int(bool(near)),"top5_recall_8A":len(near)/len(labels),
        "top5_direct_hit":int(bool(members & labels)),
        "top5_direct_recall":len(members & labels)/len(labels),
        "n_patches_top5":len(chosen),"n_top5_member_residues":len(members),
        "n_top5_centers":len(centers)}
def compute(case):
    pdb,provenance=bench.prepare(case)
    pH=float(case["pH"])
    result=analyze_all_maps(pdb_text=pdb,chain=None,pH=pH,ionic_mM=150.,temp_K=298.)
    if result["version"]!=MODEL_VERSION or result["n_maps"]!=11:raise ValueError("Wrong model version")
    if result["method"]["surface_gate"]!="total_residue_RSA":raise ValueError("Wrong surface gate")
    profiles=result["residue_profiles"]
    labels=bench.labels_for_case(case,profiles)
    coords=result["diagnostics"]["surface_ca_xyz_A"]
    # Labels are experimental regions, some buried. Coordinates for the metric
    # must include all annotated residues, not merely exposed candidates.
    from geometry import extract_ca_nodes
    coords_all={str(node["key"]):np.asarray(node["coord"],dtype=float)
        for node in extract_ca_nodes(pdb)}
    unmatched=labels-set(coords_all)
    if unmatched:raise ValueError(f"Experimental labels without coordinates: {sorted(unmatched)[:7]}")
    rows=[];details=[]
    for channel in result["map_order"]:
        patches=result["maps"][channel]["patches"]
        sorted_chem=sorted(patches,key=lambda p:(-float(p["chemistry_support"]),str(p["center_key"])))
        sorted_trsa=sorted(patches,key=lambda p:(-float(p["mean_accessibility"]),str(p["center_key"])))
        for method,candidates in [("pareto",patches),("chemistry_only",sorted_chem),
                                  ("accessibility_only",sorted_trsa)]:
            metrics=measure(candidates,labels,coords_all)
            rows.append({"case":case["id"],"protein":case["pdb_id"],
                "analysis_set":case["analysis_set"],"resolution":case["evidence_resolution"],
                "channel":channel,"method":method,"n_evidence_residues":len(labels),
                "n_total_patches":len(patches),**metrics})
        details.append({"case":case["id"],"channel":channel,
            "n_patches":len(patches),
            "pareto_top5":";".join(str(p["center_key"]) for p in patches[:5]),
            "chemistry_top5":";".join(str(p["center_key"]) for p in sorted_chem[:5]),
            "accessibility_top5":";".join(str(p["center_key"]) for p in sorted_trsa[:5])})
    provenance.update({"case":case["id"],"model":MODEL_VERSION,"pH":pH,
                       "n_evidence_residues":len(labels),"n_maps":len(result["map_order"]),
                       "evidence_resolution":case["evidence_resolution"],
                       "status":"PASS"})
    return rows,details,provenance
def shard(i,n):
    rows=[];details=[];provenance=[];errors=[]
    for j,case in enumerate(MANIFEST["cases"]):
        if j%n!=i:continue
        try:
            a,b,c=compute(case);rows+=a;details+=b;provenance.append(c)
            print("FINAL_LITERATURE_PASS",json.dumps(c),flush=True)
        except Exception as e:
            problem={"case":case["id"],"error":str(e),"trace":traceback.format_exc()[-1700:]}
            errors.append(problem);print("FINAL_LITERATURE_FAIL",json.dumps(problem),flush=True)
    put(OUT/f"results_{i}.csv",rows);put(OUT/f"ranks_{i}.csv",details)
    (OUT/f"provenance_{i}.json").write_text(json.dumps(provenance,indent=2))
    (OUT/f"errors_{i}.json").write_text(json.dumps(errors,indent=2))
    if errors:raise SystemExit(2)
def aggregate(path):
    import pandas as pd
    p=sorted(path.glob("results_[0-9].csv"))
    if len(p)!=3:raise ValueError("Missing result shard")
    df=pd.concat([pd.read_csv(x) for x in p],ignore_index=True)
    if len(df)!=9*11*3 or df.case.nunique()!=9 or df.channel.nunique()!=11:
        raise ValueError("Incomplete 9x11x3 benchmark")
    df.to_csv(path/"all_297_condition_channel_method.csv",index=False)
    rank=pd.concat([pd.read_csv(x) for x in sorted(path.glob("ranks_[0-9].csv"))],ignore_index=True)
    rank.to_csv(path/"top5_by_case_channel_method.csv",index=False)
    means=df.groupby(["analysis_set","case","protein","method"],as_index=False)[list(METRICS)].mean()
    means.to_csv(path/"per_case_all_map_means.csv",index=False)
    long=[]
    for (grp,case,protein),part in means.groupby(["analysis_set","case","protein"]):
        assert set(part.method)==set(METHODS)
        base=part.set_index("method")
        for method in METHODS[1:]:
            for metric in METRICS:
                long.append({"analysis_set":grp,"case":case,"protein":protein,
                    "baseline":method,"metric":metric,
                    "pareto_value":float(base.loc["pareto",metric]),
                    "baseline_value":float(base.loc[method,metric]),
                    "paired_delta":float(base.loc["pareto",metric]-base.loc[method,metric])})
    deltas=pd.DataFrame(long)
    deltas.to_csv(path/"per_case_paired_deltas.csv",index=False)
    # Condition means are descriptive: two carbonic anhydrase cases are dependent.
    # Protein-cluster bootstrap resamples the biological protein, not 11 maps.
    rng=np.random.default_rng(20261123)
    bootstrap=[]
    for (group,baseline,metric),part in deltas.groupby(["analysis_set","baseline","metric"]):
        byprot=part.groupby("protein").paired_delta.mean().to_numpy()
        draws=np.mean(rng.choice(byprot,(4000,len(byprot)),replace=True),axis=1)
        bootstrap.append({"analysis_set":group,"baseline":baseline,"metric":metric,
            "n_proteins":len(byprot),"mean_paired_delta":float(byprot.mean()),
            "protein_bootstrap_ci_low":float(np.quantile(draws,.025)),
            "protein_bootstrap_ci_high":float(np.quantile(draws,.975))})
    pd.DataFrame(bootstrap).to_csv(path/"protein_cluster_bootstrap.csv",index=False)
    summary={"model":MODEL_VERSION,"status":"PASS","conditions":9,"unique_PDBs":df.protein.nunique(),
        "n_comparisons":len(df),"primary_n_cases":7,"secondary_n_cases":2,
        "all_11_map_condition_averages":means.groupby(["analysis_set","method"])[list(METRICS)].mean().reset_index().to_dict("records"),
        "limits":"Retrospective, non-held-out exploratory localization; broad evidence regions differ from exact contact labels; identical candidate pool means baseline is rank ablation, not full-method independent control; map channels correlated."}
    (path/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("FINAL_LITERATURE_SUMMARY",json.dumps(summary),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--shard",type=int);p.add_argument("--shards",type=int,default=3)
    p.add_argument("--aggregate",type=Path);a=p.parse_args()
    if a.aggregate:aggregate(a.aggregate)
    elif a.shard is not None:shard(a.shard,a.shards)
    else:p.error("shard or aggregate")
