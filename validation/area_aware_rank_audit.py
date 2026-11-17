"""InterfaceScout v2 exploratory area-aware ranking study.
No tuning to experimental labels; the original v1 patches and structural coordinates remain frozen.
Primary matching criterion: Top-k surface SIDECHAIN SASA union / total exposed sidechain SASA, tolerance ±10%.
Secondary matching criterion: same Top-k residue count ±10%. pH as source benchmark.
"""
from pathlib import Path
import csv,json,sys
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"validation"/"area_aware_diagnostic";OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/"backend"))
from prepare import prepare_pdb_text
from interface_engine import analyze_all_maps
from geometry import extract_ca_nodes
import core
CASES=json.loads((ROOT/"validation"/"benchmark_manifest.json").read_text())["cases"]
SURFACE={"silica":"anionic","paa_fe3o4":"anionic","citrate_au":"anionic","mpa_au":"anionic","alumoh":"cationic","polystyrene":"hydrophobic","graphene_oxide":"anionic"}
GO={"G0":[9,46,66],
"G30":[6,7,10,13,15,17,18,22,23,27,29,34,35,36,43,45,46,50,60,64,67,68],
"G80":[2,6,9,10,15,34,43,44,45,46,50,59,64,67,68,70],
"G100":[2,3,5,6,7,8,9,10,13,15,27,34,35,43,45,46,50,55,59,60,64,66,67,68,70]}
for code,numbers in GO.items():
 CASES.append({"id":"ubiquitin_GO_"+code,"analysis_set":"same_protein_new_GO_study","protein":"Human ubiquitin","pdb_id":"1UBQ","chain":"A","surface":"graphene_oxide","pH":6.0,"gt_kind":"exact","gt_exact":numbers,"evidence_resolution":"NMR_relaxation_at_least_75pct_T2_decrease"})
def csv_save(file,rows):
 if not rows:return
 with (OUT/file).open("w",newline="",encoding="utf8") as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def truth_indices(c,nodes):
 return [i for i,n in enumerate(nodes) if n["chain"] in c["chain"].split(",") and
 ((c["gt_kind"]=="exact" and n["res_seq"] in c.get("gt_exact",[])) or
 (c["gt_kind"] in ("ranges","equivalent_chain_ranges") and any(a<=n["res_seq"]<=b for a,b in c.get("gt_ranges",[]))))]
def recall(union,gt,nodes,cutoff,eq):
 if not union:return 0.
 xyz=np.asarray([n["coord"] for n in nodes],float)
 idx={n["key"]:i for i,n in enumerate(nodes)}
 t=np.asarray(gt,dtype=int);a=np.asarray([idx[x] for x in union if x in idx],dtype=int)
 if not len(a):return 0.
 if cutoff==0:hit=np.isin(t,a)
 else:hit=(np.linalg.norm(xyz[t,None,:]-xyz[None,a,:],axis=2)<=cutoff+1e-9).any(axis=1)
 if eq:
  chains=np.asarray([n["chain"] for n in nodes])
  return max(float(hit[chains[t]==ch].mean()) for ch in set(chains[t]))
 return float(hit.mean())
def results(ordered,k):
 p=ordered[:k]
 return set(x for z in p for x in z["members"])
def diversity_rank(patches):
 """Non-tuned diversity-aware rule, applied to fixed patches without experimental labels.
 At each step rank density × coherence × orientation, penalized by fractional
 residue reuse. This is an exploratory proposal, not a previously validated estimator.
 """
 rem=list(patches);ans=[];union=set()
 while rem:
  pos=max(range(len(rem)),key=lambda i:(
    rem[i]["score_composite"]*(1.-len(set(rem[i]["members"])&union)/max(1,len(rem[i]["members"]))),
    rem[i]["score_composite"],-rem[i]["display_rank"]))
  chosen=rem.pop(pos);ans.append(chosen);union.update(chosen["members"])
 return ans
def constrained_swap(original_patches,k,field,surface_rows):
 """Exploratory deterministic single-swap local search under both ±10% budgets.
 Start from v1 Top-k (a feasible set), replace one candidate with a higher
 intrinsic area score only when both union residue count and exposed SASA
 stay close to v1. NO ground-truth labels, no model retraining.
 Results are SETS chosen separately for each k, not a nested Top-k ordering.
 """
 n=len(original_patches)
 top=[original_patches.index(q) for q in original_patches[:k]]
 mems=[set(q["members"]) for q in original_patches]
 sas={r:surface_rows[r]["sidechain_sasa"] for q in original_patches for r in q["members"]}
 def union(which):return set().union(*(mems[j] for j in which))
 tgt=union(top);budget_n=len(tgt);budget_a=sum(sas[x] for x in tgt)
 scores=[float(q[field]) for q in original_patches]
 def score(which):return sum(scores[j] for j in which)
 old=score(top)
 for round_ix in range(12):
  best=None
  for ii in range(k):
   for new in range(n):
    if new in top:continue
    test=top[:];test[ii]=new
    used=union(test)
    if abs(len(used)-budget_n)/max(budget_n,1)>.10:continue
    area=sum(sas[x] for x in used)
    if abs(area-budget_a)/max(budget_a,1)>.10:continue
    value=score(test)
    if value<=old+1e-10:continue
    if best is None or (value,-ii,-new)>(best[0],-best[1],-best[2]):best=(value,ii,new)
  if best is None:break
  old,ii,new=best;top[ii]=new
 selected=[original_patches[j] for j in top]
 return selected+[q for j,q in enumerate(original_patches) if j not in set(top)]

rows=[];qc=[];details=[]
for case in CASES:
 cid=case["id"];pid=case["pdb_id"];ph=case["pH"]
 files=sorted((ROOT/"examples"/pid).rglob("*_Bfactor.pdb"))
 if pid=="2CBA":files=[f for f in files if "pH_"+str(ph) in str(f)]
 if len(files)!=1:raise RuntimeError(f"{cid}: frozen structure not uniquely located")
 prepared,_=prepare_pdb_text(files[0].read_text(),chain=case["chain"])
 nodes=extract_ca_nodes(prepared)
 gt=truth_indices(case,nodes)
 if not gt:raise RuntimeError(cid+" no ground truth mapped")
 raw=core.analyze(core.AnalyzeRequest(pdb_text=prepared,env=core.EnvParams(pH=ph,ionic=150,temp=298)))
 exposed={r["key"]:r for r in raw["surface_residues"]}
 # use direct prepared structure; on exactly the same coordinates as ranked patches
 m=analyze_all_maps(pdb_text=prepared,chain=None,pH=ph)["maps"][SURFACE[case["surface"]]]
 p=m["patches"]
 chem=core.CHEMISTRIES[SURFACE[case["surface"]]]["favorable"]
 if not p:raise ValueError(cid+" empty patch pool")
 total_surface_area=sum(r["sidechain_sasa"] for r in exposed.values())
 out=[]
 for item in p:
  mm=item["members"];ss=item["seed_members"]
  if not all(k in exposed for k in mm):raise ValueError("missing patch member from exposed SASA")
  A=sum(exposed[k]["sidechain_sasa"] for k in mm)
  Ac=sum(exposed[k]["sidechain_sasa"]*core.state_availability(exposed[k]["res_name"],chem[exposed[k]["res_name"]][2],ph) for k in ss)
  D=Ac/A if A else 0
  rr=dict(item)
  rr["total_area_A2"]=A;rr["compatible_area_A2"]=Ac;rr["compatible_area_fraction"]=D
  rr["score_composite"]=D*float(item["patch_coherence"])*float(item["orientation_coherence"])
  out.append(rr)
 # Fixed-pool rankings, no retraining and no new patches
 rankings={
 "InterfaceScout_v1":out,
 "Area_total":sorted(out,key=lambda z:(-z["compatible_area_A2"],z["center_key"])),
 "Area_fraction":sorted(out,key=lambda z:(-z["compatible_area_fraction"],z["center_key"])),
 "Area_fraction_x_coherence_x_orientation":sorted(out,key=lambda z:(-z["score_composite"],z["center_key"])),
 "Area_fraction_x_coherence_x_orientation_diversity":diversity_rank(out)
 }
 for q in out:
  details.append({"case":cid,"pdb_id":pid,"chemistry":SURFACE[case["surface"]],"display_rank":q["display_rank"],"center_key":q["center_key"],
   "n_members":q["n_members"],"total_area_A2":q["total_area_A2"],"compatible_area_A2":q["compatible_area_A2"],
   "compatible_area_fraction":q["compatible_area_fraction"],"composite_score":q["score_composite"],
   "pareto_front":q["pareto_front"]})
 for k in (1,3,5):
  original=results(out,k)
  budget_n=len(original);budget_area=sum(exposed[r]["sidechain_sasa"] for r in original)
  for name,field in [("BudgetMatched_Area_total","compatible_area_A2"),("BudgetMatched_Area_fraction","compatible_area_fraction"),("BudgetMatched_Area_x_coherence_x_orientation","score_composite")]:
   rankings[name]=constrained_swap(out,k,field,exposed)
  for method,ranked in rankings.items():
   selected=results(ranked,k)
   used_area=sum(exposed[r]["sidechain_sasa"] for r in selected)
   err_a=abs(used_area-budget_area)/budget_area if budget_area else 0
   err_n=abs(len(selected)-budget_n)/budget_n if budget_n else 0
   row={"case":cid,"set":case["analysis_set"],"pdb_id":pid,"method":method,
     "k":k,"n_ground_truth":len(gt),"n_patches_pool":len(out),"selected_residue_count":len(selected),
     "selected_sidechain_SASA_A2":used_area,"surface_SASA_A2":total_surface_area,
     "sidechain_area_fraction":used_area/total_surface_area,
     "reference_sidechain_area_fraction":budget_area/total_surface_area,
     "relative_SASA_coverage_error":err_a,"relative_residue_coverage_error":err_n,
     "SASA_area_matched_10pct":int(err_a<=.1),"residue_count_matched_10pct":int(err_n<=.1),
     "both_coverage_matched_10pct":int(err_a<=.1 and err_n<=.1),
     "top1_same_original":int(ranked[0]["center_key"]==out[0]["center_key"]),
     "selected_centers_changed":sum(q["center_key"] not in {v["center_key"] for v in out[:k]} for q in ranked[:k])}
   for d in (0,5,8):
    row["recall_direct" if d==0 else "recall_"+str(d)+"A"]=recall(selected,gt,nodes,d,case["gt_kind"]=="equivalent_chain_ranges")
   rows.append(row)
 qc.append({"case":cid,"n_patches":len(out),"n_truth":len(gt),"n_surface_residues":len(exposed),
  "surface_sidechain_SASA_A2":total_surface_area,"n_ranking_rules":len(rankings),
  "all_share_same_patch_pool":1})
 print("AREA_BENCHMARK_CASE",cid,len(out),flush=True)
 csv_save("casewise_area_ranking_results.csv",rows)
csv_save("casewise_area_ranking_results.csv",rows)
csv_save("patch_area_descriptors_13_cases.csv",details);csv_save("case_qc.csv",qc)
summ=[]
for grp in ("primary","secondary","same_protein_new_GO_study"):
 for k in (1,3,5):
  for method in sorted(set(x["method"] for x in rows)):
   sub=[x for x in rows if x["set"]==grp and x["k"]==k and x["method"]==method]
   matched=[x for x in sub if x["both_coverage_matched_10pct"]]
   summ.append({"set":grp,"k":k,"method":method,"n_cases":len(sub),"n_unique_proteins":len({x["pdb_id"] for x in sub}),
    "n_both_budget_matched":len(matched),
    "mean_direct_ALL_UNMATCHED":float(np.mean([x["recall_direct"] for x in sub])) if sub else None,
    "mean_direct_MATCHED_ONLY":float(np.mean([x["recall_direct"] for x in matched])) if matched else None,
    "mean_near5_MATCHED_ONLY":float(np.mean([x["recall_5A"] for x in matched])) if matched else None,
    "mean_near8_MATCHED_ONLY":float(np.mean([x["recall_8A"] for x in matched])) if matched else None})
csv_save("group_summary_coverage_caution.csv",summ)
report={"n_conditions":len(CASES),"n_unique_proteins":len({x["pdb_id"] for x in CASES}),
 "n_records":len(rows),"n_patch_descriptors":len(details),
 "top5_budget_matched_counts":{method:sum(x["both_coverage_matched_10pct"] for x in rows if x["k"]==5 and x["method"]==method) for method in sorted(set(x["method"] for x in rows))},
 "scope":"Exploratory performance audit with labels already examined. Not independent confirmation or optimization.",
 "fixed_patch_pool":True,"fixed_chemistry_mapping":True,
 "GO_limits":"All GO oxidation levels assigned the same anionic proxy. Four conditions are one protein and one study.",
 "area_definition":"native exposed side-chain area sum, never nanoparticle-buried interface area",
 "primary_matched_definition":"both exposed sidechain area and residue count matched to InterfaceScout v1 TopK within ±10%",
 "budgeted_variant":"deterministic single-swap intrinsic-score optimization initialized at v1 Top-k, no labels; returns nonnested separately optimized Top-k SETS, not a single ranking; potential v1 initialization bias; transparent count of changed patches",
 "do_not_claim":"Do not compare unmatched method group means as equal-budget predictive performance."}
(OUT/"summary.json").write_text(json.dumps(report,indent=2))
print("AREA_BENCHMARK_COMPLETE",json.dumps(report),flush=True)
