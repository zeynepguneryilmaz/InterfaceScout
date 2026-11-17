"""Fixed-patch-pool Pareto component and rank sensitivity audit (13 conditions).
Scientific status: retrospective exploratory analysis; no held-out validation.
Freeze the patches: every ranking uses exactly the same InterfaceScout-generated candidates
and member sets for a given protein, pH, and chemistry channel.
"""
from __future__ import annotations
import csv,json,random,sys,math
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
from interface_engine import analyze_all_maps
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes
from coarse_patch import assign_pareto_fronts
OUT=ROOT/"validation"/"pareto_component_results";OUT.mkdir(parents=True,exist_ok=True)
CASES=json.loads((ROOT/"validation"/"benchmark_manifest.json").read_text())["cases"]
SURFACE={"silica":"anionic","paa_fe3o4":"anionic","citrate_au":"anionic","mpa_au":"anionic","alumoh":"cationic","polystyrene":"hydrophobic","graphene_oxide":"anionic"}
GO={
"G0":[9,46,66],
"G30":[6,7,10,13,15,17,18,22,23,27,29,34,35,36,43,45,46,50,60,64,67,68],
"G80":[2,6,9,10,15,34,43,44,45,46,50,59,64,67,68,70],
"G100":[2,3,5,6,7,8,9,10,13,15,27,34,35,43,45,46,50,55,59,60,64,66,67,68,70]}
for name,nums in GO.items():
 CASES.append({"id":"ubiquitin_GO_"+name,"analysis_set":"GO_new_same_protein_and_study","protein":"Human ubiquitin","pdb_id":"1UBQ","chain":"A","surface":"graphene_oxide","pH":6.0,"gt_kind":"exact","gt_exact":nums,"evidence_resolution":"NMR T2 perturbation, not atomic contact","doi":"10.1039/C9NA00807A"})
F=("chemistry_support","mean_accessibility","patch_coherence","orientation_coherence")
SEED=20261117;N_RANDOM=2000;KS=(1,3,5);DIST=(0,5,8)

def save(name,rows):
 if not rows:return
 with (OUT/name).open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def front_rank(patches,fields):
 """Nondominated sorting, using the exact public method's default secondary ordering."""
 remaining=set(range(len(patches)));order=[];front=1
 while remaining:
  current=[i for i in sorted(remaining) if not any(
    all(float(patches[j][f])>=float(patches[i][f]) for f in fields) and
    any(float(patches[j][f])>float(patches[i][f]) for f in fields)
    for j in remaining if i!=j)]
  if not current:raise RuntimeError("Empty nondominated front")
  current.sort(key=lambda i:(-float(patches[i]["patch_coherence"]),-float(patches[i]["chemistry_support"]),str(patches[i].get("center_key",""))))
  order+=current;remaining-=set(current);front+=1
 return order
def rank_candidates(patches):
 baseline=list(range(len(patches)))
 full=front_rank(patches,F)
 if baseline!=full:
  raise AssertionError("Our independent Pareto front implementation fails to reproduce frozen public ordering")
 alternatives={"InterfaceScout_full_Pareto":full}
 for drop in F:
  alternatives["Pareto_without_"+drop]=front_rank(patches,[x for x in F if x!=drop])
 for key in F:
  alternatives["Rank_by_"+key]=sorted(baseline,key=lambda i:(-float(patches[i][key]),str(patches[i].get("center_key",""))))
 alternatives["Equal_weight_four_descriptors"]=sorted(baseline,key=lambda i:(
  -sum(float(patches[i][key]) for key in F)/4.0,
  str(patches[i].get("center_key",""))))
 return alternatives
def gt_indices(c,nodes):
 ans=[]
 for i,r in enumerate(nodes):
  if r["chain"] not in c["chain"].split(","):continue
  num=int(r["res_seq"])
  if c["gt_kind"]=="exact" and num in c.get("gt_exact",[]):ans.append(i)
  if c["gt_kind"] in ("ranges","equivalent_chain_ranges") and any(lo<=num<=hi for lo,hi in c.get("gt_ranges",[])):ans.append(i)
 return ans
def recalled(members,gt,xyz,ch,equiv,cutoff):
 if not members:return 0.
 t=np.array(gt,dtype=int);m=np.array(sorted(members),dtype=int)
 if cutoff==0:
  hit=np.isin(t,m)
 else:
  hit=(np.linalg.norm(xyz[t,None,:]-xyz[None,m,:],axis=2)<=cutoff+1e-9).any(axis=1)
 if equiv:return max(float(hit[ch[t]==k].mean()) for k in set(ch[t]))
 return float(hit.mean())
def union(patches,members_by_patch,ordering,k):
 return set().union(*(members_by_patch[i] for i in ordering[:k])) if ordering else set()
def cachekey(c):return(c["pdb_id"],c["chain"],c["pH"],SURFACE[c["surface"]])
cache={}
rows=[];qc=[];null=[];rank_changes=[]
for case in CASES:
 cid=case["id"];key=cachekey(case)
 if key not in cache:
  pdb,chain,pH,surface=key
  files=sorted((ROOT/"examples"/pdb).rglob("*_Bfactor.pdb"))
  if pdb=="2CBA":files=[p for p in files if "pH_"+str(pH) in str(p)]
  if len(files)!=1:raise ValueError(f"{cid} structure ambiguous: {len(files)}")
  prepared,_=prepare_pdb_text(files[0].read_text(),chain=chain)
  nodes=extract_ca_nodes(prepared)
  out=analyze_all_maps(pdb_text=prepared,chain=None,pH=pH)
  patches=out["maps"][surface]["patches"]
  ind={n["key"]:i for i,n in enumerate(nodes)}
  members_by_patch=[{ind[k] for k in p["members"] if k in ind} for p in patches]
  raw_surface_keys=set(out["diagnostics"]["surface_residue_keys"])
  surface_set={ind[k] for k in raw_surface_keys if k in ind}
  if not surface_set:raise RuntimeError("No surface")
  cache[key]=(nodes,patches,members_by_patch,surface_set,rank_candidates(patches))
 nodes,patches,pool,surface,ranks=cache[key]
 gt=gt_indices(case,nodes)
 if not gt:raise ValueError(f"No experimental labels mapped for {cid}")
 xyz=np.array([n["coord"] for n in nodes],float);chain_labels=np.array([n["chain"] for n in nodes])
 equiv=case["gt_kind"]=="equivalent_chain_ranges"
 base=ranks["InterfaceScout_full_Pareto"]
 rng=random.Random(SEED+sum((i+1)*ord(x) for i,x in enumerate(cid)))
 for k in KS:
  target=union(patches,pool,base,k)
  if not target:continue
  cov_target=len(target)/len(surface)
  for model,ordered in ranks.items():
   u=union(patches,pool,ordered,k)
   relerr=abs(len(u)-len(target))/max(len(target),1)
   rec={d:recalled(u,gt,xyz,chain_labels,equiv,d) for d in DIST}
   rows.append({"case":cid,"group":case["analysis_set"],"pdb_id":case["pdb_id"],
    "model":model,"top_k":k,"n_candidate_patches":len(pool),"n_gt":len(gt),
    "selected_residues":len(u),"surface_residues":len(surface),"coverage":len(u)/len(surface),
    "target_coverage":cov_target,"relative_coverage_error":relerr,"coverage_matched_10pct":int(relerr<=.1),
    "rank1_same_as_full":int(ordered[0]==base[0]),"recall_direct":rec[0],"recall_5A":rec[5],"recall_8A":rec[8]})
  # Random ranking of exact same fixed patch pool; reject by coverage if >10% differs.
  # Sampling a new random list is not a new patch discovery and cannot use labels.
  null_dist={d:[] for d in DIST};accepted=0;all_trials=0
  for _ in range(N_RANDOM*8):
   if accepted>=N_RANDOM:break
   ix=rng.sample(range(len(pool)),min(k,len(pool)))
   u=union(patches,pool,ix,k)
   all_trials+=1
   if abs(len(u)-len(target))/max(len(target),1)>.1:continue
   accepted+=1
   for d in DIST:null_dist[d].append(recalled(u,gt,xyz,chain_labels,equiv,d))
  for d in DIST:
   vals=np.array(null_dist[d],float)
   original=recalled(target,gt,xyz,chain_labels,equiv,d)
   null.append({"case":cid,"group":case["analysis_set"],"top_k":k,"distance_A":d,
    "accepted":len(vals),"attempts":all_trials,"target_surface_coverage":cov_target,
    "model_recall":original,"null_mean":float(vals.mean()) if len(vals)>0 else None,
    "null_q025":float(np.quantile(vals,.025)) if len(vals)>0 else None,
    "null_q975":float(np.quantile(vals,.975)) if len(vals)>0 else None,
    "empirical_p_upper":float((1+(vals>=original-1e-12).sum())/(1+len(vals))) if len(vals)>0 else None})
 qc.append({"case":cid,"pdb":case["pdb_id"],"group":case["analysis_set"],"n_gt":len(gt),
  "n_patches":len(pool),"n_surface_residues":len(surface),"n_rank_methods":len(ranks),
  "all_same_candidates":1,"full_Pareto_reconstructed_exactly":1})
 print("PARETO_AUDIT_CASE",cid,"n_patches",len(pool),flush=True)
 save("all_models_at_k.csv",rows);save("matched_random_ranking_nulls.csv",null)
save("all_models_at_k.csv",rows);save("matched_random_ranking_nulls.csv",null);save("quality_controls.csv",qc)
# Descriptive independent-protein-cluster summary. Do NOT treat 13 conditions as 13 independent experimental replicates.
summary={}
for group in ("primary","GO_new_same_protein_and_study","secondary"):
 for k in (1,3,5):
  rr=[r for r in rows if r["group"]==group and r["top_k"]==k]
  models=sorted({r["model"] for r in rr})
  per=[]
  for model in models:
   rsub=[r for r in rr if r["model"]==model]
   if not rsub:continue
   per.append({"group":group,"top_k":k,"model":model,"n_cases":len(rsub),
    "n_unique_pdb":len(set(r["pdb_id"] for r in rsub)),
    "mean_direct":float(np.mean([r["recall_direct"] for r in rsub])),
    "mean_5A":float(np.mean([r["recall_5A"] for r in rsub])),
    "mean_8A":float(np.mean([r["recall_8A"] for r in rsub])),
    "mean_coverage":float(np.mean([r["coverage"] for r in rsub])),
    "n_coverage_matched_10pct":sum(r["coverage_matched_10pct"] for r in rsub),
    "rank1_agreement_with_full":float(np.mean([r["rank1_same_as_full"] for r in rsub]))})
  summary[group+"_"+str(k)]=per
(OUT/"group_descriptives.json").write_text(json.dumps(summary,indent=2))
print("PARETO_AUDIT_COMPLETE",json.dumps({"cases":len(CASES),"cached_structural_systems":len(cache),"rows":len(rows),"random_null_rows":len(null),"all_protein_cases_matched":all(q["all_same_candidates"] for q in qc)}),flush=True)
