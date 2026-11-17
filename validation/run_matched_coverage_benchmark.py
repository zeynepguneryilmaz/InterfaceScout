"""Matched-top-k interface benchmark. Reuses frozen public model and archived annotations.
WARNING: coverage and patch count matching are audited, not assumed. Study is exploratory.
"""
from __future__ import annotations
import csv,json,math,random,sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import wilcoxon
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
from prepare import prepare_pdb_text
from interface_engine import analyze_all_maps
from geometry import extract_ca_nodes,build_surface_geometry,ca_distance,same_face
import core
OUT=ROOT/"validation"/"matched_coverage_results";OUT.mkdir(parents=True,exist_ok=True)
CASES=json.loads((ROOT/"validation"/"benchmark_manifest.json").read_text())["cases"]
CHEM={"silica":"anionic","paa_fe3o4":"anionic","citrate_au":"anionic","mpa_au":"anionic","alumoh":"cationic","polystyrene":"hydrophobic"}
GO={
"G0":[9,46,66],
"G30":[6,7,10,13,15,17,18,22,23,27,29,34,35,36,43,45,46,50,60,64,67,68],
"G80":[2,6,9,10,15,34,43,44,45,46,50,59,64,67,68,70],
"G100":[2,3,5,6,7,8,9,10,13,15,27,34,35,43,45,46,50,55,59,60,64,66,67,68,70],
}
for code,nums in GO.items():
 CASES.append({"id":"ubiquitin_GO_"+code,"analysis_set":"external_same_protein_multiple_materials","protein":"Human ubiquitin","pdb_id":"1UBQ","chain":"A","surface":"graphene_oxide","pH":6.0,"gt_kind":"exact","gt_exact":nums,
 "evidence_resolution":"NMR_relaxation_at_least_75pct_T2_decrease_NOT_proven_direct_contact","doi":"10.1039/C9NA00807A","source_url":"https://pubs.rsc.org/en/content/articlehtml/2020/na/c9na00807a",
 "material_mapping":"anionic_proxy_FIXED_FOR_ALL_FOUR; oxidation percentage NOT represented"})
CHEM["graphene_oxide"]="anionic"
KS=(1,3,5); THRESHOLDS=(0,5,8); ROUNDS=2000;SEED=20261117
ALL20=set("ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
def write_csv(name,rows):
 if not rows:return
 with open(OUT/name,"w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def gt_ids(case,nodes):
 out=[]
 for i,r in enumerate(nodes):
  if r["chain"] not in case["chain"].split(","):continue
  n=int(r["res_seq"])
  if case["gt_kind"]=="exact" and n in case.get("gt_exact",[]):out.append(i)
  if case["gt_kind"] in ("ranges","equivalent_chain_ranges") and any(a<=n<=b for a,b in case.get("gt_ranges",[])):out.append(i)
 return out
def metric(mask,truth,xyz,chain_labels,equiv,d):
 if not mask or not truth:return 0.
 jj=np.array(sorted(mask));t=np.array(truth)
 if d==0:hit=np.isin(t,jj)
 else:hit=(np.linalg.norm(xyz[t,None,:]-xyz[None,jj,:],axis=-1)<=d+1.e-9).any(axis=1)
 if equiv:return max(float(np.mean(hit[chain_labels[t]==ch])) for ch in set(chain_labels[t]))
 return float(hit.mean())
def candidate_patches(keys,nodes,geometry,index):
 return [frozenset(index[k] for k in keys if ca_distance(center,k,geometry)<=8.0 and same_face(center,k,geometry)) for center in keys]
def matched_ranked(candidates,ranks,k,target,seed):
 """Fixed k patches, choose candidate sequence closest to target coverage from ranked priority.
 No experimental labels are used. No patch membership artificially altered.
 """
 if not candidates or not k:return set(),[],False
 # Stable restriction: prefer higher model-independent priority, while matching union size.
 ordered=sorted(range(len(candidates)),key=lambda i:(-ranks[i],i))
 selected=[];mask=set();used=set()
 for t in range(k):
  choices=ordered[:min(100,len(ordered))]
  valid=[j for j in choices if j not in used]
  if not valid:break
  # Penalize both priority rank and coverage deviation, no label-informed terms.
  target_partial=target*(t+1)/k
  pick=min(valid,key=lambda j:(abs(len(mask|candidates[j])-target_partial)/max(1,target),
    ordered.index(j)/max(1,len(ordered)),j))
  selected.append(pick);used.add(pick);mask.update(candidates[pick])
 return mask,selected,abs(len(mask)-target)/max(1,target)<=.10
def rand_matched(candidates,k,target,rng,max_attempts=ROUNDS):
 """Random geometric centers; filter to <=10% relative union-coverage difference."""
 valid=[];count=0
 for attempt in range(max_attempts*3):
  idx=rng.sample(range(len(candidates)),min(k,len(candidates)))
  union=set().union(*(candidates[i] for i in idx))
  if abs(len(union)-target)/max(1,target)<=.10:valid.append((union,idx))
  if len(valid)>=max_attempts:break
 return valid,attempt+1
all_rows=[];qc=[];null_rows=[]
for case in CASES:
 pid=case["pdb_id"];cid=case["id"];ch=CHEM[case["surface"]]
 files=sorted((ROOT/"examples"/pid).rglob("*_Bfactor.pdb"))
 if pid=="2CBA":files=[f for f in files if "pH_"+str(case["pH"]) in str(f)]
 if len(files)!=1:raise RuntimeError(f"{cid}: missing unique frozen structure")
 prepared,_=prepare_pdb_text(files[0].read_text(),chain=case["chain"])
 nodes=extract_ca_nodes(prepared)
 xyz=np.asarray([n["coord"] for n in nodes],float);chains=np.asarray([n["chain"] for n in nodes])
 index={n["key"]:i for i,n in enumerate(nodes)}
 gt=gt_ids(case,nodes)
 if not gt:raise RuntimeError(f"{cid}: no matching truth residues")
 out=analyze_all_maps(pdb_text=prepared,chain=None,pH=case["pH"])
 patches=out["maps"][ch]["patches"]
 model=[set(index[k] for k in p["members"] if k in index) for p in patches]
 v1=core.analyze(core.AnalyzeRequest(pdb_text=prepared,env=core.EnvParams(pH=case["pH"],ionic=150,temp=298)))
 geom=build_surface_geometry(v1,nodes)
 keys=sorted(geom["coords"])
 candidates=candidate_patches(keys,nodes,geom,index)
 surface={index[k] for k in keys}
 if not candidates:raise RuntimeError(cid+" no geometric patch candidates")
 channel=v1["chemistries"][ch]
 local={r["key"]:r["local_score"] for r in channel["residues"]}
 sasa=np.asarray([geom["scrsa"][k] for k in keys])
 lscore=np.asarray([float(local.get(k,0.)) for k in keys])
 # Method 1: solvent accessibility + geometry, no chemistry.
 score_sasa=np.asarray([np.mean(sasa[[keys.index(nodes[j]["key"]) for j in patch if nodes[j]["key"] in keys]]) if patch else 0 for patch in candidates])
 # Method 2: chemical compatibility only at patch center, without multiscale density or Pareto.
 score_chem=lscore
 # Method 3: chemistry-free 6/9 density ranked geometric patch; full model's patch radius and same-face preserved.
 dst=np.linalg.norm(np.asarray([geom["coords"][k] for k in keys])[:,None,:]-np.asarray([geom["coords"][k] for k in keys])[None,:,:],axis=-1)
 d6=(dst<=6)@sasa;d9=(dst<=9)@sasa
 n6=d6/d6.max() if d6.max()>0 else d6
 n9=d9/d9.max() if d9.max()>0 else d9
 score_geom=np.minimum(n6,n9)
 rng=random.Random(SEED+sum(map(ord,cid)))
 for k in KS:
  original=set().union(*(model[:k])) if model else set()
  target=len(original)
  if target==0:continue
  pred_sets={"InterfaceScout":(original,list(range(min(k,len(model)))),True)}
  for name,score in [("SASA_geometry",score_sasa),("chemistry_only",score_chem),("no_chem_multiscale",score_geom)]:
   pred_sets[name]=matched_ranked(candidates,score,k,target,SEED)
  rand,attempts=rand_matched(candidates,k,target,rng)
  qc.append({"case":cid,"set":case["analysis_set"],"pdb_id":pid,"k":k,"n_gt":len(gt),
   "n_surface":len(surface),"target_n_residues":target,"target_coverage":target/len(surface),
   "random_accepted":len(rand),"random_attempts":attempts,
   "sample_valid":int(len(rand)>=100)})
  for modelname,(chosen,centers,within) in pred_sets.items():
   for dist in THRESHOLDS:
    all_rows.append({"case":cid,"pdb_id":pid,"analysis_set":case["analysis_set"],"chemistry":ch,
      "model":modelname,"top_k":k,"distance_A":dist,"n_truth":len(gt),
      "selected_patches":len(centers),"selected_residues":len(chosen),"coverage":len(chosen)/len(surface),
      "coverage_error_fraction":abs(len(chosen)-target)/target,
      "within_10pct_coverage":int(within),"recall":metric(chosen,gt,xyz,chains,case["gt_kind"]=="equivalent_chain_ranges",dist)})
  for dist in THRESHOLDS:
   values=np.asarray([metric(s,gt,xyz,chains,case["gt_kind"]=="equivalent_chain_ranges",dist) for s,_ in rand])
   val=metric(original,gt,xyz,chains,case["gt_kind"]=="equivalent_chain_ranges",dist)
   null_rows.append({"case":cid,"analysis_set":case["analysis_set"],"top_k":k,"distance_A":dist,"random_n":len(values),
     "model_recall":val,"random_mean":float(values.mean()) if len(values) else float("nan"),
     "random_q025":float(np.quantile(values,.025)) if len(values) else float("nan"),
     "random_q975":float(np.quantile(values,.975)) if len(values) else float("nan"),
     "empirical_p":float((1+np.sum(values>=val))/(1+len(values))) if len(values) else float("nan")})
 print("MATCHED_CASE",cid,"n_truth",len(gt),"n_patches",len(model),flush=True)
 write_csv("case_results.csv",all_rows)
 write_csv("matched_random_nulls.csv",null_rows)
 write_csv("coverage_and_sampling_qc.csv",qc)
write_csv("case_results.csv",all_rows);write_csv("matched_random_nulls.csv",null_rows);write_csv("coverage_and_sampling_qc.csv",qc)
summary={"cases":len(CASES),"distinct_proteins":len(set(x["pdb_id"] for x in CASES)),"primary_cases":sum(x["analysis_set"]=="primary" for x in CASES),
"new_GO_conditions":4,"new_GO_distinct_proteins":1,
"comparators":["SASA_geometry","chemistry_only","no_chem_multiscale","matched_random_geometry"],
"coverage_acceptance":.10,"random_trials_per_case":ROUNDS,"seed":SEED,"no_posthoc_best_chemistry":True,
"known_limits":["same 1UBQ repeated for 4 GO coatings; one study-family cluster",
"GO oxidation fraction is not represented; anionic channel identical for all four coatings",
"Chemistry-only is a center-priority baseline, not a molecular physical adsorption simulation",
"Excessively broad regional truth sets are not exact contact atoms",
"Comparators with >10% union-coverage mismatch must NOT be used for primary matched comparisons",
"Random nulls with fewer than 100 accepted draws must not support an inference",
"All nine earlier cases were viewed in previous model development, so comparison is retrospective, not prospective"] }
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("MATCHED_BENCHMARK_COMPLETE",json.dumps(summary),flush=True)
