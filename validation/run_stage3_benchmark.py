"""Stage-3 independent benchmark audit; no refitting and no claims of independent experimental replication."""
from __future__ import annotations
import csv,json,sys,hashlib,random
from pathlib import Path
from collections import Counter
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from prepare import prepare_pdb_text
from interface_engine import analyze_all_maps
from geometry import extract_ca_nodes,build_surface_geometry,ca_distance,same_face
import core
BASE=Path(__file__).resolve().parents[1]
OUT=BASE/"validation"/"stage3_results";OUT.mkdir(parents=True,exist_ok=True)
CASES=json.loads((BASE/"validation"/"benchmark_manifest.json").read_text())["cases"]
FROZEN={r["case_id"]:r for r in csv.DictReader((BASE/"validation"/"benchmark_case_results.csv").open(encoding="utf-8-sig"))}
SURFACE={"silica":"anionic","paa_fe3o4":"anionic","citrate_au":"anionic","mpa_au":"anionic","alumoh":"cationic","polystyrene":"hydrophobic"}
N=int(sys.argv[1]) if len(sys.argv)>1 else 2000
SEED=20261116
ks=(1,3,5)
metrics=(("direct",0.0),("near5",5.0),("near8",8.0))
def rowsave(name,rows):
 if not rows:return
 with (OUT/name).open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def truth_idx(case,nodes):
 out=[]
 for i,n in enumerate(nodes):
  if n["chain"] not in case["chain"].split(","):continue
  seq=n["res_seq"]
  if case["gt_kind"]=="exact" and seq in case.get("gt_exact",[]):out.append(i)
  if case["gt_kind"] in ("ranges","equivalent_chain_ranges") and any(a<=seq<=b for a,b in case.get("gt_ranges",[])):out.append(i)
 return out
def score(members,truth,xyz,nodes,equiv,threshold):
 if not members:return 0.
 txyz=xyz[truth]
 mxyz=xyz[sorted(members)]
 if threshold==0:
  hits=np.isin(truth,list(members))
 else:
  dist=np.linalg.norm(txyz[:,None,:]-mxyz[None,:,:],axis=-1)
  hits=(dist<=threshold+1e-9).any(axis=1)
 if equiv:
  chains=np.array([nodes[i]["chain"] for i in truth])
  return max(float(hits[chains==ch].mean()) for ch in set(chains))
 return float(hits.mean())
def pvalue(null,observed):
 return (1+sum(x>=observed-1e-12 for x in null))/(1+len(null))
obs_rows=[];null_rows=[];case_qc=[]
for caseidx,case in enumerate(CASES):
 cid=case["id"];pdb=case["pdb_id"];ph=case["pH"]
 candidates=sorted((BASE/"examples"/pdb).rglob("*_Bfactor.pdb"))
 if pdb=="2CBA":candidates=[p for p in candidates if "pH_"+str(ph) in str(p)]
 if len(candidates)!=1:raise RuntimeError(f"{cid}: unique frozen example not found: {len(candidates)}")
 prepared,_=prepare_pdb_text(candidates[0].read_text(),chain=case["chain"])
 nodes=extract_ca_nodes(prepared);truth=truth_idx(case,nodes)
 if not truth:raise RuntimeError(f"{cid}: no ground-truth residues mapped")
 out=analyze_all_maps(pdb_text=prepared,chain=None,pH=ph)
 patches=out["maps"][SURFACE[case["surface"]]]["patches"]
 index={r["key"]:i for i,r in enumerate(nodes)}
 xyz=np.array([r["coord"] for r in nodes])
 equiv=case["gt_kind"]=="equivalent_chain_ranges"
 observed_patch_members=[{index[m] for m in p["members"] if m in index} for p in patches]
 union=set().union(*observed_patch_members[:5]) if observed_patch_members else set()
 val=score(union,truth,xyz,nodes,equiv,8.)
 saved=float(FROZEN[cid]["top5_near_8A_recall"])
 if abs(saved-val)>0.015:raise RuntimeError(f"{cid}: frozen benchmark mismatch {saved} vs {val}")
 v1=core.analyze(core.AnalyzeRequest(pdb_text=prepared,env=core.EnvParams(pH=ph,ionic=150,temp=298)))
 geom=build_surface_geometry(v1,nodes);keys=list(geom["coords"])
 # Candidate geometric neighborhoods have exactly the same radius and face condition as ranked model patches.
 geometry_pool=[{index[k] for k in keys if ca_distance(center,k,geom)<=8 and same_face(center,k,geom) and k in index} for center in keys]
 size=np.array([len(x) for x in geometry_pool])
 target=[len(x) for x in observed_patch_members[:5]]
 if len(target)<5:print("WARN_FEWER_THAN_FIVE",cid,len(target),flush=True)
 rng=random.Random(SEED+caseidx*1009)
 nulls={(model,k,metric):[] for model in ("rank_random","geometry_matched") for k in ks for metric,_ in metrics}
 eligible=range(len(observed_patch_members))
 for rep in range(N):
  shuffled=rng.sample(list(eligible),len(observed_patch_members))
  reordered=[observed_patch_members[i] for i in shuffled[:5]]
  used=set();random_geo=[]
  for t in target:
   best=min(abs(int(x)-t) for idx,x in enumerate(size) if idx not in used)
   ids=[idx for idx,x in enumerate(size) if idx not in used and abs(int(x)-t)==best]
   pick=rng.choice(ids);used.add(pick);random_geo.append(geometry_pool[pick])
  for model,arr in (("rank_random",reordered),("geometry_matched",random_geo)):
   for k in ks:
    union=set().union(*arr[:k]) if arr else set()
    for metric,threshold in metrics:
     nulls[(model,k,metric)].append(score(union,truth,xyz,nodes,equiv,threshold))
 for k in ks:
  m=set().union(*observed_patch_members[:k]) if observed_patch_members else set()
  for metric,threshold in metrics:
   value=score(m,truth,xyz,nodes,equiv,threshold)
   obs_rows.append({"case":cid,"set":case["analysis_set"],"k":k,"metric":metric,"observed":value,"evidence":case["evidence_resolution"],"n_ground_truth":len(truth)})
   for model in ("rank_random","geometry_matched"):
    arr=np.array(nulls[(model,k,metric)])
    null_rows.append({"case":cid,"set":case["analysis_set"],"k":k,"metric":metric,"null_model":model,"n":N,"observed":value,"null_mean":float(np.mean(arr)),"null_median":float(np.median(arr)),"null_q05":float(np.quantile(arr,.05)),"null_q95":float(np.quantile(arr,.95)),"p_greater_equal":pvalue(arr,value),"effect_vs_null_mean":value-float(np.mean(arr)),"random_seed":SEED+caseidx*1009})
 case_qc.append({"case":cid,"pdb":pdb,"n_total_patches":len(patches),"n_surface_residues":len(keys),"n_gt":len(truth),"frozen_near8":saved,"computed_near8":val,"n_geometric_candidates":len(geometry_pool),"status":"MATCH"})
 print("CASE_COMPLETE",cid,"patches",len(patches),"Top5near8",round(val,5),flush=True)
rowsave("observed_metrics.csv",obs_rows)
rowsave("null_results.csv",null_rows)
rowsave("case_qc.csv",case_qc)
# Posthoc descriptive summary, NOT an independent-sample significance test. HCAII repeated pH handled as paired case.
primary=[x for x in null_rows if x["set"]=="primary"]
summary={"model":"2.4.0-selected-parameters","n_permutations_per_null":N,"cases":len(CASES),"frozen_match":all(x["status"]=="MATCH" for x in case_qc),"n_metric_null_rows":len(null_rows),"no_unadjusted_p_below_0_05":sum(x["p_greater_equal"]<0.05 for x in primary),"scope":"exploratory casewise null tests; no pooled independent protein p-value; geometry matching by patch member count only and surface face not full accessibility matching","note":"No claim of chemistry-specific predictive superiority follows from geometry or ranking null alone. Previous chemistry-label null requires methodological redesign."}
(OUT/"audit_summary.json").write_text(json.dumps(summary,indent=2))
print("STAGE3_COMPLETE",json.dumps(summary),flush=True)
