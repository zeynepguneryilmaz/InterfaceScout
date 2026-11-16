"""Chemistry-off ablation against frozen predictions; nonparametric paired descriptive comparison."""
import csv,json,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
import core
from prepare import prepare_pdb_text
from interface_engine import analyze_all_maps
from coarse_patch import build_coarse_patches
from geometry import extract_ca_nodes
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"validation"/"chemistry_ablation_results";OUT.mkdir(parents=True,exist_ok=True)
MANIFEST=json.loads((ROOT/"validation"/"benchmark_manifest.json").read_text())["cases"]
FROZEN={r["case_id"]:r for r in csv.DictReader((ROOT/"validation"/"benchmark_case_results.csv").open(encoding="utf-8-sig"))}
SURFACE={"silica":"anionic","paa_fe3o4":"anionic","citrate_au":"anionic","mpa_au":"anionic","alumoh":"cationic","polystyrene":"hydrophobic"}
# Exposure-only ablation: all 20 standard amino acids receive identical chemistry membership.
# Do not alter input coordinates, scRSA, 6/9 A density, patch radius or Pareto ranking rules.
def all_chemistry_control(v1,nodes,ph,original_key):
  key="_exposure_only_audit"
  original=core.CHEMISTRIES[original_key]
  fav={rn:(0.0,"exposure-only chemistry ablation",None) for rn in (
    "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")}
  core.CHEMISTRIES[key]={"label":"Exposure only, all residues compatible","surface_group":"none","description":"All standard amino acids compatible, no pH state chemistry","favorable":fav,"repulsive":{}}
  try:
   channel=core.chemistry_map(v1["surface_residues"],core.build_distances(v1["surface_residues"]),key,ph)
   x=dict(v1,chemistries={key:channel})
   return build_coarse_patches(v1_result=x,chemistry=key,ca_nodes=nodes)
  finally:
   del core.CHEMISTRIES[key]
def truth_indices(case,nodes):
 out=[]
 for i,n in enumerate(nodes):
  if n["chain"] not in case["chain"].split(","):continue
  k=int(n["res_seq"]); typ=case["gt_kind"]
  if typ=="exact" and k in case.get("gt_exact",[]):out.append(i)
  elif typ in ("ranges","equivalent_chain_ranges") and any(a<=k<=b for a,b in case.get("gt_ranges",[])):out.append(i)
 return out
def recall(patches,index,truth,xyz,nodes,equiv,k,t):
 sel={index[v] for p in patches[:k] for v in p["members"] if v in index}
 if not sel:return 0.
 if t==0: hit=np.isin(truth,list(sel))
 else:
  ds=np.linalg.norm(xyz[np.array(truth),None,:]-xyz[np.array(sorted(sel))[None,:,:],:],axis=-1)
  hit=np.min(ds,axis=1)<=t+1.e-9
 if equiv:
  chains=np.array([nodes[i]["chain"] for i in truth])
  return max(float(hit[chains==ch].mean()) for ch in set(chains))
 return float(hit.mean())
rows=[]
for c in MANIFEST:
 cid=c["id"];ph=c["pH"];pdb=c["pdb_id"]
 files=sorted((ROOT/"examples"/pdb).rglob("*_Bfactor.pdb"))
 if pdb=="2CBA":files=[p for p in files if "pH_"+str(ph) in str(p)]
 if len(files)!=1:raise RuntimeError(f"{cid} missing frozen PDB example")
 prepared,_=prepare_pdb_text(files[0].read_text(),chain=c["chain"])
 nodes=extract_ca_nodes(prepared)
 xyz=np.array([n["coord"] for n in nodes]);index={n["key"]:i for i,n in enumerate(nodes)}
 truth=truth_indices(c,nodes)
 if not truth:raise RuntimeError(cid+" no mapped evidence")
 v1=core.analyze(core.AnalyzeRequest(pdb_text=prepared,env=core.EnvParams(pH=ph,ionic=150,temp=298)))
 true_map=analyze_all_maps(pdb_text=prepared,chain=None,pH=ph)["maps"][SURFACE[c["surface"]]]["patches"]
 nul_map=all_chemistry_control(v1,nodes,ph,SURFACE[c["surface"]])
 equiv=c["gt_kind"]=="equivalent_chain_ranges"
 observed=recall(true_map,index,truth,xyz,nodes,equiv,5,8)
 assert abs(observed-float(FROZEN[cid]["top5_near_8A_recall"]))<0.015,(cid,observed)
 for k in (1,3,5):
  for name,threshold in (("direct",0),("near5",5),("near8",8)):
   a=recall(true_map,index,truth,xyz,nodes,equiv,k,threshold)
   b=recall(nul_map,index,truth,xyz,nodes,equiv,k,threshold)
   rows.append({"case":cid,"set":c["analysis_set"],"protein":c["protein"],"chemistry":SURFACE[c["surface"]],"top_k":k,"metric":name,"observed":a,"exposure_only":b,"delta":a-b,"original_n_patches":len(true_map),"exposure_only_n_patches":len(nul_map),"n_experimental_residues":len(truth)})
 print("ABLATION_CASE",cid,"Top5 near8",round(observed,5),"exposure_only",round(recall(nul_map,index,truth,xyz,nodes,equiv,5,8),5),flush=True)
with (OUT/"chemistry_off_ablation.csv").open("w",newline="") as f:
 w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
primary=[r for r in rows if r["set"]=="primary"]
summary={"n_conditions":len(MANIFEST),"primary_conditions":sum(c["analysis_set"]=="primary" for c in MANIFEST),"n_comparisons":len(rows),"frozen_benchmark_match":True,"nonindependent_hcaii_pH":True,"analysis":"deterministic one-time chemistry-off ablation; not a permutation test"}
for k in (1,3,5):
 for m in ("direct","near5","near8"):
  rr=[r for r in primary if r["top_k"]==k and r["metric"]==m]
  summary[f"top{k}_{m}"]={"mean_chemistry":float(np.mean([r["observed"] for r in rr])),"mean_exposure_only":float(np.mean([r["exposure_only"] for r in rr])),"mean_delta":float(np.mean([r["delta"] for r in rr])),"positive":sum(r["delta"]>1e-12 for r in rr),"negative":sum(r["delta"]< -1e-12 for r in rr),"ties":sum(abs(r["delta"])<=1e-12 for r in rr)}
(OUT/"ablation_summary.json").write_text(json.dumps(summary,indent=2))
print("ABLATION_COMPLETE",json.dumps(summary),flush=True)
