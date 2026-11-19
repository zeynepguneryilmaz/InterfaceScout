#!/usr/bin/env python3
"""Independent post-freeze structure-only numerical/scale robustness check.

16 unique PDBs drawn from the 44-condition application panel and not among the
8 adsorption-label-free parameter development structures. No label is used;
model parameters remain 200 SASA points and 6/9 Å radii.
Compare 200 vs 1000-point scRSA and 6/9 vs nine alternate radius pairs.
"""
import sys, csv, json, gzip, time, urllib.request, traceback
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
import core
from prepare import prepare_pdb_text
from model_settings import RADIUS_SELECTION_CANDIDATES_A
ROOT=Path("is44_independent_parameter_holdout");ROOT.mkdir(exist_ok=True)
PDBS=["1LYZ","1MBN","1HRC","1UBQ","2GB1","2OED","1JNJ","2PTN","2CBA","3ECA","1A5R","2Y0G","2HHB","4F5S","2HAV","7KWI"]
PAIRS=list(RADIUS_SELECTION_CANDIDATES_A)
REFPAIR=(6.0,9.0)
def jacc(a,b):
 a,b=set(a),set(b)
 return len(a&b)/len(a|b) if a|b else None
def hotspot_dict(rows,pairs):
 exposed=[r for r in rows if r["surface_exposed"]]
 if not exposed:return {}
 D=core.build_distances(exposed)
 nm={float(radius):(D<=float(radius)).astype(np.float32) for pair in pairs for radius in pair}
 results={}
 for channel,cfg in core.CHEMISTRIES.items():
  favored=cfg["favorable"]
  local=np.asarray([r["scrsa"]*core.state_availability(r["res_name"],favored[r["res_name"]][2],7.4) if r["res_name"] in favored else 0.0 for r in exposed],dtype=np.float32)
  if np.max(local)<=0:continue
  for inner,outer in pairs:
   v1=nm[inner]@local;v2=nm[outer]@local
   n1=v1/(max(v1) if max(v1)>0 else 1);n2=v2/(max(v2) if max(v2)>0 else 1)
   h=100*np.minimum(n1,n2);g=100*np.sqrt(n1*n2)
   pos=[i for i in range(len(exposed)) if h[i]>0]
   pos=sorted(pos,key=lambda i:(-round(float(h[i]),3),-round(float(g[i]),3),str(exposed[i]["key"])))
   results[(channel,inner,outer)]=[exposed[i]["key"] for i in pos[:10]]
 return results
manifest=[];sasa=[];pairrows=[]
for pid in PDBS:
 tic=time.time()
 try:
  req=urllib.request.Request(f"https://files.rcsb.org/download/{pid}.pdb",headers={"User-Agent":"InterfaceScout independent numerical validation"})
  with urllib.request.urlopen(req,timeout=60) as response: raw=response.read().decode("utf-8","replace")
  prepared,prep=prepare_pdb_text(raw,chain="A")
  core.SASA_POINTS=200;rows200=core._surface_residues(prepared,pH=7.4)
  core.SASA_POINTS=1000;rows1000=core._surface_residues(prepared,pH=7.4)
  core.SASA_POINTS=200
  a={r["key"]:r for r in rows200};b={r["key"]:r for r in rows1000}
  keys=sorted(set(a)&set(b))
  if len(keys)<2:raise Exception("too few matched residues")
  mae=float(np.mean([abs(float(a[k]["scrsa_raw"])-float(b[k]["scrsa_raw"])) for k in keys]))
  surface200={k for k in a if a[k]["surface_exposed"]}
  surface1000={k for k in b if b[k]["surface_exposed"]}
  h200=hotspot_dict(rows200,PAIRS);h1000=hotspot_dict(rows1000,[REFPAIR])
  overlap=[]
  for channel in core.CHEMISTRIES:
   key=(channel,*REFPAIR)
   if key in h200 and key in h1000:
    v=jacc(h200[key],h1000[key])
    if v is not None:overlap.append(v)
  sasa.append({"PDB":pid,"n_residues":len(keys),"surface200":len(surface200),"surface1000":len(surface1000),"scRSA_MAE":round(mae,6),
    "surface_jaccard":round(jacc(surface200,surface1000),6),"median_top10_Jaccard":round(float(np.median(overlap)),6) if overlap else None,"n_chem":len(overlap)})
  perpair={}
  for pair in PAIRS:
   x=[]
   for channel in core.CHEMISTRIES:
    ka=(channel,*REFPAIR);kb=(channel,*pair)
    if ka in h200 and kb in h200:
     v=jacc(h200[ka],h200[kb])
     if v is not None:x.append(v)
   if x:
    perpair[pair]=float(np.mean(x))
    pairrows.append({"PDB":pid,"candidate_inner_A":pair[0],"candidate_outer_A":pair[1],"mean_top10_Jaccard_to_frozen_6_9":round(float(np.mean(x)),6),
        "median_top10_Jaccard_to_frozen_6_9":round(float(np.median(x)),6),"n_channels":len(x)})
  manifest.append({"PDB":pid,"status":"ok","n_residues":len(keys),"elapsed_s":round(time.time()-tic,2)})
  print(pid,"ok",len(keys),"MAE",mae,"SASA Jacc",jacc(surface200,surface1000),"top10",np.median(overlap) if overlap else None,flush=True)
 except Exception as e:
  manifest.append({"PDB":pid,"status":"error","error":str(e)[:300],"elapsed_s":round(time.time()-tic,2)})
  traceback.print_exc()
for name,rows in (("sasa_results.csv",sasa),("radius_altpair_results.csv",pairrows),("run_manifest.csv",manifest)):
 if rows:
  with (ROOT/name).open("w",newline="") as f:
   wr=csv.DictWriter(f,fieldnames=sorted({k for r in rows for k in r}));wr.writeheader();wr.writerows(rows)
summary={"n_requested":len(PDBS),"n_completed":len(sasa),"n_failed":len(PDBS)-len(sasa),"development_PDBs":["1CRN","1R69","1SHG","2PPN","4AKE","1OMP","1TIM","5CSC"],
 "mean_scrsa_MAE":round(float(np.mean([r["scRSA_MAE"] for r in sasa])),6) if sasa else None,
 "median_scrsa_MAE":round(float(np.median([r["scRSA_MAE"] for r in sasa])),6) if sasa else None,
 "median_exposure_jaccard":round(float(np.median([r["surface_jaccard"] for r in sasa])),6) if sasa else None,
 "worst_exposure_jaccard":round(float(np.min([r["surface_jaccard"] for r in sasa])),6) if sasa else None,
 "median_top10_jaccard":round(float(np.median([r["median_top10_Jaccard"] for r in sasa if r["median_top10_Jaccard"] is not None])),6) if sasa else None,
 "qualifier":"Holdout structural numerical stability, not additional experimental adsorption benchmark or optimization of hyperparameters."}
(ROOT/"summary.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2),flush=True)
