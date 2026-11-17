"""Locked external MD contact benchmark: 2025 native hen egg-white lysozyme + SWCNT.
Original paper: Sci Rep 15,11593, DOI 10.1038/s41598-025-96435-3,
System-1 native lysozyme/SWCNT, PDB 1LYZ, 1 microsecond, 0.5nm contact cutoff.
All 14 residues are explicitly reported in its full text. No labels are used by predictions.
This is *agreement with a single MD study*, not experimental confirmation.
"""
from pathlib import Path
import csv,json,sys,urllib.request,random
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"validation"/"independent_md_functional"/"lysozyme_swc";OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/"backend"))
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes
from interface_engine import analyze_all_maps
import core
PDB="1LYZ";CHAIN="A";PH=7.0 # paper does not report pH here; reference reference-only pH
DOI="10.1038/s41598-025-96435-3"
source="https://www.nature.com/articles/s41598-025-96435-3"
residues={61:"ARG",62:"TRP",63:"TRP",68:"ARG",69:"THR",70:"PRO",73:"ARG",
103:"ASN",106:"ASN",107:"ALA",111:"TRP",112:"ARG",113:"ASN",116:"LYS"}
sites={"BS1":[61,62,63],"BS2":[68,69,70,73],"BS3":[103,106,107],"BS4":[111,112,113,116]}
data=urllib.request.urlopen("https://files.rcsb.org/download/1LYZ.pdb",timeout=45).read().decode("utf8")
prepared,meta=prepare_pdb_text(data,chain=CHAIN);(OUT/"1LYZ_prepared.pdb").write_text(prepared)
nodes=extract_ca_nodes(prepared)
lookup={int(x["res_seq"]):x for x in nodes if x["chain"]==CHAIN}
for n,aa in residues.items():
 if n not in lookup or lookup[n]["res_name"]!=aa:raise ValueError(f"MD article numbering does not match PDB 1LYZ at {n}, expected {aa}, got {lookup.get(n)}")
xyz=np.asarray([x["coord"] for x in nodes])
nindex={r["key"]:i for i,r in enumerate(nodes)}
gt=np.asarray([nindex[lookup[i]["key"]] for i in residues])
v1=core.analyze(core.AnalyzeRequest(pdb_text=prepared,env=core.EnvParams(pH=PH,ionic=150,temp=298)))
surface={r["key"]:r for r in v1["surface_residues"]};sa=sum(x["sidechain_sasa"] for x in surface.values())
maps=analyze_all_maps(pdb_text=prepared,chain=None,pH=PH)["maps"]
def score(sel,gt,d):
 ids=np.array([nindex[x] for x in sel if x in nindex],dtype=int)
 if not len(ids):return 0.
 if d==0:return float(np.isin(gt,ids).mean())
 return float((np.linalg.norm(xyz[gt,None,:]-xyz[None,ids,:],axis=2)<=d+1e-9).any(axis=1).mean())
def save(fn,rows):
 if not rows:return
 with (OUT/fn).open("w",newline="",encoding="utf8") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
# Freeze the relevant chemistry channels independent of outcome:
# carbon graphitic -> pi_carbon; alternative nonpolar channel -> hydrophobic.
# Evaluate both, separately, not selecting the better channel post hoc.
primary=["pi_carbon","hydrophobic"];secondary=[x for x in maps if x not in primary]
rows=[];random_rows=[];site_rows=[]
for name in primary+sorted(secondary):
 patches=maps[name]["patches"]
 if not patches:continue
 m=[set(p["members"]) for p in patches]
 for k in (1,3,5):
  pred=set().union(*m[:k]);area=sum(surface[x]["sidechain_sasa"] for x in pred if x in surface)
  base={"chemistry":name,"analysis_tier":"primary_prespecified" if name in primary else "exploratory_secondary",
   "k":k,"n_patches":len(m),"n_selected_residues":len(pred),"sidechain_area_A2":area,
   "area_fraction":area/sa,"n_MD_contact_residues":14,"pH_model_reference":PH,
   "pH_source":"NOT_REPORTED_REFERENCE_DEFAULT","doi":DOI}
  for dist in (0,5,8):
   row=dict(base,distance_A=dist,recall=score(pred,gt,dist))
   rows.append(row)
  for label,ns in sites.items():
   idx=np.asarray([nindex[lookup[i]["key"]] for i in ns])
   site_rows.append(dict(base,site=label,n_site_residues=len(ns),direct_recall=score(pred,idx,0),near5_recall=score(pred,idx,5)))
 # rank-random baseline: same actual candidate patches and exact k.
 # Conditioning on selected area and selected count ±10%; explicit acceptance.
 rng=random.Random(20261117+sum(ord(c) for c in name))
 for k in (1,3,5):
  pred=set().union(*m[:k])
  target_area=sum(surface[x]["sidechain_sasa"] for x in pred if x in surface)
  target_n=len(pred)
  accepted=[];attempts=0
  for a in range(10000):
   attempts+=1
   inds=rng.sample(range(len(m)),min(k,len(m)));sel=set().union(*(m[i] for i in inds))
   area=sum(surface[x]["sidechain_sasa"] for x in sel if x in surface)
   if abs(len(sel)-target_n)/max(1,target_n)>.1:continue
   if abs(area-target_area)/max(1,target_area)>.1:continue
   accepted.append(sel)
   if len(accepted)>=2000:break
  for dist in (0,5,8):
   obs=score(pred,gt,dist)
   vals=np.array([score(x,gt,dist) for x in accepted])
   random_rows.append({"chemistry":name,"k":k,"distance_A":dist,"n_matched_random":len(vals),"attempts":attempts,
    "observed_recall":obs,"null_mean":float(vals.mean()) if len(vals) else None,
    "null_q025":float(np.quantile(vals,.025)) if len(vals) else None,
    "null_q975":float(np.quantile(vals,.975)) if len(vals) else None,
    "one_sided_empirical_p":float((1+sum(vals>=obs-1e-12))/(len(vals)+1)) if len(vals) else None,
    "inferential_status":"EXPLORATORY_ONLY" if len(vals)>=100 else "INSUFFICIENT_MATCHED_RANDOM_DRAWS"})
 print("CHANNEL",name,"patches",len(patches),flush=True)
save("chemistry_topk_MD_recall.csv",rows);save("four_site_coverage.csv",site_rows)
save("matched_rank_random_nulls.csv",random_rows)
save("source_MD_contacts.csv",[{"resnum":k,"residue":v,"site":next(z for z,ns in sites.items() if k in ns),"doi":DOI,"evidence":"MD_native_system1_contact_under_0.5nm","source_url":source} for k,v in residues.items()])
summary={"source":source,"doi":DOI,"pdb":"1LYZ","n_ground_truth":len(residues),"all_numbered_residues_identity_match":True,
 "prespecified_channels":primary,"channels_analyzed":len({r["chemistry"] for r in rows}),
 "pH_warning":"Paper pH not established from text; using 7.0 as reference value is exploratory, must rerun if source conditions found.",
 "md_tier":"external computational corroboration only, NOT independent experimental validation",
 "no_contact_trajectory":"Only author-reported 14 contact residue labels, not a time-resolved MD trajectory",
 "potential_overfitting":"Model and candidate selection have been iteratively studied on other systems; this is an exploratory new-case assessment"}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("MD_LYSOZYME_COMPLETE",json.dumps(summary),flush=True)
