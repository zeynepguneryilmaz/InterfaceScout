"""Direct provenance from RNase-A adsorption paper 10.1021/la503854a.
This is a **functional application / AAL-MS structural-response** assessment,
not a ground-truth direct adsorption contact benchmark.
Frozen native 7RSA structure: 124 aa; source paper did not supply adsorbed 3-D coordinates.
Surface channels predeclared: glass=anionic, HDPE=hydrophobic;
PMMA=hbond_acceptor CARBONYL PROXY (surface is more complex, exploratory).
No activity values are fed to InterfaceScout or chemistry channel selection.
"""
from pathlib import Path
import csv,io,json,urllib.request
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"validation/independent_md_functional/rnase_activity";OUT.mkdir(parents=True,exist_ok=True)
import sys
sys.path.insert(0,str(ROOT/"backend"))
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes
from interface_engine import analyze_all_maps
import core
DOI="10.1021/la503854a";SRC="https://pmc.ncbi.nlm.nih.gov/articles/PMC4270395/"
# These numeric values and ordinal structural-response label sets are transcribed
# from authenticated Europe-PMC fulltext Tables 3, 4 and 5.
activity=[
 ("glass",.03,38,8,.08),("HDPE",.03,43,6,.10),("PMMA",.03,33,5,.08),
 ("glass",1.0,39,9,.16),("HDPE",1.0,35,8,.17),("PMMA",1.0,45,9,.16)]
# 'strong protected' profile <= -0.72 ; 'moderate protected' -0.72 < profile < -0.12
p04={
 ("glass",.03):(["K1","R39","R85","E111"],["K7","R10","K31","D38","D53","K98"]),
 ("HDPE",.03):(["D53","K61","K66","Y76","K91","K98"],["K104","R39","K41","E49","Y73","Y92"]),
 ("PMMA",.03):(["E111"],["D38","R39","K66","Y76","K91","Y92","K98"]),
 ("glass",1.0):([],["K1","K7","E9","D38","K61","K66","Y76","K91","K104","H105","H119"]),
 ("HDPE",1.0):(["K1","K7","K31","K61","K66","K91","K98"],["K37","D38","Y76","K104","H105","H119"]),
 ("PMMA",1.0):(["K1","K66","K31","R85","K7"],["R33","K37","D38","R39","K41","D53","K61","Y76","K91","Y92","E111","H105","H119"])}
CHEM={"glass":"anionic","HDPE":"hydrophobic","PMMA":"hbond_acceptor"}
PDB="7RSA";CHAIN="A";PH=7.4
raw=urllib.request.urlopen("https://files.rcsb.org/download/7RSA.pdb",timeout=50).read().decode("utf8")
prepared,_=prepare_pdb_text(raw,chain=CHAIN);(OUT/"7RSA_prepared.pdb").write_text(prepared)
nodes=extract_ca_nodes(prepared);resseq={int(x["res_seq"]):x for x in nodes}
assert len(nodes)==124,(len(nodes),[x["res_seq"] for x in nodes[:3]])
def indices(labels):
 inds=[]
 for q in labels:
  res=int(q[1:]);name={"K":"LYS","R":"ARG","E":"GLU","D":"ASP","Y":"TYR","H":"HIS"}.get(q[0])
  assert res in resseq and resseq[res]["res_name"]==name,("AAL-MS mismatch",q,resseq.get(res))
  inds.append(res)
 return inds
for a,b in p04.values():indices(a);indices(b)
index={x["key"]:i for i,x in enumerate(nodes)}
xyz=np.asarray([x["coord"] for x in nodes])
native=core.analyze(core.AnalyzeRequest(pdb_text=prepared,env=core.EnvParams(pH=PH,ionic=150,temp=298)))
surf={x["key"]:x for x in native["surface_residues"]}
maps=analyze_all_maps(pdb_text=prepared,chain=None,pH=PH)["maps"]
active=[12,41,119] # active-site catalytic + phosphate anchoring, reference RCSB/RNase A
for n in active:assert n in resseq
def min_d(sel,site):
 if not sel:return None
 aa=np.asarray([index[x] for x in sel if x in index])
 jj=np.asarray([index[resseq[n]["key"]] for n in site])
 return float(np.linalg.norm(xyz[aa,None,:]-xyz[None,jj,:],axis=-1).min())
def recall(sel,truth,d):
 if not sel or not truth:return 0.
 a=np.asarray([index[x] for x in sel if x in index]);b=np.asarray([index[resseq[n]["key"]] for n in truth])
 if d==0:return float(np.isin(b,a).mean())
 return float((np.linalg.norm(xyz[b,None,:]-xyz[None,a,:],axis=-1)<=d+1e-9).any(axis=1).mean())
records=[];residues=[]
for material,conc,val,ci,coverage in activity:
 ch=CHEM[material];patches=maps[ch]["patches"]
 strong,mod=p04[(material,conc)];gt_strong=indices(strong);gt_mod=indices(mod);gt_all=sorted(set(gt_strong+gt_mod))
 for k in (1,3,5):
  sel=set(x for p in patches[:k] for x in p["members"])
  area=sum(surf[x]["sidechain_sasa"] for x in sel if x in surf)
  records.append({"protein":"Bovine pancreatic RNase A","pdb_id":PDB,"surface":material,
   "model_channel":ch,"channel_warning":"PMMA_IS_EXPLORATORY_PROXY" if material=="PMMA" else "SURFACE_MAPPING_APPROXIMATE",
   "solution_concentration_mg_ml":conc,"activity_percent":val,"activity_CI95_reported_halfwidth":ci,
   "surface_coverage_ug_cm2":coverage,"k":k,"n_selected":len(sel),"selected_SASA_A2":area,
   "n_protected_strong":len(gt_strong),"n_protected_all":len(gt_all),
   "protected_strong_direct_recall":recall(sel,gt_strong,0),
   "protected_all_direct_recall":recall(sel,gt_all,0),
   "protected_all_near5_recall":recall(sel,gt_all,5),
   "protected_all_near8_recall":recall(sel,gt_all,8),
   "active_site_min_CA_A":min_d(sel,active),
   "active_site_direct_overlap":int(any(n in {resseq[x]["res_seq"] for x in []} for n in [])) if False else int(any(resseq[n]["key"] in sel for n in active)),
   "DOI":DOI,"source_table":"3,4,5",
   "interpretation":"AAL-MS protected surface region: contact-OR-structural response, NOT confirmed physical contacts"})
 for key in sorted(set(gt_all+active)):
  residues.append({"surface":material,"concentration_mg_ml":conc,"resnum":key,"resname":resseq[key]["res_name"],
    "strongly_protected":int(key in gt_strong),"moderately_protected":int(key in gt_mod),
    "catalytic_or_site_anchor":int(key in active),"source_doi":DOI})
def save(name,rows):
 with (OUT/name).open("w",newline="",encoding="utf8") as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
save("rnase_activity_functional_site_screen.csv",records)
save("rnase_AALMS_residue_ground_truth.csv",residues)
# Explicit check: IS predicted ranking is concentration-INVARIANT, although
# measured bioactivity surface ranking changes with exposure concentration.
source_rank={}
for conc in (.03,1.0):
 candidates=[(s,v) for s,c,v,_,_ in activity if c==conc]
 source_rank[str(conc)]=sorted(candidates,key=lambda x:-x[1])
summary={"doi":DOI,"protein":"RNase A","pdb":"7RSA","n_conditions":len(activity),
 "n_unique_proteins":1,"n_surfaces":3,"n_concentrations":2,
 "source_activity_ranks":source_rank,"sources":"Source Table 3 activity and Tables 4-5 AAL-MS residue protection",
 "accuracy_caveat":"No direct interface residues proven; surface protection can be protein-protein effects",
 "concentration_caveat":"v1 lacks concentration input and produces the SAME structural ranking at both concentrations, while activity ordering changes",
 "surface_caveat":"Anionic glass / hydrophobic HDPE / hbond acceptor PMMA are simplified proxies. No calibrated physical model of silica silanols or PMMA groups",
 "outcome":"Quantitative functional screen completed; no valid claim of activity prediction on these six observations"}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("RNASE_FUNCTIONAL_COMPLETE",json.dumps(summary),flush=True)
