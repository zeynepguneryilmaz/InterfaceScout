"""Frozen InterfaceScout v1 benchmark against external published MD contact lists.
Source: Sun et al., PLOS ONE 2014, DOI 10.1371/journal.pone.0107696.
3 proteins, 9 reported material conditions, ONE independent study family.
The contact lists are directly transcribed from article body paragraphs around Fig 5,6,7,14.
No source-derived labels are used to select ranking or chemistry.

Comparators:
 - Original frozen v1, same patches
 - Mean scRSA (surface-accessibility-only)
 - Compatibility fraction (chemistry-only)
 - Multiscale patch coherence (chemical-density/persistence ablation)
 - Random ranking conditioned on SAME residue-count & sidechain-SASA Top-k budget (+/-10%)
Reported "unmatched" rankings are descriptive and NOT equal-budget head-to-head comparisons.
Same patch candidate pool for deterministic baselines evaluates RANKING, not patch-generation.
"""
from __future__ import annotations
import csv,json,sys,urllib.request,random
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
O=ROOT/"validation"/"frozen_v1_external_silica_md";O.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/"backend"))
from prepare import prepare_pdb_text
from interface_engine import analyze_all_maps
from geometry import extract_ca_nodes
import core
DOI="10.1371/journal.pone.0107696"
URL="https://journals.plos.org/plosone/article?id="+DOI
# label (type,residue_number), input model and mapping predeclared BEFORE prediction.
# Coated-surface mapping is one-component PROXY: -CH3->hydrophobic, -COOH->anionic,
# -OH->hbond_acceptor; flat/uncoated silica->anionic at pH 7.4 reference.
# NH2 group from article lacks explicitly enumerated contacts -> excluded.
CASES=[
 ("cytc_silica4","3NWV","Cytochrome c",4,"bare_silica","anionic","PHE36 GLY37 ARG38 LYS39 PRO44 LYS55 ILE58 GLU61 ASP62 ASN103"),
 ("cytc_silica11","3NWV","Cytochrome c",11,"bare_silica","anionic","ALA92 ILE95 ALA96 LYS99 LYS100"),
 ("rnase_silica4","7RSA","RNase A",4,"bare_silica","anionic","SER15 SER16 THR17 SER18 SER21 SER50 LEU51 ALA52 ASP53 GLN55"),
 ("rnase_silica11","7RSA","RNase A",11,"bare_silica","anionic","THR3 GLU9 MET13 SER15 SER16 LEU51 ALA52 GLN55 PRO114"),
 ("lys_silica4","3LN2","Lysozyme",4,"bare_silica","anionic","CYS6 GLU7 ARG10 ARG14 LEU15 CYS128 GLY129 VAL130"),
 ("lys_silica11","3LN2","Lysozyme",11,"bare_silica","anionic","GLU7 ARG10 THR11 ARG14 LEU15 TYR20 HIS78 GLN86 ASN88 ALA90 ASP91 ALA94 LYS97"),
 ("cytc_methyl2","3NWV","Cytochrome c",2,"methyl","hydrophobic","GLY1 ASP2 VAL3 GLU4 ALA96 LYS99 LYS100"),
 ("cytc_carboxyl2","3NWV","Cytochrome c",2,"carboxyl","anionic","GLY1 ASP2 GLU4 GLU61 LYS99 LYS100 ASN103"),
 ("cytc_hydroxyl2","3NWV","Cytochrome c",2,"hydroxyl","hbond_acceptor","GLY1 ASP2 GLU4 ARG38 GLU61 LYS99 LYS100 ASN103")]
PH=7.4
AA={"ALA":"ALA","ARG":"ARG","ASN":"ASN","ASP":"ASP","CYS":"CYS","GLN":"GLN","GLU":"GLU","GLY":"GLY","HIS":"HIS","ILE":"ILE","LEU":"LEU","LYS":"LYS","MET":"MET","PHE":"PHE","PRO":"PRO","SER":"SER","THR":"THR","TRP":"TRP","TYR":"TYR","VAL":"VAL"}
def parse(text):
 result=[]
 for token in text.split():
  name=token[:3];num=int(token[3:])
  if name not in AA:raise ValueError("Unsupported residue label: "+token)
  result.append((name,num))
 return result
def load(pdb):
 data=urllib.request.urlopen("https://files.rcsb.org/download/"+pdb+".pdb",timeout=50).read().decode("utf8")
 prepared,_=prepare_pdb_text(data,chain="A")
 (O/(pdb+"_prepared.pdb")).write_text(prepared)
 ns=extract_ca_nodes(prepared)
 return prepared,ns
def recall(sel,truth,xyz,key_index,dist):
 if not sel:return 0.
 q=np.array([key_index[v] for v in sel if v in key_index],int)
 t=np.array(truth,int)
 if not len(q):return 0.
 if dist==0:return float(np.isin(t,q).mean())
 return float((np.linalg.norm(xyz[t,None,:]-xyz[None,q,:],axis=-1)<=dist+1e-9).any(axis=1).mean())
def write(name,rows):
 if not rows:return
 with (O/name).open("w",newline="",encoding="utf8") as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
records=[];qc=[];nulls=[];labels=[];case_details=[]
cache={}
for cid,pdb,protein,diameter,material,channel,contact in CASES:
 try:
  if pdb not in cache:cache[pdb]=load(pdb)
  prepared,nodes=cache[pdb]
  by_num={int(r["res_seq"]):r for r in nodes if r["chain"]=="A"}
  listed=parse(contact)
  wrong=[]
  for aa,n in listed:
   rr=by_num.get(n)
   if rr is None or rr["res_name"]!=aa:wrong.append((aa,n,rr["res_name"] if rr else None))
  if wrong:
   qc.append({"case":cid,"pdb":pdb,"status":"EXCLUDED_RESIDUE_SEQUENCE_MISMATCH","details":str(wrong)})
   print("EXCLUDED_SEQUENCE_MISMATCH",cid,wrong,flush=True);continue
  idx={r["key"]:i for i,r in enumerate(nodes)}
  truth=[idx[by_num[n]["key"]] for aa,n in listed]
  xyz=np.asarray([r["coord"] for r in nodes],float)
  # Every surface chemical mapping is declared above; no best-of-11 channel reporting.
  v1=core.analyze(core.AnalyzeRequest(pdb_text=prepared,env=core.EnvParams(pH=PH,ionic=150,temp=298)))
  exposed={r["key"]:r for r in v1["surface_residues"]}
  mp=analyze_all_maps(pdb_text=prepared,chain=None,pH=PH)["maps"].get(channel)
  patches=mp["patches"] if mp else []
  if not patches:qc.append({"case":cid,"pdb":pdb,"status":"EXCLUDED_EMPTY_PATCH_CHANNEL","details":channel});continue
  # Frozen v1 all patch order and membership unchanged.
  ordered={
   "IS_v1":list(patches),
   "Accessibility_only":sorted(patches,key=lambda p:(-p["mean_accessibility"],str(p["center_key"]))),
   "Chemistry_only":sorted(patches,key=lambda p:(-p["chemistry_support"],str(p["center_key"]))),
   "Coherence_only":sorted(patches,key=lambda p:(-p["patch_coherence"],str(p["center_key"])))
  }
  side={k:float(r["sidechain_sasa"]) for k,r in exposed.items()}
  def members(arr,k):return set(x for patch in arr[:k] for x in patch["members"])
  def coverage(s):return sum(side.get(x,0) for x in s)
  for aa,n in listed:
   labels.append({"case":cid,"protein":protein,"pdb":pdb,"surface":material,"diameter_nm":diameter,
      "channel":channel,"residue_num":n,"residue_name":aa,"DOI":DOI,
      "evidence_type":"SIMULATED_CONTACT_LIST_IN_PUBLISHER_ARTICLE",
      "article_url":URL,"experimentally_validated":0})
  for k in (1,3,5):
   original=members(patches,k);area0=coverage(original);n0=len(original)
   for method,ranking in ordered.items():
    s=members(ranking,k);area=coverage(s)
    erra=abs(area-area0)/max(area0,1e-12);errn=abs(len(s)-n0)/max(n0,1)
    for dist in (0,5,8):
     records.append({"case":cid,"family":"Sun2014_ONE_STUDY","protein":protein,"pdb":pdb,"diameter_nm":diameter,"surface":material,
     "chemistry_channel":channel,"mapping_quality":"COARSE_PROXY_NO_SIZE_OR_CURVATURE",
     "k":k,"method":method,"distance_A":dist,"recall":recall(s,truth,xyz,idx,dist),
     "n_contact_labels":len(truth),"n_selected":len(s),"sidechain_sasa_A2":area,
     "IS_v1_sasa_A2":area0,"IS_v1_n_selected":n0,
     "match_A2_10pct":int(erra<=.1),"match_count_10pct":int(errn<=.1),
     "both_matched_10pct":int(erra<=.1 and errn<=.1),
     "DOI":DOI})
   # Rank-random within SAME actual patches with both area AND count matched; failed draws flagged.
   rng=random.Random(20261117+sum(ord(t) for t in cid)*127+k)
   pool=[set(p["members"]) for p in patches]
   accepted=[];attempts=0
   for _ in range(25000):
    if len(accepted)>=2000:break
    attempts+=1
    pids=rng.sample(range(len(pool)),min(k,len(pool)))
    s=set().union(*(pool[z] for z in pids))
    a=coverage(s)
    if abs(len(s)-n0)/max(n0,1)>.1 or abs(a-area0)/max(area0,1e-12)>.1:continue
    accepted.append(s)
   for dist in (0,5,8):
    scores=np.array([recall(s,truth,xyz,idx,dist) for s in accepted],float)
    observed=recall(original,truth,xyz,idx,dist)
    nulls.append({"case":cid,"protein":protein,"k":k,"distance_A":dist,"attempts":attempts,
    "accepted":len(accepted),"observed":observed,
    "mean_null":float(scores.mean()) if len(scores) else None,
    "null_q025":float(np.quantile(scores,.025)) if len(scores) else None,
    "null_q975":float(np.quantile(scores,.975)) if len(scores) else None,
    "empirical_upper_tail_p":float((1+sum(scores>=observed-1e-12))/(1+len(scores))) if len(scores) else None,
    "interpretation":"NO_INFERENCE_FEWER_THAN_100_MATCHED_DRAWS" if len(scores)<100 else "DESCRIPTIVE_NOT_INDEPENDENT_REPLICATES"})
  qc.append({"case":cid,"pdb":pdb,"status":"COMPUTED","details":f"{len(listed)} verified labels, {len(patches)} patches, channel={channel}"})
  print("COMPUTED_CASE",cid,"labels",len(listed),"patches",len(patches),flush=True)
 except Exception as e:
  qc.append({"case":cid,"pdb":pdb,"status":"EXCLUDED_RUNTIME_OR_MAPPING_ERROR","details":str(e)})
  print("EXCLUDED_ERROR",cid,repr(e),flush=True)
write("casewise_recall.csv",records);write("source_contact_residues.csv",labels);write("case_audit.csv",qc);write("matched_rank_random_nulls.csv",nulls)
summary={"DOI":DOI,"source":URL,"n_predeclared_conditions":len(CASES),
 "n_computed_conditions":sum(x["status"]=="COMPUTED" for x in qc),
 "n_excluded":sum(x["status"]!="COMPUTED" for x in qc),
 "n_unique_proteins":len({x["protein"] for x in records}),
 "n_independent_study_families":1,
 "v1_version":"Unchanged public GitHub main backend",
 "prior_exposure":"Article selected after literature review; independent of scoring parameters but not prospectively preregistered",
 "primary_caveats":[
 "Simulated contact labels are NOT independent experiments",
 "Single PLOS article is one evidence family even though 9 conditions",
 "Bare silica, coated groups, and size are approximated by fixed chemistry channels",
 "4nm and 11nm silica use same predicted patch ranking; model has no curvature/diameter input",
 "Within-pool ranking baselines share v1-generated candidate patches, so they do not test relative candidate generation",
 "Secondary channel selection was NOT optimized from performance",
 "Only protein-side residue lists explicitly enumerated in article are labels; nearby helices not silently added"
 ],
 "mapping_source":"DOI publication body Figure 5/6/7/14",
 "control_coverage_rule":"TopK within ±10% both unique exposed sidechain SASA and unique residue count"}
(O/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print("EXTERNAL_SILICA_MD_COMPLETE",json.dumps(summary),flush=True)
