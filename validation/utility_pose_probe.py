"""Frozen v1 *operational* MD start-orientation selection + GB1 mutation triage.
This is a geometric, retrospective SCREEN, NOT trajectory MD / prospective utility proof.
No interaction energies, relaxation, equilibrium probability, or dynamics assumed.
All methods nominate k starting orientations from the *same* frozen IS patch candidates.
Primary endpoint: is an author-reported MD contact-face direction within 45 degrees
of at least one nominated candidate normal (protein centroid -> patch center)?
Secondary: 30/60 degrees, best angle, first-success rank, randomized orientation null.
One article family contains 9 silica conditions; never treat these as independent replicates.
"""
from __future__ import annotations
import csv,json,sys,urllib.request,random,ast
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"validation"/"utility_pose_probe"
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/"backend"))
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes
from interface_engine import analyze_all_maps
import core
# Import frozen source case labels without importing or executing the benchmark.
tree=ast.parse((ROOT/"validation"/"external_silica_md_9_conditions.py").read_text())
CASES=next(ast.literal_eval(node.value) for node in tree.body if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=="CASES" for t in node.targets))
CASES.append(("lys_SWCNT_native","1LYZ","Lysozyme",0,"SWCNT","pi_carbon",
              "ARG61 TRP62 TRP63 ARG68 THR69 PRO70 ARG73 ASN103 ASN106 ALA107 TRP111 ARG112 ASN113 LYS116"))
# Secondary sensitivity: SWCNT hydrophobic map alternative (same experimental evidence, NOT new case).
CASES.append(("lys_SWCNT_native_hydrophobic_sensitivity","1LYZ","Lysozyme",0,"SWCNT","hydrophobic",
              "ARG61 TRP62 TRP63 ARG68 THR69 PRO70 ARG73 ASN103 ASN106 ALA107 TRP111 ARG112 ASN113 LYS116"))
PLOS_DOI="10.1371/journal.pone.0107696";CNT_DOI="10.1038/s41598-025-96435-3"
def unit(v):
 a=np.linalg.norm(v)
 return v/a if a>1e-12 else np.zeros(3)
def angle(u,v):
 return float(np.rad2deg(np.arccos(np.clip(np.dot(u,v),-1,1))))
def save(fn,data):
 if not data:return
 with (OUT/fn).open("w",newline="",encoding="utf8") as f:
  fields=list(dict.fromkeys(key for row in data for key in row))
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
def make(pdb):
 data=urllib.request.urlopen("https://files.rcsb.org/download/"+pdb+".pdb",timeout=50).read().decode("utf8")
 prepared,_=prepare_pdb_text(data,chain="A")
 (OUT/(pdb+"_prepared.pdb")).write_text(prepared)
 ns=extract_ca_nodes(prepared)
 return prepared,ns
cache={};results=[];qc=[];matched_residues=[];rankings_all=[]
for cid,pdb,protein,diameter,material,channel,contacts in CASES:
 try:
  if pdb not in cache:cache[pdb]=make(pdb)
  prepared,nodes=cache[pdb]
  centroid=np.mean(np.array([r["coord"] for r in nodes]),axis=0)
  lookup={r["res_seq"]:r for r in nodes if r["chain"]=="A"}
  labels=[(x[:3],int(x[3:])) for x in contacts.split()]
  mismatch=[(a,n,lookup.get(n,{}).get("res_name")) for a,n in labels if lookup.get(n,{}).get("res_name")!=a]
  if mismatch:raise ValueError("Identity mismatch "+str(mismatch))
  # Source contact normal is mean of UNIT radial CA vectors; no interaction force field involved.
  reference_vectors=[unit(lookup[n]["coord"]-centroid) for a,n in labels]
  ref_resultant=np.mean(reference_vectors,axis=0)
  face_coherence=float(np.linalg.norm(ref_resultant))
  ref=unit(ref_resultant)
  maps=analyze_all_maps(pdb_text=prepared,chain=None,pH=7.4)["maps"]
  p=maps[channel]["patches"]
  coords={r["key"]:r["coord"] for r in nodes}
  if not p:raise ValueError("Empty predicted patch pool")
  normals=[unit(coords[row["center_key"]]-centroid) for row in p]
  angles=np.array([angle(u,ref) for u in normals])
  # Critical: random null is UNIFORM PATCH INDEX sampling, not uniform over 3D orientations;
  # controls use the same actual candidate patch pool and same k trials.
  ranks={
   "IS_v1":list(range(len(p))),
   "Accessibility_only":sorted(range(len(p)),key=lambda i:(-p[i]["mean_accessibility"],str(p[i]["center_key"]))),
   "Chemistry_only":sorted(range(len(p)),key=lambda i:(-p[i]["chemistry_support"],str(p[i]["center_key"]))),
   "Coherence_only":sorted(range(len(p)),key=lambda i:(-p[i]["patch_coherence"],str(p[i]["center_key"])))
  }
  # Other baselines are NOT independent patch generation. They test prioritization within frozen pool.
  for m,order in ranks.items():
   for i,patchid in enumerate(order[:min(len(order),5)]):
    rankings_all.append({"case":cid,"method":m,"rank":i+1,"center":p[patchid]["center_key"],
       "n_patch_residues":p[patchid]["n_members"],"angle_to_published_MD_contact_face_degrees":float(angles[patchid])})
   for k in (1,3,5):
    picks=order[:k]
    for threshold in (30,45,60):
     results.append({"case":cid,"source_family":"Sun2014_PLOS" if "silica" in cid or "cytc_" in cid else "LysozymeSWCNT2025",
      "protein":protein,"pdb":pdb,"surface":material,"diameter_nm":diameter,"channel":channel,
      "chemistry_mapping":"fixed_proxy_no_material_atomistic_model","method":m,"k":k,
      "threshold_degrees":threshold,
      "n_candidates_generated":len(p),"n_source_MD_contacts":len(labels),
      "reference_face_resultant_length":face_coherence,
      "face_ambiguity_flag":int(face_coherence<.5),
      "best_angle_degrees":float(min(angles[picks])),
      "orientation_hit":int(np.min(angles[picks])<=threshold),
      "first_success_rank":next((z+1 for z,idx in enumerate(order) if angles[idx]<=threshold),None),
      "doi":PLOS_DOI if material!="SWCNT" else CNT_DOI,
      "analysis_tier":"POSTHOC_GEOMETRIC_UTILITY_PROXY_NOT_MD_RUN"})
  # random equal k patch orientations, 5000 independent Monte Carlo draws
  rng=random.Random(20261117+sum((i+1)*ord(ch) for i,ch in enumerate(cid)))
  for k in (1,3,5):
   nmc=5000
   sample=[rng.sample(range(len(p)),min(k,len(p))) for _ in range(nmc)]
   for threshold in (30,45,60):
    hit=np.array([np.min(angles[pick])<=threshold for pick in sample],dtype=bool)
    obs=int(np.min(angles[ranks["IS_v1"][:k]])<=threshold)
    results.append({"case":cid,"source_family":"Sun2014_PLOS" if material!="SWCNT" else "LysozymeSWCNT2025",
     "protein":protein,"pdb":pdb,"surface":material,"diameter_nm":diameter,"channel":channel,
     "chemistry_mapping":"fixed_proxy_no_material_atomistic_model","method":"Random_same_patch_pool",
     "k":k,"threshold_degrees":threshold,"n_candidates_generated":len(p),
     "n_source_MD_contacts":len(labels),"reference_face_resultant_length":face_coherence,
     "face_ambiguity_flag":int(face_coherence<.5),
     "best_angle_degrees":float(np.mean([np.min(angles[pick]) for pick in sample])),
     "orientation_hit":float(hit.mean()),"first_success_rank":None,
     "empirical_upper_tail_p_for_IS":float((1+int(sum(hit>=obs)))/(1+nmc)),
     "doi":PLOS_DOI if material!="SWCNT" else CNT_DOI,
     "analysis_tier":"PATCH_INDEX_RANDOM_NULL_5000_NO_MD_RUN"})
  for aa,n in labels:
   matched_residues.append({"case":cid,"pdb":pdb,"residue_num":n,"residue_name":aa,"doi":PLOS_DOI if material!="SWCNT" else CNT_DOI,
    "model_sequence_identity_verified":1})
  qc.append({"case":cid,"status":"COMPUTED","n_source_contacts":len(labels),"n_candidate_patches":len(p),"source_face_resultant_length":face_coherence})
  print("POSE_UTILITY_CASE",cid,"MD_face_coherence",round(face_coherence,3),
    "top1_degrees",round(float(angles[0]),1),"top5_min_degrees",round(float(min(angles[:5])),1),flush=True)
 except Exception as e:
  qc.append({"case":cid,"status":"EXCLUDED_ERROR","details":str(e)})
  print("POSE_UTILITY_EXCLUSION",cid,repr(e),flush=True)
save("MD_orientation_screen_all_methods.csv",results)
save("MD_pose_ranked_candidate_angles.csv",rankings_all)
save("published_MD_residue_ground_truth.csv",matched_residues);save("case_qc.csv",qc)
# GB1 mutation PRIORITIZATION decision test: rank variants by distance to IS Top-k
# and compare to the same rank using geometry-only/accessibility and coherence; 
# endpoint HOTSPOT RETRIEVAL for D46 and T53 and observed kinetic effect.
# These seven variants are retrospective from ONE publication (not independent intervention).
try:
 prepared,nodes=make("2GB1")
 xy={r["key"]:r["coord"] for r in nodes};bys={r["res_seq"]:r for r in nodes}
 model=analyze_all_maps(pdb_text=prepared,chain=None,pH=7.4)["maps"]["anionic"]["patches"]
 variants={"D46A":2.87,"T53A":7.65,"F30L":30.91,"T25A":12.02,"L7A":18.13,"T11A":18.02,"T16A":21.69}
 wildtype=17.87
 orders={"IS_v1":model,
 "Accessibility_only":sorted(model,key=lambda x:(-x["mean_accessibility"],str(x["center_key"]))),
 "Chemistry_only":sorted(model,key=lambda x:(-x["chemistry_support"],str(x["center_key"]))),
 "Coherence_only":sorted(model,key=lambda x:(-x["patch_coherence"],str(x["center_key"])))}
 mutation_rows=[];mutation_scores=[]
 for_name='not used'
 for meth,rank in orders.items():
  for k in (1,3,5):
   sel=set(x for q in rank[:k] for x in q["members"])
   def proxy_distance(variant):
    n=int(variant[1:-1]);res=bys[n];center=xy[res["key"]]
    return min(np.linalg.norm(center-xy[x]) for x in sel)
   decision=sorted(variants,key=lambda z:(proxy_distance(z),-wildtype/variants[z],z))
   # IMPORTANT: sort ties by residue number, NOT any experimental label/kinetic effect
   decision=sorted(variants,key=lambda z:(proxy_distance(z),z))
   top2=set(decision[:2])
   for i,v in enumerate(decision):
    mutation_rows.append({"model":meth,"k":k,"rank_mutation":i+1,"variant":v,
      "ca_distance_A":float(proxy_distance(v)),"experimental_kon_Minv_sinv":variants[v],
      "experimental_fold_deceleration":wildtype/variants[v],
      "hotspot_published":int(v in ("D46A","T53A")),"selected_for_top2":int(v in top2)})
   mutation_scores.append({"model":meth,"k":k,"top2_precision_hotspot":len(top2&{"D46A","T53A"})/2,
    "top2_hotspot_recall":len(top2&{"D46A","T53A"})/2,
    "top2_mean_experimental_fold_deceleration":float(np.mean([wildtype/variants[v] for v in top2])),
    "chance_two_from_seven_both_hotspots":1/21,
    "source_doi":"10.3390/biom16010114","interpretation":"retrospective test, NOT prospectively validated experimental savings"})
 save("GB1_mutant_prioritization_all_methods.csv",mutation_rows)
 save("GB1_mutant_decision_scores.csv",mutation_scores)
 print("GB1_MUTATION_UTILITY",json.dumps(mutation_scores),flush=True)
except Exception as e:
 print("GB1_MUTANT_SCREEN_FAILED",str(e),flush=True)
 (OUT/"GB1_ERROR.txt").write_text(str(e))
summary={"n_conditions":len(CASES),"n_computed":sum(z["status"]=="COMPUTED" for z in qc),
 "n_unique_pdbs":len({x["pdb"] for x in results}),"single_face_coherence_threshold":.5,
 "primary_angular_hit_degrees":45,"secondary_thresholds_degrees":[30,60],
 "not_executed":"No actual material force field MD, adsorption free energies, relaxation simulations, or wet lab performed.",
 "interpretation":"Geometric decision-utility proxy: ability to nominate k orientations near published MD contact-face direction. NOT true MD acceleration or confirmed savings.",
 "chemistry_proxy":"PLOS silica sizes not represented; SWCNT chemistry pi_carbon proxy; original chemistry ranking frozen",
 "prospective_validation_status":"Retrospective; independent experimental and physically explicit MD prospective protocol remains outstanding",
 "source_dois":[PLOS_DOI,CNT_DOI,"10.3390/biom16010114"]}
(OUT/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
print("UTILITY_POSE_PROBE_COMPLETE",json.dumps(summary),flush=True)
