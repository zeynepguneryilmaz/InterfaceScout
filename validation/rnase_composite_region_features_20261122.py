"""Exploratory, UNFITTED composite patch descriptors for RNase A validation.

Purely additive to production InterfaceScout. The existing 11 maps, SASA,
8 A patch scale and 6/9 A multiscale settings remain unchanged.
The experimental labels are used only AFTER feature/rank generation.
Results are not an independent holdout because RNase A labels informed this study.
"""
import csv,json,hashlib,math,random,itertools,urllib.request,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from interface_engine import analyze_all_maps
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes,build_surface_geometry
import core
import rnase_material_discrimination_20261122 as old

OUT=Path("rnase_feature_audit"); OUT.mkdir(exist_ok=True)
IDS=["7RSA","1KF5"]
BUDGETS=[5,10,15]
def save(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with (OUT/name).open("w",newline="",encoding="utf8") as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
  for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
def descriptor_rows(analysis,geom,selected):
 S=[k for k in analysis["diagnostics"]["surface_residue_keys"] if k in geom["coords"]]
 assert len(S)>5
 props={m:{str(r["key"]):float(r.get("local_score") or 0) for r in analysis["maps"][m]["residues"]} for m in selected}
 out=[]
 for seed in S:
  nei=[v for v in S if float(np.linalg.norm(geom["coords"][seed]-geom["coords"][v]))<=8.00001 and
       (seed==v or float(np.dot(geom["normals"][seed],geom["normals"][v]))>0)]
  assert nei
  d={m:sum(props[m].get(k,0.0) for k in nei)/len(nei) for m in selected}
  f={m:sum(props[m].get(k,0.0)>0 for k in nei)/len(nei) for m in selected}
  out.append({"center_key":seed,"neighborhood":nei,"n_neighbors":len(nei),"density":d,
              "fraction":f,"mean_scrsa":sum(geom["scrsa"][k] for k in nei)/len(nei)})
 # Normalize per-channel by its own largest observed support. Fixed data normalization,
 # not an experimental label–fitted weighting.
 max_support={m:max(v["density"][m] for v in out) for m in selected}
 for r in out:
  r["normalized_densities"]={m:(r["density"][m]/max_support[m] if max_support[m]>0 else 0) for m in selected}
  r["min_normalized_density"]=min(r["normalized_densities"].values())
  r["min_fraction"]=min(r["fraction"].values())
  r["all_maps_present"]=all(r["fraction"][m]>0 for m in selected)
  pairs=[]
  if len(selected)==2:
   a,b=selected
   for x in r["neighborhood"]:
    if props[a].get(x,0)<=0:continue
    for y in r["neighborhood"]:
     if x==y or props[b].get(y,0)<=0:continue
     dist=float(np.linalg.norm(geom["coords"][x]-geom["coords"][y]))
     if dist<=6.00001:
      pairs.append((x,y,round(dist,3)))
  r["distinct_residue_pairs_6A"]=len(pairs)
  r["nearest_pairs_preview"]=pairs[:8]
  r["chemistry_signature"]={m:round(v,5) for m,v in r["normalized_densities"].items()}
 # Percentile is relative to all SAME-protein neighborhoods, not conditions.
 for m in selected:
  z=np.asarray([r["density"][m] for r in out])
  for r in out:r.setdefault("within_protein_percentile",{})[m]=round(float(np.mean(z<=r["density"][m])),5)
 # Deterministic Pareto-like ordering: prioritize balanced support and local coverage
 # (unvalidated exploratory ranking, evaluated against SASA / chemistry baseline).
 out.sort(key=lambda r:(-r["min_normalized_density"],-r["min_fraction"],
                          -r["distinct_residue_pairs_6A"],-r["mean_scrsa"],r["center_key"]))
 return out
def rank_assayed(rows,U,analysis,geom,selected):
 prop={m:{r["key"]:float(r.get("local_score") or 0) for r in analysis["maps"][m]["residues"]} for m in selected}
 output=[];seen=set()
 for row in rows:
  n=sorted((k for k in row["neighborhood"] if k in U),
           key=lambda k:(-sum(prop[m].get(k,0)>0 for m in selected),
                         -sum(prop[m].get(k,0) for m in selected),k))
  for k in n:
   if k not in seen:seen.add(k);output.append(k)
 for k in sorted(U):
  if k not in seen:output.append(k)
 assert len(output)==len(U)
 return output
def main():
 diagnostics=[];positions=[];comp=[];permutation=[]
 for pid in IDS:
  with urllib.request.urlopen("https://files.rcsb.org/download/"+pid+".pdb",timeout=50) as r:raw=r.read()
  (OUT/(pid+".pdb")).write_bytes(raw)
  a=analyze_all_maps(pdb_text=raw.decode(errors="replace"),chain="A",pH=7.4,ionic_mM=10,temp_K=298.15)
  prof={x["key"]:x for x in a["residue_profiles"] if x["chain"]=="A"}
  U={k for k,v in prof.items() if v.get("surface_exposed") and int(v["res_seq"]) in old.ASSAYD}
  prep,_=prepare_pdb_text(raw.decode(errors="replace"),chain="A")
  v1=core.analyze(core.AnalyzeRequest(pdb_text=prep,env=core.EnvParams(pH=7.4,ionic=10,temp=298.15)))
  geom=build_surface_geometry(v1,extract_ca_nodes(prep))
  assert U.issubset(geom["coords"])
  bymat={}; descriptions={}
  for material,selected in old.CHEM.items():
   d=descriptor_rows(a,geom,selected);byrank=rank_assayed(d,U,a,geom,selected)
   bymat[material]=byrank
   descriptions[material]=d
   for rank,r in enumerate(d,1):
    positions.append({"pdb":pid,"material":material,"center_key":r["center_key"],
     "center_aa":f'{prof[r["center_key"]]["res_name"]}{prof[r["center_key"]]["res_seq"]}' if r["center_key"] in prof else "",
     "region_rank":rank,"n_neighbors":r["n_neighbors"],"neighborhood":r["neighborhood"],
     "within_protein_percentile":r["within_protein_percentile"],
     "n_distinct_pairs_6A":r["distinct_residue_pairs_6A"],
     "min_fraction":round(r["min_fraction"],5),"min_normalized_density":round(r["min_normalized_density"],5),
     "all_maps_present":r["all_maps_present"],"chemistry_signature":r["chemistry_signature"],
     "near_pairs":r["nearest_pairs_preview"]})
   diagnostics.append({"pdb":pid,"material":material,"selected":selected,"n_regions":len(d),
    "n_all_maps_present":sum(r["all_maps_present"] for r in d),
    "n_balance_atleast_0p5":sum(r["min_normalized_density"]>=.5 for r in d),
    "n_balance_atleast_0p75":sum(r["min_normalized_density"]>=.75 for r in d),
    "n_with_distinct_close_pairs":sum(r["distinct_residue_pairs_6A"]>0 for r in d),
    "top5_centers":[r["center_key"] for r in d[:5]],
    "top10_eligible_assayed":byrank[:10]})
  for conc in old.GT:
   for K in BUDGETS:
    matrix=[]
    for prediction_mat in old.MATS:
     row=[]
     for observed_mat in old.MATS:
      exp={key for key in U if old.G[conc][observed_mat].get(old.idx(key),99)<=1}
      hits=len(set(bymat[prediction_mat][:K])&exp)
      recall=hits/len(exp) if exp else None
      matrixcell={"pdb":pid,"concentration_mg_ml":conc,"budget":K,"prediction_material":prediction_mat,
        "observed_material":observed_mat,"n_eligible":len(U),"n_experimental_eligible":len(exp),
        "hits":hits,"recall":recall,"hit_residues":sorted(old.ASSAYD[old.idx(k)] for k in set(bymat[prediction_mat][:K])&exp)}
      comp.append(matrixcell);row.append(recall or 0)
     matrix.append(row)
    obs=old.diag_score(matrix)
    null=[old.diag_score([[matrix[i][perm[j]] for j in range(3)] for i in range(3)])
          for perm in itertools.permutations(range(3))]
    permutation.append({"pdb":pid,"conc_mg_ml":conc,"budget":K,"method":"balanced_patch_descriptors_exploratory",
        "diagonal_minus_cross_recall":obs,"diagonal":np.mean([matrix[i][i] for i in range(3)]),
        "cross":np.mean([matrix[i][j] for i in range(3) for j in range(3) if i!=j]),
        "six_label_permutation_p":sum(x>=obs-1e-10 for x in null)/6,"recall_matrix":matrix})
  print("RUN_DONE",pid,"assayed exposed",len(U),flush=True)
 save("regional_features.csv",positions)
 save("regional_summary.csv",diagnostics)
 save("cross_condition_3x3.csv",comp)
 save("exploratory_specificity.csv",permutation)
 (OUT/"protocol.json").write_text(json.dumps({
   "status":"exploratory_posthoc_development_data_NO_INDEPENDENT_VALIDATION",
   "feature_definition":"same 8A CA geometry + same outward hemisphere; all canonical map local_score values; balance=min of within-map max normalized local density",
   "ranking":"min normalized local density then min channel fraction then close distinct pairs then mean scRSA",
   "feature_gate":"No experimental labels used to construct features or order patch regions",
   "test_labels":"RNase A 2014 DOI 10.1021/la503854a Tables 4/5, accessibility loss not direct adsorption contact",
   "scope":"same PDB same pH 7.4; ionic=10 is not physically modeled",
   "independent_test_needed":True,
   "comparison_warning":"This is exploratory on the data that motivated the change; not an independent performance validation",
   "selected_material_maps":old.CHEM,"PDB":IDS},indent=2))
 print("SUMMARY",json.dumps({"diagnostics":diagnostics,"specificity":permutation}),flush=True)
if __name__=="__main__":main()
