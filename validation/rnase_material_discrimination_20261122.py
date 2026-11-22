"""InterfaceScout RNase A literature-based material discrimination, unfit pilot.
DOI 10.1021/la503854a, original Tables 4 and 5. Observations are accessibility
changes (NOT direct material contacts). Same PDB/pH for all material conditions.
"""
import csv,hashlib,json,random,sys,urllib.request,itertools
from pathlib import Path
import numpy as np
from scipy.stats import hypergeom
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from interface_engine import analyze_all_maps
from research_workflows import combine_maps
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes,build_surface_geometry,ca_distance,same_face
import core
O=Path("rnase_cross_condition_results");O.mkdir(exist_ok=True)
MATS=["Glass","HDPE","PMMA"]
CHEM={"Glass":["anionic","hbond_donor"],"HDPE":["hydrophobic"],"PMMA":["hydrophobic","hbond_acceptor"]}
DOI="10.1021/la503854a"
# Strongly and moderately REDUCED AAL/MS label-accessibility bins.
GT={
"0.03":{"Glass":["K1 R39 R85 E111","K7 R10 K31 D38 D53 K98"],
"HDPE":["D53 K61 K66 Y76 K91 K98","K104 R39 K41 E49 Y73 Y92"],
"PMMA":["E111","D38 R39 K66 Y76 K91 Y92 K98"]},
"1.00":{"Glass":["","K1 K7 E9 D38 K61 K66 Y76 K91 K104 H105 H119"],
"HDPE":["K1 K7 K31 K61 K66 K91 K98","K37 D38 Y76 K104 H105 H119"],
"PMMA":["K1 K66 K31 R85 K7","R33 K37 D38 R39 K41 D53 K61 Y76 K91 Y92 E111 H105 H119"]}}
# Measured, NOT model-selected, residues, from Table 4/5 all five bins.
ASSAY="K1 E2 K7 E9 R10 H12 D14 H15 Y25 K31 R33 K37 D38 R39 K41 H48 E49 D53 K61 K66 Y73 Y76 D83 R85 E86 K91 Y92 Y97 K98 K104 H105 E111 Y115 H119 D121".split()
ASSAYD={int(s[1:]):s for s in ASSAY}
AA=dict(zip("ACDEFGHIKLMNPQRSTVWY","ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()))
R3={v:k for k,v in AA.items()}
N_RANDOM=3000
BUDGETS=[5,10,15]
def clean(s):return s.split()
G={c:{m:{int(v[1:]):cat for cat,raw in enumerate(cats) for v in clean(raw)} for m,cats in groups.items()} for c,groups in GT.items()}
def save(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with (O/name).open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
  for row in rows:w.writerow({k:json.dumps(row[k]) if isinstance(row.get(k),(list,dict)) else row.get(k) for k in keys})
def idx(k):return int(k.split(":")[1])
def rank_from_patches(m,U):
 rows={r["key"]:r for r in m["residues"]};out=[];seen=set()
 def add(k):
  if k in U and k not in seen:seen.add(k);out.append(k)
 for p in m["patches"]:
  keys=sorted(p["members"],key=lambda k:(-float(rows.get(k,{}).get("propensity") or 0),k))
  for k in keys:add(k)
 for r in sorted(m["residues"],key=lambda r:(-float(r.get("propensity") or 0),r["key"])):add(r["key"])
 return out
def colocation_rank(analysis,channels,U,perchan,rsa):
 if len(channels)==1:return perchan[channels[0]],None,{}
 combo=combine_maps(analysis,channels,max_regions=None)
 chem={m:{r["key"]:float(r.get("local_score") or 0) for r in analysis["maps"][m]["residues"]} for m in channels}
 out=[];seen=set();sources={}
 def add(k,source):
  if k in U and k not in seen:out.append(k);seen.add(k);sources[k]=source
 for region in combo["regions"]:
  if not region["all_maps_colocated"]:continue
  keys=sorted(region["member_keys"],key=lambda k:(-sum(chem[m].get(k,0)>0 for m in channels),k))
  for k in keys:add(k,"colocated")
 ncolo=len(out)
 for m in channels:
  for k in perchan[m]:add(k,"single_channel_fallback")
 for k in sorted(U,key=lambda k:(-rsa[k],k)):add(k,"rsa_fallback")
 return out,combo,{"n_colocated_unique_exposed_assayed":ncolo,"selection_source":sources}
def metric(picked,ground,U):
 hits=len(picked&ground);n=len(ground);k=len(picked);N=len(U)
 return {"hit":hits,"recall":hits/n if n else None,"precision":hits/k if k else None,
 "expected_random_hits":k*n/N if N else 0,
 "hypergeom_one_sided_p":float(hypergeom.sf(hits-1,N,n,k)) if k and n else None}
def diag_score(mat):
 return sum(mat[i][i] for i in range(3))/3-sum(mat[i][j] for i in range(3) for j in range(3) if j!=i)/6
def main():
 allmetrics=[];deltas=[];regions=[];ranked=[];mappings=[];controls=[];provenance=[]
 for pid in ["7RSA","1KF5"]:
  print("RUN",pid,flush=True)
  u=f"https://files.rcsb.org/download/{pid}.pdb"
  with urllib.request.urlopen(u,timeout=60) as f: raw=f.read()
  (O/(pid+".pdb")).write_bytes(raw)
  a=analyze_all_maps(pdb_text=raw.decode(errors="replace"),chain="A",pH=7.4,ionic_mM=10,temp_K=298.15)
  prof={r["key"]:r for r in a["residue_profiles"] if r["chain"]=="A"}
  if not prof:raise RuntimeError("No chain A profiles "+pid)
  seqmap={int(r["res_seq"]):r for r in prof.values()}
  bad=[(n,s) for n,s in ASSAYD.items() if n not in seqmap or R3.get(seqmap[n]["res_name"])!=s[0]]
  if bad:raise RuntimeError("NUMBERING/IDENTITY_FAIL "+pid+" "+str(bad))
  U={k for k,v in prof.items() if v.get("surface_exposed") and int(v["res_seq"]) in ASSAYD}
  assert len(U)>=15
  rsa={k:float(prof[k]["scrsa_raw"]) for k in U}
  maps=a["maps"]
  perchan={m:rank_from_patches(maps[m],U) for m in maps}
  ranks={};method_ranks={}
  for mat,chans in CHEM.items():
   rc,combo,meta=colocation_rank(a,chans,U,perchan,rsa)
   if len(rc)<len(U):
    for k in sorted(U,key=lambda k:(-rsa[k],k)):
     if k not in rc:rc.append(k)
   best=perchan[chans[0]]+[k for k in sorted(U) if k not in perchan[chans[0]]]
   chem_rows={name:{r["key"]:float(r.get("propensity") or 0) for r in maps[name]["residues"]} for name in chans}
   chemistry=sorted(U,key=lambda k:(-max(chem_rows[name].get(k,0) for name in chans),k))
   sasa=sorted(U,key=lambda k:(-rsa[k],k))
   method_ranks[mat]={"combined":rc,"best_single_patch":best,"chemistry_only":chemistry,"sasa_only":sasa}
   if combo:
    regions.append({"pdb":pid,"material":mat,"selected_maps":"+".join(chans),"total_regions":combo["regions_total"],
      "all_map_colocated_regions":sum(r["all_maps_colocated"] for r in combo["regions"]),
      "exact_same_residue_overlap":combo["n_exact_overlap"],"exposed_assayed_colocated":meta["n_colocated_unique_exposed_assayed"],
      "top10_colocated":sum(meta["selection_source"].get(x)=="colocated" for x in rc[:10]),
      "top10_keys":" ".join(rc[:10])})
   for method,order in method_ranks[mat].items():
    for rank,key in enumerate(order,1):
     ranked.append({"pdb":pid,"material":mat,"method":method,"rank":rank,"residue":ASSAYD[idx(key)],
       "key":key,"rsa":rsa[key]})
  for num,expected in ASSAYD.items():
   p=seqmap[num];mappings.append({"pdb":pid,"residue":expected,"observed":R3[p["res_name"]]+str(num),
      "is_exposed":bool(p["surface_exposed"]),"rsa":float(p["scrsa_raw"])})
  # Spatial control same 8A, same outward-hemisphere rule as the canonical patch geometry.
  prep,_=prepare_pdb_text(raw.decode(errors="replace"),chain="A")
  v1=core.analyze(core.AnalyzeRequest(pdb_text=prep,env=core.EnvParams(pH=7.4,ionic=10,temp=298.15)))
  geom=build_surface_geometry(v1,extract_ca_nodes(prep))
  assert U.issubset(geom["coords"])
  rng=random.Random(20261122+sum(map(ord,pid)))
  random_prefix={K:[] for K in BUDGETS if K<=len(U)}
  for j in range(N_RANDOM):
   seeds=list(U);rng.shuffle(seeds);seq=[];seen=set()
   for seed in seeds:
    group=[k for k in U if ca_distance(seed,k,geom)<=8 and same_face(seed,k,geom)]
    rng.shuffle(group)
    for k in group:
     if k not in seen:seen.add(k);seq.append(k)
   for K in random_prefix:random_prefix[K].append(set(seq[:K]))
  for conc in GT:
   for K in random_prefix:
    for method in ["combined","best_single_patch","chemistry_only","sasa_only"]:
     values=[]
     for pm in MATS:
      row=[]
      for gm in MATS:
       gt={k for k in U if G[conc][gm].get(idx(k),99)<=1}
       q=metric(set(method_ranks[pm][method][:K]),gt,U)
       allmetrics.append({"pdb":pid,"concentration_mg_ml":conc,"prediction_material":pm,"experimental_material":gm,
          "method":method,"budget":K,"n_eligible":len(U),"n_gt_eligible":len(gt),
          **q,"hit_residues":" ".join(sorted(ASSAYD[idx(k)] for k in set(method_ranks[pm][method][:K])&gt))})
       row.append(q["recall"] or 0)
      values.append(row)
     observed=diag_score(values)
     perms=[diag_score([[values[i][perm[j]] for j in range(3)] for i in range(3)]) for perm in itertools.permutations(range(3))]
     deltas.append({"pdb":pid,"concentration_mg_ml":conc,"budget":K,"method":method,
       "diagonal_minus_cross_recall":observed,
       "diagonal_mean_recall":sum(values[i][i] for i in range(3))/3,
       "offdiag_mean_recall":sum(values[i][j] for i in range(3) for j in range(3) if i!=j)/6,
       "six_label_permutation_p":sum(v>=observed-1e-9 for v in perms)/6,
       "recall_matrix":values})
    obs=next(x for x in deltas if x["pdb"]==pid and x["concentration_mg_ml"]==conc and x["budget"]==K and x["method"]=="combined")
    null_deltas=[]
    for t in range(N_RANDOM//3):
     mat=[]
     for i in range(3):
      sample=random_prefix[K][3*t+i]
      mat.append([metric(sample,{k for k in U if G[conc][gm].get(idx(k),99)<=1},U)["recall"] or 0 for gm in MATS])
     null_deltas.append(diag_score(mat))
    controls.append({"pdb":pid,"concentration_mg_ml":conc,"budget":K,"null":"random 8A same-face patch",
       "n_draws":len(null_deltas),"null_mean":float(np.mean(null_deltas)),
       "null_95pct_lower":float(np.percentile(null_deltas,2.5)),
       "null_95pct_upper":float(np.percentile(null_deltas,97.5)),
       "combined_observed":obs["diagonal_minus_cross_recall"],
       "p_upper_tail":sum(x>=obs["diagonal_minus_cross_recall"]-1e-9 for x in null_deltas)/len(null_deltas)})
  provenance.append({"pdb":pid,"file_sha256":hashlib.sha256(raw).hexdigest(),"residues":len(prof),
    "surface_exposed":sum(v["surface_exposed"] for v in prof.values()),
    "assayed":len(ASSAYD),"assayed_and_exposed":len(U),
    "assayed_buried":sorted(ASSAYD[idx(k)] for k,v in prof.items() if not v["surface_exposed"] and idx(k) in ASSAYD),
    "model_version":a.get("version"),"pdb_url":u})
  print("SUCCESS",pid,"eligible",len(U),flush=True)
 for name,rows in [("cross_condition_matrix.csv",allmetrics),("specificity.csv",deltas),
   ("combined_region_diagnostics.csv",regions),("per_residue_ranks.csv",ranked),
   ("assayed_residue_mapping.csv",mappings),("spatial_null.csv",controls)]:
  save(name,rows)
 (O/"protocol.json").write_text(json.dumps({"experimental_source":DOI,"experimental_tables":{"0.03":"Table 4","1.00":"Table 5"},
   "raw_decreased_bins":GT,"measured_residues":ASSAYD,
   "prespecified_material_map_combinations":CHEM,
   "analysis_pH":7.4,"buffer":"10 mM potassium phosphate pH7.4, 25C",
   "ionic_note":"ionic_mM=10 is a recorded input; the code does not physically model ionic strength",
   "static_map_caveat":"No protein adsorption, unfolding or protein-protein interaction prediction",
   "selection":"exact original canonical map patch order, combined co-localized neighborhoods first; deterministic fallback",
   "mapping_scope":"only chemically labelled and native structure solvent-exposed residues",
   "patch_random_null":"random seeds + 8A C-alpha and outward hemisphere neighborhoods; 1000 draws/condition/K",
   "PDB":provenance},indent=2))
 print("RESULT",json.dumps({"specificity":deltas,"colocation":regions,"provenance":provenance},default=str),flush=True)
if __name__=="__main__":main()
