"""Independent function-aware screening, no measured enzyme activity claims."""
from pathlib import Path
import csv, io, json, sys, urllib.request
import numpy as np
from Bio.PDB import PDBParser
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from interface_engine import analyze_all_maps
from prepare import prepare_pdb_text
OUT=Path("validation/stage2_results");OUT.mkdir(parents=True,exist_ok=True)
CASES=[
 {"name":"HCAII_pH6.3","pdb":"2CBA","chain":"A","pH":6.3,"core":[64,94,96,119],"extended":[64,94,96,119,106,199],"reference":"https://pmc.ncbi.nlm.nih.gov/articles/PMC5755581/"},
 {"name":"HCAII_pH8.3","pdb":"2CBA","chain":"A","pH":8.3,"core":[64,94,96,119],"extended":[64,94,96,119,106,199],"reference":"https://pmc.ncbi.nlm.nih.gov/articles/PMC5755581/"},
 {"name":"Lysozyme_pH7.4","pdb":"1LYZ","chain":"A","pH":7.4,"core":[35,52],"extended":[35,52],"reference":"https://pdbj.org/mine/functional_details/1lyz"}
]
def download(p):
 u="https://files.rcsb.org/download/"+p+".pdb"
 with urllib.request.urlopen(urllib.request.Request(u,headers={"User-Agent":"InterfaceScout-stage2/1.0"}),timeout=60) as r:return r.read().decode("utf8","replace")
def parse(text,pdb):
 return next(iter(PDBParser(QUIET=True).get_structure(pdb,io.StringIO(text))))
def heavy_coords(residue):
 return np.array([a.coord for a in residue if a.element not in ("H","D")],dtype=float)
def atom_distance(patch_res,active_res):
 aa=[heavy_coords(r) for r in patch_res if len(heavy_coords(r))];bb=[heavy_coords(r) for r in active_res if len(heavy_coords(r))]
 if not aa or not bb:return None
 a=np.vstack(aa);b=np.vstack(bb);return float(np.min(np.linalg.norm(a[:,None,:]-b[None,:,:],axis=-1)))
def ca_distance(patch_res,active_res):
 a=np.array([r["CA"].coord for r in patch_res if "CA" in r]);b=np.array([r["CA"].coord for r in active_res if "CA" in r])
 return float(np.min(np.linalg.norm(a[:,None,:]-b[None,:,:],axis=-1))) if len(a) and len(b) else None
rows=[];detail=[];meta=[]
for case in CASES:
 raw=download(case["pdb"]);(OUT/(case["pdb"]+"_original.pdb")).write_text(raw)
 prepared, prep_report=prepare_pdb_text(raw,chain=case["chain"])
 (OUT/(case["name"]+"_prepared.pdb")).write_text(prepared)
 model=parse(prepared,case["pdb"])
 chain=model[case["chain"]]
 residues={r.id[1]:r for r in chain if r.id[0]==" " and "CA" in r}
 for typ in ("core","extended"):
  missing=[n for n in case[typ] if n not in residues]
  if missing:raise ValueError(case["name"]+": missing functional residues "+repr(missing))
 out=analyze_all_maps(pdb_text=prepared,chain=case["chain"],pH=case["pH"])
 if out["n_maps"]!=11:raise AssertionError("not 11 maps")
 meta.append({k:v for k,v in case.items() if k not in ("core","extended")} | {"core_residues":",".join(map(str,case["core"])),"extended_residues":",".join(map(str,case["extended"])),"n_surface_residues":out["diagnostics"]["n_surface_residues"]})
 for mapname,mp in out["maps"].items():
  patches=mp["patches"]
  for k in (1,3,5):
   selected=patches[:k]
   patch_rs=set()
   for p in selected:
    for m in p["member_residues"]:
     if str(m.get("chain"))==case["chain"]:
      rid=int(m["res_seq"])
      if rid in residues:patch_rs.add(rid)
   member=[residues[i] for i in sorted(patch_rs)]
   row={"case":case["name"],"pdb_id":case["pdb"],"pH":case["pH"],"chemistry":mapname,"top_k":k,"n_patches":len(patches),"union_residues":len(member),"residue_ids":";".join(map(str,sorted(patch_rs)))}
   for typ in ("core","extended"):
    active=[residues[i] for i in case[typ]]
    row[typ+"_CA_min_A"]=ca_distance(member,active)
    row[typ+"_heavy_min_A"]=atom_distance(member,active)
    row[typ+"_direct_overlap"]=len(patch_rs.intersection(case[typ]))
    row[typ+"_within_5A_heavy"]=int(row[typ+"_heavy_min_A"]<=5) if row[typ+"_heavy_min_A"] is not None else None
    row[typ+"_within_8A_heavy"]=int(row[typ+"_heavy_min_A"]<=8) if row[typ+"_heavy_min_A"] is not None else None
   rows.append(row)
  for p in patches:
   members=[residues[int(m["res_seq"])] for m in p["member_residues"] if str(m.get("chain"))==case["chain"] and int(m["res_seq"]) in residues]
   z={"case":case["name"],"chemistry":mapname,"rank":p["display_rank"],"pareto_front":p["pareto_front"],"center_key":p["center_key"],"n_members":len(members),"chemistry_support":p["chemistry_support"],"mean_accessibility":p["mean_accessibility"],"patch_coherence":p["patch_coherence"],"orientation_coherence":p["orientation_coherence"]}
   for typ in ("core","extended"):
    a=[residues[i] for i in case[typ]]
    z[typ+"_CA_min_A"]=ca_distance(members,a)
    z[typ+"_heavy_min_A"]=atom_distance(members,a)
   detail.append(z)
 print("CASE",case["name"],"maps",len(out["maps"]),"surface",out["diagnostics"]["n_surface_residues"],flush=True)
def csvout(fname,data):
 if not data:return
 with (OUT/fname).open("w",newline="") as f:
  writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
csvout("functional_site_top1_top3_top5.csv",rows)
csvout("patch_level_functional_distances.csv",detail)
csvout("input_provenance.csv",meta)
(OUT/"methods.json").write_text(json.dumps({"scope":"exploratory structural screening, not prediction of enzyme activity or adsorption orientation","cases":CASES,"patch_atom_distance":"minimum nonhydrogen atom-atom Euclidean distance between cumulative ranked patch member residues and selected functional site residues; a zero distance may denote shared active site residue","patch_CA_distance":"minimum C-alpha distance","excluded_cofactor":"Structure preparation drops hetero groups including HCAII zinc; does not model holo active-site electrostatics","pH_effect":"pH enters residue state weighting only; no pH-induced structural rearrangement","n_maps":11},indent=2))
try:
 import matplotlib
 matplotlib.use("Agg")
 import matplotlib.pyplot as plt
 chemistry=list(dict.fromkeys(r["chemistry"] for r in rows))
 fig,axs=plt.subplots(1,3,figsize=(15,7),sharey=True)
 for ax,case in zip(axs,CASES):
  rr={r["chemistry"]:r for r in rows if r["case"]==case["name"] and r["top_k"]==5}
  vals=[rr[c]["core_heavy_min_A"] if rr[c]["core_heavy_min_A"] is not None else np.nan for c in chemistry]
  ax.barh(range(len(chemistry)),vals)
  ax.set_yticks(range(len(chemistry)),chemistry)
  ax.invert_yaxis();ax.set_title(case["name"].replace("_"," "))
  ax.set_xlabel("Minimum heavy-atom distance (Å)")
 fig.suptitle("Function-aware screening: Top-5 candidate regions vs catalytic residues")
 fig.tight_layout();fig.savefig(OUT/"Fig_stage2_functional_distance.png",dpi=300,bbox_inches="tight")
except Exception as e:print("PLOT_ERROR",repr(e),flush=True)
print("STAGE2_COMPLETE",json.dumps({"cases":len(CASES),"maps":len(rows)//3,"topk_rows":len(rows),"patch_rows":len(detail),"min_distance_examples":[{"case":c["name"],"closest_top5":min((r["core_heavy_min_A"] for r in rows if r["case"]==c["name"] and r["top_k"]==5 and r["core_heavy_min_A"] is not None),default=None)} for c in CASES]}),flush=True)
