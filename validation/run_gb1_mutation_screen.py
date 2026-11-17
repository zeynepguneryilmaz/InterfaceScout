"""Independent experimental GB1 sulfate-latex mutant pilot.
Published 2026 Biomolecules DOI 10.3390/biom16010114 Table 2:
k_on for seven alanine/core variants, measured in PBS.
The source mechanism identifies docking residues D46 and T53.
This measures mutation kinetic effects, not a direct contact mask.
"""
import csv,json,sys,urllib.request
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"validation/independent_md_functional/gb1_mutagenesis";OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/"backend"))
from prepare import prepare_pdb_text
from geometry import extract_ca_nodes
from interface_engine import analyze_all_maps
import core
source="https://www.mdpi.com/2218-273X/16/1/114";DOI="10.3390/biom16010114"
variants={"D46A":2.87,"T53A":7.65,"F30L":30.91,"T25A":12.02,"L7A":18.13,"T11A":18.02,"T16A":21.69}
WT=17.87
raw=urllib.request.urlopen("https://files.rcsb.org/download/2GB1.pdb",timeout=40).read().decode()
prepared,_=prepare_pdb_text(raw,chain="A")
(OUT/"2GB1_prepared.pdb").write_text(prepared)
nodes=extract_ca_nodes(prepared)
nr={int(x["res_seq"]):x for x in nodes};name_to_3={"D":"ASP","T":"THR","F":"PHE","L":"LEU"}
for variant in variants:
 n=int(variant[1:-1])
 if n not in nr or nr[n]["res_name"]!=name_to_3[variant[0]]:raise ValueError(f"2GB1 numbering or identity mismatch {variant} {nr.get(n)}")
env=core.EnvParams(pH=7.4,ionic=150,temp=298)
results=analyze_all_maps(pdb_text=prepared,chain=None,pH=env.pH)["maps"]
names=["anionic","hydrophobic","hbond_acceptor","cationic"]
xyz={x["key"]:x["coord"] for x in nodes}
rows=[]
for v,kon in variants.items():
 n=int(v[1:-1]);key=nr[n]["key"]
 for ch in names:
  p=results[ch]["patches"]
  for k in (1,3,5):
   selected=set(x for patch in p[:k] for x in patch["members"])
   dist=min([float(np.linalg.norm(xyz[key]-xyz[m])) for m in selected],default=float("nan"))
   rows.append({"variant":v,"residue":n,"residue_name":nr[n]["res_name"],"kon_mutant_Minv_sinv":kon,"kon_WT_Minv_sinv":WT,
    "kon_WT_div_mutant":WT/kon,"delta_ln_kon":float(np.log(WT/kon)),
    "reported_transition_docking_hotspot":int(v in ("D46A","T53A")),
    "PDB":"2GB1","channel":ch,"top_k":k,
    "ranked_patch_contains_residue":int(key in selected),"distance_CA_to_selected_A":dist,
    "near5_A":int(dist<=5),"pH_model_reference":7.4,"pH_source":"SOURCE_PBS_BUFFER_PH_NOT_CONFIRMED",
    "surface_mapping_note":"negative sulfate latex is mapped to anionic. Other channels exploratory counterfactual; nitrate/sulfate chemistry is not explicitly modeled.",
    "DOI":DOI,"evidence":"site-directed mutation kinetic effect, not direct atomic contact"})
with (OUT/"gb1_mutation_vs_IS.csv").open("w",newline="",encoding="utf8") as f:
 w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
summary={"source":source,"doi":DOI,"structure":"2GB1","n_mutants":7,"n_independent_proteins":1,
 "pdb_sequence_mapping_all_seven":True,
 "source":"Biomolecules 2026 Table 2 k_on and residue-level docking interpretation",
 "predefined_surface":"anionic, negatively charged sulfate latex proxy; non-native surface representation",
 "critical_limitation":"Measured mutation effects can include adsorption-induced conformational changes and need not imply direct contact. GB1 internal Phe30 is functionally linked but buried in native state.",
 "no_stated_predictive_claim":"Descriptive pilot only, 1 protein and previously inspected study"}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2))
print("GB1_MUTAGENESIS_COMPLETE",json.dumps(summary),flush=True)
