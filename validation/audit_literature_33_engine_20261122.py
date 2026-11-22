"""Real InterfaceScout 11-map compatibility audit for 33 historical literature PDBs.

TECHNICAL ONLY. Ligands and waters stripped in software's conventional
protein-only preparation; this may change experimental state. NMR first
model used for technical compatibility only. No adsorption prediction.
"""
import argparse,csv,json,sys,traceback,hashlib,subprocess,string
from pathlib import Path
import urllib.request
import gemmi
ROOT=Path(__file__).resolve().parents[1]
OUT=Path("literature33_engine");OUT.mkdir(exist_ok=True)
IDS="4F5S 2VUF 2CBA 1UBQ 1JNJ 3ECA 2HAV 1TGU 3NWV 7RSA 3LN2 1LYZ 1HRC 1MBN 2HHB 2PTN 1A5R 2Y0G 256L 1AML 2OED 2GB1 7KWI 1BEB 1F6R 1H3X 1HZH 2VH7 3GHG 3STN 5J17 6O5E 8ONG".split()
CANON=set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
def read(pid):
 url=f"https://files.rcsb.org/download/{pid}-assembly1.cif"
 with urllib.request.urlopen(url,timeout=75) as r:data=r.read()
 tmp=OUT/f"{pid}-assembly1.cif";tmp.write_bytes(data)
 return tmp,hashlib.sha256(data).hexdigest()
def prepare(pid):
 path,sha=read(pid);st=gemmi.read_structure(str(path))
 st.remove_alternative_conformations();st.remove_hydrogens()
 clean=gemmi.Structure();clean.name=pid;clean.cell=st.cell;clean.spacegroup_hm=st.spacegroup_hm
 model=gemmi.Model(1);mapping=[]
 avail=string.ascii_uppercase+string.ascii_lowercase+string.digits
 for chain in st[0]:
  residues=[r for r in chain.get_polymer() if r.name in CANON]
  if not residues:continue
  if len(mapping)>=len(avail):raise RuntimeError("too_many_polymer_chains_for_PDB")
  ch=gemmi.Chain(avail[len(mapping)])
  for res in residues:ch.add_residue(res.clone())
  model.add_chain(ch);mapping.append({"original_chain":chain.name,"prepared_chain":ch.name,
                                      "n_amino_acids":len(residues)})
 if not mapping:raise ValueError("no canonical protein residues")
 clean.add_model(model)
 text=clean.make_pdb_string()
 if not any(line.startswith("ATOM") for line in text.splitlines()):raise ValueError("no PDB ATOM records")
 (OUT/f"{pid}_protein_only.pdb").write_text(text,encoding="utf8")
 return {"pdb":pid,"original_sha256":sha,"protein_only_sha256":hashlib.sha256(text.encode()).hexdigest(),
         "pdb_bytes":len(text),"n_chains":len(mapping),"chain_mapping":mapping,
         "preparation_impact":"Removed water, ligands, ions, noncanonical partners; first structural model; altloc collapse (NOT structural validation)"}
def work(pid):
 sys.path.insert(0,str(ROOT/"backend"))
 from interface_engine import analyze_all_maps
 from research_workflows import combine_maps
 t=(OUT/f"{pid}_protein_only.pdb").read_text()
 a=analyze_all_maps(pdb_text=t,chain=None,pH=7.4,ionic_mM=150,temp_K=298)
 order=list(a.get("map_order",[]))
 if len(order)!=11 or len(set(order))!=11:raise RuntimeError("11 canonical maps missing")
 profiles=a["residue_profiles"];keys=[r["key"] for r in profiles]
 if len(set(keys))!=len(keys):raise RuntimeError("duplicate keys")
 mapped={m:len(a["maps"][m].get("patches",[])) for m in order}
 combinations={}
 pairs=[
  ("anionic","hbond_donor"),("hydrophobic","pi_carbon"),
  ("hydrophobic","hbond_acceptor"),("cationic","hbond_acceptor")
 ]
 for pair in pairs:
  z=combine_maps(a,pair,max_regions=None)
  assert z["regional_profile_method"]["ranking_unchanged"]
  combinations["+".join(pair)]={
    "n_regions":z["regions_total"],
    "n_exact_overlap":z["n_exact_overlap"],
    "n_regions_co_localized":sum(r["all_maps_colocated"] for r in z["regions"]),
    "n_with_pairwise_CA_distance":sum(r["regional_profile"]["min_distinct_intermap_CA_distance_A"] is not None for r in z["regions"])
  }
 return {"n_maps":len(order),"n_residues":len(profiles),"n_surface":a["diagnostics"]["n_surface_residues"],
         "n_patches_per_map":mapped,"combined_checks":combinations}
def main():
 args=argparse.ArgumentParser();args.add_argument("--worker");args.add_argument("--shard",type=int,default=0);args.add_argument("--shards",type=int,default=6);a=args.parse_args()
 if a.worker:
  try:print("RESULT_JSON "+json.dumps(work(a.worker)),flush=True)
  except BaseException:traceback.print_exc();sys.exit(1)
  return
 results=[]
 for pid in IDS[a.shard::a.shards]:
  r={"pdb":pid}
  try:
   r.update(prepare(pid))
   p=subprocess.run([sys.executable,__file__,"--worker",pid],capture_output=True,text=True,timeout=720)
   if p.returncode!=0 or "RESULT_JSON " not in p.stdout:raise RuntimeError(p.stderr[-2300:] or p.stdout[-1000:] or "process failure")
   r.update(json.loads(p.stdout.split("RESULT_JSON ",1)[1].strip().splitlines()[0]))
   r["status"]="PASS_TECHNICAL_ONLY"
  except Exception as ex:
   r["status"]="FAIL_ENGINE_OR_PREP";r["error"]=f"{type(ex).__name__}: {ex}"
  results.append(r)
  print("ENGINE_AUDIT",pid,r["status"],r.get("n_residues"),str(r.get("error",""))[:170],flush=True)
 folder=OUT/f"shard_{a.shard}";folder.mkdir(exist_ok=True)
 (folder/"engine_records.json").write_text(json.dumps(results,indent=2))
 with (folder/"engine_records.csv").open("w",newline="") as f:
  keys=["pdb","status","n_maps","n_residues","n_surface","n_chains","n_patches_per_map","combined_checks","chain_mapping","original_sha256","protein_only_sha256","preparation_impact","error"]
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
  for r in results:w.writerow({k:json.dumps(r.get(k),ensure_ascii=False) if isinstance(r.get(k),(dict,list)) else r.get(k,"") for k in keys})
 if len(results)!=len(IDS[a.shard::a.shards]):raise RuntimeError("cohort incomplete")
if __name__=="__main__":main()
