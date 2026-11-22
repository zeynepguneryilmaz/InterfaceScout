"""Independent residue-level sequence and assembly audit of clean 48 PDBs."""
import csv,hashlib,json,re,requests,time
from pathlib import Path
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
import gemmi
ROOT=Path(__file__).resolve().parents[1]
P=Path("strict_apo_final_qa");P.mkdir(exist_ok=True)
ROWS=list(csv.DictReader(open(ROOT/"validation"/"strict_apo_48_selected_candidates_v2_final.csv",encoding="utf-8")))
assert len(ROWS)==48 and len(set(x["pdb_id"] for x in ROWS))==48
STD={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
SIDE={"ALA":"CB","ARG":"CB CG CD NE CZ NH1 NH2","ASN":"CB CG OD1 ND2","ASP":"CB CG OD1 OD2","CYS":"CB SG","GLN":"CB CG CD OE1 NE2","GLU":"CB CG CD OE1 OE2","GLY":"","HIS":"CB CG ND1 CD2 CE1 NE2","ILE":"CB CG1 CG2 CD1","LEU":"CB CG CD1 CD2","LYS":"CB CG CD CE NZ","MET":"CB CG SD CE","PHE":"CB CG CD1 CD2 CE1 CE2 CZ","PRO":"CB CG CD","SER":"CB OG","THR":"CB OG1 CG2","TRP":"CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2","TYR":"CB CG CD1 CD2 CE1 CE2 CZ OH","VAL":"CB CG1 CG2"}
def fetch(url):
 for k in range(4):
  try:
   r=requests.get(url,timeout=35);r.raise_for_status();return r.content
  except Exception:
   if k==3:raise
   time.sleep(k+1)
def one(row):
 pid=row["pdb_id"]
 result={"pdb":pid,"flags":[],"chains":[]}
 try:
  content=fetch(f"https://files.rcsb.org/download/{pid}-assembly1.cif")
  result["sha_match"]=hashlib.sha256(content).hexdigest()==row["assembly_sha256"]
  if not result["sha_match"]:result["flags"].append("ASSEMBLY_HASH_CHANGED")
  cif=P/(pid+"-assembly1.cif");cif.write_bytes(content)
  st=gemmi.read_structure(str(cif))
  e=json.loads(fetch(f"https://data.rcsb.org/rest/v1/core/entry/{pid}"))
  entities=[]
  for ent in e["rcsb_entry_container_identifiers"]["polymer_entity_ids"]:
   entity=json.loads(fetch(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{ent}"))
   ep=entity.get("entity_poly",{})
   if "polypeptide" not in str(ep.get("type","")).lower():continue
   seq=re.sub("[^A-Z]","",ep.get("pdbx_seq_one_letter_code_can","").upper())
   auth=entity.get("rcsb_polymer_entity_container_identifiers",{}).get("auth_asym_ids") or []
   entities.append({"seq":seq,"auth":auth})
  for ch in st[0]:
   residues=[r for r in ch.get_polymer() if r.name in STD]
   if not residues:continue
   original=ch.name.split("-")[0]
   refs=[x for x in entities if original in x["auth"]]
   if len(refs)!=1:result["flags"].append("CHAIN_ENTITY_MAPPING");continue
   seq=refs[0]["seq"]
   labels=[int(r.label_seq) for r in residues]
   if any(i<=0 for i in labels):
    result["flags"].append("NO_LABEL_SEQ_"+ch.name)
    continue
   internal_missing=sum(max(0,b-a-1) for a,b in zip(labels,labels[1:]))
   mismatches=[]
   for r,label in zip(residues,labels):
    if label<=len(seq):
     letter=gemmi.find_tabulated_residue(r.name).one_letter_code
     if letter!=seq[label-1]:mismatches.append({"label_seq":label,"observed":letter,"expected":seq[label-1]})
    else:mismatches.append({"label_seq":label,"observed":r.name,"expected":"OUT_OF_RANGE"})
   side=0
   alt=0
   for r in residues:
    atoms={a.name for a in r if a.element.name not in ("H","D")}
    side+=bool(set(SIDE[r.name].split())-atoms)
    alt+=any(a.altloc not in (chr(0)," ",".","?") for a in r)
   c={"chain":ch.name,"declared_len":len(seq),"observed":len(residues),"terminal_missing_N":max(0,labels[0]-1),
      "terminal_missing_C":max(0,len(seq)-labels[-1]),"internal_missing_positions":internal_missing,
      "sequence_mismatch_count":len(mismatches),"sequence_mismatch_examples":mismatches[:5],
      "missing_sidechain_residues":side,"altloc_residues":alt}
   result["chains"].append(c)
   if internal_missing:result["flags"].append("INTERNAL_MISSING_"+ch.name)
   if mismatches:result["flags"].append("SEQ_MISMATCH_"+ch.name)
   if side:result["flags"].append("SIDECHAIN_MISSING_"+ch.name)
   if alt:result["flags"].append("ALTLOC_"+ch.name)
  if not result["chains"]:result["flags"].append("NO_CHAINS")
 except Exception as err:result["flags"].append("FETCH_PARSE_FAILED");result["error"]=str(err)
 return result
with ThreadPoolExecutor(max_workers=12) as ex:records=list(ex.map(one,ROWS))
families=Counter()
topos=Counter()
for row in ROWS:
 for f in set(row["CATH_H"].split(";")):
  families[f]+=1
  topos[".".join(f.split(".")[:3])]+=1
assert all(v==1 for v in families.values()),"CATH H nonredundancy violated"
assert all(v==1 for v in topos.values()),"CATH T nonredundancy violated"
for x in records:print("RESIDUE_QA "+json.dumps(x,separators=(",",":")),flush=True)
summary={"n":len(records),"nonredundant_H":len(families),"nonredundant_T":len(topos),
 "n_hash_match":sum(x.get("sha_match",False) for x in records),
 "n_without_flags":sum(not x["flags"] for x in records),
 "flagged":[{"pdb":x["pdb"],"flags":x["flags"]} for x in records if x["flags"]]}
(P/"residue_48.json").write_text(json.dumps(records,indent=2))
(P/"summary.json").write_text(json.dumps(summary,indent=2))
with (P/"residue_48.csv").open("w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=["pdb","sha_match","n_chains","internal_missing","sequence_mismatches","sidechain_missing","altloc","flags"])
 w.writeheader()
 for r in records:
  chains=r["chains"]
  w.writerow({"pdb":r["pdb"],"sha_match":r.get("sha_match"),
  "n_chains":len(chains),"internal_missing":sum(x["internal_missing_positions"] for x in chains),
  "sequence_mismatches":sum(x["sequence_mismatch_count"] for x in chains),
  "sidechain_missing":sum(x["missing_sidechain_residues"] for x in chains),
  "altloc":sum(x["altloc_residues"] for x in chains),"flags":";".join(r["flags"])})
print("RESIDUE_QA_SUMMARY "+json.dumps(summary,separators=(",",":")),flush=True)
if summary["n_without_flags"]!=48:raise SystemExit("Some strict dataset entries require structural review")
