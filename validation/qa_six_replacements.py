"""Targeted six-candidate physical-coordinate QA. Screening only, no automatic scientific acceptance."""
import csv, json, requests, re, hashlib, sys
from pathlib import Path
from collections import defaultdict
import gemmi

IDS=["1C5E","1PDB","1A1X","1P32","1STN","2CBA","3T9L","1QYS"]
STANDARD={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
SIDE={
 "ALA":"CB","ARG":"CB CG CD NE CZ NH1 NH2","ASN":"CB CG OD1 ND2","ASP":"CB CG OD1 OD2",
 "CYS":"CB SG","GLN":"CB CG CD OE1 NE2","GLU":"CB CG CD OE1 OE2","GLY":"",
 "HIS":"CB CG ND1 CD2 CE1 NE2","ILE":"CB CG1 CG2 CD1","LEU":"CB CG CD1 CD2",
 "LYS":"CB CG CD CE NZ","MET":"CB CG SD CE","PHE":"CB CG CD1 CD2 CE1 CE2 CZ",
 "PRO":"CB CG CD","SER":"CB OG","THR":"CB OG1 CG2","TRP":"CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
 "TYR":"CB CG CD1 CD2 CE1 CE2 CZ OH","VAL":"CB CG1 CG2"
}
base="https://data.rcsb.org/rest/v1/core"
s=requests.Session();s.headers["User-Agent"]="InterfaceScout-targeted-six-structure-QA/1.0"
out=Path("six_candidate_qa");out.mkdir(exist_ok=True)
def get(path):
 r=s.get(path,timeout=35);r.raise_for_status();return r
def metadata(pid):
 e=get(f"{base}/entry/{pid}").json()
 a=get(f"{base}/assembly/{pid}/1").json()
 entities={}
 for ent in e["rcsb_entry_container_identifiers"]["polymer_entity_ids"]:
  p=get(f"{base}/polymer_entity/{pid}/{ent}").json()
  poly=p.get("entity_poly",{})
  if "polypeptide" not in poly.get("type","").lower():continue
  seq=re.sub("[^A-Z]","",poly.get("pdbx_seq_one_letter_code_can","").upper())
  auth=p.get("rcsb_polymer_entity_container_identifiers",{}).get("auth_asym_ids") or []
  entities[str(ent)]={"seq_len":len(seq),"auth_chains":auth,"description":p.get("rcsb_polymer_entity",{}).get("pdbx_description","")}
 return e,a,entities
def parse(pid,assembly):
 url=f"https://files.rcsb.org/download/{pid}{'-assembly1' if assembly else ''}.cif"
 payload=get(url).content
 sha=hashlib.sha256(payload).hexdigest()
 path=out/(pid+("-assembly1" if assembly else "")+".cif")
 path.write_bytes(payload)
 st=gemmi.read_structure(str(path))
 model=st[0]
 results=[]
 for c in model:
  residues=[r for r in c.get_polymer() if r.name in STANDARD]
  if not residues:continue
  ca=sum("CA" in [a.name for a in r] for r in residues)
  miss=[]
  altloc_count=0
  for r in residues:
   atoms={a.name for a in r if a.element.name not in ("H","D")}
   required=set(SIDE[r.name].split())
   if required - atoms:miss.append({"residue":str(r.seqid),"name":r.name,"missing":sorted(required-atoms)})
   if any(a.altloc not in (chr(0)," ",".","?") for a in r):altloc_count+=1
  results.append({"chain":c.name,"residues":len(residues),"ca_present":ca,
   "missing_sidechain_residues":len(miss),"missing_sidechain_examples":miss[:8],
   "altloc_residues":altloc_count,
   "first_residue":str(residues[0].seqid) if residues else "",
   "last_residue":str(residues[-1].seqid) if residues else ""})
 return {"sha256":sha,"polymer_chain_count":len(results),"protein_residues":sum(r["residues"] for r in results),"chains":results}
rows=[]
for pid in IDS:
 try:
  e,a,entities=metadata(pid)
  ass=parse(pid,True);asym=parse(pid,False)
  refine=e.get("refine") or []
  rfree=next((r.get("ls_R_factor_R_free") for r in refine if r.get("ls_R_factor_R_free") is not None),None)
  rr=e.get("rcsb_entry_info",{}).get("resolution_combined",[])
  item={"pdb":pid,"method":[x.get("method") for x in e.get("exptl",[])],"resolution_A":rr[0] if rr else None,
    "r_free":rfree,"assembly_oligomer":a.get("pdbx_struct_assembly",{}).get("oligomeric_details"),
    "entities":entities,"assembly":ass,"asymmetric_unit":asym}
  rows.append(item)
  print("SIX_QA "+json.dumps(item,separators=(",",":")),flush=True)
 except Exception as ex:
  item={"pdb":pid,"error":str(ex)};rows.append(item);print("SIX_ERROR "+json.dumps(item),flush=True)
(out/"six_results.json").write_text(json.dumps(rows,indent=2))
if any("error" in i for i in rows):sys.exit("At least one QA entry failed")
