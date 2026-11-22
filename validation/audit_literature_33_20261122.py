"""33 literature PDB candidate audit: actual archive records + actual coordinates.
Structural, experimental evidence and holdout-independence status are DISTINCT.
Output is not an experimental adsorption validation result.
"""
from __future__ import annotations
import csv, json, hashlib, gzip, re, sys, urllib.request, time, math, traceback, argparse
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import gemmi
from Bio.Align import PairwiseAligner

ROOT=Path(__file__).resolve().parents[1]
OUT=Path("literature33_audit"); OUT.mkdir(exist_ok=True)
IDS="4F5S 2VUF 2CBA 1UBQ 1JNJ 3ECA 2HAV 1TGU 3NWV 7RSA 3LN2 1LYZ 1HRC 1MBN 2HHB 2PTN 1A5R 2Y0G 256L 1AML 2OED 2GB1 7KWI 1BEB 1F6R 1H3X 1HZH 2VH7 3GHG 3STN 5J17 6O5E 8ONG".split()
REF=list(csv.DictReader(open(ROOT/"validation"/"strict_apo_48_selected_candidates_v2_final.csv")))
DEV_IDS=[r["pdb_id"] for r in REF]
assert len(IDS)==33 and len(set(IDS))==33 and len(DEV_IDS)==48 and not (set(IDS)&set(DEV_IDS))
CANON=set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
SIDE={a:tuple(v.split()) for a,v in {
"ALA":"CB","ARG":"CB CG CD NE CZ NH1 NH2","ASN":"CB CG OD1 ND2","ASP":"CB CG OD1 OD2","CYS":"CB SG","GLN":"CB CG CD OE1 NE2","GLU":"CB CG CD OE1 OE2","GLY":"","HIS":"CB CG ND1 CD2 CE1 NE2","ILE":"CB CG1 CG2 CD1","LEU":"CB CG CD1 CD2","LYS":"CB CG CD CE NZ","MET":"CB CG SD CE","PHE":"CB CG CD1 CD2 CE1 CE2 CZ","PRO":"CB CG CD","SER":"CB OG","THR":"CB OG1 CG2","TRP":"CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2","TYR":"CB CG CD1 CD2 CE1 CE2 CZ OH","VAL":"CB CG1 CG2"}.items()}
METALS={"ZN","MN","MG","CA","FE","FE2","CU","CO","CD","NI","MO","W","K","NA","SR","HG"}
NUCLEIC={"DA","DT","DC","DG","DI","DU","A","C","G","U","I"}
EXPERIMENT_ANCHORS={"2VUF":[373,389,403,410,414,428],"2CBA":[10,37],
 "1UBQ":[2,3,15,17,18],"1JNJ":[3,30,31,33,36,37,55,56,57,61],
 "7RSA":[1,7,10,31,38,39,53,61,66,76,85,91,92,98,104,111],
 "1HRC":[22,25,27,79,86],"6O5E":[162,207,255,347]}
# Annotated anchors from previously compiled literature curation.
# Ranges use endpoints here; chemical identity and construct matching are NOT certified.
HEAD={"User-Agent":"InterfaceScout-33Candidate-ScientificAudit/2026.11"}
def fetch(url,tries=3,timeout=45):
 err=None
 for i in range(tries):
  try:
   req=urllib.request.Request(url,headers=HEAD)
   with urllib.request.urlopen(req,timeout=timeout) as f:return f.read()
  except Exception as e:
   err=e;time.sleep(1+i)
 raise RuntimeError(f"{url}: {type(err).__name__}: {err}")
def js(url):return json.loads(fetch(url))
def num(x):
 try:
  z=float(x);return z if math.isfinite(z) else None
 except (ValueError,TypeError):return None
def metadata(pid):
 e=js(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
 ents=[]
 for ent in e.get("rcsb_entry_container_identifiers",{}).get("polymer_entity_ids",[]):
  x=js(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{ent}")
  poly=x.get("entity_poly",{})
  if "polypeptide" not in str(poly.get("type","")).lower():continue
  seq=re.sub("[^ACDEFGHIKLMNPQRSTVWY]","",poly.get("pdbx_seq_one_letter_code_can","").upper())
  desc=x.get("rcsb_polymer_entity",{})
  ents.append({"entity_id":str(ent),"seq":seq,"length":len(seq),"description":desc.get("pdbx_description"),
    "mutation_annotation":desc.get("pdbx_mutation"),"auth_chains":x.get("rcsb_polymer_entity_container_identifiers",{}).get("auth_asym_ids") or []})
 info=e.get("rcsb_entry_info",{})
 return {"title":e.get("struct",{}).get("title",""),"methods":[a.get("method") for a in e.get("exptl",[])],
  "resolution_A":(info.get("resolution_combined") or [None])[0],
  "entity_count":len(ents),"entities":ents,
  "deposition_date":e.get("rcsb_accession_info",{}).get("deposit_date"),
  "release_date":e.get("rcsb_accession_info",{}).get("initial_release_date")}
def wwpdb(pid):
 url=f"https://files.rcsb.org/pub/pdb/validation_reports/{pid[1:3].lower()}/{pid.lower()}/{pid.lower()}_validation.xml.gz"
 try:
  raw=fetch(url,tries=2,timeout=25)
  if raw[:2]==b"\x1f\x8b":raw=gzip.decompress(raw)
  import xml.etree.ElementTree as ET
  root=ET.fromstring(raw)
  entry=next((x for x in root.iter() if x.tag.split("}")[-1].lower()=="entry"),None)
  if entry is None:raise ValueError("Entry validation XML missing")
  d=entry.attrib
  metrics={name:num(d.get(name)) for name in
    ["clashscore","percent-rama-outliers","percent-rota-outliers","percent-RSRZ-outliers","DCC_Rfree","PDB-Rfree"]}
  return {"available":True,"source":url,"sha256":hashlib.sha256(raw).hexdigest(),"metrics":metrics}
 except Exception as ex:return {"available":False,"source":url,"error":str(ex)}
def cath_load():
 import requests
 from collections import defaultdict
 url="https://download.cathdb.info/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt"
 try:
  raw=fetch(url,timeout=120)
  cats=defaultdict(set)
  for row in raw.decode("utf8","replace").splitlines():
   fields=row.split()
   if len(fields)>=5 and len(fields[0])>=7 and fields[0][0].isdigit():
    cats[fields[0][:4].upper()].add(".".join(fields[1:5]))
  return cats,{"source":url,"sha256":hashlib.sha256(raw).hexdigest(),"n_pdb":len(cats)}
 except Exception as ex:return defaultdict(set),{"error":str(ex),"source":url}
def coord(pid,entities):
 url=f"https://files.rcsb.org/download/{pid}-assembly1.cif"
 raw=fetch(url,tries=3,timeout=75)
 temp=OUT/(pid+"-assembly1.cif");temp.write_bytes(raw)
 st=gemmi.read_structure(str(temp))
 if not len(st):raise ValueError("missing model")
 chains=[];ligands=Counter();num_res=0;aa_names=set();present_num=set();present_name={}
 for ch in st[0]:
  seq=[r for r in ch.get_polymer() if r.name in CANON]
  if not seq:continue
  auth=ch.name.split("-")[0]
  refs=[e for e in entities if auth in e["auth_chains"]]
  declared=refs[0]["length"] if len(refs)==1 else None
  nca=sum("CA" in {a.name for a in r} for r in seq)
  side=0;alt=0;nonca=0;num_res+=len(seq)
  for residue in seq:
   atoms={a.name for a in residue if a.element.name not in ("H","D")}
   side+=int(bool(set(SIDE[residue.name])-atoms))
   alt+=int(any(a.altloc not in (chr(0)," ",".","?") for a in residue))
   nonca+=int("CA" not in atoms)
   present_num.add(residue.seqid.num)
   present_name[residue.seqid.num]=residue.name
   aa_names.add(residue.name)
  chains.append({"chain":ch.name,"auth_chain":auth,"mapped_entity":refs[0]["entity_id"] if len(refs)==1 else None,
   "declared_length":declared,"observed_canonical":len(seq),"n_CA":nca,
   "CA_coverage":round(nca/declared,4) if declared else None,
   "sidechain_missing_residues":side,"altloc_residues":alt,"CA_missing_residues":nonca})
 for ch in st[0]:
  for r in ch:
   if r.name in CANON or r.name in {"HOH","DOD","WAT"}:continue
   if len(r)==0:continue
   # Count only noncanonical NONpolymer species; don't count DNA as chemical ligand.
   if r.name in NUCLEIC:ligands["NUCLEIC:"+r.name]+=1
   elif r.name in METALS:ligands["METAL:"+r.name]+=1
   else:
    heavy=sum(a.element.name not in ("H","D") for a in r)
    if heavy>=6:ligands["NONTRIVIAL:"+r.name]+=1
 return {"sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw),"chains":chains,"n_chains":len(chains),
   "n_observed_residues":num_res,"ligand_records":dict(ligands),
   "observed_auth_residue_numbers":sorted(present_num),
   "anchor_presence":{str(x):x in present_num for x in EXPERIMENT_ANCHORS.get(pid,[])}}
def seqsimilarity(a,b,aligner):
 top=None
 for e1 in a["entities"]:
  for e2 in b["entities"]:
   x,y=e1["seq"],e2["seq"]
   if len(x)<25 or len(y)<25:continue
   try:al=aligner.align(x,y)[0]
   except Exception:continue
   v=[(int(i),int(j)) for i,j in zip(*al.indices) if i>=0 and j>=0]
   n=len(v)
   if not n:continue
   identity=sum(x[i]==y[j] for i,j in v)/n;cov=n/min(len(x),len(y))
   z={"identity":round(identity,5),"shorter_coverage":round(cov,5),"aligned":n,
      "reference_entity":e2["entity_id"],"query_entity":e1["entity_id"]}
   if top is None or identity*cov>top["identity"]*top["shorter_coverage"]:top=z
 return top
def row(pid,allmeta,cath):
 x={"pdb":pid}
 try:
  x.update(allmeta[pid]);x["archive_metadata"]="PASS"
  x["geometry"]=coord(pid,x["entities"]);x["coordinate_status"]="PASS"
  x["wwpdb"]=wwpdb(pid)
  x["CATH_H"]=sorted(cath.get(pid,[]))
  x["CATH_T"]=sorted({".".join(q.split(".")[:3]) for q in x["CATH_H"]})
  chains=x["geometry"]["chains"];sides=sum(q["sidechain_missing_residues"] for q in chains)
  covs=[q["CA_coverage"] for q in chains if q["CA_coverage"] is not None]
  cov=min(covs) if covs else None
  muts=[(e["entity_id"],e["mutation_annotation"]) for e in x["entities"] if e.get("mutation_annotation") not in (None,"?","",0)]
  ligs=list(x["geometry"]["ligand_records"])
  v=x["wwpdb"].get("metrics",{});res=num(x.get("resolution_A"))
  x["strict_Xray_development_gate"]="NOT_APPLICABLE_NON_XRAY" if x["methods"]!=["X-RAY DIFFRACTION"] else "FAIL_OR_INCOMPLETE"
  blockers=[]
  if x["methods"]==["X-RAY DIFFRACTION"]:
   if res is None or res>2.5:blockers.append("resolution_gt2.5_or_missing")
   if cov is None or cov<.98:blockers.append("CA_coverage_lt98_or_unknown")
   if sides>0:blockers.append("missing_sidechains")
   if sum(q["altloc_residues"] for q in chains):blockers.append("alternate_coordinates")
   if muts:blockers.append("mutation_annotation")
   if ligs:blockers.append("nonwater_ligand_or_hetero_partner")
   if len(chains)>4:blockers.append("assembly_gt4_chains")
   if sum(q["observed_canonical"] for q in chains)>1500:blockers.append("assembly_gt1500_aa")
   if not x["wwpdb"].get("available"):blockers.append("wwpdb_unavailable")
   if v:
    lim={"clashscore":15,"percent-rama-outliers":1,"percent-rota-outliers":5,
      "percent-RSRZ-outliers":5}
    for key,threshold in lim.items():
     if v.get(key) is None:blockers.append("wwpdb_metric_unknown_"+key)
     elif v[key]>threshold:blockers.append("wwpdb_metric_high_"+key)
    rfree=v.get("DCC_Rfree") if v.get("DCC_Rfree") is not None else v.get("PDB-Rfree")
    if rfree is None or rfree>.28:blockers.append("Rfree_gt0.28_or_missing")
   x["strict_Xray_development_gate"]="PASS" if not blockers else "FAIL_OR_INCOMPLETE"
  if not chains:blockers.append("no_canonical_protein_chains")
  x["strict_blockers"]=blockers
  if x["methods"]!=["X-RAY DIFFRACTION"]:
   x["application_structural_status"]="CONDITIONAL_NONXRAY_REQUIRES_METHOD_QA"
  elif not chains or cov is None or cov<.90:
   x["application_structural_status"]="HOLD_COVERAGE_OR_CHAIN"
  elif sides>0 or ligs or muts:
   x["application_structural_status"]="CONDITIONAL_CONSTRUCT_OR_CHEMISTRY"
  else:x["application_structural_status"]="STRUCTURALLY_USABLE_WITH_CONTEXT_REVIEW"
  x["sidechain_missing_count"]=sides;x["min_CA_coverage"]=cov
  x["mutation_annotations"]=muts
 except Exception as ex:
  x["archive_metadata"]="PASS" if pid in allmeta else "FAIL"
  x["coordinate_status"]="FAIL"
  x["error"]=f"{type(ex).__name__}: {ex}"
  x["application_structural_status"]="INCOMPLETE_FETCH_OR_PARSE"
 return x
def main():
 args=argparse.ArgumentParser();args.add_argument("--shard",type=int,default=0)
 args.add_argument("--shards",type=int,default=6);args= args.parse_args()
 selected=IDS[args.shard::args.shards]
 cats,catinfo=cath_load()
 # Fetch reference and panel metadata; source-metadata is immutable within this CI run.
 pool_ids=list(dict.fromkeys(DEV_IDS+IDS))
 allmeta={};errors={}
 with ThreadPoolExecutor(max_workers=12) as pool:
  futures={pool.submit(metadata,pid):pid for pid in pool_ids}
  for future in as_completed(futures):
   pid=futures[future]
   try:allmeta[pid]=future.result()
   except Exception as e:errors[pid]=str(e)
 aligner=PairwiseAligner(mode="local",match_score=2,mismatch_score=-1,open_gap_score=-8,extend_gap_score=-1)
 output=[]
 for pid in selected:
  if pid in allmeta:r=row(pid,allmeta,cats)
  else:r={"pdb":pid,"archive_metadata":"FAIL","coordinate_status":"NOT_RUN","error":errors.get(pid,"metadata_missing"),"application_structural_status":"INCOMPLETE_FETCH_OR_PARSE"}
  conflicts=[];unresolved=[];top=[]
  if pid in allmeta:
   for other in DEV_IDS:
    if other not in allmeta:
     unresolved.append(other);continue
    matched=seqsimilarity(allmeta[pid],allmeta[other],aligner)
    if matched:
     if matched["identity"]>=.30 and matched["shorter_coverage"]>=.70 and matched["aligned"]>=40:
      conflicts.append({"pdb":other,**matched})
     top.append({"pdb":other,**matched})
  top.sort(key=lambda z:(z["identity"]*z["shorter_coverage"],z["aligned"]),reverse=True)
  ch=set(r.get("CATH_H",[])); ct=set(r.get("CATH_T",[]))
  overlap_h=[{"pdb":q["pdb_id"],"CATH_H":sorted(ch & set(q["CATH_H"].split(";")))} for q in REF
    if ch & set(q["CATH_H"].split(";"))]
  overlap_t=[{"pdb":q["pdb_id"],"CATH_T":sorted(ct & {".".join(z.split(".")[:3]) for z in q["CATH_H"].split(";") if z})} for q in REF
    if ct & {".".join(z.split(".")[:3]) for z in q["CATH_H"].split(";") if z}]
  r["sequence_reference_complete"]=not unresolved and pid in allmeta
  r["dev48_seq_conflicts"]=conflicts
  r["closest_dev48_sequence"]=top[:3]
  r["dev48_CATH_H_conflicts"]=overlap_h
  r["dev48_CATH_T_conflicts"]=overlap_t
  r["dev48_independence_screen"]="NOT_VERIFIED" if not r["sequence_reference_complete"] or not ch else ("POTENTIAL_HOMOLOGY" if conflicts or overlap_h else "NO_SCREENED_CONFLICT")
  r["exact_PDB_overlap_dev48"]=pid in DEV_IDS
  r["original_23_or_new_10"]="original23" if pid in IDS[:23] else "new10"
  output.append(r)
  print("AUDIT",pid,r["archive_metadata"],r["coordinate_status"],r.get("strict_Xray_development_gate"),
    r["application_structural_status"],r["dev48_independence_screen"],flush=True)
 folder=OUT/f"shard_{args.shard}";folder.mkdir(exist_ok=True)
 with (folder/"records.json").open("w") as f:json.dump(output,f,indent=2,ensure_ascii=False)
 fields=["pdb","original_23_or_new_10","title","methods","resolution_A","entity_count",
  "archive_metadata","coordinate_status","strict_Xray_development_gate","application_structural_status",
  "min_CA_coverage","sidechain_missing_count","mutation_annotations","strict_blockers",
  "dev48_independence_screen","dev48_seq_conflicts","dev48_CATH_H_conflicts","dev48_CATH_T_conflicts",
  "closest_dev48_sequence","exact_PDB_overlap_dev48","CATH_H","error"]
 with (folder/"records.csv").open("w",newline="",encoding="utf8") as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for r in output:
   w.writerow({k:json.dumps(r.get(k),ensure_ascii=False) if isinstance(r.get(k),(dict,list)) else r.get(k,"") for k in fields})
 (folder/"methods.json").write_text(json.dumps({"source_ids":IDS,"reference_dev48":DEV_IDS,
  "screen_30pct":">=30% aligned identity, >=70% coverage of shorter, >=40 aligned aa",
  "CATH_info":catinfo,"n_references_unavailable":len([k for k in DEV_IDS if k not in allmeta]),
  "failed_metadata":errors,"strict_gate":"preexisting 48-PDB development x-ray criteria; not compulsory for application",
  "no_repair":True,"no_experimental_ground_truth_from_MD":True},indent=2))
 if len(output)!=len(selected):raise RuntimeError("incomplete shard")
 if sum(r["coordinate_status"]=="PASS" for r in output)==0:raise RuntimeError("No coordinates audited; infrastructure failure")
if __name__=="__main__":main()
