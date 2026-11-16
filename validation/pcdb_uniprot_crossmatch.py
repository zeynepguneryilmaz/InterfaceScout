import csv,io,json,urllib.request,urllib.parse,urllib.error,time,zipfile,collections,hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
BASE=Path(__file__).resolve().parents[1];OUT=BASE/"validation/pcdb_crossmatch";OUT.mkdir(parents=True,exist_ok=True)
G="https://raw.githubusercontent.com/choulab210/Meta-Analysis-of-Nanoparticle-and-Protein-Corona-Interactions-in-Biological-Media/main/Datasets/"
def get(url):
 for i in range(4):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"InterfaceScout academic evidence audit","Accept":"application/json"})
   with urllib.request.urlopen(req,timeout=60) as r:return r.read()
  except Exception:
   if i==3:raise
   time.sleep(2*(i+1))
def lookup(row):
 pid,entity=row["pdb_id"],row["entity_id"]
 try:
  obj=json.loads(get("https://data.rcsb.org/rest/v1/core/polymer_entity/"+pid+"/"+entity))
  ent=obj.get("rcsb_polymer_entity",{}) or {}
  acc=(obj.get("rcsb_polymer_entity_container_identifiers",{}) or {}).get("uniprot_ids",[]) or []
  if isinstance(acc,str):acc=[acc]
  return {"pdb_id":pid,"chain":row["chain"],"entity_id":entity,"protein_name":ent.get("pdbx_description",""),"uniprot_ids":";".join(sorted(set(x.strip().upper() for x in acc if x))),"metadata_status":"OK"}
 except Exception as e:
  return {"pdb_id":pid,"chain":row["chain"],"entity_id":entity,"protein_name":"","uniprot_ids":"","metadata_status":"ERROR:"+str(e)[:150]}
def save(name,rows):
 if not rows:return
 with (OUT/name).open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
with zipfile.ZipFile(BASE/"validation/input_300_proteins.zip") as z:
 rows=[r for r in csv.DictReader(io.StringIO(z.read("manifest.csv").decode("utf-8-sig"))) if r["pdb_id"]!="3IYK"]
meta=[]
with ThreadPoolExecutor(max_workers=4) as pool:
 for future in as_completed([pool.submit(lookup,r) for r in rows]):
  meta.append(future.result())
  if len(meta)%50==0:print("RCSB_META_PROGRESS",len(meta),flush=True)
meta.sort(key=lambda x:x["pdb_id"]);save("pdb_uniprot_crosswalk.csv",meta)
proteins_csv=get(G+"proteinIDs.csv").decode("utf-8-sig")
protein_ids={row.strip().upper() for row in proteins_csv.splitlines()[1:] if row.strip()}
# The full PC-DB has a 2,497-ID accession roster and a 1,513-protein wide
# human-only matrix in which columns are protein NAMES, NOT UniProt IDs.
# Both independent pieces of evidence must match; never treat an ID-list hit
# alone as evidence of a sample with positive adsorbed abundance.
data=get(G+"Dataset-for-protein-freq-analysis-v2.csv").decode("utf-8-sig")
ds=list(csv.DictReader(io.StringIO(data)));header=list(ds[0]) if ds else []
pcdb_meta=list(csv.DictReader(io.StringIO(get(G+"PC-DB_for_Meta-Analysis.csv").decode("utf-8-sig"))))
sample_metadata={x.get("NP Entry ID"):x for x in pcdb_meta}
def norm_name(s):
 s=str(s or "").lower().split(",")[0]
 s=s.replace("–","-").replace("—","-").replace("−","-")
 s=s.replace("protein,","protein")
 return " ".join(re.sub(r"[^a-z0-9]+"," ",s).split())
header_map={}
for h in header[2:]:
 header_map.setdefault(norm_name(h),[]).append(h)
summaryrows=[];sample_rows=[]
for r in meta:
 accessions=[x for x in r["uniprot_ids"].split(";") if x]
 matches=[x for x in accessions if x in protein_ids]
 # Only exact normalized name matches; no fuzzy matches or unverified aliases.
 columns=header_map.get(norm_name(r["protein_name"]),[]) if matches else []
 nonzero=0
 for column in columns:
  for datum in ds:
   raw=str(datum.get(column,"")).strip()
   try:v=float(raw)
   except (TypeError,ValueError):continue
   if not (v>0):continue
   nonzero+=1
   sample=sample_metadata.get(datum.get("NP Entry ID"),{})
   sample_rows.append({"pdb_id":r["pdb_id"],"chain":r["chain"],"uniprot_accession":";".join(matches),
     "matched_protein_name":column,"corona_measurement":v,
     "np_entry_id":datum.get("NP Entry ID",""),
     "study_id":datum.get("Study ID",""),"doi":sample.get("DOI",""),
     "material":sample.get("Core NP Material",""),
     "surface_modification":sample.get("Surface Modification",""),
     "zeta_mV":sample.get("Zeta Potential (mV)",""),
     "species":sample.get("Species",""),
     "source_dataset":"choulab210/Meta-Analysis-of-Nanoparticle-and-Protein-Corona-Interactions-in-Biological-Media",
     "evidence_type":"POTENTIAL_CORONA_PROTEOMICS_ABUNDANCE_NOT_CONTACT_SITE",
     "review_requirement":"Confirm UniProt accession to name mapping and species in PC-DB; broad protein name is not unique"})
 summaryrows.append({**r,"uniprot_ID_matched":int(bool(matches)),"matched_uniprot":";".join(matches),
   "present_as_proteomics_column":int(bool(columns)),"matched_column_names":";".join(columns),
   "measured_nonzero_samples":nonzero,
   "evidence_status":"PROVISIONAL_CORONA_ABUNDANCE_MATCH" if nonzero else ("UNIPROT_ID_PRESENT_ONLY" if matches else "NO_UNIPROT_MATCH_IN_PCDB"),
   "contact_residue_ground_truth":0})
save("pcdb_exact_uniprot_299.csv",summaryrows)
save("pcdb_matched_samples.csv",sample_rows or [{"pdb_id":"","uniprot_accession":"","np_entry_id":"","study_id":"","doi":"","material":"","corona_measurement":"","evidence_type":"NO_POSITIVE_MATCH"}])
counts={"atlas_proteins":len(meta),"rcsb_metadata_ok":sum(x["metadata_status"]=="OK" for x in meta),"atlas_with_any_uniprot":sum(bool(x["uniprot_ids"]) for x in meta),
"protein_ids_in_dataset":len(protein_ids),"dataset_rows":len(ds),"protein_name_columns":len(header)-2,
"exact_uniprot_matches":sum(x["uniprot_ID_matched"] for x in summaryrows),"name_columns_matched":sum(x["present_as_proteomics_column"] for x in summaryrows),
"positive_proteomics_matches":sum(x["measured_nonzero_samples"]>0 for x in summaryrows),
"positive_sample_rows":len(sample_rows),"source":"https://github.com/choulab210/Meta-Analysis-of-Nanoparticle-and-Protein-Corona-Interactions-in-Biological-Media",
"scope":"Protein-level proteomic detection, never residue-level adsorption orientation or contact site"}
(OUT/"summary.json").write_text(json.dumps(counts,indent=2))
print("PCDB_CROSSMATCH_COMPLETE",json.dumps(counts),flush=True)
