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
data=get(G+"Dataset-for-protein-freq-analysis-v2.csv").decode("utf-8-sig")
ds=list(csv.DictReader(io.StringIO(data)));header=list(ds[0]) if ds else []
header_map={h.strip().upper():h for h in header}
# A few protein identifiers may have spaces in headers; exact normalized ID only.
summaryrows=[];sample_rows=[]
for r in meta:
 accessions=[x for x in r["uniprot_ids"].split(";") if x]
 matches=[x for x in accessions if x in protein_ids]
 columns=[(x,header_map[x]) for x in matches if x in header_map]
 nonzero=0
 for accession,column in columns:
  for datum in ds:
   raw=str(datum.get(column,"")).strip()
   try:v=float(raw)
   except (TypeError,ValueError):continue
   if v<=0:continue
   nonzero+=1
   sample_rows.append({"pdb_id":r["pdb_id"],"chain":r["chain"],"uniprot_accession":accession,
     "corona_measurement":v,
     "np_entry_id":datum.get("NP Entry ID",datum.get("NP Entry No","")),
     "study_id":datum.get("Study ID",""),"doi":datum.get("DOI",""),
     "material":datum.get("Core NP Material",""),
     "surface_modification":datum.get("Surface Modification",""),
     "zeta_mV":datum.get("Zeta Potential (mV)",""),
     "species":datum.get("Species",""),
     "source_dataset":"choulab210/Meta-Analysis-of-Nanoparticle-and-Protein-Corona-Interactions-in-Biological-Media",
     "evidence_type":"CORONA_PROTEOMICS_ABUNDANCE_NOT_CONTACT_SITE"})
 summaryrows.append({**r,"uniprot_ID_matched":int(bool(matches)),"matched_uniprot":";".join(matches),
   "present_as_proteomics_column":int(bool(columns)),"measured_nonzero_samples":nonzero,
   "evidence_status":"PROTEIN_CORONA_ABUNDANCE_MATCH" if nonzero else ("UNIPROT_ID_PRESENT_WITHOUT_POSITIVE_SAMPLE" if matches else "NO_UNIPROT_MATCH_IN_PCDB"),
   "contact_residue_ground_truth":0})
save("pcdb_exact_uniprot_299.csv",summaryrows)
save("pcdb_matched_samples.csv",sample_rows or [{"pdb_id":"","uniprot_accession":"","np_entry_id":"","study_id":"","doi":"","material":"","corona_measurement":"","evidence_type":"NO_POSITIVE_MATCH"}])
counts={"atlas_proteins":len(meta),"rcsb_metadata_ok":sum(x["metadata_status"]=="OK" for x in meta),"atlas_with_any_uniprot":sum(bool(x["uniprot_ids"]) for x in meta),
"protein_ids_in_dataset":len(protein_ids),"dataset_rows":len(ds),"protein_id_columns":sum(k in protein_ids for k in header_map),
"exact_uniprot_matches":sum(x["uniprot_ID_matched"] for x in summaryrows),
"positive_proteomics_matches":sum(x["measured_nonzero_samples"]>0 for x in summaryrows),
"positive_sample_rows":len(sample_rows),"source":"https://github.com/choulab210/Meta-Analysis-of-Nanoparticle-and-Protein-Corona-Interactions-in-Biological-Media",
"scope":"Protein-level proteomic detection, never residue-level adsorption orientation or contact site"}
(OUT/"summary.json").write_text(json.dumps(counts,indent=2))
print("PCDB_CROSSMATCH_COMPLETE",json.dumps(counts),flush=True)
