"""Auditable corpus-screen of all 299 frozen InterfaceScout PDB entries.
Important: a literature hit is a CANDIDATE requiring protein and surface verification, not ground truth.
"""
import csv,io,json,re,sys,time,random,urllib.request,urllib.parse,urllib.error,zipfile,collections,html
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"validation"/"literature_299";OUT.mkdir(parents=True,exist_ok=True)
ARCHIVE=ROOT/"validation/input_300_proteins.zip"
HEAD={"User-Agent":"InterfaceScout-literature-audit/1.0 (academic reproducibility)","Accept":"application/json"}
DATE="2026-11-17"
def api(url,tries=5):
 for attempt in range(tries):
  try:
   req=urllib.request.Request(url,headers=HEAD)
   with urllib.request.urlopen(req,timeout=38) as r:return json.load(r)
  except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError) as e:
   if isinstance(e,urllib.error.HTTPError) and e.code==404:return {}
   if isinstance(e,urllib.error.HTTPError) and e.code in (400,401,403):raise
   if attempt==tries-1:raise
   time.sleep(min(3*(attempt+1),15))
def clean(s):
 return re.sub(r"\s+"," ",html.unescape(re.sub("<[^>]+>"," ",s or ""))).strip()
def names_for(entry):
 pid,entity=entry["pdb_id"],entry["entity_id"]
 obj=api("https://data.rcsb.org/rest/v1/core/polymer_entity/"+pid+"/"+entity)
 desc=(obj.get("rcsb_polymer_entity",{}) or {}).get("pdbx_description","")
 # Additional entity names may be provided as comma-separated PDB molecule descriptions.
 if not desc:desc="Unknown / lookup incomplete"
 acc=(obj.get("rcsb_polymer_entity",{}) or {}).get("rcsb_uniprot_accession",[]) or []
 if not isinstance(acc,list):acc=[acc]
 uniprot=[str(s) for s in acc if s]
 return {"protein_name":clean(desc),"uniprot_accession":";".join(uniprot),"metadata_resolved":bool(obj)}
material_pat=re.compile(r"\b(silica|silicon\s+dioxide|glass|quartz|gold|silver|platinum|graphen\w*|carbon\s+nanotube\w*|carbon\s+black|polystyrene|polymer\w*|polyethylene|polyurethane|polymeric|hydroxyapatite|alumina|aluminum\s+hydroxide|magnetite|iron\s+oxide|metal\s+oxide|solid\s+surface|material\s+surface|solid[-\s]liquid\s+interface|silicon\s+nitride|titanium\s+dioxide|titania|nanoparticl\w*|nanomaterial\w*|biochar|nanotube\w*|substrate|biomaterial\w*|nanofib\w*|nanocarrier\w*|self[-\s]assembled\s+monolayer|surface[-\s]bound|immobiliz\w*|immobilis\w*|electrode\w*|plasmonic)\b",re.I)
ads_pat=re.compile(r"\b(adsorp\w*|immobiliz\w*|immobilis\w*|surface[-\s]bound|surface[-\s]induced|protein\s+corona|nanoparticle[-\s]protein|protein[-\s]nanoparticle|contact\s+residue\w*|orientation|interfacial\s+interactions?)\b",re.I)
md_pat=re.compile(r"\b(molecular\s+dynamics|atomistic\s+simulat\w*|coarse[-\s]grain\w*|metadynamics|umbrella\s+sampling|monte\s+carlo|steered\s+molecular\s+dynamics|enhanced\s+sampling|molecular\s+simulation\w*|computational\s+adsorption)\b",re.I)
contact_pat=re.compile(r"\b(contact\s+residu\w*|binding\s+orientat\w*|residue[-\s]level|NMR\s+chemical\s+shift|HDX[-\s]MS|hydrogen[-\s]deuterium|cross[-\s]link\w*|interface\s+mapping|adsorption\s+site|adsorption\s+orientat\w*|protein\s+orientation)\b",re.I)
material_cats={"silica":r"silica|silicon dioxide|glass|quartz","gold":r"gold|AuNP|gold nanoparticle","graphitic":r"graphen|carbon nanotube|graphitic|carbon black|biochar","polymer":r"polystyrene|polymer|polyethylene|polyurethane|PMMA|PVC|nanofib","metal_oxide":r"oxide|titania|magnetite|Fe3O4|alumina","calcium_phosphate":r"hydroxyapatite|calcium phosphate","generic_surface":r"nanoparticl|solid surface|material surface|substrate|immobili"}
def material_type(s):
 for k,pattern in material_cats.items():
  if re.search(pattern,s,re.I):return k
 return "not_identified"
def epmc(query,page=10):
 url="https://www.ebi.ac.uk/europepmc/webservices/rest/search?"+urllib.parse.urlencode({"query":query,"format":"json","resultType":"core","pageSize":page})
 obj=api(url)
 return int(obj.get("hitCount",0)),obj.get("resultList",{}).get("result",[])
def qphrase(s):
 s=re.sub(r"[^a-zA-Z0-9 /-]"," ",s)
 s=clean(s);return s[:75]
def cands_for(entry,meta):
 pid,name=entry["pdb_id"],meta["protein_name"]
 # Precise query for existing PDB identifier (usually low recall).
 query1='TITLE_ABS:"'+pid+'" AND (adsorp* OR immobiliz* OR nanoparticl* OR surface* OR "molecular dynamics")'
 # Protein name query is intentionally BOUNDED, not generic "protein".
 invalid=not name or name.lower() in {"unknown / lookup incomplete","uncharacterized protein","hypothetical protein","protein","unknown"} or len(name)<6
 query2=('TITLE_ABS:"'+qphrase(name)+'" AND (adsorp* OR immobiliz* OR nanoparticl* OR "molecular dynamics" OR nanomaterial* OR surface*)') if not invalid else None
 searches=[]
 for lab,q in [("PDB_ID",query1),("PROTEIN_NAME",query2)]:
  if not q:continue
  try:
   count,results=epmc(q,page=8)
   searches.append((lab,q,count,results,None))
  except Exception as e:searches.append((lab,q,0,[],repr(e)))
 return searches
def work(row):
 try:
  meta=names_for(row)
  searches=cands_for(row,meta)
  return row,meta,searches,None
 except Exception as e:return row,{"protein_name":"UNRESOLVED","uniprot_accession":"","metadata_resolved":False},[],repr(e)
def csvsave(name,rows):
 if not rows:return
 with (OUT/name).open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
with zipfile.ZipFile(ARCHIVE) as z:
 original=list(csv.DictReader(io.StringIO(z.read("manifest.csv").decode("utf-8-sig"))))
 rows=[r for r in original if r["pdb_id"]!="3IYK"]
if len(rows)!=299:raise RuntimeError("Not 299 frozen structures")
master=[];article_records=[];query_records=[];exceptions=[];done=0
# modest parallelism, with retries for throttled requests
with ThreadPoolExecutor(max_workers=4) as pool:
 futures=[pool.submit(work,row) for row in rows]
 for fut in as_completed(futures):
  row,meta,searches,error=fut.result()
  pid=row["pdb_id"];seen=set();unique=[];candidate_count=0;direct_contact_count=0
  for method,q,hitcount,results,query_error in searches:
   query_records.append({"pdb_id":pid,"protein_name":meta["protein_name"],"method":method,"query":q,"hit_count":hitcount,"returned":len(results),"error":query_error or ""})
   if query_error:exceptions.append({"pdb_id":pid,"reason":"QUERY_FAILED","detail":str(query_error)[:240]})
   for a in results:
    identifier=a.get("doi") or a.get("id") or a.get("pmid")
    if not identifier or identifier in seen:continue
    seen.add(identifier)
    title=clean(a.get("title",""));abstract=clean(a.get("abstractText",""))
    material=bool(material_pat.search(title+" "+abstract))
    ads=bool(ads_pat.search(title+" "+abstract))
    md=bool(md_pat.search(title+" "+abstract))
    contact=bool(contact_pat.search(title+" "+abstract))
    type_hint="MD_or_simulation" if md else ("adsorption_experiment_candidate" if material and ads else "general_material_article")
    # Candidate status ONLY; no experimental localization can be verified from title/abstract alone.
    candidate=material and (ads or md)
    if candidate:candidate_count+=1
    if candidate and contact:direct_contact_count+=1
    evidence_label="POTENTIAL_CONTACT_METHOD_MENTIONED" if candidate and contact else ("POTENTIAL_SURFACE_STUDY" if candidate else "NOT_SURFACE_RELEVANT_FROM_ABSTRACT")
    article_records.append({
      "pdb_id":pid,"chain":row["chain"],"protein_name":meta["protein_name"],"uniprot_accession":meta["uniprot_accession"],
      "doi":a.get("doi",""),"pmid":a.get("pmid",""),"pmcid":a.get("pmcid",""),
      "article_title":title,"year":a.get("pubYear",""),"journal":clean(a.get("journalTitle","")),
      "search_method":method,"material_class_hint":material_type(title+" "+abstract),
      "md_terms":int(md),"adsorption_terms":int(ads),"material_terms":int(material),"contact_terms":int(contact),
      "candidate_for_review":int(candidate),"study_type_hint":type_hint,"evidence_label":evidence_label,
      "protein_identity_verified":0,"material_interface_verified":0,"contact_residues_extracted":0,
      "screening_status":"AUTOMATED_CANDIDATE_NOT_MANUALLY_VERIFIED",
      "abstract_excerpt":abstract[:650],"source_url":("https://doi.org/"+a["doi"] if a.get("doi") else ("https://europepmc.org/article/MED/"+str(a.get("pmid")) if a.get("pmid") else "")),
    })
  if error:exceptions.append({"pdb_id":pid,"reason":"METADATA_FAILED","detail":str(error)[:240]})
  master.append({"pdb_id":pid,"chain":row["chain"],"entity_id":row["entity_id"],"n_residues":row["n_residues"],
    "protein_name":meta["protein_name"],"uniprot_accession":meta["uniprot_accession"],"metadata_resolved":int(meta["metadata_resolved"]),
    "searches_attempted":len(searches),"retrieved_article_count":len(seen),"potential_material_adsorption_articles":candidate_count,
    "contact_related_terms_in_abstract":direct_contact_count,
    "manual_protein_verification":"PENDING","manual_material_verification":"PENDING",
    "reliable_interface_ground_truth":"NOT_ESTABLISHED",
    "literature_screen_result":("CANDIDATES_REQUIRING_FULLTEXT_REVIEW" if candidate_count else "NO_QUALIFYING_CANDIDATE_IN_AUTOMATED_QUERIES"),
    "no_publication_claim":"No candidate found is NOT evidence that no paper exists",
    "source":"RCSB PDB Data API + Europe PMC REST search", "search_date":DATE})
  done+=1
  if done%20==0:
   print("PROGRESS",done,"/299","candidate_rows",sum(int(x["candidate_for_review"]) for x in article_records),"errors",len(exceptions),flush=True)
   csvsave("protein_screen_299.csv",sorted(master,key=lambda x:x["pdb_id"]))
   csvsave("candidate_articles_raw.csv",article_records)
csvsave("protein_screen_299.csv",sorted(master,key=lambda x:x["pdb_id"]))
csvsave("candidate_articles_raw.csv",sorted(article_records,key=lambda x:(x["pdb_id"],-x["candidate_for_review"],x["year"]),reverse=False))
csvsave("query_audit.csv",sorted(query_records,key=lambda x:x["pdb_id"]))
csvsave("api_errors.csv",exceptions or [{"pdb_id":"","reason":"NONE","detail":""}])
summary={
  "structures_screened":len(master),"structure_metadata_resolved":sum(x["metadata_resolved"] for x in master),
  "total_abstract_candidates":sum(x["candidate_for_review"] for x in article_records),
  "proteins_with_candidates":sum(x["potential_material_adsorption_articles"]>0 for x in master),
  "proteins_without_automated_candidates":sum(x["potential_material_adsorption_articles"]==0 for x in master),
  "total_retrieved_articles":len(article_records),"query_calls_attempted":len(query_records),
  "query_errors":sum(bool(x["error"]) for x in query_records),
  "confirmed_independent_protein_material_interfaces":0,
  "methodology":"RCSB polymer entity description; Europe PMC PDB-ID and name title/abstract queries; automated material + adsorption/MD keyword triage",
  "very_important_limitations":[
    "This is a title/abstract machine screening, not an exhaustive fulltext systematic review",
    "PDB codes often do not appear in articles about the same protein",
    "Molecular species identity and material interface must be verified independently",
    "MD is not independent experimental ground truth; must check whether its force field or starting structure leads to circular validation",
    "No negative evidence from absent literature results",
    "Direct contact residues are not extracted or confirmed by keyword detection",
    "Publication and protein synonyms not present in RCSB may cause false negatives",
  ]}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False))
(OUT/"README.md").write_text("""# InterfaceScout: 299-protein literature evidence inventory
The title/abstract search is a **high-recall first pass** and every candidate **requires manual verification**.
protein_screen_299.csv: one row per original atlas protein.
candidate_articles_raw.csv: article-level candidates with DOI, PMID, source URL, matched name and triage terms.
query_audit.csv: exact API queries and hit counts, including failures.
api_errors.csv: connection and metadata issues.
summary.json: honest coverage.
No automatic hits are promoted to verified experimental contact sites; a simulation of a related protein is not equivalent to experimental binding-interface data.
Follow-up: manually inspect full texts, verify identity via UniProt sequence and species, check protein–material system, and extract contact-residue definitions with simulation trajectories/experiments, pH, model, and DOI.
""")
print("LITERATURE_SCREEN_COMPLETE",json.dumps(summary,ensure_ascii=False),flush=True)
