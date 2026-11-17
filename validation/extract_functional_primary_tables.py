"""Retrieve original open-access publisher-indexed tables without interpolation or OCR.
Audit source, all table labels, rows, footnotes; no auto-claimed activity ranking.
"""
from pathlib import Path
import urllib.request,urllib.error,time,json,xml.etree.ElementTree as ET,csv,re
R=Path(__file__).resolve().parents[1];O=R/"validation"/"independent_md_functional"/"source_tables";O.mkdir(parents=True,exist_ok=True)
cases={
"PMC4270395":"RNaseA_fused_glass_HDPE_PMMA_activity_2014",
"PMC4312349":"Lysozyme_fused_glass_HDPE_PMMA_activity_2015",
"PMC12839133":"GB1_latex_mutagenesis_2026",
"PMC10843815":"R2ab_polystyrene_HDX_NMR_2024"}
def get(url):
 for t in range(6):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"InterfaceScout academic open access provenance checker","Accept":"application/xml"})
   with urllib.request.urlopen(req,timeout=60) as res:return res.read()
  except (urllib.error.URLError,TimeoutError) as e:
   if t==5:raise
   time.sleep(2*(t+1))
def node_text(node):
 return " ".join(" ".join(node.itertext()).split()) if node is not None else ""
rows=[];log=[]
for pmcid,topic in cases.items():
 url="https://www.ebi.ac.uk/europepmc/webservices/rest/"+pmcid+"/fullTextXML"
 try:
  raw=get(url);(O/(pmcid+".xml")).write_bytes(raw)
  root=ET.fromstring(raw)
  article_title=node_text(root.find(".//article-title"))
  doi=""
  for aid in root.findall(".//article-meta/article-id"):
   if aid.get("pub-id-type")=="doi":doi=node_text(aid);break
  tables=root.findall(".//table-wrap")
  for idx,t in enumerate(tables,1):
   tid=t.get("id","");label=node_text(t.find("label"));caption=node_text(t.find("caption"))
   footnotes=" | ".join(node_text(q) for q in t.findall(".//table-wrap-foot/fn"))
   tr=t.findall(".//tr")
   if not tr:tr=t.findall(".//table/row")
   for ri,row in enumerate(tr):
    cells=[node_text(x) for x in list(row) if x.tag.rsplit("}",1)[-1] in ("th","td","entry")]
    if not cells:continue
    rows.append({"pmcid":pmcid,"topic":topic,"doi":doi,"title":article_title,"table_id":tid,"table_number":idx,"table_label":label,
     "table_caption":caption,"table_footnotes":footnotes,"row_index":ri,"columns_json":json.dumps(cells,ensure_ascii=False),
     "source_url":url,"verification_status":"DIRECT_EUROPE_PMC_FULLTEXT_XML"})
  log.append({"pmcid":pmcid,"doi":doi,"title":article_title,"xml_bytes":len(raw),"n_tables":len(tables),
    "n_rows":sum(x["pmcid"]==pmcid for x in rows),"status":"OK"})
  print("ARTICLE_EXTRACTED",pmcid,"tables",len(tables),"rows",log[-1]["n_rows"],flush=True)
 except Exception as ex:
  log.append({"pmcid":pmcid,"status":"ERROR","reason":str(ex)})
  print("ARTICLE_ERROR",pmcid,str(ex),flush=True)
with (O/"all_tables_rows.csv").open("w",newline="",encoding="utf8") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ["pmcid","status"]);w.writeheader();w.writerows(rows)
(O/"download_audit.json").write_text(json.dumps(log,indent=2))
print("FULLTEXT_TABLE_AUDIT_DONE",json.dumps(log),flush=True)
