import csv,io,zipfile,urllib.request,json,re,html,sys
from pathlib import Path
from bs4 import BeautifulSoup
R=Path(__file__).resolve().parents[1];O=R/"validation/bad_crossmatch";O.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(R/"validation/input_300_proteins.zip") as z:
 panel=[r for r in csv.DictReader(io.StringIO(z.read("manifest.csv").decode("utf-8-sig"))) if r["pdb_id"]!="3IYK"]
ids={r["pdb_id"] for r in panel}
url_candidates=["https://molecularsense.com/bad-2/","https://bionanoinfo.com/bad/"]
allrows=[];logs=[]
for url in url_candidates:
 try:
  req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (InterfaceScout academic audit)"})
  with urllib.request.urlopen(req,timeout=80) as r:s=r.read().decode("utf-8","replace")
  soup=BeautifulSoup(s,"html.parser");tables=soup.find_all("table")
  counts=[]
  for t in tables:
   head=[x.get_text(" ",strip=True) for x in t.find_all("th")]
   if not any("PDB" in x for x in head):continue
   for tr in t.find_all("tr"):
    cells=tr.find_all("td")
    if not cells:continue
    vals=[c.get_text(" ",strip=True) for c in cells]
    if len(vals)!=len(head):continue
    row=dict(zip(head,vals))
    pdbval=row.get("PDB","").strip().upper()
    # direct PDB code reference (not sequence identity yet).
    allrows.append({"source_page":url,**row,"pdb_normalized":pdbval,"exact_pdb_in_atlas":int(pdbval in ids)})
   counts.append({"head":head[:6],"rows":len(t.find_all("tr"))})
  logs.append({"source":url,"html_bytes":len(s),"tables":len(tables),"matches":counts})
 except Exception as e:logs.append({"source":url,"error":repr(e)})
def csvout(p,rows):
 if rows:
  cols=list(dict.fromkeys(k for r in rows for k in r))
  with open(p,"w",newline="",encoding="utf8") as f:
   w=csv.DictWriter(f,cols);w.writeheader();w.writerows(rows)
csvout(O/"bad_extracted_records.csv",allrows)
match=[r for r in allrows if r["exact_pdb_in_atlas"]]
csvout(O/"bad_exact_pdb_matches.csv",match or [{"source_page":"","PDB":"","Protein":"","DOI":"","exact_pdb_in_atlas":0}])
summary={"db_rows_parsed":len(allrows),"exact_pdb_match_rows":len(match),"unique_atlas_proteins_exact_pdb":len({r["pdb_normalized"] for r in match}),"source_logs":logs,"limitations":"An exact PDB match is stronger identity evidence than text matching, but these records provide adsorption amount, not protein-side contact residues. Repeated BAD sources may duplicate rows."}
(O/"summary.json").write_text(json.dumps(summary,indent=2))
print("BAD_CROSSMATCH_COMPLETE",json.dumps(summary),flush=True)
