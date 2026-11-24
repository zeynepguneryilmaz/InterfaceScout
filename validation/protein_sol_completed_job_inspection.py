"""Inspect one known public Protein-sol result to verify actual data path.

No new jobs are submitted. Result URLs were created by a prior authorized
submission of public ISV76 PDBs (2OED and 1HEL).
"""
import json,re,time
from pathlib import Path
import requests
from bs4 import BeautifulSoup
BASE="https://protein-sol.manchester.ac.uk"
JOBS={"2OED":"6eee809aa274d2adaae4","1HEL":"e737fa0b17ba6f4e36b1"}
OUT=Path("external_protein_sol_inspect");OUT.mkdir(exist_ok=True)
results=[]
with requests.Session() as s:
    s.headers["User-Agent"]="InterfaceScout methods-validation/1.0"
    for pid,jobid in JOBS.items():
        url=f"{BASE}/results/patches/run-{jobid}/results.html"
        item={"PDB":pid,"url":url}
        try:
            rsp=s.get(url,timeout=(15,45))
            html=rsp.text
            (OUT/(pid+"_results.html")).write_text(html,encoding="utf-8")
            soup=BeautifulSoup(html,"html.parser")
            item.update({"http_status":rsp.status_code,"bytes":len(rsp.content),
                "title":soup.title.get_text(" ",strip=True) if soup.title else "",
                "preview":soup.get_text(" ",strip=True)[:1700],
                "references":sorted(set(re.findall(r'[\w/.\-]+\.(?:pdb|csv|txt|js)(?:\?[\w=+-]+)?',html,re.I)))[:100],
                "links":[{"href":a.get("href"),"text":a.get_text(" ",strip=True)}
                         for a in soup.find_all("a",href=True)][:35],
                "scripts":[x.get("src") for x in soup.find_all("script",src=True)],
                "pdb_mentions":len(re.findall(r"\.pdb",html,re.I)),
                "potential_label":"Potential field PDB" in html,
                "npp_label":"NPP ratio PDB" in html})
        except Exception as exc:item["error"]=str(exc)
        if pid=="2OED":
            jsurl=BASE+"/js/patches/viewStructure.js"
            try:
                js=s.get(jsurl,timeout=(15,35))
                (OUT/"viewStructure.js").write_text(js.text,encoding="utf-8")
                item["viewer_js_status"]=js.status_code
                item["viewer_js_bytes"]=len(js.content)
                item["viewer_references"]=sorted(set(re.findall(
                    r'[^"\\s]{1,120}\\.(?:pdb|txt|csv|json|gz)',js.text,re.I)))[:50]
            except Exception as e:
                item["viewer_js_error"]=str(e)
        results.append(item)
        print("PROTEIN_SOL_RESULT_CHECK "+json.dumps(item,ensure_ascii=False),flush=True)
(OUT/"RESULTS.json").write_text(json.dumps(results,indent=2),encoding="utf-8")
