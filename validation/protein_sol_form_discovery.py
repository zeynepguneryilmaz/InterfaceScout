"""Inspect Protein-sol Patches public form, endpoints, input restrictions.

No structures submitted. Avoid submitting anything until a real valid HTML
form and response contract have been identified. This is a server-only service
whose CPU/RAM cannot be measured from an external client.
"""
from __future__ import annotations
import json
from pathlib import Path
from urllib.parse import urljoin, urlparse
import requests
from bs4 import BeautifulSoup

OUT=Path("external_protein_sol")
OUT.mkdir(exist_ok=True)
BASE="https://www.protein-sol.manchester.ac.uk"
paths=["/patches","/examples/patches/1hzh/results"]
result={"service":"Protein-sol Patches","access_method":"public_web",
        "server_resource_usage_unavailable":True,"observations":[]}
with requests.Session() as s:
    s.headers.update({"User-Agent":"InterfaceScoutScientificBenchmark/1.0 (academic, 5 structures maximum)"})
    for path in paths:
        item={"url":urljoin(BASE,path)}
        try:
            r=s.get(item["url"],timeout=35)
            item.update({"http_status":r.status_code,"final_url":r.url,
                         "elapsed_s":round(r.elapsed.total_seconds(),4),
                         "content_type":r.headers.get("Content-Type"),
                         "bytes":len(r.content)})
            soup=BeautifulSoup(r.content,"html.parser")
            forms=[]
            for form in soup.find_all("form"):
                forms.append({
                    "method":form.get("method","GET").upper(),
                    "action":urljoin(r.url,form.get("action") or ""),
                    "enctype":form.get("enctype","application/x-www-form-urlencoded"),
                    "inputs":[{"name":x.get("name"),"type":x.get("type"),
                               "required":bool(x.has_attr("required"))}
                              for x in form.find_all(["input","textarea","select"])],
                })
            item["forms"]=forms
            item["anchors"]=[
                {"text":x.get_text(" ",strip=True)[:70],
                 "url":urljoin(r.url,x.get("href",""))}
                for x in soup.find_all("a",href=True)
                if any(q in x.get_text(" ",strip=True).lower()
                       for q in ["download","pdb","example","patches"])
            ][:30]
            item["external_script_urls"]=[
                urljoin(r.url,x.get("src"))
                for x in soup.find_all("script",src=True)
            ][:25]
            (OUT/(path.strip("/").replace("/","_")+".html")).write_text(
                str(soup),encoding="utf-8")
        except Exception as exc:
            item["error"]=type(exc).__name__+": "+str(exc)[:500]
        result["observations"].append(item)
        print("PROTEIN_SOL_INSPECT "+json.dumps(item),flush=True)
(OUT/"FORM_DISCOVERY.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
