"""Protein-sol Patches: real public-structure submission from GitHub Actions.

Use only the five preselected public ISV76 PDBs; record all HTTP evidence.
"HTTP 200" is not proof a scientific result was calculated. Provider output
is not directly interchangeable with InterfaceScout APBS/SASA quantities.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re,sys,time,traceback
from pathlib import Path
from urllib.parse import urljoin,urlparse
import requests
from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
from isv_full_80_qa import prepare
from prepare import prepare_pdb_text

ACTIVE={q["proposed_pdb_id"]:q for q in csv.DictReader(
    (ROOT/"validation/ISV_active_candidates_76_v3.csv").open(
        encoding="utf-8",newline=""))}
PIDS={"2OED","1HEL","7KWI","1GCT","1HC1"}
assert PIDS.issubset(ACTIVE)
OUT=Path("external_protein_sol_real").resolve()
BASE="https://protein-sol.manchester.ac.uk"
POST=BASE+"/cgi-bin/patches/patches.php"

def allowed(url):
    u=urlparse(url)
    return u.scheme=="https" and u.hostname=="protein-sol.manchester.ac.uk"

def extract_links(base,soup):
    links=[]
    for tag in soup.find_all(["a","iframe","script"],src=True)+soup.find_all("a",href=True):
        raw=tag.get("href") or tag.get("src")
        if not raw: continue
        u=urljoin(base,raw)
        if allowed(u):
            links.append({"label":tag.get_text(" ",strip=True)[:100],"url":u})
    return list({x["url"]:x for x in links}.values())

def capture_page(folder,name,url,session):
    if not allowed(url):
        return {"status":"REFUSED_NON_PROVIDER_HOST","url":url}
    start=time.monotonic()
    try:
        rsp=session.get(url,timeout=(15,45),allow_redirects=True)
        meta={"url":url,"final_url":rsp.url,"http_status":rsp.status_code,
              "duration_s":round(time.monotonic()-start,3),
              "content_type":rsp.headers.get("Content-Type",""),
              "bytes":len(rsp.content)}
        if not allowed(rsp.url):
            meta["status"]="UNEXPECTED_EXTERNAL_REDIRECT";return meta
        if len(rsp.content)>25_000_000:
            meta["status"]="TOO_LARGE";return meta
        suffix=".html" if "html" in meta["content_type"].lower() else ".bin"
        (folder/(name+suffix)).write_bytes(rsp.content)
        meta["status"]="SAVED" if rsp.ok else "HTTP_ERROR"
        if suffix==".html":
            soup=BeautifulSoup(rsp.content,"html.parser")
            meta["title"]=soup.title.get_text(" ",strip=True) if soup.title else ""
            meta["text_excerpt"]=soup.get_text(" ",strip=True)[:1000]
            meta["links"]=extract_links(rsp.url,soup)[:80]
        return meta
    except Exception as exc:
        return {"status":"FETCH_FAILED","url":url,"exception":str(exc)[:600]}

def run(pid):
    folder=OUT/pid;folder.mkdir(parents=True,exist_ok=True)
    source,original=prepare(ACTIVE[pid])
    pdb,_=prepare_pdb_text(original.read_text(encoding="utf-8"))
    payload=pdb.encode()
    (folder/(pid+"_submitted.pdb")).write_bytes(payload)
    report={"pdb_id":pid,"service":"Protein-sol Patches",
            "endpoint":POST,"input_sha256":hashlib.sha256(payload).hexdigest(),
            "input_bytes":len(payload),
            "original_assembly_sha256":source["assembly_cif_sha256"],
            "submitted_public_RCSB_PDB":True,
            "server_cpu_ram":"UNAVAILABLE_TO_EXTERNAL_CLIENT",
            "scientific_job_status":"NOT_VERIFIED",
            "pages":[]}
    t0=time.monotonic()
    try:
        with requests.Session() as sess:
            sess.headers.update({"User-Agent":"InterfaceScout-method-comparison/1.0 academic evaluation"})
            rsp=sess.post(POST,files={"uploadedFile":(pid+".pdb",payload,"chemical/x-pdb")},
                          data={"multiprediction":"Upload"},
                          timeout=(20,160),allow_redirects=True)
            report["upload"]={
                "status_code":rsp.status_code,
                "final_url":rsp.url,
                "elapsed_s":round(time.monotonic()-t0,3),
                "content_type":rsp.headers.get("Content-Type",""),
                "bytes":len(rsp.content),
                "redirects":[{"status":z.status_code,"url":z.url} for z in rsp.history]}
            if not allowed(rsp.url):
                report["scientific_job_status"]="UNEXPECTED_REDIRECT"
            elif len(rsp.content)<25_000_000:
                ext=".html" if "html" in rsp.headers.get("Content-Type","").lower() else ".bin"
                (folder/("upload_response"+ext)).write_bytes(rsp.content)
                if ext==".html":
                    soup=BeautifulSoup(rsp.content,"html.parser")
                    links=extract_links(rsp.url,soup)
                    report["upload"]["title"]=soup.title.get_text(" ",strip=True) if soup.title else ""
                    report["upload"]["text_excerpt"]=soup.get_text(" ",strip=True)[:1700]
                    report["upload"]["links"]=links[:100]
                    is_result=(re.search(r"/results?(/|$)",rsp.url,re.I) is not None or
                               any(re.search(r"/results?(/|$)",x["url"],re.I) for x in links))
                    report["scientific_job_status"]=("RESULT_LINK_DETECTED_NOT_VALIDATED"
                        if is_result else "UPLOAD_RETURNED_NO_RESULT")
                    potential=[x["url"] for x in links
                               if re.search(r"/results?(/|$)|\.pdb($|\?)|\.csv($|\?)|\.txt($|\?)",
                                            x["url"],re.I)]
                    if re.search(r"/results?(/|$)",rsp.url,re.I):
                        potential=[rsp.url]+potential
                    # Follow at most eight provider-host links; no job polling.
                    for idx,u in enumerate(dict.fromkeys(potential)):
                        if idx>=8: break
                        report["pages"].append(capture_page(
                            folder,f"provider_result_{idx:02d}",u,sess))
                else:
                    report["scientific_job_status"]="NON_HTML_RESPONSE_NOT_VALIDATED"
            if rsp.status_code>=400:
                report["scientific_job_status"]="UPLOAD_HTTP_ERROR"
    except Exception as exc:
        report["scientific_job_status"]="NETWORK_OR_SERVER_ERROR"
        report["error"]=str(exc)[:1000]
        report["traceback_tail"]=traceback.format_exc()[-2400:]
    report["local_elapsed_s"]=round(time.monotonic()-t0,3)
    report["downloaded_provider_files"]=[
        {"name":z.name,"bytes":z.stat().st_size}
        for z in sorted(folder.iterdir()) if z.is_file() and z.name.startswith("provider_result_")]
    (folder/"SUMMARY.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("PROTEIN_SOL_REAL "+json.dumps({
        "pdb_id":pid,"status":report["scientific_job_status"],
        "http":report.get("upload",{}).get("status_code"),
        "url":report.get("upload",{}).get("final_url"),
        "elapsed_s":report["local_elapsed_s"],
        "result_files":len(report["downloaded_provider_files"]),
        "excerpt":report.get("upload",{}).get("text_excerpt","")[:550],
        "error":report.get("error")},ensure_ascii=False),flush=True)
    # A submission with unverified results cannot count as a PASS.
    if report["scientific_job_status"] not in (
        "RESULT_LINK_DETECTED_NOT_VALIDATED",):
        raise SystemExit("Scientific result not verified; see SUMMARY.json")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--pdb",choices=sorted(PIDS),required=True)
    args=parser.parse_args()
    run(args.pdb)

if __name__=="__main__":
    main()
