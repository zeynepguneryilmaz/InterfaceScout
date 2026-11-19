#!/usr/bin/env python3
"""Public-source fulltext/SI retrieval, with per-file provenance and strict evidence labeling.
No residue is declared experimentally adsorbed by regex alone.
"""
import requests,re,json,os,time,subprocess,shutil,hashlib,html,urllib.parse
from pathlib import Path
from bs4 import BeautifulSoup
ROOT=Path("is44_literature_supplement_audit");ROOT.mkdir(exist_ok=True)
DOIS="""10.1002/cphc.202300505
10.1002/smll.202305684
10.1016/j.csbj.2025.06.004
10.1016/j.jbc.2026.113533
10.1016/j.jcis.2005.03.064
10.1021/ac201889j
10.1021/acs.analchem.9b00482
10.1021/acs.biomac.1c00495
10.1021/acs.jpcb.5c07826
10.1021/acs.jpcb.6b07646
10.1021/acs.jpcb.7b01637
10.1021/acs.jpcc.6b08469
10.1021/acs.langmuir.3c02609
10.1021/acsnano.9b01622
10.1021/nl101746v
10.1038/s41467-022-34749-w
10.1038/s41541-019-0115-7
10.1038/s41598-025-96435-3
10.1039/C8SM01909F
10.1039/C8TB01129J
10.1039/D2SC04434J
10.1039/D4SM00912F
10.1371/journal.pone.0107696
10.1529/biophysj.104.054809
10.3390/biom16010114
10.3390/ijms20143539""".splitlines()
s=requests.Session();s.headers.update({"User-Agent":"InterfaceScoutEvidenceAudit/1.1 (noncommercial, academic research)","Accept":"*/*"})
def get(url,timeout=24):
    try:
        r=s.get(url,timeout=timeout)
        if r.ok and len(r.content)<18000000:return r
    except Exception as e:return None
    return None
def save_text(p,t):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(str(t),encoding="utf-8")
def extract_main(xml):
    soup=BeautifulSoup(xml,"xml")
    for x in soup(["ref-list","permissions","table-wrap-foot","fn-group"]):x.decompose()
    texts=[]
    for el in soup.find_all(["p","title","table","fig"]):
        if el.find_parent(["p","table","fig"]):continue
        t=" ".join(el.get_text(" ",strip=True).split())
        if len(t)>7:texts.append(t)
    return "\n".join(texts)
def extract_hits(txt):
    pat=re.compile(r"(?:(?:[A-Z][a-z]{2})|(?:[A-Z]))\s?[-]?\d{1,4}\b|(?:residue|peptide|contact|lysine|surface|HDX|adsorption).{0,95}",re.I)
    lines=txt.splitlines()
    return [x[:650] for x in lines if pat.search(x)][:500]
def record_binary(doc_dir,meta,url,data,typ):
    n=len(meta["assets"])+1
    suff=".pdf" if data[:4]==b"%PDF" else (".xml" if data[:5].lower().startswith(b"<?xml") else ".bin")
    if suff==".bin":
        if data[:2]==b"PK":suff=".zip"
        elif data[:5].lower()==b"<html":return
    path=doc_dir/f"asset_{n:02d}{suff}"
    path.write_bytes(data)
    entry={"url":url,"file":str(path),"sha256":hashlib.sha256(data).hexdigest(),"bytes":len(data),"type":typ}
    if suff==".pdf":
        try:
            subprocess.run(["pdftotext","-layout",str(path),str(path)+".txt"],timeout=18,capture_output=True)
            t=Path(str(path)+".txt").read_text(encoding="utf-8",errors="replace")
            entry["text_length"]=len(t)
            save_text(doc_dir/f"asset_{n:02d}_residue_candidate_lines.txt","\n".join(extract_hits(t)))
        except Exception as e:entry["extract_error"]=str(e)[:120]
    meta["assets"].append(entry)
    return entry
metadata=[]
for i,doi in enumerate(DOIS,1):
    safe=re.sub("[^a-zA-Z0-9]+","_",doi)
    folder=ROOT/f"{i:02d}_{safe}";folder.mkdir(exist_ok=True)
    meta={"doi":doi,"doi_link":"https://doi.org/"+doi,"pmcid":"","europepmc_status":"","assets":[],"source_urls":[],"errors":[]}
    # EuropePMC: identify fulltext and PMCID
    u="https://www.ebi.ac.uk/europepmc/webservices/rest/search?query="+urllib.parse.quote('DOI:"'+doi+'"')+"&format=json&pageSize=3"
    r=get(u)
    if r:
        try:
            articles=r.json().get("resultList",{}).get("result",[])
            if articles:
                a=articles[0];meta["pmcid"]=a.get("pmcid","");meta["title"]=a.get("title","")
        except Exception as e:meta["errors"].append("EuropePMC parse: "+str(e)[:150])
    if meta["pmcid"]:
        pmc=meta["pmcid"]
        u=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML"
        r=get(u)
        if r and "<article" in r.text[:400]:
            meta["europepmc_status"]="fulltext";meta["source_urls"].append(u)
            save_text(folder/"fulltext.xml",r.text)
            plain=extract_main(r.text);save_text(folder/"fulltext.txt",plain)
            save_text(folder/"fulltext_residue_candidate_lines.txt","\n".join(extract_hits(plain)))
            soup=BeautifulSoup(r.text,"xml")
            si=[]
            for el in soup.find_all(["supplementary-material","media","ext-link"]):
                href=el.get("href") or el.get("xlink:href") or el.get("rid")
                if href and any(z in (str(href)+el.get_text(" ",strip=True)).lower() for z in (".pdf",".xlsx",".csv",".docx","suppl","support")):
                    si.append(str(href))
            meta["supplementary_listed"]=list(dict.fromkeys(si))[:40]
            # direct supplementary bin endpoint
            for href in meta["supplementary_listed"][:9]:
                if not href.lower().endswith((".pdf",".xlsx",".csv",".zip",".docx")):continue
                urls=[]
                if href.startswith("http"):urls.append(href)
                else:
                    urls+=[f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/bin/{href}",
                           f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles"]
                for cand in urls[:2]:
                    rr=get(cand,35)
                    if rr and (rr.content[:4]==b"%PDF" or rr.content[:2]==b"PK"):
                        record_binary(folder,meta,cand,rr.content,"pmc-supplementary");break
        else: meta["europepmc_status"]="no-fulltextXML"
    # open-access RSC supplement exact DOI suffix, where available
    m=re.match(r"10\.1039/([CD]\d)([A-Z]{2})(\w+)",doi,re.I)
    if m:
        code=doi.split("/")[-1].lower()
        url=f"https://www.rsc.org/suppdata/{m[1].lower()}/{m[2].lower()}/{code}/{code}1.pdf"
        rr=get(url,35)
        if rr and rr.content[:4]==b"%PDF":
            record_binary(folder,meta,url,rr.content,"RSC-ESI")
        else:meta["source_urls"].append(url+" [not fetched]")
    # ACS Figshare metadata: 1) public supplementary page 2) api articles filtered resource DOI
    if "10.1021/" in doi:
        page="https://acs.figshare.com/doi/suppl/"+doi
        rr=get(page)
        if rr:
            text=rr.text
            meta["source_urls"].append(page)
            arlinks=sorted(set(re.findall(r'https?://[^\\" ]+(?:_si_0\d+|download|supplementary)[^\\" ]*',text,re.I)))
            arlinks=[html.unescape(z).replace("&amp;","&") for z in arlinks][:8]
            for url in arlinks:
                if "acs.figshare.com" not in url and "ndownloader.figshare.com" not in url:continue
                result=get(url,35)
                if result and (result.content[:4]==b"%PDF" or result.content[:2]==b"PK"):
                    record_binary(folder,meta,url,result.content,"ACS-figshare")
        api="https://api.figshare.com/v2/articles?resource_doi="+urllib.parse.quote(doi)
        rr=get(api)
        if rr:
            try:
                data=rr.json()
                meta["figshare_hits"]=[{"id":x.get("id"),"title":x.get("title")} for x in data[:8]] if isinstance(data,list) else []
            except:pass
    metadata.append(meta)
    save_text(folder/"metadata.json",json.dumps(meta,indent=2))
    print(i,doi,"pmc",meta["pmcid"],"OA",meta["europepmc_status"],"assets",len(meta["assets"]),flush=True)
save_text(ROOT/"audit_manifest.json",json.dumps(metadata,indent=2))
summary={"total_doi":len(metadata),"pmc_fulltext":sum(x["europepmc_status"]=="fulltext" for x in metadata),"asset_count":sum(len(x["assets"]) for x in metadata),"unretrieved_fulltext":[x["doi"] for x in metadata if x["europepmc_status"]!="fulltext"],"assets":[{"doi":x["doi"],"urls":[z["url"] for z in x["assets"]]} for x in metadata]}
save_text(ROOT/"summary.json",json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2),flush=True)
