import requests,json,time,csv,html
from pathlib import Path
p=Path("is44_publication_metadata");p.mkdir(exist_ok=True)
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
s=requests.Session();s.headers.update({"User-Agent":"IS44-manuscript-references/1.0 (mailto:research@openai.com)"})
out=[]
for doi in DOIS:
 d={"DOI":doi,"URL":"https://doi.org/"+doi}
 try:
  r=s.get("https://api.crossref.org/works/"+doi,timeout=13);r.raise_for_status();m=r.json()["message"]
  d.update({"title":m.get("title",[""])[0],"journal":m.get("container-title",[""])[0],"volume":m.get("volume",""),"issue":m.get("issue",""),"page":m.get("page",m.get("article-number","")),"year":(m.get("published",m.get("issued",{})).get("date-parts",[[None]])[0][0]),"authors":", ".join([((a.get("family","")+", "+a.get("given","")).strip(", ")) for a in m.get("author",[])[:12]]),"publisher":m.get("publisher",""),"status":"crossref-confirmed"})
 except Exception as e:d["status"]="not-retrieved";d["error"]=str(e)[:160]
 out.append(d);print(doi,d.get("status"),d.get("title","")[:95],flush=True)
 time.sleep(.2)
(p/"metadata.json").write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding="utf-8")
with (p/"metadata.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=["DOI","URL","authors","title","journal","year","volume","issue","page","publisher","status","error"],extrasaction="ignore");w.writeheader();w.writerows(out)
print("retrieved",sum(x["status"]=="crossref-confirmed" for x in out),"of",len(out))
