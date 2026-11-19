#!/usr/bin/env python3
"""Fetch public ACS Figshare SI by exact article ID, Wiley and RSC supplements.
Track every file URL, expected bytes, HTTP, SHA, and extracted text.
"""
import requests, re, json, hashlib, zipfile, io, os, subprocess
from pathlib import Path
root=Path("is44_SI_direct");root.mkdir(exist_ok=True)
s=requests.Session();s.headers.update({"User-Agent":"IS44 evidence audit noncommercial science (mailto:research@example.org)"})
LINKS=[
("10.1021/ac201889j",2613988),
("10.1021/acs.analchem.9b00482",7813913),
("10.1021/acs.biomac.1c00495",15087873),
("10.1021/acs.jpcb.5c07826",31119064),
("10.1021/acs.jpcb.7b01637",4829560),
("10.1021/acs.jpcb.7b01637",4829557),
("10.1021/acs.jpcc.6b08469",4033365),
("10.1021/acs.langmuir.3c02609",24942759),
("10.1021/acsnano.9b01622",8153465),
("10.1021/acsnano.9b01622",8153462),
("10.1021/nl101746v",2743963)]
MAX_BYTES=20_000_000
records=[]
def req(url,timeout=32):
 try:return s.get(url,timeout=timeout,allow_redirects=True)
 except Exception as e:print("ERROR",url,str(e)[:130],flush=True);return None
def save_file(path,body,sourceurl,doi,typ):
 if len(body)>MAX_BYTES or len(body)<200:return None
 typ=typ.lower()
 if body[:4]==b'%PDF':typ='pdf'
 elif body[:2]==b'PK':typ='zip'
 elif body[:4].lower() in [b'<htm',b'<!do']:return None
 path.parent.mkdir(parents=True,exist_ok=True)
 path.write_bytes(body)
 entry={"doi":doi,"file":str(path),"type":typ,"source_url":sourceurl,"bytes":len(body),"sha256":hashlib.sha256(body).hexdigest()}
 records.append(entry)
 if typ=="pdf":
  output=path.with_suffix(".txt")
  try:
   x=subprocess.run(["pdftotext","-layout",str(path),str(output)],capture_output=True,timeout=24)
   if x.returncode:return entry
   t=output.read_text(encoding='utf-8',errors='replace')
   lines=t.splitlines()
   relevant=[]
   patt=re.compile(r"\b([A-Z][a-z]{2}|[A-Z])\s?[-]?\d{1,4}\b|hdx|deuter|methylation|proteolysis|contact|residues|table\s*S?\d|figure\s*S?\d|peptide",re.I)
   for idx,line in enumerate(lines):
    if patt.search(line):
     relevant.append({"line":idx+1,"text":' '.join(lines[max(0,idx-1):min(len(lines),idx+2)]).strip()[:600]})
   path.with_name(path.stem+"_candidate_lines.json").write_text(json.dumps(relevant[:1500],ensure_ascii=False,indent=2),encoding="utf-8")
  except Exception as e:entry["extract_error"]=str(e)[:120]
 elif typ=="zip":
  try:
   z=zipfile.ZipFile(io.BytesIO(body))
   entry["zip_members"]=[{"name":x.filename,"bytes":x.file_size} for x in z.infolist()][:100]
   for k,m in enumerate(z.infolist()[:12]):
    if m.file_size<MAX_BYTES and m.filename.lower().endswith(".pdf"):
     b=z.read(m)
     save_file(path.with_name(path.stem+f"_attachment{k}.pdf"),b,sourceurl+"#"+m.filename,doi,"pdf")
  except Exception as e:entry["zip_error"]=str(e)[:120]
 return entry
for doi,aid in LINKS:
 print("ARTICLE",doi,aid,flush=True)
 url=f"https://api.figshare.com/v2/articles/{aid}"
 r=req(url)
 if r is None:continue
 if not r.ok:
  print('ARTICLE_HTTP',r.status_code,aid,flush=True);continue
 try:j=r.json()
 except:continue
 d=root/re.sub("[^a-z0-9]+","_",doi.lower());d.mkdir(exist_ok=True)
 (d/f"article_{aid}.json").write_text(json.dumps(j,indent=2),encoding='utf-8')
 for item in j.get('files',[])[:16]:
  id=item.get('id');name=item.get('name',f"{id}.bin")
  link=item.get('download_url') or f"https://ndownloader.figshare.com/files/{id}"
  print("FILE",aid,id,name,item.get('size'),link,flush=True)
  if item.get('size',0)>MAX_BYTES:continue
  r=req(link,60)
  if r is None:continue
  print("FILE_HTTP",r.status_code,len(r.content),flush=True)
  if r.ok:
   p=d/f"{aid}_{id}_{re.sub('[^a-zA-Z0-9._-]+','_',name)[:90]}"
   save_file(p,r.content,link,doi,p.suffix[1:])
extra=[
 ("10.1002/smll.202305684","https://onlinelibrary.wiley.com/doi/supinfo/10.1002/smll.202305684","smll202305684-SuppMat"),
 ("10.1002/smll.202305684","https://onlinelibrary.wiley.com/doi/suppl/10.1002/smll.202305684/supinfo/smll202305684-sup-0001-SuppMat.pdf","smll202305684-SuppMat"),
 ("10.1002/smll.202305684","https://advanced.onlinelibrary.wiley.com/doi/supinfo/10.1002/smll.202305684","smll202305684-SuppMat"),
 ("10.1039/D4SM00912F","https://www.rsc.org/suppdata/d4/sm/d4sm00912f/d4sm00912f1.pdf","d4sm00912f1"),
 ("10.1039/C8SM01909F","https://www.rsc.org/suppdata/c8/sm/c8sm01909f/c8sm01909f1.pdf","c8sm01909f1"),
 ("10.1039/C8TB01129J","https://www.rsc.org/suppdata/c8/tb/c8tb01129j/c8tb01129j1.pdf","c8tb01129j1"),
 ("10.1039/D2SC04434J","https://www.rsc.org/suppdata/d2/sc/d2sc04434j/d2sc04434j1.pdf","d2sc04434j1")]
for doi,link,name in extra:
 print("EXTRA",doi,link,flush=True)
 r=req(link,42)
 if r is None:continue
 print("EXTRA_HTTP",r.status_code,len(r.content),r.headers.get("content-type"),flush=True)
 if r.ok and r.content.startswith(b'%PDF'):
  d=root/re.sub("[^a-z0-9]+","_",doi.lower())
  save_file(d/(name+".pdf"),r.content,link,doi,"pdf")
(root/"retrieval_manifest.json").write_text(json.dumps(records,indent=2),encoding="utf-8")
print(json.dumps({"downloaded":len(records),"bytes":sum(z["bytes"] for z in records),"doi_with_files":sorted(set(z["doi"] for z in records))},indent=2),flush=True)
