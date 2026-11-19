import requests,hashlib,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
root=Path("is44_other_supplements");root.mkdir(exist_ok=True)
targets=[
("10.1038/s41541-019-0115-7","asparaginase","https://pmc.ncbi.nlm.nih.gov/articles/PMC6538755/bin/41541_2019_115_MOESM1_ESM.pdf"),
("10.1038/s41541-019-0115-7","asparaginase","https://static-content.springer.com/esm/art%3A10.1038%2Fs41541-019-0115-7/MediaObjects/41541_2019_115_MOESM1_ESM.pdf"),
("10.1002/smll.202305684","small2024","https://pmc.ncbi.nlm.nih.gov/articles/PMC11209821/bin/smll202305684-sup-0001-SuppMat.pdf"),
("10.1002/smll.202305684","small2024","https://pmc.ncbi.nlm.nih.gov/articles/PMC10350082/bin/smll202305684-sup-0001-SuppMat.pdf"),
("10.1002/smll.202305684","small2024","https://pmc.ncbi.nlm.nih.gov/articles/PMC11209821/bin/nihms1958479-supplement-sup-0001.pdf"),
("10.1002/smll.202305684","small2024","https://onlinelibrary.wiley.com/doi/suppl/10.1002/smll.202305684/supinfo/smll202305684-sup-0001-SuppMat.pdf"),
("10.1039/D4SM00912F","rsc2025","https://www.rsc.org/suppdata/d4/sm/d4sm00912f/d4sm00912f1.pdf"),
("10.1039/D4SM00912F","rsc2025","https://pubs.rsc.org/en/content/articlepdf/2025/sm/d4sm00912f"),
("10.1039/D4SM00912F","rsc2025","https://www.rsc.org/suppdata/d4/sm/d4sm00912f/d4sm00912f1.pdf"),
("10.1038/s41598-025-96435-3","swcnt2025","https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-025-96435-3/MediaObjects/41598_2025_96435_MOESM1_ESM.pdf")
]
s=requests.Session();s.headers['User-Agent']='Mozilla/5.0'
def get(arg):
 doi,tag,url=arg
 try:
  r=s.get(url,timeout=18)
  valid=r.status_code==200 and r.content.startswith(b'%PDF')
  res={"doi":doi,"tag":tag,"url":url,"status":r.status_code,"size":len(r.content),"valid_pdf":valid}
  if valid:
   h=hashlib.sha256(r.content).hexdigest();p=root/f"{tag}_{h[:12]}.pdf";p.write_bytes(r.content);res["path"]=str(p);res["sha256"]=h
  return res
 except Exception as e: return {"doi":doi,"tag":tag,"url":url,"error":str(e)[:100]}
with ThreadPoolExecutor(max_workers=6) as pool:out=list(pool.map(get,targets))
(root/"manifest.json").write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
