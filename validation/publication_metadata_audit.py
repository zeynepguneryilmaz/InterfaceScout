"""Fetch primary deposited citations and verify exact assembly context."""
import csv,json,hashlib,concurrent.futures,time
from pathlib import Path
import requests,gemmi
ROOT=Path(__file__).resolve().parents[1];OUT=Path('publication_metadata');OUT.mkdir(exist_ok=True)
def get(url):
 for i in range(4):
  try:
   r=requests.get(url,timeout=30);r.raise_for_status();return r
  except Exception:
   if i==3:raise
   time.sleep(i+1)
def rows(p):return list(csv.DictReader((ROOT/'validation'/p).open()))
a=rows('strict_apo_48_selected_candidates_v2_final.csv');b=rows('ISV_active_candidates_76_v3.csv')
ids=[(x['pdb_id'],'TV48') for x in a]+[(x['proposed_pdb_id'],'ISV76') for x in b];assert len(ids)==124
def citation(item):
 pid,panel=item
 try:
  j=get('https://data.rcsb.org/rest/v1/core/entry/'+pid).json()
  (OUT/f'{pid}_entry.json').write_text(json.dumps(j,indent=2))
  c=j.get('rcsb_primary_citation') or next((c for c in j.get('citation',[]) if c.get('id')=='primary'),{})
  doi=c.get('pdbx_database_id_doi') or ''
  authors=c.get('rcsb_authors') or c.get('authors') or []
  title=c.get('title') or ''
  journal=c.get('journal_abbrev') or ''
  year=c.get('year') or ''
  pages='-'.join(str(v) for v in [c.get('page_first'),c.get('page_last')] if v)
  vol=c.get('journal_volume') or ''
  pmid=c.get('pdbx_database_id_pub_med') or ''
  status='deposited_article_metadata' if journal and journal.lower() not in ('to be published','unpublished') and title and authors and year else 'deposition_only_or_unresolved'
  check=''
  if doi:
   try:
    x=get('https://api.crossref.org/works/'+doi).json()['message'];(OUT/f'{pid}_crossref.json').write_text(json.dumps(x,indent=2))
    check='DOI_resolves_crossref'
   except Exception as e:check='DOI_lookup_failed:'+str(e)
  if not doi and pmid:
   try:
    x=get('https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:'+str(pmid)+'%20AND%20SRC:MED&format=json').json()
    (OUT/f'{pid}_europepmc.json').write_text(json.dumps(x,indent=2))
    hits=x['resultList']['result']
    if len(hits)==1:doi=hits[0].get('doi','');check='PMID_resolved'
   except Exception as e:check='PMID_lookup_failed:'+str(e)
  return dict(panel=panel,pdb_id=pid,data_doi=j.get('rcsb_accession_info',{}).get('deposit_date','') and '10.2210/pdb'+pid.lower()+'/pdb',title=title,authors='; '.join(authors),journal=journal,year=year,volume=vol,pages=pages,article_doi=doi,pmid=pmid,status=status,verification=check,entry_title=j.get('struct',{}).get('title',''),deposition_authors='; '.join(x.get('name','') for x in j.get('audit_author',[])))
 except Exception as e:return dict(panel=panel,pdb_id=pid,status='fetch_failed',verification=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:result=list(ex.map(citation,ids))
fields=list(dict.fromkeys(k for r in result for k in r))
with (OUT/'124_citations.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(result)
STD=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
for pid in ['3IPJ','4LE4','4LE3']:
 for scope,suffix in [('assembly1','-assembly1.cif'),('asu','.cif')]:
  try:
   p=OUT/(pid+suffix);p.write_bytes(get('https://files.rcsb.org/download/'+pid+suffix).content);st=gemmi.read_structure(str(p));m=st[0];protein=[(c.name,str(r.seqid),r.name,x) for c in m for r in c if r.name in STD for x in r if x.element.name not in ['H','D']]
   lig=[]
   for c in m:
    for r in c:
     if r.name in STD or r.name in ('HOH','DOD','WAT'):continue
     atoms=[x for x in r if x.element.name not in ['H','D']]
     if not atoms:continue
     hits={}
     for x in atoms:
      for ch,num,name,y in protein:
       d=x.pos.dist(y.pos)
       if d<=5:
        key=f'{ch}:{num}:{name}';hits[key]=min(hits.get(key,999),d)
     lig.append(dict(chain=c.name,residue=str(r.seqid),name=r.name,n_heavy=len(atoms),elements=sorted(set(x.element.name for x in atoms)),contacts_A={k:round(v,3) for k,v in sorted(hits.items(),key=lambda kv:kv[1])}))
   sha=hashlib.sha256(p.read_bytes()).hexdigest();frozen=next((x['assembly_sha256'] for x in a if x['pdb_id']==pid),None)
   q=dict(pdb=pid,scope=scope,sha256=sha,frozen_sha256=frozen,hash_matches=(sha==frozen if frozen and scope=='assembly1' else None),chains=[dict(chain=c.name,standard_residues=sum(r.name in STD for r in c)) for c in m],ligands=lig)
   (OUT/f'{pid}_{scope}_audit.json').write_text(json.dumps(q,indent=2))
   print('STRUCTURE',json.dumps(q),flush=True)
  except Exception as e:print('STRUCTURE_ERROR',pid,scope,str(e),flush=True)
from collections import Counter
print('CITATIONS',Counter(x['status'] for x in result),flush=True)
assert len(result)==124
