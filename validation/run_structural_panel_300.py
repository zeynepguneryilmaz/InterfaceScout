import csv,io,random,json,re,sys,urllib.request,hashlib,itertools,statistics
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict,Counter
from Bio.PDB import PDBParser
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'backend'))
from interface_engine import analyze_all_maps
root=Path('validation/structural_panel_output');root.mkdir(parents=True,exist_ok=True)
url='https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt'
seed=20261116;rng=random.Random(seed)
def get(u,timeout=30):
 req=urllib.request.Request(u,headers={'User-Agent':'InterfaceScout academic validation'})
 with urllib.request.urlopen(req,timeout=timeout) as r:return r.read()
cluster_bytes=get(url,90);clusters=cluster_bytes.decode().splitlines()
ids=list(range(len(clusters)));rng.shuffle(ids)
bins=[(60,119,60),(120,199,75),(200,349,90),(350,499,45),(500,800,30)]
selected=[[] for _ in bins];seen=set();fails=Counter();attempts=0
def inspect(item):
 j,pdb,entity=item
 try:
  text=get('https://files.rcsb.org/download/'+pdb+'.pdb').decode('utf8','replace')
  method=next((x[10:].strip() for x in text.splitlines() if x.startswith('EXPDTA')),'')
  if not any(x in method for x in ('X-RAY DIFFRACTION','ELECTRON MICROSCOPY','SOLUTION NMR')):return None,'method'
  compound=' '.join(x[10:].strip() for x in text.splitlines() if x.startswith('COMPND'))
  matches=re.findall(r'MOL_ID:\\s*(\\d+)\\s*;(?:(?!MOL_ID:).)*?CHAIN:\\s*([^;]+);',compound)
  chains=[v for k,v in matches if k==entity]
  if not chains:return None,'entity_unresolved'
  chain=chains[0].split(',')[0].strip()
  structure=PDBParser(QUIET=True).get_structure(pdb,io.StringIO(text));model=next(iter(structure))
  if chain not in model:return None,'chain_missing'
  allowed=set('ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
  length=sum(1 for r in model[chain] if r.id[0]==' ' and r.resname in allowed and 'CA' in r)
  if not 60<=length<=800:return None,'length'
  return dict(pdb_id=pdb,entity_id=entity,cluster_index=j,chain=chain,n_residues=length,method=method,sha256=hashlib.sha256(text.encode()).hexdigest(),pdb_text=text),'ok'
 except Exception as e:return None,type(e).__name__
pool=ThreadPoolExecutor(max_workers=8)
for batchstart in range(0,min(18000,len(ids)),64):
 if all(len(group)>=quota for group,(_,_,quota) in zip(selected,bins)):break
 items=[]
 for j in ids[batchstart:batchstart+64]:
  members=clusters[j].split()
  if not members:continue
  token=rng.choice(members)
  if re.fullmatch(r'[A-Za-z0-9]{4}_[0-9]+',token):
   pdb,entity=token.split('_');items.append((j,pdb.upper(),entity))
 for r,reason in pool.map(inspect,items):
  attempts+=1
  if r is None:fails[reason]+=1;continue
  if r['pdb_id'] in seen:continue
  index=next((i for i,(lo,hi,quota) in enumerate(bins) if lo<=r['n_residues']<=hi and len(selected[i])<quota),None)
  if index is None:continue
  seen.add(r['pdb_id']);selected[index].append(r)
 if batchstart%512==0:print('SELECTION_PROGRESS',attempts,sum(map(len,selected)),list(map(len,selected)),flush=True)
pool.shutdown(wait=True)
panel=[r for group in selected for r in group]
(root/'pdb').mkdir(exist_ok=True)
for r in panel:(root/'pdb'/(r['pdb_id']+'_'+r['chain']+'.pdb')).write_text(r.pop('pdb_text'))
with (root/'manifest.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(panel[0]));w.writeheader();w.writerows(panel)
protocol={'seed':seed,'rcsb_url':url,'rcsb_cluster_file_sha256':hashlib.sha256(cluster_bytes).hexdigest(),'bins':bins,'selected':len(panel),'attempts':attempts,'failures':dict(fails),'prediction_condition':{'pH':7.4,'ionic_mM':150,'temp_K':298},'limitations':['Cluster nonredundancy is RCSB-defined with coverage requirements','No adsorption ground truth','PDB entity to PDB COMPND MOL_ID mapping verified by parsed chain']}
(root/'protocol.json').write_text(json.dumps(protocol,indent=2))
print('SELECTION_COMPLETE',len(panel),list(map(len,selected)),flush=True)
if len(panel)!=300:raise RuntimeError('Unable to curate 300 eligible structures')
results=[];pair_rows=[];fail_rows=[]
for i,r in enumerate(panel,1):
 try:
  a=analyze_all_maps(pdb_text=(root/'pdb'/(r['pdb_id']+'_'+r['chain']+'.pdb')).read_text(),chain=r['chain'])
  tops={}
  for k,v in a['maps'].items():
   u=set(x for patch in v['patches'][:5] for x in patch['members']);tops[k]=u
   results.append(dict(pdb_id=r['pdb_id'],chain=r['chain'],n_residues=r['n_residues'],map=k,n_patches=len(v['patches']),top5_n_residues=len(u),n_surface=a['diagnostics']['n_surface_residues']))
  for k,l in itertools.combinations(tops,2):
   uni=tops[k]|tops[l]
   pair_rows.append(dict(pdb_id=r['pdb_id'],map_1=k,map_2=l,jaccard=(len(tops[k]&tops[l])/len(uni) if uni else ''),both_nonempty=int(bool(tops[k]) and bool(tops[l]))))
  print('ANALYZED',i,'/',len(panel),r['pdb_id'],flush=True)
 except Exception as e:
  fail_rows.append({'pdb_id':r['pdb_id'],'error':repr(e)})
  print('ANALYSIS_ERROR',r['pdb_id'],repr(e),flush=True)
for filename,rows in [('map_metrics.csv',results),('pairwise_map_jaccard.csv',pair_rows),('analysis_errors.csv',fail_rows)]:
 if rows:
  with (root/filename).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
pairs=defaultdict(list)
for r in pair_rows:
 if r['both_nonempty'] and r['jaccard']!='':pairs[(r['map_1'],r['map_2'])].append(r['jaccard'])
summary={'selected':len(panel),'analyzed':len(set(x['pdb_id'] for x in results)),'failures':len(fail_rows),'chemistry_map_count':11,'total_maps':len(results),'median_jaccard':{'|'.join(p):round(statistics.median(v),5) for p,v in pairs.items()},'scope':'structural behavior only; no adsorption-site predictive validation'}
(root/'summary.json').write_text(json.dumps(summary,indent=2))
print('ANALYSIS_COMPLETE',json.dumps({'selected':summary['selected'],'analyzed':summary['analyzed'],'failures':summary['failures'],'maps':summary['total_maps']}),flush=True)
