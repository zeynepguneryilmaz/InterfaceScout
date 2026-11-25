"""Controlled robustness of descriptive region definitions; no adsorption labels."""
import csv,json,hashlib,sys,copy,itertools,traceback
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'backend'),str(ROOT/'validation')]
import core,coarse_patch as cp,interface_engine as ie,strict_apo_compat_48 as prep
from model_settings import MODEL_VERSION
shard,nshards,reference,out=sys.argv[1:];shard,nshards=int(shard),int(nshards)
ref=Path(reference);out=Path(out);out.mkdir(exist_ok=True)
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def unique(n):
 x=list(ref.rglob(n));assert len(x)==1,(n,x);return x[0]
panel=read(ROOT/'validation/strict_apo_48_selected_candidates_v2_final.csv')
old={(r['protein'],r['channel']):r for r in read(unique('all_528_maps.csv'))}
provenance={r['protein']:r for r in json.loads(unique('provenance.json').read_text())}
def sha(x):return hashlib.sha256(x).hexdigest()
def union(p,k,field):return set(q for a in p[:k] for q in (a[field] if field!='center_key' else [a[field]]))
def jac(a,b):return len(a&b)/len(a|b) if a|b else None
def metric(pid,ch,variant,p0,p1,group):
 r=dict(pdb=pid,channel=ch,group=group,variant=variant,n_reference=len(p0),n_variant=len(p1))
 for k in (1,3,5):
  for f in ('members','center_key','seed_members'):
   a,b=union(p0,k,f),union(p1,k,f)
   r[f'jaccard_{f}_top{k}']=jac(a,b)
   r[f'n_reference_{f}_top{k}']=len(a);r[f'n_variant_{f}_top{k}']=len(b)
 a={p['center_key'] for p in p0 if p['pareto_front']==1};b={p['center_key'] for p in p1 if p['pareto_front']==1}
 r['front1_jaccard']=jac(a,b)
 return r
def run(row):
 pid=row['pdb_id']
 prep.prepare(pid);cif=prep.OUT/f'{pid}-assembly1.cif';pdb=prep.OUT/f'{pid}_prepared.pdb'
 assert sha(cif.read_bytes())==row['assembly_sha256']==provenance[pid]['source_sha256']
 assert sha(pdb.read_text().encode())==provenance[pid]['prepared_sha256']
 ctx=ie._prepare_shared_context(pH=7.4,ionic_mM=150.,temp_K=298.,pdb_id=None,pdb_text=pdb.read_text(),chain=None)
 data=ctx['core_result'];surface=data['surface_residues'];dist=core.build_distances(surface);ca=ctx['ca_nodes']
 assert MODEL_VERSION=='2.5.0-dual-accessibility'
 rows=[]
 for public,key,*_ in ie.PUBLIC_MAP_SPECS:
  original=data['chemistries'][key];fav=core.CHEMISTRIES[key]['favorable']
  def build():return cp.build_coarse_patches(v1_result=data,chemistry=key,ca_nodes=ca)
  baseline=build();o=old[(pid,public)]
  assert len(baseline)==int(o['n_patches']) and union(baseline,5,'members')==set(filter(None,o['top5_members'].split(';'))) and union(baseline,5,'center_key')==set(filter(None,o['top5_centers'].split(';'))),(pid,public,'baseline mismatch')
  # No labels, no fitting; exact production pipeline with one rule omitted.
  for aa in sorted(fav):
   core.CHEMISTRIES[key]['favorable']={a:v for a,v in fav.items() if a!=aa}
   try:
    data['chemistries'][key]=core.chemistry_map(surface,dist,key,7.4)
    rows.append(metric(pid,public,'omit_'+aa,baseline,build(),'class_omission'))
   finally:
    core.CHEMISTRIES[key]['favorable']=fav;data['chemistries'][key]=original
  # Explicit positive-only hydropathy alternative, same exposed geometry.
  if key=='hydrophobic':
   kd={'ILE':4.5,'VAL':4.2,'LEU':3.8,'PHE':2.8,'CYS':2.5,'MET':1.9,'ALA':1.8}
   core.CHEMISTRIES[key]['favorable']={a:(0.,'positive KD alternative',None) for a in kd}
   alt=[dict(r,scrsa=r['scrsa']*kd.get(r['res_name'],0.)/4.5) for r in surface]
   try:
    data['chemistries'][key]=core.chemistry_map(alt,dist,key,7.4)
    rows.append(metric(pid,public,'positive_KD',baseline,build(),'hydropathy_alternative'))
   finally:
    core.CHEMISTRIES[key]['favorable']=fav;data['chemistries'][key]=original
  # One change only: remove hemisphere condition from center suppression and membership.
  sf=cp.same_face
  try:
   cp.same_face=lambda *a:True
   rows.append(metric(pid,public,'no_same_face',baseline,build(),'geometry'))
  finally:cp.same_face=sf
  # Preserve actual Pareto fronts; only change lexicographic priority within front.
  fields=('chemistry_support','mean_accessibility','patch_coherence','orientation_coherence')
  for first in fields:
   second=tuple(f for f in fields if f!=first)
   ranked=sorted(baseline,key=lambda p:(p['pareto_front'],-float(p[first]),*[ -float(p[f]) for f in second],p['center_key']))
   rows.append(metric(pid,public,'first_'+first,baseline,ranked,'within_front_order'))
 return rows
rows=[];fail=[]
for i,r in enumerate(panel):
 if i%nshards!=shard:continue
 try:rows.extend(run(r));print('PASS',r['pdb_id'],flush=True)
 except Exception as e:fail.append({'pdb':r['pdb_id'],'error':str(e),'trace':traceback.format_exc()});print('FAIL',r['pdb_id'],str(e),flush=True)
 if rows:
  with (out/f'robustness_{shard}.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 (out/f'failures_{shard}.json').write_text(json.dumps(fail,indent=2))
sys.exit(bool(fail))
