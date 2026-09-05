from __future__ import annotations
import json, sys, tempfile, urllib.request
from pathlib import Path
from collections import deque
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from Bio.PDB import PDBParser

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import backend.main as ism

ROOT=Path('validation/packages/all15')
REG=pd.read_csv(ROOT/'condition_registry.csv')
OUT=ROOT/'results_full_capacity'; OUT.mkdir(parents=True,exist_ok=True)
PDBDIR=OUT/'pdb_cache'; PDBDIR.mkdir(exist_ok=True)
SASA_POINTS=120
GNM_CUTOFF=10.0
RIN_CUTOFF=8.0
IONIC_MM=150.0
TEMP_K=298.0

# Preserve benchmark pH when explicitly established; otherwise screening assumption.
PH_EXPLICIT={'LYZ_SWCNT':7.4}
def cond_ph(cid): return PH_EXPLICIT.get(cid,7.0)
def ph_source(cid): return 'benchmark-explicit' if cid in PH_EXPLICIT else 'screening-assumption-not-source-frozen'

CHANNEL_EXPANSION={
 'anionic':['anionic'],'cationic':['cationic'],'hydrophobic':['hydrophobic'],'pi_carbon':['pi_carbon'],'oxide':['oxide'],
 'oxide_or_hbond':['oxide','hbond_donor','hbond_acceptor'],
 'pi_carbon_or_hydrophobic':['pi_carbon','hydrophobic'],
 'hbond_donor_acceptor':['hbond_donor','hbond_acceptor'],
 'cationic_or_hbond':['cationic','hbond_donor','hbond_acceptor'],
}

def ensure_pdb(pid):
 p=PDBDIR/f'{pid}.pdb'
 if not p.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',p)
 return p

def pct(x):
 x=np.asarray(x,float)
 if len(x)<=1:return np.ones_like(x)
 return (rankdata(x,method='average')-1)/(len(x)-1)

def key(ch,r): return f"{ch.id}:{int(r.id[1])}:{str(r.id[2]).strip()}"

def structural(pid):
 p=ensure_pdb(pid); s=PDBParser(QUIET=True).get_structure(pid,str(p)); model=next(s.get_models())
 keys=[]; X=[]
 for ch in model:
  for r in ch:
   if 'CA' in r: keys.append(key(ch,r)); X.append(np.asarray(r['CA'].coord,float))
 X=np.vstack(X); n=len(X); centroid=X.mean(0)
 # GNM
 K=np.zeros((n,n),float); dist=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2)
 for i in range(n):
  js=np.where((dist[i]<=GNM_CUTOFF)&(dist[i]>0))[0]
  for j in js:
   if j>i: K[i,j]=K[j,i]=-1.; K[i,i]+=1.; K[j,j]+=1.
 vals,vecs=np.linalg.eigh(K); pos=vals>1e-8
 cov=(vecs[:,pos]*(1/vals[pos]))@vecs[:,pos].T if pos.any() else np.zeros((n,n))
 msf=np.diag(cov); mobility=(msf-msf.min())/(msf.max()-msf.min()) if msf.max()>msf.min() else np.zeros(n)
 radial=np.linalg.norm(X-centroid,axis=1)
 # RIN C-alpha unweighted
 A=((dist<=RIN_CUTOFF)&(dist>0)).astype(int); degree=A.sum(1).astype(float)
 INF=1e9; SP=np.where(A>0,1.,INF); np.fill_diagonal(SP,0.)
 for k0 in range(n): SP=np.minimum(SP,SP[:,k0,None]+SP[None,k0,:])
 clos=np.zeros(n)
 for i in range(n):
  f=(SP[i]<INF)&(np.arange(n)!=i); clos[i]=f.sum()/SP[i,f].sum() if f.any() and SP[i,f].sum()>0 else 0
 bet=np.zeros(n)
 for ss in range(n):
  S=[]; P=[[] for _ in range(n)]; sig=np.zeros(n); sig[ss]=1; dd=-np.ones(n,int); dd[ss]=0; Q=deque([ss])
  while Q:
   v=Q.popleft(); S.append(v)
   for w in np.where(A[v]>0)[0]:
    if dd[w]<0: Q.append(w); dd[w]=dd[v]+1
    if dd[w]==dd[v]+1: sig[w]+=sig[v]; P[w].append(v)
  delta=np.zeros(n)
  while S:
   w=S.pop()
   if sig[w]>0:
    for v in P[w]: delta[v]+=(sig[v]/sig[w])*(1+delta[w])
   if w!=ss: bet[w]+=delta[w]
 bet*=.5
 if n>2: bet/=((n-1)*(n-2)/2)
 return {keys[i]:{'gnm_mobility':float(mobility[i]),'radial_prominence_A':float(radial[i]),'rin_degree':float(degree[i]),'rin_closeness':float(clos[i]),'rin_betweenness':float(bet[i])} for i in range(n)}

def condition_features(pid,chem,ph):
 p=ensure_pdb(pid); old=ism.SASA_POINTS; ism.SASA_POINTS=SASA_POINTS
 try:
  _,res,atoms,_=ism.build_surface_residues(p,ph)
  env=ism.EnvParams(pH=ph,ionic=IONIC_MM,temp=TEMP_K)
  with tempfile.TemporaryDirectory(prefix='is_full15_') as td:
   apbs_status=ism.attach_apbs_auxiliary(p,res,atoms,env,Path(td))
 finally: ism.SASA_POINTS=old
 surf=[r for r in res if r['surface_exposed']]; D=ism.build_distances(surf); mp=ism.chemistry_map(surf,D,chem,ph)
 core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']}
 expected=ism.CHEMISTRIES[chem].get('expected_phi_sign')
 rows=[]
 for r in surf:
  phi=r.get('phi'); comp=None
  if expected=='positive' and phi is not None: comp=max(float(phi),0.)
  elif expected=='negative' and phi is not None: comp=max(-float(phi),0.)
  rows.append({'key':r['key'],'res_name':r['res_name'],'res_seq':r['res_seq'],'chain':r['chain'],'scrsa':r['scrsa'],'core_patch':core.get(r['key'],0.),'phi':phi,'apbs_compat':comp})
 return pd.DataFrame(rows),apbs_status,expected

struct_cache={}; feature_cache={}; allrows=[]; condsummary=[]
for _,c in REG.iterrows():
 pid=c['pdb']; cid=c['condition_id']; ph=cond_ph(cid)
 if pid not in struct_cache: struct_cache[pid]=structural(pid)
 chans=CHANNEL_EXPANSION.get(str(c['chemistry_channel']),[])
 if not chans:
  condsummary.append({'condition_id':cid,'protein':c['protein'],'pdb':pid,'channel_variant':'UNMAPPED','status':'not-computed-unmapped-chemistry','ground_truth_status':c['ground_truth_status'],'package_status':c['package_status']}); continue
 for chem in chans:
  cachekey=(pid,chem,ph)
  if cachekey not in feature_cache: feature_cache[cachekey]=condition_features(pid,chem,ph)
  df,apbs_status,expected=feature_cache[cachekey]
  st=struct_cache[pid]
  for nm in ['gnm_mobility','radial_prominence_A','rin_degree','rin_closeness','rin_betweenness']:
   df[nm]=[st.get(k,{}).get(nm,np.nan) for k in df['key']]
  core_r=pct(df['core_patch'].to_numpy(float)); gnm_r=pct(df['gnm_mobility'].fillna(0).to_numpy(float)); radial_r=pct(df['radial_prominence_A'].fillna(0).to_numpy(float))
  if expected in ('positive','negative') and df['apbs_compat'].notna().any(): apbs_r=pct(df['apbs_compat'].fillna(0).to_numpy(float)); apbs_used=True
  else: apbs_r=np.ones(len(df)); apbs_used=False
  raw=core_r*gnm_r*radial_r*apbs_r
  score=100*raw/raw.max() if len(raw) and raw.max()>0 else np.zeros(len(raw))
  df['interfacescout_score']=score
  df['rank']=pd.Series(score).rank(ascending=False,method='min').astype(int)
  df['condition_id']=cid; df['protein']=c['protein']; df['pdb']=pid; df['surface_or_condition']=c['surface_or_condition']; df['channel_variant']=chem
  df['tier']=c['tier']; df['ground_truth_type']=c['ground_truth_type']; df['ground_truth_status']=c['ground_truth_status']; df['analysis_role']=c['analysis_role']; df['package_status']=c['package_status']; df['pH']=ph; df['pH_source']=ph_source(cid); df['ionic_mM']=IONIC_MM; df['apbs_status']=apbs_status; df['apbs_used_in_score']=apbs_used
  allrows.append(df)
  top=df.sort_values(['interfacescout_score','key'],ascending=[False,True]).head(10)
  condsummary.append({'condition_id':cid,'protein':c['protein'],'pdb':pid,'surface_or_condition':c['surface_or_condition'],'channel_variant':chem,'tier':c['tier'],'ground_truth_status':c['ground_truth_status'],'package_status':c['package_status'],'n_surface_res':len(df),'pH':ph,'pH_source':ph_source(cid),'apbs_status':apbs_status,'apbs_used_in_score':apbs_used,'top1':top.iloc[0]['key'] if len(top) else None,'top1_score':float(top.iloc[0]['interfacescout_score']) if len(top) else None,'top10_keys':';'.join(top['key'].astype(str))})

ALL=pd.concat(allrows,ignore_index=True) if allrows else pd.DataFrame(); SUM=pd.DataFrame(condsummary)
ALL.to_csv(OUT/'all15_all_conditions_residue_scores.csv',index=False); SUM.to_csv(OUT/'all15_condition_score_summary.csv',index=False)
# Coverage by protein and status; a true full-capacity run means every PDB has generated scores for every mappable registry condition.
coverage=SUM.groupby(['protein','pdb']).agg(n_condition_variants=('condition_id','size'),n_unique_conditions=('condition_id','nunique'),n_success=('status',lambda x:0) if 'status' in SUM.columns else ('condition_id','size')).reset_index()
coverage.to_csv(OUT/'all15_protein_coverage.csv',index=False)
meta={'n_registry_conditions':int(REG.condition_id.nunique()),'n_proteins':int(REG.pdb.nunique()),'n_computed_condition_variants':int(len(SUM[SUM.channel_variant!='UNMAPPED'])),'n_residue_score_rows':int(len(ALL)),'score_definition':'percentile(core patch persistence) * percentile(GNM mobility) * conditional percentile(APBS compatibility) * percentile(radial prominence), normalized to 0-100 within condition/channel','RIN_role':'computed and exported as auxiliary candidate descriptors; not included in primary IS score','pH_policy':'LYZ_SWCNT 7.4 from existing benchmark; all other conditions pH 7.0 screening assumption until source-condition audit freezes values','ionic_strength_mM':IONIC_MM,'sasa_points_exploratory':SASA_POINTS}
(OUT/'all15_full_capacity_metadata.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2)); print(SUM[['condition_id','protein','channel_variant','n_surface_res','top1','apbs_status']].to_string(index=False))
