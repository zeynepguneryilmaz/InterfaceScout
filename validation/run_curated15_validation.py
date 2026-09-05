from __future__ import annotations
import json, sys, tempfile, urllib.request
from pathlib import Path
from collections import deque
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import backend.main as ism

ROOT=Path('validation/packages/all15')
OUT=ROOT/'results_curated_validation'; OUT.mkdir(parents=True,exist_ok=True)
PDBDIR=OUT/'pdb_cache'; PDBDIR.mkdir(exist_ok=True)
REG=pd.read_csv(ROOT/'condition_registry.csv')
ASM=pd.read_csv(ROOT/'curation'/'assembly_map.csv')
SRC=pd.read_csv(ROOT/'curation'/'source_condition_audit.csv')

# Canonical structural settings for validation rerun.
SASA_POINTS=200; GNM_CUTOFF=10.0; RIN_CUTOFF=8.0; IONIC_MM_DEFAULT=150.0; TEMP_K_DEFAULT=298.0

# Curated chain and residue filters. WNV source studies only domain III.
CHAIN_RULES={r.pdb:set(str(r.curated_polymer_chains).split(';')) for _,r in ASM.iterrows()}
RESIDUE_RANGE={'2HG0':(8,109)}

# Source-matched pH where verified; otherwise development assumption retained but explicitly tagged.
PH_RULES={'LYZ_SWCNT':7.4,'ACHE_COOH_SAM':7.0,'ACHE_NH2_SAM':7.0}

def cond_ph(cid): return PH_RULES.get(cid,7.0)
def ph_status(cid): return 'source-matched' if cid in PH_RULES else 'development-assumption-not-source-frozen'

CHANNEL_EXPANSION={'anionic':['anionic'],'cationic':['cationic'],'hydrophobic':['hydrophobic'],'pi_carbon':['pi_carbon'],'oxide':['oxide'],'oxide_or_hbond':['oxide','hbond_donor','hbond_acceptor'],'pi_carbon_or_hydrophobic':['pi_carbon','hydrophobic'],'hbond_donor_acceptor':['hbond_donor','hbond_acceptor'],'cationic_or_hbond':['cationic','hbond_donor','hbond_acceptor']}

# Quantitative residue-level ground truth. Existing 10 positive conditions + WNV frequent-contact regions.
# WNV paper original E numbering 346-352 and 397-400 maps to 2HG0 PDB numbering 55-61 and 106-109.
BENCH=[
 {'id':'FN_COO','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'anionic','pH':7.0,'truth':[1469],'truth_definition':'final anchor residue'},
 {'id':'FN_NH3','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'cationic','pH':7.0,'truth':[1312,1509],'truth_definition':'final anchor residues'},
 {'id':'FN_CH3_HEAD','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1454,1455,1456,1457,1478,1479,1480,1481,1509],'truth_definition':'orientation-specific explicit contacts'},
 {'id':'FN_CH3_SIDE','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1275,1276,1355,1376,1378,1454,1455,1456,1457],'truth_definition':'orientation-specific explicit contacts'},
 {'id':'FN_CH3_BETA','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1431,1446,1457,1458,1459,1461,1463,1464,1466,1475,1476],'truth_definition':'orientation-specific explicit contacts'},
 {'id':'CYTC_CH3','protein':'Cytochrome c','pdb':'3NWV','chem':'hydrophobic','pH':7.0,'truth':[1,2,3,4,96,99,100],'truth_definition':'explicit residue contacts'},
 {'id':'CYTC_COOH','protein':'Cytochrome c','pdb':'3NWV','chem':'anionic','pH':7.0,'truth':[1,2,4,61,99,100,103],'truth_definition':'explicit residue contacts'},
 {'id':'LYZ_SWCNT','protein':'Lysozyme','pdb':'1LYZ','chem':'pi_carbon','pH':7.4,'truth':[61,62,63,68,69,70,73,103,106,107,111,112,113,116],'truth_definition':'14 explicit MD contact residues'},
 {'id':'CHT_CNT_HYDRO','protein':'Alpha-chymotrypsin','pdb':'4CHA','chem':'hydrophobic','pH':7.0,'truth':[3,5,6,7,8,10,77,114,115,116],'truth_definition':'favorable binding/contact residues'},
 {'id':'MB_CIT_AUNP','protein':'Myoglobin','pdb':'1MBN','chem':'anionic','pH':7.0,'truth':[43,45,96,97,98,99,146,147,148,149,150,151,152,153],'truth_definition':'persistent CG-MD binding residues'},
 {'id':'WNV_GRAPHENE','protein':'WNV E protein domain III','pdb':'2HG0','chem':'pi_carbon','pH':7.0,'truth':[55,56,57,58,59,60,61,106,107,108,109],'truth_definition':'source-reported >35% frequent-contact regions; original E residues 346-352 and 397-400 mapped to 2HG0 domain-III numbering'},
]
NEG={'id':'CHT_CNT_PI_NEG','protein':'Alpha-chymotrypsin','pdb':'4CHA','chem':'pi_carbon','pH':7.0,'truth':[3,5,6,7,8,10,77,114,115,116],'truth_definition':'intentional chemistry-mismatch negative control'}

class CuratedSelect(Select):
 def __init__(self,pid,modelid): self.pid=pid; self.modelid=modelid
 def accept_model(self,m): return int(m.id==self.modelid)
 def accept_chain(self,c): return int(str(c.id) in CHAIN_RULES[self.pid])
 def accept_residue(self,r):
  if not is_aa(r,standard=True): return 0
  if self.pid in RESIDUE_RANGE:
   lo,hi=RESIDUE_RANGE[self.pid]; return int(lo<=int(r.id[1])<=hi)
  return 1

def getpdb(pid):
 raw=PDBDIR/f'{pid}_raw.pdb'; clean=PDBDIR/f'{pid}_curated.pdb'
 if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',raw)
 if not clean.exists():
  s=PDBParser(QUIET=True).get_structure(pid,str(raw)); m=next(s.get_models()); io=PDBIO(); io.set_structure(s); io.save(str(clean),CuratedSelect(pid,m.id))
 return clean

def pct(x):
 x=np.asarray(x,float)
 if len(x)<=1:return np.ones_like(x)
 return (rankdata(x,method='average')-1)/(len(x)-1)

def auc(y,s):
 p=s[y==1]; n=s[y==0]
 if not len(p) or not len(n): return np.nan
 return float(np.mean([(a>b)+0.5*(a==b) for a in p for b in n]))
def ap(y,s):
 if not y.sum():return np.nan
 o=np.argsort(-s,kind='stable'); yy=y[o]; pr=np.cumsum(yy)/(np.arange(len(yy))+1); return float((pr*yy).sum()/yy.sum())
def rec(y,s,k):
 if not y.sum():return np.nan
 o=np.argsort(-s,kind='stable')[:min(k,len(s))]; return float(y[o].sum()/y.sum())
def spatial(truthkeys,topkeys,coords,R):
 vals=[]
 for t in truthkeys:
  if t not in coords:continue
  ds=[np.linalg.norm(coords[t]-coords[q]) for q in topkeys if q in coords]
  vals.append(bool(ds and min(ds)<=R))
 return float(np.mean(vals)) if vals else np.nan
def nearest(truthkeys,topkeys,coords):
 vals=[]
 for t in truthkeys:
  if t not in coords:continue
  ds=[np.linalg.norm(coords[t]-coords[q]) for q in topkeys if q in coords]
  if ds:vals.append(min(ds))
 return float(np.median(vals)) if vals else np.nan

def structural(pid):
 p=getpdb(pid); s=PDBParser(QUIET=True).get_structure(pid,str(p)); m=next(s.get_models()); keys=[]; X=[]
 for ch in m:
  for r in ch:
   if 'CA' in r: keys.append(f'{ch.id}:{int(r.id[1])}:{str(r.id[2]).strip()}'); X.append(np.asarray(r['CA'].coord,float))
 X=np.vstack(X); n=len(X); ctr=X.mean(0); D=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2)
 K=np.zeros((n,n))
 for i in range(n):
  for j in np.where((D[i]<=GNM_CUTOFF)&(D[i]>0))[0]:
   if j>i:K[i,j]=K[j,i]=-1.;K[i,i]+=1.;K[j,j]+=1.
 vals,V=np.linalg.eigh(K); pos=vals>1e-8; C=(V[:,pos]*(1/vals[pos]))@V[:,pos].T if pos.any() else np.zeros((n,n)); msf=np.diag(C); mob=(msf-msf.min())/(msf.max()-msf.min()) if msf.max()>msf.min() else np.zeros(n)
 radial=np.linalg.norm(X-ctr,axis=1)
 A=((D<=RIN_CUTOFF)&(D>0)).astype(int); degree=A.sum(1).astype(float)
 INF=1e9; SP=np.where(A>0,1.,INF); np.fill_diagonal(SP,0.)
 for k in range(n):SP=np.minimum(SP,SP[:,k,None]+SP[None,k,:])
 clos=np.zeros(n)
 for i in range(n):
  f=(SP[i]<INF)&(np.arange(n)!=i);clos[i]=f.sum()/SP[i,f].sum() if f.any() and SP[i,f].sum()>0 else 0
 bet=np.zeros(n)
 for ss in range(n):
  S=[];P=[[] for _ in range(n)];sig=np.zeros(n);sig[ss]=1;dd=-np.ones(n,int);dd[ss]=0;Q=deque([ss])
  while Q:
   v=Q.popleft();S.append(v)
   for w in np.where(A[v]>0)[0]:
    if dd[w]<0:Q.append(w);dd[w]=dd[v]+1
    if dd[w]==dd[v]+1:sig[w]+=sig[v];P[w].append(v)
  delta=np.zeros(n)
  while S:
   w=S.pop()
   if sig[w]>0:
    for v in P[w]:delta[v]+=(sig[v]/sig[w])*(1+delta[w])
   if w!=ss:bet[w]+=delta[w]
 bet*=.5
 if n>2:bet/=((n-1)*(n-2)/2)
 return {keys[i]:{'gnm':float(mob[i]),'radial':float(radial[i]),'degree':float(degree[i]),'closeness':float(clos[i]),'betweenness':float(bet[i])} for i in range(n)}

def environment(pid,ph,ionic=IONIC_MM_DEFAULT,temp=TEMP_K_DEFAULT):
 p=getpdb(pid);old=ism.SASA_POINTS;ism.SASA_POINTS=SASA_POINTS
 try:
  _,res,atoms,_=ism.build_surface_residues(p,ph);env=ism.EnvParams(pH=ph,ionic=ionic,temp=temp)
  with tempfile.TemporaryDirectory(prefix='is_curated_') as td:status=ism.attach_apbs_auxiliary(p,res,atoms,env,Path(td))
 finally:ism.SASA_POINTS=old
 surf=[r for r in res if r['surface_exposed']];return surf,ism.build_distances(surf),status

SCACHE={};ECACHE={}
def score_condition(pid,chem,ph):
 if pid not in SCACHE:SCACHE[pid]=structural(pid)
 ek=(pid,ph)
 if ek not in ECACHE:ECACHE[ek]=environment(pid,ph)
 surf,D,status=ECACHE[ek];mp=ism.chemistry_map(surf,D,chem,ph);core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']};st=SCACHE[pid];expected=ism.CHEMISTRIES[chem].get('expected_phi_sign')
 rows=[]
 for r in surf:
  phi=r.get('phi');comp=None
  if expected=='positive' and phi is not None:comp=max(float(phi),0.)
  elif expected=='negative' and phi is not None:comp=max(-float(phi),0.)
  z={'key':r['key'],'res_name':r['res_name'],'res_seq':int(r['res_seq']),'chain':r['chain'],'scrsa':float(r['scrsa']),'core':float(core.get(r['key'],0.)),'phi':phi,'apbs_compat':comp,'x':r['x'],'y':r['y'],'z':r['z']};z.update(st.get(r['key'],{}));rows.append(z)
 df=pd.DataFrame(rows);cr=pct(df.core);gr=pct(df.gnm);rr=pct(df.radial)
 if expected in ('positive','negative') and df.apbs_compat.notna().any():ar=pct(df.apbs_compat.fillna(0));apbs_used=True
 else:ar=np.ones(len(df));apbs_used=False
 dr=1-pct(df.degree); clr=1-pct(df.closeness); br=1-pct(df.betweenness)
 methods={'Core':cr,'Core_GNM':cr*gr,'Core_APBS':cr*ar,'Core_GNM_APBS':cr*gr*ar,'IS_candidate':cr*gr*ar*rr,'IS_RIN_degree':cr*gr*ar*rr*dr,'IS_RIN_closeness':cr*gr*ar*rr*clr,'IS_RIN_betweenness':cr*gr*ar*rr*br}
 for name,v in methods.items():df[name]=100*v/v.max() if len(v) and np.max(v)>0 else np.zeros(len(v))
 return df,status,apbs_used

# Curated full 15-protein screening across registry.
full=[];fullsum=[]
for _,c in REG.iterrows():
 pid=c.pdb;cid=c.condition_id;ph=cond_ph(cid);chans=CHANNEL_EXPANSION.get(str(c.chemistry_channel),[])
 for chem in chans:
  df,status,apbs_used=score_condition(pid,chem,ph)
  for col,val in {'condition_id':cid,'protein':c.protein,'pdb':pid,'surface_or_condition':c.surface_or_condition,'channel_variant':chem,'pH':ph,'pH_status':ph_status(cid),'ionic_mM':IONIC_MM_DEFAULT,'apbs_status':status,'apbs_used':apbs_used,'ground_truth_status':c.ground_truth_status,'analysis_role':c.analysis_role}.items():df[col]=val
  full.append(df);top=df.sort_values('IS_candidate',ascending=False).head(10)
  fullsum.append({'condition_id':cid,'protein':c.protein,'pdb':pid,'channel_variant':chem,'n_surface':len(df),'pH':ph,'pH_status':ph_status(cid),'apbs_status':status,'apbs_used':apbs_used,'top10':';'.join(top.key.astype(str))})
FULL=pd.concat(full,ignore_index=True);FULL.to_csv(OUT/'curated_all15_residue_scores.csv',index=False);pd.DataFrame(fullsum).to_csv(OUT/'curated_all15_condition_summary.csv',index=False)

# Quantitative ablation on curated explicit residue-level set.
metrics=[];qdetail=[]
for b in BENCH+[NEG]:
 df,status,apbs_used=score_condition(b['pdb'],b['chem'],b['pH']);truthkeys=[k for k,r in zip(df.key,df.res_seq) if int(r) in set(b['truth'])];y=df.key.isin(truthkeys).astype(int).to_numpy();coords={r.key:np.array([r.x,r.y,r.z],float) for _,r in df.iterrows()}
 for method in ['Core','Core_GNM','Core_APBS','Core_GNM_APBS','IS_candidate','IS_RIN_degree','IS_RIN_closeness','IS_RIN_betweenness']:
  s=df[method].to_numpy(float);order=np.argsort(-s,kind='stable');top=df.key.iloc[order[:10]].tolist();metrics.append({'condition_id':b['id'],'protein':b['protein'],'pdb':b['pdb'],'method':method,'n_surface':len(df),'n_truth_requested':len(b['truth']),'n_truth_mapped':int(y.sum()),'AUROC':auc(y,s),'AP':ap(y,s),'Recall@5':rec(y,s,5),'Recall@10':rec(y,s,10),'Recall@20':rec(y,s,20),'Spatial@5A_top10':spatial(truthkeys,top,coords,5),'Spatial@8A_top10':spatial(truthkeys,top,coords,8),'Spatial@10A_top10':spatial(truthkeys,top,coords,10),'median_nearest_top10_A':nearest(truthkeys,top,coords),'truth_definition':b['truth_definition'],'apbs_status':status,'apbs_used':apbs_used,'is_negative_control':b['id']==NEG['id']})
 df['truth']=y;df['condition_id']=b['id'];df['protein']=b['protein'];df['truth_definition']=b['truth_definition'];qdetail.append(df)
MET=pd.DataFrame(metrics);MET.to_csv(OUT/'curated_quantitative_ablation.csv',index=False);pd.concat(qdetail,ignore_index=True).to_csv(OUT/'curated_quantitative_residue_scores.csv',index=False)

# Primary aggregate excludes intentional negative control.
PRI=MET[~MET.is_negative_control].copy();mc=['AUROC','AP','Recall@5','Recall@10','Recall@20','Spatial@5A_top10','Spatial@8A_top10','Spatial@10A_top10','median_nearest_top10_A']
agg=[]
for method,g in PRI.groupby('method'):
 row={'method':method,'n_conditions':g.condition_id.nunique(),'n_proteins':g.protein.nunique()}
 for m in mc:row['condition_median_'+m]=g[m].median();row['condition_mean_'+m]=g[m].mean();row['protein_macro_median_'+m]=g.groupby('protein')[m].median().median();row['protein_macro_mean_'+m]=g.groupby('protein')[m].median().mean()
 agg.append(row)
AGG=pd.DataFrame(agg);AGG.to_csv(OUT/'curated_quantitative_aggregate.csv',index=False)

# Protein-clustered bootstrap: sample proteins with replacement, each sampled protein contributes its condition median.
rng=np.random.default_rng(20260905); proteins=sorted(PRI.protein.unique());B=10000;boot=[]
for method,g in PRI.groupby('method'):
 pm=g.groupby('protein')[mc].median()
 vals={m:[] for m in mc}
 for _ in range(B):
  samp=rng.choice(proteins,size=len(proteins),replace=True)
  for m in mc:vals[m].append(float(np.median([pm.loc[p,m] for p in samp])))
 for m in mc:
  a=np.asarray(vals[m]);boot.append({'method':method,'metric':m,'bootstrap_protein_macro_median':float(np.median(a)),'ci2.5':float(np.quantile(a,.025)),'ci97.5':float(np.quantile(a,.975)),'n_proteins':len(proteins),'B':B})
pd.DataFrame(boot).to_csv(OUT/'curated_clustered_bootstrap.csv',index=False)

# Paired protein-level deltas vs Core.
base=PRI[PRI.method=='Core'].groupby('protein')[mc].median();drows=[]
for method in [x for x in PRI.method.unique() if x!='Core']:
 cur=PRI[PRI.method==method].groupby('protein')[mc].median()
 for m in mc:
  d=(cur[m]-base[m]).dropna();arr=[]
  for _ in range(B):arr.append(float(np.median(rng.choice(d.to_numpy(),size=len(d),replace=True))))
  arr=np.asarray(arr);drows.append({'method':method,'metric':m,'observed_median_delta_vs_Core':float(np.median(d)),'proteins_improved':int((d>0).sum()),'proteins_worsened':int((d<0).sum()),'proteins_tied':int((d==0).sum()),'bootstrap_ci2.5':float(np.quantile(arr,.025)),'bootstrap_ci97.5':float(np.quantile(arr,.975))})
pd.DataFrame(drows).to_csv(OUT/'curated_component_deltas_by_protein.csv',index=False)

# LOO reporting (not fitted CV): report left-out protein performance and aggregate of remaining proteins.
loo=[]
for left in proteins:
 for method,g in PRI.groupby('method'):
  hold=g[g.protein==left];train=g[g.protein!=left]
  for m in mc:loo.append({'left_out_protein':left,'method':method,'metric':m,'holdout_condition_median':hold[m].median(),'remaining_protein_macro_median':train.groupby('protein')[m].median().median(),'note':'No coefficients are fitted; this is leave-one-protein-out robustness reporting, not training CV.'})
pd.DataFrame(loo).to_csv(OUT/'curated_leave_one_protein_out.csv',index=False)

meta={'date':'2026-09-05','canonical_sasa_points':SASA_POINTS,'gnm_cutoff_A':GNM_CUTOFF,'rin_cutoff_A':RIN_CUTOFF,'default_ionic_mM_for_unspecified_sources':IONIC_MM_DEFAULT,'n_registry_conditions':int(REG.condition_id.nunique()),'n_full_proteins':int(REG.pdb.nunique()),'n_quantitative_primary_conditions':len(BENCH),'n_quantitative_primary_proteins':len(set(b['protein'] for b in BENCH)),'negative_controls':1,'WNV_mapping':'2HG0 PDB residues 8-109 = original E residues 299-400; frequent-contact GT PDB 55-61 and 106-109 = original 346-352 and 397-400','important':'This is model-development evidence because benchmark outcomes were inspected before feature freeze. Conditions with non-source-frozen pH/ionic strength remain developmental, not final independent validation.'}
(OUT/'CURATED_VALIDATION_METADATA.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2));print(AGG[['method','n_conditions','n_proteins','protein_macro_median_AUROC','protein_macro_median_AP','protein_macro_median_Recall@10','protein_macro_median_Recall@20']].to_string(index=False))
