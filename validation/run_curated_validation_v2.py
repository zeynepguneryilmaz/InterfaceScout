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
OUT=ROOT/'results_curated_validation_v2'; OUT.mkdir(parents=True,exist_ok=True)
PDBDIR=OUT/'pdb_cache'; PDBDIR.mkdir(exist_ok=True)
SASA_POINTS=200; GNM_CUTOFF=10.0; RIN_CUTOFF=8.0

# Authoritative structure curation for quantitative/stress systems.
CHAIN_RULES={'1FNF':{'A'},'3NWV':{'A'},'1LYZ':{'A'},'4CHA':{'A','B','C'},'1MBN':{'A'},'2HG0':{'A'},'7RSA':{'A'}}
RANGE_RULES={'1FNF':(1236,1509),'2HG0':(8,109)}

# Primary residue-level development set. No repeated condition is treated as a new protein.
PRIMARY=[
 {'id':'FN_COO','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'anionic','pH':7.0,'ionic':1000.0,'temp':300.0,'truth':[1469],'source':'Liamas et al. IJMS 2018 doi:10.3390/ijms19113321','condition_note':'pH is development assumption; 1 M NaCl source charged-surface trajectory'},
 {'id':'FN_NH3','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'cationic','pH':7.0,'ionic':800.0,'temp':300.0,'truth':[1312,1509],'source':'Liamas et al. IJMS 2018 doi:10.3390/ijms19113321','condition_note':'pH is development assumption; 0.8 M NaCl source charged-surface trajectory'},
 {'id':'FN_CH3_HEAD','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'ionic':50.0,'temp':300.0,'truth':[1454,1455,1456,1457,1478,1479,1480,1481,1509],'source':'Liamas et al. IJMS 2018 doi:10.3390/ijms19113321','condition_note':'source construct 1236-1509; 0.05 M base solvent setup'},
 {'id':'FN_CH3_SIDE','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'ionic':50.0,'temp':300.0,'truth':[1275,1276,1355,1376,1378,1454,1455,1456,1457],'source':'Liamas et al. IJMS 2018 doi:10.3390/ijms19113321','condition_note':'source construct 1236-1509; 0.05 M base solvent setup'},
 {'id':'FN_CH3_BETA','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'ionic':50.0,'temp':300.0,'truth':[1431,1446,1457,1458,1459,1461,1463,1464,1466,1475,1476],'source':'Liamas et al. IJMS 2018 doi:10.3390/ijms19113321','condition_note':'source construct 1236-1509; 0.05 M base solvent setup'},
 {'id':'CYTC_CH3','protein':'Cytochrome c','pdb':'3NWV','chem':'hydrophobic','pH':7.4,'ionic':150.0,'temp':298.0,'truth':[1,2,3,4,96,99,100],'source':'Sun et al. PLOS ONE 2014 doi:10.1371/journal.pone.0107696','condition_note':'pH source-aligned physiological; ionic value is continuum reference, not source bulk molarity'},
 {'id':'CYTC_COOH','protein':'Cytochrome c','pdb':'3NWV','chem':'anionic','pH':7.4,'ionic':150.0,'temp':298.0,'truth':[1,2,4,61,99,100,103],'source':'Sun et al. PLOS ONE 2014 doi:10.1371/journal.pone.0107696','condition_note':'pH source-aligned; APBS ionic reference tested separately for sensitivity'},
 {'id':'LYZ_SWCNT','protein':'Lysozyme','pdb':'1LYZ','chem':'pi_carbon','pH':7.4,'ionic':150.0,'temp':298.0,'truth':[61,62,63,68,69,70,73,103,106,107,111,112,113,116],'source':'Scientific Reports 2025 doi:10.1038/s41598-025-96435-3','condition_note':'APBS is neutral for pi-carbon channel'},
 {'id':'CHT_CNT_HYDRO','protein':'Alpha-chymotrypsin','pdb':'4CHA','chem':'hydrophobic','pH':7.0,'ionic':150.0,'temp':300.0,'truth':[3,5,6,7,8,10,77,114,115,116],'source':'Scientific Reports 2015 doi:10.1038/srep09297','condition_note':'chains A+B+C form one enzyme; pH not source-frozen'},
 {'id':'MB_CIT_AUNP','protein':'Myoglobin','pdb':'1MBN','chem':'anionic','pH':7.0,'ionic':150.0,'temp':298.0,'truth':[43,45,96,97,98,99,146,147,148,149,150,151,152,153],'source':'Int J Mol Sci 2019 doi:10.3390/ijms20143539','condition_note':'citrate represented by generic anionic proxy; ionic/pH not source-frozen'},
 {'id':'WNV_GRAPHENE','protein':'WNV E protein domain III','pdb':'2HG0','chem':'pi_carbon','pH':7.0,'ionic':150.0,'temp':300.0,'truth':[55,56,57,58,59,60,61,106,107,108,109],'source':'Nanoscale Advances 2026 doi:10.1039/D5NA00988J','condition_note':'2HG0 PDB residues 8-109 correspond to original E 299-400; GT is >35% frequent-contact regions'},
 {'id':'RNASE_SILICA_4NM','protein':'RNase A','pdb':'7RSA','chem':'oxide','pH':7.4,'ionic':150.0,'temp':298.0,'truth':[15,16,17,18,21,50,51,52,53,55],'source':'Sun et al. PLOS ONE 2014 doi:10.1371/journal.pone.0107696','condition_note':'exact key adsorption residues; radius not represented explicitly, but this condition is used as one residue-level development case'},
]

STRESS=[
 {'id':'CYTC_SILICA_4NM','protein':'Cytochrome c','pdb':'3NWV','chem':'oxide','pH':7.4,'ionic':150.0,'temp':298.0,'truth':[36,37,38,39,44,55,58,61,62,103],'source':'Sun et al. PLOS ONE 2014','note':'curvature stress; exact reported 4-nm contact residues'},
 {'id':'CYTC_SILICA_11NM','protein':'Cytochrome c','pdb':'3NWV','chem':'oxide','pH':7.4,'ionic':150.0,'temp':298.0,'truth':[92,95,96,99,100],'source':'Sun et al. PLOS ONE 2014','note':'curvature stress; only strong-contact residue subset encoded, broader 1-17 near-surface region omitted'},
 {'id':'RNASE_SILICA_11NM','protein':'RNase A','pdb':'7RSA','chem':'oxide','pH':7.4,'ionic':150.0,'temp':298.0,'truth':[3,9,13,15,16,51,52,55,114],'source':'Sun et al. PLOS ONE 2014','note':'curvature stress; exact directly interacting residues'},
]
NEG={'id':'CHT_CNT_PI_NEG','protein':'Alpha-chymotrypsin','pdb':'4CHA','chem':'pi_carbon','pH':7.0,'ionic':150.0,'temp':300.0,'truth':[3,5,6,7,8,10,77,114,115,116],'source':'derived from CHT_CNT_HYDRO truth','note':'intentional wrong chemistry negative control'}

class Sel(Select):
 def __init__(self,pid,mid):self.pid=pid;self.mid=mid
 def accept_model(self,m):return int(m.id==self.mid)
 def accept_chain(self,c):return int(str(c.id) in CHAIN_RULES[self.pid])
 def accept_residue(self,r):
  if not is_aa(r,standard=True):return 0
  if self.pid in RANGE_RULES:
   lo,hi=RANGE_RULES[self.pid];return int(lo<=int(r.id[1])<=hi)
  return 1

def getpdb(pid):
 raw=PDBDIR/f'{pid}_raw.pdb';clean=PDBDIR/f'{pid}_curated.pdb'
 if not raw.exists():urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',raw)
 if not clean.exists():
  s=PDBParser(QUIET=True).get_structure(pid,str(raw));m=next(s.get_models());io=PDBIO();io.set_structure(s);io.save(str(clean),Sel(pid,m.id))
 return clean

def pct(x):
 x=np.asarray(x,float)
 if len(x)<=1:return np.ones_like(x)
 return (rankdata(x,method='average')-1)/(len(x)-1)
def auc(y,s):
 p=s[y==1];n=s[y==0]
 if len(p)==0 or len(n)==0:return np.nan
 return float(np.mean([(a>b)+.5*(a==b) for a in p for b in n]))
def ap(y,s):
 if y.sum()==0:return np.nan
 o=np.argsort(-s,kind='stable');yy=y[o];pr=np.cumsum(yy)/(np.arange(len(yy))+1);return float((pr*yy).sum()/yy.sum())
def recall(y,s,k):
 if y.sum()==0:return np.nan
 o=np.argsort(-s,kind='stable')[:min(k,len(s))];return float(y[o].sum()/y.sum())
def spatial(truth,top,coords,R):
 z=[]
 for t in truth:
  if t not in coords:continue
  d=[np.linalg.norm(coords[t]-coords[q]) for q in top if q in coords];z.append(bool(d and min(d)<=R))
 return float(np.mean(z)) if z else np.nan
def nearest(truth,top,coords):
 z=[]
 for t in truth:
  if t not in coords:continue
  d=[np.linalg.norm(coords[t]-coords[q]) for q in top if q in coords]
  if d:z.append(min(d))
 return float(np.median(z)) if z else np.nan

def structure(pid):
 s=PDBParser(QUIET=True).get_structure(pid,str(getpdb(pid)));m=next(s.get_models());keys=[];X=[]
 for ch in m:
  for r in ch:
   if 'CA' in r:keys.append(f'{ch.id}:{int(r.id[1])}:{str(r.id[2]).strip()}');X.append(np.asarray(r['CA'].coord,float))
 X=np.vstack(X);n=len(X);D=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2);ctr=X.mean(0)
 K=np.zeros((n,n))
 for i in range(n):
  js=np.where((D[i]<=GNM_CUTOFF)&(D[i]>0))[0]
  for j in js:
   if j>i:K[i,j]=K[j,i]=-1.;K[i,i]+=1.;K[j,j]+=1.
 ev,V=np.linalg.eigh(K);pos=ev>1e-8;C=(V[:,pos]*(1/ev[pos]))@V[:,pos].T if pos.any() else np.zeros((n,n));msf=np.diag(C);mob=(msf-msf.min())/(msf.max()-msf.min()) if msf.max()>msf.min() else np.zeros(n)
 radial=np.linalg.norm(X-ctr,axis=1)
 A=((D<=RIN_CUTOFF)&(D>0)).astype(int);degree=A.sum(1).astype(float)
 # sparse BFS/Brandes rather than dense Floyd-Warshall; quantitative structures are small but this is scalable.
 clos=np.zeros(n);bet=np.zeros(n)
 for ss in range(n):
  S=[];P=[[] for _ in range(n)];sig=np.zeros(n);sig[ss]=1.;dist=-np.ones(n,int);dist[ss]=0;Q=deque([ss])
  while Q:
   v=Q.popleft();S.append(v)
   for w in np.where(A[v]>0)[0]:
    if dist[w]<0:dist[w]=dist[v]+1;Q.append(w)
    if dist[w]==dist[v]+1:sig[w]+=sig[v];P[w].append(v)
  f=dist>0;closeness_den=dist[f].sum();clos[ss]=f.sum()/closeness_den if closeness_den>0 else 0
  delta=np.zeros(n)
  while S:
   w=S.pop()
   if sig[w]>0:
    for v in P[w]:delta[v]+=(sig[v]/sig[w])*(1+delta[w])
   if w!=ss:bet[w]+=delta[w]
 bet*=.5
 if n>2:bet/=((n-1)*(n-2)/2)
 return {keys[i]:{'gnm':float(mob[i]),'radial':float(radial[i]),'rin_degree':float(degree[i]),'rin_closeness':float(clos[i]),'rin_betweenness':float(bet[i])} for i in range(n)},n

SC={};EC={}
def env(pid,ph,ionic,temp):
 key=(pid,ph,ionic,temp)
 if key in EC:return EC[key]
 old=ism.SASA_POINTS;ism.SASA_POINTS=SASA_POINTS
 try:
  _,res,atoms,_=ism.build_surface_residues(getpdb(pid),ph);e=ism.EnvParams(pH=ph,ionic=ionic,temp=temp)
  with tempfile.TemporaryDirectory(prefix='is_v2_') as td:status=ism.attach_apbs_auxiliary(getpdb(pid),res,atoms,e,Path(td))
 finally:ism.SASA_POINTS=old
 surf=[r for r in res if r['surface_exposed']];EC[key]=(surf,ism.build_distances(surf),status);return EC[key]

def score(case,ionic_override=None):
 pid=case['pdb'];ph=case['pH'];ionic=case['ionic'] if ionic_override is None else ionic_override;temp=case['temp']
 if pid not in SC:SC[pid]=structure(pid)[0]
 surf,D,status=env(pid,ph,ionic,temp);mp=ism.chemistry_map(surf,D,case['chem'],ph);core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']};expected=ism.CHEMISTRIES[case['chem']].get('expected_phi_sign');rows=[]
 for r in surf:
  phi=r.get('phi');comp=None
  if expected=='positive' and phi is not None:comp=max(float(phi),0.)
  elif expected=='negative' and phi is not None:comp=max(-float(phi),0.)
  z={'key':r['key'],'res_name':r['res_name'],'res_seq':int(r['res_seq']),'chain':r['chain'],'scrsa':float(r['scrsa']),'core_raw':float(core.get(r['key'],0.)),'phi':phi,'apbs_compat':comp,'x':r['x'],'y':r['y'],'z':r['z']};z.update(SC[pid].get(r['key'],{}));rows.append(z)
 df=pd.DataFrame(rows);cr=pct(df.core_raw);gr=pct(df.gnm);rr=pct(df.radial)
 if expected in ('positive','negative') and df.apbs_compat.notna().any():ar=pct(df.apbs_compat.fillna(0));apbs_used=True
 else:ar=np.ones(len(df));apbs_used=False
 dr=1-pct(df.rin_degree);clr=1-pct(df.rin_closeness);br=1-pct(df.rin_betweenness)
 vals={'Core':cr,'Core_GNM':cr*gr,'Core_APBS':cr*ar,'Core_GNM_APBS':cr*gr*ar,'IS_candidate':cr*gr*ar*rr,'IS_RIN_degree':cr*gr*ar*rr*dr,'IS_RIN_closeness':cr*gr*ar*rr*clr,'IS_RIN_betweenness':cr*gr*ar*rr*br}
 for k,v in vals.items():df[k]=100*v/v.max() if len(v) and np.max(v)>0 else np.zeros(len(v))
 return df,status,apbs_used,ionic

def evaluate(case,role,ionic_override=None):
 df,status,apbs_used,ionic=score(case,ionic_override);truth=set(case['truth']);truthkeys=df[df.res_seq.isin(truth)].key.tolist();y=df.key.isin(truthkeys).astype(int).to_numpy();coords={r.key:np.array([r.x,r.y,r.z],float) for _,r in df.iterrows()};rows=[]
 for method in ['Core','Core_GNM','Core_APBS','Core_GNM_APBS','IS_candidate','IS_RIN_degree','IS_RIN_closeness','IS_RIN_betweenness']:
  s=df[method].to_numpy(float);o=np.argsort(-s,kind='stable');top=df.key.iloc[o[:10]].tolist();rows.append({'condition_id':case['id'],'protein':case['protein'],'pdb':case['pdb'],'role':role,'method':method,'ionic_mM':ionic,'n_surface':len(df),'n_truth_requested':len(truth),'n_truth_mapped':int(y.sum()),'AUROC':auc(y,s),'AP':ap(y,s),'Recall@5':recall(y,s,5),'Recall@10':recall(y,s,10),'Recall@20':recall(y,s,20),'Spatial@5A_top10':spatial(truthkeys,top,coords,5),'Spatial@8A_top10':spatial(truthkeys,top,coords,8),'Spatial@10A_top10':spatial(truthkeys,top,coords,10),'median_nearest_top10_A':nearest(truthkeys,top,coords),'apbs_status':status,'apbs_used':apbs_used,'source':case['source'],'note':case.get('condition_note',case.get('note',''))})
 df['condition_id']=case['id'];df['protein']=case['protein'];df['role']=role;df['truth']=y;df['ionic_mM']=ionic;return rows,df

# Structure QC before metrics.
qc=[]
for pid in CHAIN_RULES:
 st,n=structure(pid);p=getpdb(pid);s=PDBParser(QUIET=True).get_structure(pid,str(p));m=next(s.get_models());chains=sorted({c.id for c in m});seqs=[int(r.id[1]) for c in m for r in c if is_aa(r,standard=True)]
 qc.append({'pdb':pid,'chains':';'.join(chains),'n_CA':n,'min_resseq':min(seqs) if seqs else None,'max_resseq':max(seqs) if seqs else None,'range_rule':str(RANGE_RULES.get(pid,'all'))})
pd.DataFrame(qc).to_csv(OUT/'structure_qc.csv',index=False)

allm=[];alld=[]
for c in PRIMARY:
 r,d=evaluate(c,'PRIMARY');allm+=r;alld.append(d)
for c in STRESS:
 r,d=evaluate(c,'STRESS');allm+=r;alld.append(d)
r,d=evaluate(NEG,'NEGATIVE_CONTROL');allm+=r;alld.append(d)
MET=pd.DataFrame(allm);MET.to_csv(OUT/'all_metrics.csv',index=False);pd.concat(alld,ignore_index=True).to_csv(OUT/'residue_scores.csv',index=False)
PRI=MET[MET.role=='PRIMARY'].copy();metrics=['AUROC','AP','Recall@5','Recall@10','Recall@20','Spatial@5A_top10','Spatial@8A_top10','Spatial@10A_top10','median_nearest_top10_A']
agg=[]
for method,g in PRI.groupby('method'):
 pm=g.groupby('protein')[metrics].median();row={'method':method,'n_conditions':g.condition_id.nunique(),'n_proteins':g.protein.nunique()}
 for x in metrics:
  row[f'condition_median_{x}']=g[x].median();row[f'condition_mean_{x}']=g[x].mean();row[f'protein_macro_median_{x}']=pm[x].median();row[f'protein_macro_mean_{x}']=pm[x].mean()
 agg.append(row)
AGG=pd.DataFrame(agg);AGG.to_csv(OUT/'primary_aggregate.csv',index=False)

# protein-clustered bootstrap
rng=np.random.default_rng(20260905);proteins=sorted(PRI.protein.unique());B=10000;boots=[]
for method,g in PRI.groupby('method'):
 pm=g.groupby('protein')[metrics].median()
 for x in metrics:
  arr=np.empty(B)
  for b in range(B):
   samp=rng.choice(proteins,size=len(proteins),replace=True);arr[b]=np.median([pm.loc[p,x] for p in samp])
  boots.append({'method':method,'metric':x,'protein_macro_median':float(pm[x].median()),'bootstrap_median':float(np.median(arr)),'ci2.5':float(np.quantile(arr,.025)),'ci97.5':float(np.quantile(arr,.975)),'B':B,'n_proteins':len(proteins)})
pd.DataFrame(boots).to_csv(OUT/'protein_clustered_bootstrap.csv',index=False)

# paired protein-level component deltas vs Core
base=PRI[PRI.method=='Core'].groupby('protein')[metrics].median();deltas=[]
for method in [m for m in PRI.method.unique() if m!='Core']:
 cur=PRI[PRI.method==method].groupby('protein')[metrics].median()
 for x in metrics:
  d=(cur[x]-base[x]).dropna();arr=np.empty(B)
  for b in range(B):arr[b]=np.median(rng.choice(d.to_numpy(),size=len(d),replace=True))
  deltas.append({'method':method,'metric':x,'observed_median_delta_vs_Core':float(np.median(d)),'proteins_improved':int((d>0).sum()),'proteins_worsened':int((d<0).sum()),'proteins_tied':int((d==0).sum()),'bootstrap_ci2.5':float(np.quantile(arr,.025)),'bootstrap_ci97.5':float(np.quantile(arr,.975))})
pd.DataFrame(deltas).to_csv(OUT/'component_deltas_by_protein.csv',index=False)

# leave-one-protein-out robustness reporting, no training/fitting
loo=[]
for left in proteins:
 for method,g in PRI.groupby('method'):
  hold=g[g.protein==left];rest=g[g.protein!=left]
  for x in metrics:loo.append({'left_out_protein':left,'method':method,'metric':x,'holdout_condition_median':hold[x].median(),'remaining_protein_macro_median':rest.groupby('protein')[x].median().median()})
pd.DataFrame(loo).to_csv(OUT/'leave_one_protein_out.csv',index=False)

# APBS ionic-strength sensitivity only for charged cases whose source molarity is not frozen.
sens=[]
for cid in ['CYTC_COOH','MB_CIT_AUNP']:
 c=next(x for x in PRIMARY if x['id']==cid)
 for ionic in [0.0,50.0,150.0,300.0]:
  rr,_=evaluate(c,'IONIC_SENSITIVITY',ionic)
  sens += [x for x in rr if x['method'] in ['Core','Core_APBS','Core_GNM_APBS','IS_candidate']]
pd.DataFrame(sens).to_csv(OUT/'apbs_ionic_sensitivity.csv',index=False)

meta={'version':'curated-validation-v2','date':'2026-09-05','primary_conditions':len(PRIMARY),'primary_unique_proteins':len(set(c['protein'] for c in PRIMARY)),'stress_conditions':len(STRESS),'negative_controls':1,'SASA_POINTS':SASA_POINTS,'GNM_cutoff_A':GNM_CUTOFF,'RIN_cutoff_A':RIN_CUTOFF,'construct_corrections':{'1FNF':'chain A residues 1236-1509','2HG0':'chain A PDB residues 8-109'},'FN_source':'doi:10.3390/ijms19113321','FN_ionic_mM':{'COO':1000,'NH3':800,'CH3':50},'PLOS_pH':7.4,'WNV_truth':'PDB 55-61,106-109 corresponding original E 346-352,397-400 >35% contact regions','RNase4_truth':'15,16,17,18,21,50,51,52,53,55','interpretation':'Model-development evidence, not prospective independent validation. No fitted coefficients. Protein-clustered summaries treat repeated conditions within one protein as non-independent.'}
(OUT/'METADATA.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2));print(AGG[['method','n_conditions','n_proteins','protein_macro_median_AUROC','protein_macro_median_AP','protein_macro_median_Recall@10','protein_macro_median_Recall@20']].to_string(index=False))
