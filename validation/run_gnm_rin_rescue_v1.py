from pathlib import Path
from collections import deque
import sys, urllib.request, json
import numpy as np, pandas as pd
from Bio.PDB import PDBParser,PDBIO,Select
from Bio.PDB.Polypeptide import is_aa
sys.path.insert(0,str(Path(__file__).resolve().parents[1])); import backend.main as ism

OUT=Path('validation/prospective_holdout/results_gnm_rin_rescue_v1'); OUT.mkdir(parents=True,exist_ok=True)
PDBDIR=OUT/'pdb_cache'; PDBDIR.mkdir(exist_ok=True)
SYSTEMS=[
 {'id':'GB3_CIT_AUNP_MUT','protein':'GB3','pdb':'2OED','chain':'A','pH':6.4,'chem':'anionic','truth':set([4,13,50]),'gt_type':'residue'},
 {'id':'B2M_CIT_AUNP_WT','protein':'beta-2-microglobulin','pdb':'1JNJ','chain':'A','pH':7.7,'chem':'anionic','truth':set([2,3,26,28,29,30,33,55,56,58,59]),'gt_type':'residue'},
 {'id':'HSA_PAA_FE3O4_XLMS','protein':'HSA','pdb':'2VUF','chain':'A','pH':7.0,'chem':'anionic','truth':set(list(range(373,390))+list(range(403,411))+list(range(414,429))),'gt_type':'peptide'},
]
VARIANTS=['Core','GNM_rescue','RIN_rescue','GNM_RIN_rescue']
GNM_CUTOFF=10.0; RIN_CUTOFF=8.0; FLOOR=0.05
class Sel(Select):
 def __init__(self,ch): self.ch=ch
 def accept_model(self,m): return int(m.id==0)
 def accept_chain(self,c): return int(c.id==self.ch)
 def accept_residue(self,r): return int(is_aa(r,standard=True))

def ensure_pdb(pid,ch):
 raw=PDBDIR/f'{pid}_raw.pdb'; clean=PDBDIR/f'{pid}_{ch}.pdb'
 if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',raw)
 if not clean.exists():
  s=PDBParser(QUIET=True).get_structure(pid,str(raw)); io=PDBIO(); io.set_structure(s); io.save(str(clean),Sel(ch))
 return clean

def structural(pdb):
 s=PDBParser(QUIET=True).get_structure('x',str(pdb)); model=next(s.get_models()); keys=[]; X=[]
 for ch in model:
  for r in ch:
   if is_aa(r,standard=True) and 'CA' in r:
    keys.append(f"{ch.id}:{int(r.id[1])}:{str(r.id[2]).strip()}"); X.append(np.asarray(r['CA'].coord,float))
 X=np.vstack(X); n=len(X); D=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2)
 K=np.zeros((n,n),float)
 for i in range(n):
  for j in np.where((D[i]<=GNM_CUTOFF)&(D[i]>0))[0]:
   if j>i: K[i,j]=K[j,i]=-1.; K[i,i]+=1.; K[j,j]+=1.
 vals,vecs=np.linalg.eigh(K); pos=vals>1e-8; cov=(vecs[:,pos]*(1/vals[pos]))@vecs[:,pos].T if pos.any() else np.zeros((n,n)); msf=np.diag(cov)
 mob=(msf-msf.min())/(msf.max()-msf.min()) if msf.max()>msf.min() else np.zeros(n)
 A=((D<=RIN_CUTOFF)&(D>0)).astype(int); degree=A.sum(1).astype(float)
 return {k:{'gnm':float(mob[i]),'degree':float(degree[i])} for i,k in enumerate(keys)}

def auc(y,s):
 p=s[y==1]; n=s[y==0]
 if len(p)==0 or len(n)==0:return np.nan
 return float(np.mean([(a>b)+0.5*(a==b) for a in p for b in n]))
def ap(y,s):
 if y.sum()==0:return np.nan
 o=np.argsort(-s,kind='stable'); yy=y[o]; pr=np.cumsum(yy)/(np.arange(len(yy))+1); return float((pr*yy).sum()/yy.sum())
def rec(y,s,k):
 if y.sum()==0:return np.nan
 return float(y[np.argsort(-s,kind='stable')[:min(k,len(s))]].sum()/y.sum())

def score_variant(allres,struct,variant,chem,pH):
 gvals=np.array([struct[r['key']]['gnm'] for r in allres]); dvals=np.array([struct[r['key']]['degree'] for r in allres])
 gq=float(np.quantile(gvals,0.75)); dq=float(np.quantile(dvals,0.25))
 selected=[]
 for r in allres:
  exposed=bool(r['surface_exposed']); g=struct[r['key']]['gnm']; d=struct[r['key']]['degree']
  rescue=(variant=='GNM_rescue' and g>=gq) or (variant=='RIN_rescue' and d<=dq) or (variant=='GNM_RIN_rescue' and g>=gq and d<=dq)
  if variant=='Core': rescue=False
  if exposed or rescue:
   rr=dict(r); rr['scrsa']=max(float(rr['scrsa']),FLOOR) if rescue and not exposed else float(rr['scrsa']); rr['rescued']=bool(rescue and not exposed); selected.append(rr)
 D=ism.build_distances(selected); mp=ism.chemistry_map(selected,D,chem,pH); core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']}
 keys=[r['key'] for r in selected]; seq=np.array([int(r['res_seq']) for r in selected]); score=np.array([core.get(k,0.) for k in keys]); rescued=np.array([bool(r.get('rescued',False)) for r in selected])
 return selected,D,keys,seq,score,rescued,gq,dq

metrics=[]; diagnostics=[]; residues=[]
old=ism.SASA_POINTS; ism.SASA_POINTS=200
try:
 for S in SYSTEMS:
  pdb=ensure_pdb(S['pdb'],S['chain']); _,allres,_,_=ism.build_surface_residues(pdb,S['pH']); st=structural(pdb)
  for r in allres:
   x=st.get(r['key'],{}); diagnostics.append({'system':S['id'],'key':r['key'],'res_seq':r['res_seq'],'res_name':r['res_name'],'scrsa':r['scrsa'],'surface_exposed':r['surface_exposed'],'gnm_mobility':x.get('gnm'),'rin_degree':x.get('degree'),'ground_truth':int(int(r['res_seq']) in S['truth'])})
  for V in VARIANTS:
   selected,D,keys,seq,score,rescued,gq,dq=score_variant(allres,st,V,S['chem'],S['pH']); y=np.isin(seq,list(S['truth'])).astype(int); top=np.argsort(-score,kind='stable')[:10]
   present_truth=sum(int(r['res_seq']) in S['truth'] for r in allres); mapped_truth=int(y.sum())
   nearest=[]
   for ti in np.where(y==1)[0]: nearest.append(float(np.min(D[ti,top])) if len(top) else np.nan)
   row={'system':S['id'],'protein':S['protein'],'variant':V,'n_candidates':len(selected),'truth_present_structure':present_truth,'truth_mapped_candidate':mapped_truth,'truth_rescued':int(sum(rescued & (y==1))),'gnm_q75':gq,'rin_degree_q25':dq,'AUROC':auc(y,score),'AP':ap(y,score),'Recall@5':rec(y,score,5),'Recall@10':rec(y,score,10),'Recall@20':rec(y,score,20),'Spatial@5A_top10':float(np.mean(np.array(nearest)<=5)) if nearest else np.nan,'Spatial@8A_top10':float(np.mean(np.array(nearest)<=8)) if nearest else np.nan,'Spatial@10A_top10':float(np.mean(np.array(nearest)<=10)) if nearest else np.nan,'median_nearest_top10_A':float(np.nanmedian(nearest)) if nearest else np.nan,'top10':';'.join(keys[i] for i in top)}; metrics.append(row)
   order=np.argsort(-score,kind='stable'); ranks=np.empty(len(order),int); ranks[order]=np.arange(1,len(order)+1)
   for i,r in enumerate(selected): residues.append({'system':S['id'],'variant':V,'key':r['key'],'res_seq':r['res_seq'],'res_name':r['res_name'],'scrsa':r['scrsa'],'rescued':bool(rescued[i]),'gnm_mobility':st[r['key']]['gnm'],'rin_degree':st[r['key']]['degree'],'score':float(score[i]),'rank':int(ranks[i]),'ground_truth':int(int(r['res_seq']) in S['truth'])})
finally: ism.SASA_POINTS=old
pd.DataFrame(metrics).to_csv(OUT/'metrics.csv',index=False); pd.DataFrame(diagnostics).to_csv(OUT/'all_residue_diagnostics.csv',index=False); pd.DataFrame(residues).to_csv(OUT/'residue_scores.csv',index=False)
(OUT/'METADATA.json').write_text(json.dumps({'protocol':'GNM_RIN_RESCUE_PROTOCOL_V1','model_status':'exploratory; frozen IS-v1-Core-primary unchanged','SASA_POINTS':200,'GNM_cutoff_A':GNM_CUTOFF,'RIN_cutoff_A':RIN_CUTOFF,'rescue_floor_scrsa':FLOOR,'rules':'top quartile GNM; bottom quartile RIN degree; conjunction for combined rescue'},indent=2))
print(pd.DataFrame(metrics).to_string(index=False))
