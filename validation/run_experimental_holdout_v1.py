from pathlib import Path
import sys, urllib.request, json
import numpy as np, pandas as pd
from Bio.PDB import PDBParser,PDBIO,Select
from Bio.PDB.Polypeptide import is_aa
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));import backend.main as ism
ROOT=Path('validation/prospective_holdout');OUT=ROOT/'results_experimental_v1';OUT.mkdir(parents=True,exist_ok=True);PDBDIR=OUT/'pdb_cache';PDBDIR.mkdir(exist_ok=True)
CASES=[
 {'id':'UBQ_AUNP_NMR','protein':'human ubiquitin','pdb':'1UBQ','chains':{'A'},'chem':'anionic','ph':7.7,'truth':[2,3,15,16,17,18],'gt':'residue'},
 {'id':'B2M_CIT_AUNP_WT','protein':'beta-2-microglobulin','pdb':'1JNJ','chains':{'A'},'chem':'anionic','ph':7.7,'truth':[2,3,26,28,29,30,33,55,56,58,59],'gt':'residue'},
 {'id':'HSA_PAA_FE3O4_XLMS','protein':'human serum albumin','pdb':'2VUF','chains':{'A'},'chem':'anionic','ph':7.0,'regions':[(373,389),(403,410),(414,428)],'gt':'region'},
]
class S(Select):
 def __init__(self,c,mid):self.c=c;self.mid=mid
 def accept_model(self,m):return int(m.id==self.mid)
 def accept_chain(self,ch):return int(ch.id in self.c['chains'])
 def accept_residue(self,r):return int(is_aa(r,standard=True))
def getpdb(c):
 raw=PDBDIR/f"{c['pdb']}_raw.pdb";clean=PDBDIR/f"{c['pdb']}_curated.pdb"
 if not raw.exists():urllib.request.urlretrieve(f"https://files.rcsb.org/download/{c['pdb']}.pdb",raw)
 if not clean.exists():
  s=PDBParser(QUIET=True).get_structure(c['pdb'],str(raw));m=next(s.get_models());io=PDBIO();io.set_structure(s);io.save(str(clean),S(c,m.id))
 return clean
def auc(y,s):
 p=s[y==1];n=s[y==0]
 return float(np.mean([(a>b)+.5*(a==b) for a in p for b in n])) if len(p) and len(n) else np.nan
def ap(y,s):
 if y.sum()==0:return np.nan
 o=np.argsort(-s,kind='stable');yy=y[o];pr=np.cumsum(yy)/(np.arange(len(yy))+1);return float((pr*yy).sum()/yy.sum())
def rec(y,s,k):return float(y[np.argsort(-s,kind='stable')[:min(k,len(s))]].sum()/y.sum()) if y.sum() else np.nan
def spatial(truthkeys,topkeys,coords,R):
 if not truthkeys:return np.nan
 hit=0
 for t in truthkeys:
  if any(np.linalg.norm(coords[t]-coords[k])<=R for k in topkeys):hit+=1
 return hit/len(truthkeys)
def nearest(truthkeys,topkeys,coords):
 if not truthkeys or not topkeys:return np.nan
 ds=[min(np.linalg.norm(coords[t]-coords[k]) for k in topkeys) for t in truthkeys]
 return float(np.median(ds))
rows=[];resrows=[];old=ism.SASA_POINTS;ism.SASA_POINTS=200
try:
 for c in CASES:
  _,res,atoms,_=ism.build_surface_residues(getpdb(c),c['ph']);surf=[r for r in res if r['surface_exposed']];D=ism.build_distances(surf);mp=ism.chemistry_map(surf,D,c['chem'],c['ph']);core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']};keys=[r['key'] for r in surf];seq=np.array([int(r['res_seq']) for r in surf]);score=np.array([core.get(k,0.) for k in keys],float);coords={r['key']:np.array([r['x'],r['y'],r['z']],float) for r in surf};order=np.argsort(-score,kind='stable');topkeys=[keys[i] for i in order[:10]]
  if c['gt']=='residue': truth=set(c['truth']);y=np.isin(seq,list(truth)).astype(int);truthkeys=[keys[i] for i,v in enumerate(y) if v]
  else:
   mask=np.zeros(len(seq),dtype=bool)
   for a,b in c['regions']:mask|=((seq>=a)&(seq<=b))
   y=mask.astype(int);truthkeys=[keys[i] for i,v in enumerate(y) if v]
  rows.append({'candidate_id':c['id'],'protein':c['protein'],'pdb':c['pdb'],'pH':c['ph'],'ground_truth_type':c['gt'],'n_surface':len(surf),'n_truth_requested':len(c.get('truth',[])) if c['gt']=='residue' else sum(b-a+1 for a,b in c['regions']),'n_truth_mapped':int(y.sum()),'AUROC':auc(y,score) if c['gt']=='residue' else np.nan,'AP':ap(y,score) if c['gt']=='residue' else np.nan,'Recall@5':rec(y,score,5) if c['gt']=='residue' else np.nan,'Recall@10':rec(y,score,10) if c['gt']=='residue' else np.nan,'Recall@20':rec(y,score,20) if c['gt']=='residue' else np.nan,'Spatial@5A_top10':spatial(truthkeys,topkeys,coords,5),'Spatial@8A_top10':spatial(truthkeys,topkeys,coords,8),'Spatial@10A_top10':spatial(truthkeys,topkeys,coords,10),'median_nearest_top10_A':nearest(truthkeys,topkeys,coords),'top10':';'.join(topkeys)})
  for i,r in enumerate(surf):resrows.append({'candidate_id':c['id'],'key':keys[i],'res_name':r['res_name'],'res_seq':int(r['res_seq']),'chain':r['chain'],'scrsa':float(r['scrsa']),'core_is_raw':float(score[i]),'rank':int(np.where(order==i)[0][0])+1,'ground_truth':int(y[i])})
finally:ism.SASA_POINTS=old
R=pd.DataFrame(rows);R.to_csv(OUT/'experimental_metrics.csv',index=False);pd.DataFrame(resrows).to_csv(OUT/'experimental_residue_scores.csv',index=False)
meta={'model':'IS-v1-Core-primary','frozen_before_scoring':True,'SASA_POINTS':200,'chemistry':'anionic for all three cases','note':'No model, cutoff, ontology, or score rule was changed after experimental GT freeze.'};(OUT/'METADATA.json').write_text(json.dumps(meta,indent=2))
print(R.to_string(index=False));print(json.dumps(meta,indent=2))
