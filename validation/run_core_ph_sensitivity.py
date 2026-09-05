from pathlib import Path
import sys, urllib.request
import numpy as np, pandas as pd
from scipy.stats import rankdata
from Bio.PDB import PDBParser,PDBIO,Select
from Bio.PDB.Polypeptide import is_aa
sys.path.insert(0,str(Path(__file__).resolve().parents[1]));import backend.main as ism
ROOT=Path('validation/packages/all15');OUT=ROOT/'results_core_ph_sensitivity';OUT.mkdir(parents=True,exist_ok=True);PDBDIR=OUT/'pdb_cache';PDBDIR.mkdir(exist_ok=True)
CASES=[
 {'id':'CHT_CNT_HYDRO','protein':'Alpha-chymotrypsin','pdb':'4CHA','chains':{'A','B','C'},'chem':'hydrophobic','truth':[3,5,6,7,8,10,77,114,115,116]},
 {'id':'MB_CIT_AUNP','protein':'Myoglobin','pdb':'1MBN','chains':{'A'},'chem':'anionic','truth':[43,45,96,97,98,99,146,147,148,149,150,151,152,153]},
 {'id':'WNV_GRAPHENE','protein':'WNV E protein domain III','pdb':'2HG0','chains':{'A'},'range':(8,109),'chem':'pi_carbon','truth':[55,56,57,58,59,60,61,106,107,108,109]},
]
class S(Select):
 def __init__(self,c,mid):self.c=c;self.mid=mid
 def accept_model(self,m):return int(m.id==self.mid)
 def accept_chain(self,ch):return int(ch.id in self.c['chains'])
 def accept_residue(self,r):
  if not is_aa(r,standard=True):return 0
  if 'range' in self.c:return int(self.c['range'][0]<=int(r.id[1])<=self.c['range'][1])
  return 1
def pdb(c):
 raw=PDBDIR/f"{c['pdb']}_raw.pdb";clean=PDBDIR/f"{c['pdb']}_curated.pdb"
 if not raw.exists():urllib.request.urlretrieve(f"https://files.rcsb.org/download/{c['pdb']}.pdb",raw)
 if not clean.exists():
  s=PDBParser(QUIET=True).get_structure(c['pdb'],str(raw));m=next(s.get_models());io=PDBIO();io.set_structure(s);io.save(str(clean),S(c,m.id))
 return clean
def auc(y,s):
 p=s[y==1];n=s[y==0];return float(np.mean([(a>b)+.5*(a==b) for a in p for b in n]))
def ap(y,s):
 o=np.argsort(-s,kind='stable');yy=y[o];pr=np.cumsum(yy)/(np.arange(len(yy))+1);return float((pr*yy).sum()/yy.sum())
def rec(y,s,k):return float(y[np.argsort(-s,kind='stable')[:min(k,len(s))]].sum()/y.sum())
rows=[]
old=ism.SASA_POINTS;ism.SASA_POINTS=200
try:
 for c in CASES:
  for ph in [6.5,7.0,7.5]:
   _,res,atoms,_=ism.build_surface_residues(pdb(c),ph);surf=[r for r in res if r['surface_exposed']];D=ism.build_distances(surf);mp=ism.chemistry_map(surf,D,c['chem'],ph);core={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']};keys=[r['key'] for r in surf];seq=np.array([int(r['res_seq']) for r in surf]);score=np.array([core.get(k,0.) for k in keys]);y=np.isin(seq,c['truth']).astype(int)
   top=np.argsort(-score,kind='stable')[:10]
   rows.append({'condition_id':c['id'],'protein':c['protein'],'pH':ph,'n_surface':len(surf),'n_truth_mapped':int(y.sum()),'AUROC':auc(y,score),'AP':ap(y,score),'Recall@5':rec(y,score,5),'Recall@10':rec(y,score,10),'Recall@20':rec(y,score,20),'top10':';'.join(keys[i] for i in top)})
finally:ism.SASA_POINTS=old
R=pd.DataFrame(rows);R.to_csv(OUT/'core_ph_sensitivity.csv',index=False)
base=R[R.pH==7.0].set_index('condition_id');d=[]
for _,r in R[R.pH!=7.0].iterrows():
 b=base.loc[r.condition_id];d.append({'condition_id':r.condition_id,'pH':r.pH,'delta_AUROC_vs_pH7':r.AUROC-b.AUROC,'delta_AP_vs_pH7':r.AP-b.AP,'delta_R10_vs_pH7':r['Recall@10']-b['Recall@10'],'top10_identical_to_pH7':r.top10==b.top10})
pd.DataFrame(d).to_csv(OUT/'core_ph_sensitivity_deltas.csv',index=False)
print(R.to_string(index=False));print(pd.DataFrame(d).to_string(index=False))
