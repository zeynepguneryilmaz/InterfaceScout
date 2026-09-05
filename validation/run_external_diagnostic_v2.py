from __future__ import annotations
import json, urllib.request, sys
from pathlib import Path
import numpy as np, pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import backend.main as ism
OUT=Path('validation/external_results_v2'); OUT.mkdir(parents=True,exist_ok=True)

BENCH=[
 {'id':'LYZ_SWCNT','pdb':'1LYZ','chains':['A'],'chem':'pi_carbon','pH':7.4,'truth':[61,62,63,68,69,70,73,103,106,107,111,112,113,116]},
 {'id':'CHT_CNT_hydrophobic','pdb':'4CHA','chains':['A','B','C'],'chem':'hydrophobic','pH':7.0,'truth':[3,5,6,7,8,10,77,114,115,116]},
 {'id':'CHT_CNT_pi','pdb':'4CHA','chains':['A','B','C'],'chem':'pi_carbon','pH':7.0,'truth':[3,5,6,7,8,10,77,114,115,116]},
 {'id':'FN_NH3','pdb':'1FNF','chains':['A'],'chem':'cationic','pH':7.0,'truth':[1312,1509]},
 {'id':'FN_COO','pdb':'1FNF','chains':['A'],'chem':'anionic','pH':7.0,'truth':[1469]},
 {'id':'FN10_HAP_Ca','pdb':'1TTF','chains':['A'],'chem':'hydroxyapatite','pH':7.0,'truth':[6,7,9,91,92,93]},
 {'id':'FN10_HAP_phosphate','pdb':'1TTF','chains':['A'],'chem':'phosphate','pH':7.0,'truth':[6,7,9,91,92,93]},
]
CHAINSETS={b['pdb']:set(b['chains']) for b in BENCH}
class ProteinSelect(Select):
 def __init__(self,pid): self.pid=pid
 def accept_chain(self,c): return 1 if str(c.id) in CHAINSETS[self.pid] else 0
 def accept_residue(self,r): return 1 if is_aa(r,standard=True) else 0

def getpdb(pid):
 raw=OUT/f'{pid}_raw.pdb'; clean=OUT/f'{pid}.pdb'
 if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',raw)
 s=PDBParser(QUIET=True).get_structure(pid,str(raw)); io=PDBIO();io.set_structure(s);io.save(str(clean),ProteinSelect(pid)); return clean

def auc(y,s):
 p=s[y==1];n=s[y==0]
 return float(np.mean([(a>b)+.5*(a==b) for a in p for b in n])) if len(p) and len(n) else np.nan
def ap(y,s):
 if not y.sum(): return np.nan
 o=np.argsort(-s,kind='stable'); q=y[o]; prec=np.cumsum(q)/(np.arange(len(q))+1); return float((prec*q).sum()/q.sum())
def recall(y,s,k):
 if not y.sum():return np.nan
 o=np.argsort(-s,kind='stable')[:min(k,len(s))];return float(y[o].sum()/y.sum())
def spatial(truth,top,coords,R):
 vals=[]
 for t in truth:
  if t not in coords:continue
  ds=[np.linalg.norm(coords[t]-coords[p]) for p in top if p in coords]
  vals.append(bool(ds and min(ds)<=R))
 return float(np.mean(vals)) if vals else np.nan
def perm(y,s,exp,n=1000,seed=20260905):
 rng=np.random.default_rng(seed);obs=float(s[y==1].mean());q=pd.qcut(exp,4,labels=False,duplicates='drop')
 need={int(b):int(((q==b)&(y==1)).sum()) for b in np.unique(q)}; null=np.empty(n)
 pools={b:np.where(q==b)[0] for b in need}
 for z in range(n):
  idx=[]
  for b,c in need.items():
   if c:idx.extend(rng.choice(pools[b],c,replace=False))
  null[z]=s[np.array(idx,dtype=int)].mean()
 return float((1+(null>=obs).sum())/(n+1)),obs,float(null.mean())

cache={}; rows=[]; residue=[]
for b in BENCH:
 ck=(b['pdb'],b['pH'])
 if ck not in cache:
  _,allres,_,_=ism.build_surface_residues(getpdb(b['pdb']),b['pH']);surf=[r for r in allres if r['surface_exposed']];cache[ck]=(surf,ism.build_distances(surf))
 surf,D=cache[ck]; mp=ism.chemistry_map(surf,D,b['chem'],b['pH'])
 sm={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']}
 keys=[r['key'] for r in surf];scores=np.array([sm.get(k,0) for k in keys]);exp=np.array([r['scrsa'] for r in surf]);coords={r['key']:np.array([r['x'],r['y'],r['z']]) for r in surf}
 truthkeys=[r['key'] for r in surf if int(r['res_seq']) in b['truth']]; y=np.array([int(k in truthkeys) for k in keys])
 order=np.argsort(-scores,kind='stable');top=[keys[i] for i in order];pval,obs,null=perm(y,scores,exp) if y.sum() else (np.nan,np.nan,np.nan)
 ranks={k:top.index(k)+1 for k in truthkeys}
 rows.append({**b,'n_surface':len(keys),'n_truth_mapped':int(y.sum()),'AUROC':auc(y,scores),'AP':ap(y,scores),'Recall@5':recall(y,scores,5),'Recall@10':recall(y,scores,10),'Recall@20':recall(y,scores,20),'Spatial@5A_top10':spatial(truthkeys,top[:10],coords,5),'Spatial@8A_top10':spatial(truthkeys,top[:10],coords,8),'Spatial@10A_top10':spatial(truthkeys,top[:10],coords,10),'permutation_p':pval,'truth_mean_M':obs,'null_mean_M':null,'best_truth_rank':min(ranks.values()) if ranks else None,'median_truth_rank':float(np.median(list(ranks.values()))) if ranks else None,'truth_ranks':json.dumps(ranks)})
 for r,k,yy,ss in zip(surf,keys,y,scores):residue.append([b['id'],k,r['res_name'],r['res_seq'],r['scrsa'],yy,ss])

df=pd.DataFrame(rows);df.to_csv(OUT/'external_diagnostic_summary.csv',index=False);pd.DataFrame(residue,columns=['benchmark','key','res_name','res_seq','scRSA','truth','M']).to_csv(OUT/'residue_scores.csv',index=False);(OUT/'external_diagnostic_summary.json').write_text(df.to_json(orient='records',indent=2))
print(df[['id','n_truth_mapped','AUROC','AP','Recall@10','Spatial@8A_top10','permutation_p','best_truth_rank','median_truth_rank']].to_string(index=False))
