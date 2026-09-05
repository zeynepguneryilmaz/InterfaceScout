from __future__ import annotations
import json, math, urllib.request
from pathlib import Path
import numpy as np, pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import backend.main as ism

OUT=Path('validation/external_results'); OUT.mkdir(parents=True,exist_ok=True)
PH=7.4

# Pre-freeze diagnostic only: all systems below were inspected in literature before final model freeze.
BENCHMARKS=[
 {'id':'LYZ_SWCNT','pdb':'1LYZ','surface':'SWCNT','chem':'pi_carbon','pH':7.4,
  'truth':[61,62,63,68,69,70,73,103,106,107,111,112,113,116],
  'truth_type':'explicit contact residues <=0.5 nm','source':'Scientific Reports 2025, doi:10.1038/s41598-025-96435-3'},
 {'id':'CHT_CNT_hydrophobic','pdb':'4CHA','surface':'CNT aqueous','chem':'hydrophobic','pH':7.0,
  'truth':[3,5,6,7,8,10,77,114,115,116],
  'truth_type':'significant favorable per-residue binding-energy contributors','source':'Scientific Reports 2015, doi:10.1038/srep09297'},
 {'id':'CHT_CNT_pi','pdb':'4CHA','surface':'CNT aqueous','chem':'pi_carbon','pH':7.0,
  'truth':[3,5,6,7,8,10,77,114,115,116],
  'truth_type':'same literature contact-energy set; pi-carbon map diagnostic','source':'Scientific Reports 2015, doi:10.1038/srep09297'},
 {'id':'FN_NH3','pdb':'1FNF','surface':'NH3+ SAM','chem':'cationic','pH':7.0,
  'truth':[1312,1509], 'truth_type':'final-stage anchoring residues','source':'Materials 2018, doi:10.3390/ma11122570'},
 {'id':'FN_COO','pdb':'1FNF','surface':'COO- SAM','chem':'anionic','pH':7.0,
  'truth':[1469], 'truth_type':'final-stage anchoring residue','source':'Materials 2018, doi:10.3390/ma11122570'},
 {'id':'FN10_HAP_Ca','pdb':'1TTF','surface':'HAp(001)','chem':'hydroxyapatite','pH':7.0,
  'truth':[6,7,9,91,92,93], 'truth_type':'residues within 0.6 nm of HAp','source':'RSC Adv 2014, doi:10.1039/C3RA47381C'},
 {'id':'FN10_HAP_phosphate','pdb':'1TTF','surface':'HAp(001)','chem':'phosphate','pH':7.0,
  'truth':[6,7,9,91,92,93], 'truth_type':'same HAp contact set; phosphate-site map diagnostic','source':'RSC Adv 2014, doi:10.1039/C3RA47381C'},
]

class Std(Select):
 def accept_residue(self,r): return 1 if is_aa(r,standard=True) else 0

def getpdb(pid):
 raw=OUT/f'{pid}_raw.pdb'; clean=OUT/f'{pid}.pdb'
 if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',raw)
 s=PDBParser(QUIET=True).get_structure(pid,str(raw)); io=PDBIO();io.set_structure(s);io.save(str(clean),Std());return clean

def auc(labels,scores):
 p=scores[labels==1];n=scores[labels==0]
 if not len(p) or not len(n): return np.nan
 return float(np.mean([(a>b)+.5*(a==b) for a in p for b in n]))
def ap(labels,scores):
 npos=int(labels.sum())
 if not npos:return np.nan
 o=np.argsort(-scores,kind='stable'); y=labels[o]
 pr=np.cumsum(y)/(np.arange(len(y))+1)
 return float(np.sum(pr*y)/npos)
def recallk(labels,scores,k):
 n=int(labels.sum());
 if not n:return np.nan
 o=np.argsort(-scores,kind='stable')[:min(k,len(scores))]
 return float(labels[o].sum()/n)
def perm_p(labels,scores,exposure,nperm=10000,seed=20260905):
 # exposure-matched: positives sampled from same scRSA quartile counts
 rng=np.random.default_rng(seed); obs=float(scores[labels==1].mean()); q=pd.qcut(exposure,4,labels=False,duplicates='drop')
 need={int(x):int(np.sum((q==x)&(labels==1))) for x in np.unique(q)}
 vals=[]
 for _ in range(nperm):
  idx=[]
  for b,c in need.items():
   pool=np.where(q==b)[0]
   if c: idx.extend(rng.choice(pool,size=c,replace=False).tolist())
  vals.append(float(scores[idx].mean()) if idx else 0.)
 vals=np.asarray(vals); return float((1+np.sum(vals>=obs))/(nperm+1)),obs,float(np.mean(vals))

def spatial_recovery(truth_keys,top_keys,coords,R):
 if not truth_keys:return np.nan
 hit=0
 for t in truth_keys:
  if t not in coords: continue
  d=[np.linalg.norm(coords[t]-coords[p]) for p in top_keys if p in coords]
  if d and min(d)<=R: hit+=1
 return hit/len(truth_keys)

rows=[]; residue_rows=[]
for b in BENCHMARKS:
 old=ism.SASA_POINTS; ism.SASA_POINTS=200
 try:
  _,allres,_,_=ism.build_surface_residues(getpdb(b['pdb']),b['pH'])
 finally: ism.SASA_POINTS=old
 surface=[r for r in allres if r['surface_exposed']]
 D=ism.build_distances(surface); mp=ism.chemistry_map(surface,D,b['chem'],b['pH'])
 # M score for every exposed residue, including chemistry-incompatible centers
 score_by={x['center_key']:x['multiscale_persistence']/100 for x in mp['patch_centers']}
 coords={r['key']:np.array([r['x'],r['y'],r['z']],float) for r in surface}
 # truth mapping by residue number; include any chains with that unique resseq
 truth_keys=[]
 for r in surface:
  if int(r['res_seq']) in b['truth']:truth_keys.append(r['key'])
 keys=[r['key'] for r in surface]; labels=np.array([1 if k in truth_keys else 0 for k in keys],int)
 scores=np.array([score_by.get(k,0.) for k in keys],float)
 expos=np.array([r['scrsa'] for r in surface],float)
 top_order=[keys[i] for i in np.argsort(-scores,kind='stable')]
 pp,obs,nullmean=perm_p(labels,scores,expos) if labels.sum() else (np.nan,np.nan,np.nan)
 anchor_ranks={}
 for t in truth_keys:
  anchor_ranks[t]=top_order.index(t)+1 if t in top_order else None
 row={
  **{k:b[k] for k in ['id','pdb','surface','chem','pH','truth_type','source']},
  'n_surface':len(keys),'n_truth_mapped':int(labels.sum()),'AUROC':auc(labels,scores),'AP':ap(labels,scores),
  'Recall@5':recallk(labels,scores,5),'Recall@10':recallk(labels,scores,10),'Recall@20':recallk(labels,scores,20),
  'Spatial@5A_top10':spatial_recovery(truth_keys,top_order[:10],coords,5.0),
  'Spatial@8A_top10':spatial_recovery(truth_keys,top_order[:10],coords,8.0),
  'Spatial@10A_top10':spatial_recovery(truth_keys,top_order[:10],coords,10.0),
  'permutation_p':pp,'observed_truth_mean_M':obs,'null_mean_M':nullmean,
  'best_truth_rank':min(anchor_ranks.values()) if anchor_ranks else None,
  'median_truth_rank':float(np.median(list(anchor_ranks.values()))) if anchor_ranks else None,
  'truth_ranks_json':json.dumps(anchor_ranks,sort_keys=True),
 }
 rows.append(row)
 for r,k,l,s in zip(surface,keys,labels,scores):residue_rows.append([b['id'],k,r['res_name'],r['res_seq'],r['scrsa'],l,s])

pd.DataFrame(rows).to_csv(OUT/'external_diagnostic_summary.csv',index=False)
pd.DataFrame(residue_rows,columns=['benchmark','key','res_name','res_seq','scRSA','truth','M']).to_csv(OUT/'external_diagnostic_residue_scores.csv',index=False)
(Path(OUT/'external_diagnostic_summary.json')).write_text(json.dumps(rows,indent=2))
print(pd.DataFrame(rows)[['id','n_truth_mapped','AUROC','AP','Recall@10','Spatial@8A_top10','permutation_p','best_truth_rank','median_truth_rank']].to_string(index=False))
