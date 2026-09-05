from __future__ import annotations

import json, math, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import backend.main as ism

OUT=Path('validation/external_results'); OUT.mkdir(parents=True, exist_ok=True)
PH_DEFAULT=7.0

class FirstChain(Select):
    def __init__(self,cid): self.cid=cid
    def accept_chain(self,c): return 1 if str(c.id)==self.cid else 0
    def accept_residue(self,r): return 1 if is_aa(r,standard=True) else 0

def download(pid):
    raw=OUT/f'{pid}_raw.pdb'; clean=OUT/f'{pid}.pdb'
    if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',raw)
    if not clean.exists():
        s=PDBParser(QUIET=True).get_structure(pid,str(raw)); model=next(s.get_models())
        cid=next(c.id for c in model if any(is_aa(r,standard=True) for r in c))
        io=PDBIO(); io.set_structure(s); io.save(str(clean),FirstChain(str(cid)))
    return clean

def residues(pid,ph):
    _,res,_,_=ism.build_surface_residues(download(pid),ph)
    return res

def dist(surface):
    if not surface:return np.empty((0,0))
    x=np.asarray([[r['x'],r['y'],r['z']] for r in surface],float)
    d=x[:,None,:]-x[None,:,:]
    return np.sqrt((d*d).sum(2))

def score_map(res,chem,ph):
    surf=[r for r in res if r['scrsa_raw']>=ism.SC_RSA_THRESHOLD]
    d=dist(surf); defs=ism.CHEMISTRIES[chem]['favorable']; n=len(surf)
    local=np.zeros(n)
    for i,r in enumerate(surf):
        if r['res_name'] in defs:
            st=defs[r['res_name']][2]
            local[i]=r['scrsa']*ism.state_availability(r['res_name'],st,ph)
    norms=[]
    for R in ism.PATCH_RADII_A:
        a=(d<=R).astype(float)@local if n else np.zeros(0)
        m=a.max() if len(a) else 0
        norms.append(a/m if m>0 else np.zeros_like(a))
    M=100*np.minimum(norms[0],norms[1]) if n else np.zeros(0)
    return {r['res_seq']:float(M[i]) for i,r in enumerate(surf)}, {r['res_seq']:np.array([r['x'],r['y'],r['z']],float) for r in surf}, [r['res_seq'] for r in res]

def auroc(labels,scores):
    y=np.asarray(labels); s=np.asarray(scores)
    pos=s[y==1]; neg=s[y==0]
    if len(pos)==0 or len(neg)==0:return np.nan
    return float(np.mean([(p>n)+0.5*(p==n) for p in pos for n in neg]))

def average_precision(labels,scores):
    order=np.argsort(-np.asarray(scores),kind='mergesort'); y=np.asarray(labels)[order]
    if y.sum()==0:return np.nan
    prec=np.cumsum(y)/(np.arange(len(y))+1)
    return float(np.sum(prec*y)/y.sum())

def recall_at(labels,scores,k):
    y=np.asarray(labels); order=np.argsort(-np.asarray(scores),kind='mergesort')
    return float(y[order[:k]].sum()/y.sum()) if y.sum() else np.nan

def spatial_recovery(gt,score_by_res,coords,R,k=10):
    pred=[r for r,_ in sorted(score_by_res.items(), key=lambda z:(-z[1],z[0]))[:k] if zval((r,_))>0]
    if not gt:return np.nan
    hit=0
    for g in gt:
        if g not in coords: continue
        if any(p in coords and np.linalg.norm(coords[g]-coords[p])<=R for p in pred): hit+=1
    denom=sum(g in coords for g in gt)
    return hit/denom if denom else np.nan

def zval(z): return z[1]

def anchor_rank(gt,score_by_res,allres):
    full={r:score_by_res.get(r,0.0) for r in allres}
    order=sorted(full,key=lambda r:(-full[r],r))
    ranks={r:i+1 for i,r in enumerate(order)}
    vals=[ranks[g] for g in gt if g in ranks]
    return min(vals) if vals else np.nan

# Ground truths manually transcribed from cited literature.
BENCH=[
 {"id":"FNIII_NH3","pdb":"1FNF","chem":"cationic","ph":7.0,"gt":[1312,1509],"source":"Liamas et al. 2018; NH3+ SAM anchors"},
 {"id":"FNIII_COO","pdb":"1FNF","chem":"anionic","ph":7.0,"gt":[1469],"source":"Liamas et al. 2018; COO- SAM anchor"},
 {"id":"FNIII_CH3_head","pdb":"1FNF","chem":"hydrophobic","ph":7.0,"gt":[1454,1455,1456,1457,1478,1479,1480,1481,1509],"source":"Liamas et al. 2018; CH3 head-on contacts"},
 {"id":"LYZ_SWCNT","pdb":"1LYZ","chem":"pi_carbon","ph":7.0,"gt":[61,62,63,68,69,70,73,103,106,107,111,112,113,116],"source":"Scientific Reports 2025; native lysozyme/SWCNT <0.5 nm"},
 {"id":"MB_citrate_AuNP","pdb":"1MBN","chem":"anionic","ph":7.0,"gt":[43,45,96,97,98,99,146,147,148,149,150,151,152,153],"source":"Menziani et al. 2019; persistent CG contacts"},
 {"id":"FNIII10_HAP","pdb":"1TTF","chem":"hydroxyapatite","ph":7.4,"gt":[6,7,9,91,92,93],"source":"RSC Adv 2014; residues within 0.6 nm of HAP"},
]

rows=[]; detail=[]
for b in BENCH:
    res=residues(b['pdb'],b['ph']); sm,coords,allres=score_map(res,b['chem'],b['ph'])
    labels=[1 if r in set(b['gt']) else 0 for r in allres]
    scores=[sm.get(r,0.0) for r in allres]
    row={
      'benchmark':b['id'],'pdb':b['pdb'],'chemistry':b['chem'],'pH':b['ph'],'n_residues':len(allres),'n_gt':len(b['gt']),
      'AUROC':auroc(labels,scores),'AP':average_precision(labels,scores),
      'Recall@5':recall_at(labels,scores,5),'Recall@10':recall_at(labels,scores,10),'Recall@20':recall_at(labels,scores,20),
      'Spatial@5A':spatial_recovery(b['gt'],sm,coords,5,10),'Spatial@8A':spatial_recovery(b['gt'],sm,coords,8,10),'Spatial@10A':spatial_recovery(b['gt'],sm,coords,10,10),
      'Best_anchor_rank':anchor_rank(b['gt'],sm,allres),'source_note':b['source']}
    rows.append(row)
    for g in b['gt']:
        detail.append({'benchmark':b['id'],'gt_residue':g,'score':sm.get(g,0.0),'rank':anchor_rank([g],sm,allres),'surface_scored':g in sm})

DF=pd.DataFrame(rows); DF.to_csv(OUT/'external_benchmark_summary.csv',index=False)
pd.DataFrame(detail).to_csv(OUT/'external_benchmark_residue_detail.csv',index=False)
summary={
 'n_benchmarks':len(DF),
 'median_AUROC':float(DF.AUROC.median()),
 'median_AP':float(DF.AP.median()),
 'median_Recall@10':float(DF['Recall@10'].median()),
 'median_Spatial@8A':float(DF['Spatial@8A'].median()),
 'median_Spatial@10A':float(DF['Spatial@10A'].median()),
 'benchmarks':DF.to_dict('records')}
(OUT/'external_benchmark_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
