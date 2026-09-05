from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import backend.main as ism

BASE=Path('validation/packages/all15/results_ready_batch')
PDBDIR=BASE/'pdb_cache'
OUT=Path('validation/packages/all15/results_gnm_exploratory')
OUT.mkdir(parents=True,exist_ok=True)
CUTOFF_A=10.0
SASA_POINTS=120

BENCH=[
{'id':'FN_COO','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'anionic','pH':7.0,'truth':[1469]},
{'id':'FN_NH3','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'cationic','pH':7.0,'truth':[1312,1509]},
{'id':'FN_CH3_HEAD','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1454,1455,1456,1457,1478,1479,1480,1481,1509]},
{'id':'FN_CH3_SIDE','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1275,1276,1355,1376,1378,1454,1455,1456,1457]},
{'id':'FN_CH3_BETA','protein':'Fibronectin III8-10','pdb':'1FNF','chem':'hydrophobic','pH':7.0,'truth':[1431,1446,1457,1458,1459,1461,1463,1464,1466,1475,1476]},
{'id':'CYTC_CH3','protein':'Cytochrome c','pdb':'3NWV','chem':'hydrophobic','pH':7.0,'truth':[1,2,3,4,96,99,100]},
{'id':'CYTC_COOH','protein':'Cytochrome c','pdb':'3NWV','chem':'anionic','pH':7.0,'truth':[1,2,4,61,99,100,103]},
{'id':'LYZ_SWCNT','protein':'Lysozyme','pdb':'1LYZ','chem':'pi_carbon','pH':7.4,'truth':[61,62,63,68,69,70,73,103,106,107,111,112,113,116]},
{'id':'CHT_CNT_HYDRO','protein':'Alpha-chymotrypsin','pdb':'4CHA','chem':'hydrophobic','pH':7.0,'truth':[3,5,6,7,8,10,77,114,115,116]},
{'id':'MB_CIT_AUNP','protein':'Myoglobin','pdb':'1MBN','chem':'anionic','pH':7.0,'truth':[43,45,96,97,98,99,146,147,148,149,150,151,152,153]},
]

def rkey(chain,res):
    return f"{chain.id}:{int(res.id[1])}:{str(res.id[2]).strip()}"

def gnm_msf(pid):
    s=PDBParser(QUIET=True).get_structure(pid,str(PDBDIR/f'{pid}.pdb'))
    model=next(s.get_models()); keys=[]; xyz=[]
    for ch in model:
        for r in ch:
            if 'CA' in r:
                keys.append(rkey(ch,r)); xyz.append(np.asarray(r['CA'].coord,float))
    X=np.vstack(xyz); n=len(X); K=np.zeros((n,n),float)
    for i in range(n-1):
        d=np.linalg.norm(X[i+1:]-X[i],axis=1)
        js=np.where(d<=CUTOFF_A)[0]+i+1
        for j in js:
            K[i,j]=K[j,i]=-1.0; K[i,i]+=1.0; K[j,j]+=1.0
    vals,vecs=np.linalg.eigh(K)
    pos=vals>1e-8
    cov=(vecs[:,pos]*(1.0/vals[pos]))@vecs[:,pos].T
    msf=np.diag(cov)
    lo,hi=float(msf.min()),float(msf.max())
    norm=(msf-lo)/(hi-lo) if hi>lo else np.zeros_like(msf)
    return dict(zip(keys,norm)),dict(zip(keys,msf))

def surface_scores(pid,chem,ph):
    old=ism.SASA_POINTS; ism.SASA_POINTS=SASA_POINTS
    try: _,allres,_,_=ism.build_surface_residues(PDBDIR/f'{pid}.pdb',ph)
    finally: ism.SASA_POINTS=old
    surf=[r for r in allres if r['surface_exposed']]
    D=ism.build_distances(surf); mp=ism.chemistry_map(surf,D,chem,ph)
    sm={x['center_key']:x['multiscale_persistence']/100.0 for x in mp['patch_centers']}
    return {r['key']:{'score':sm.get(r['key'],0.0),'seq':int(r['res_seq'])} for r in surf}

def auc(y,s):
    p=s[y==1]; n=s[y==0]
    if not len(p) or not len(n): return np.nan
    return float(np.mean([(a>b)+0.5*(a==b) for a in p for b in n]))

def ap(y,s):
    if not y.sum(): return np.nan
    o=np.argsort(-s,kind='stable'); yy=y[o]; pr=np.cumsum(yy)/(np.arange(len(yy))+1)
    return float((pr*yy).sum()/yy.sum())

def rec(y,s,k):
    if not y.sum(): return np.nan
    o=np.argsort(-s,kind='stable')[:min(k,len(s))]
    return float(y[o].sum()/y.sum())

def rank01(x):
    o=np.argsort(np.argsort(x,kind='stable'),kind='stable').astype(float)
    return o/(len(x)-1) if len(x)>1 else np.zeros_like(x,dtype=float)

gnm_cache={}; rows=[]; detail=[]
for b in BENCH:
    if b['pdb'] not in gnm_cache: gnm_cache[b['pdb']]=gnm_msf(b['pdb'])
    gnorm,graw=gnm_cache[b['pdb']]
    ref=surface_scores(b['pdb'],b['chem'],b['pH']); keys=list(ref)
    static=np.asarray([ref[k]['score'] for k in keys],float)
    mobility=np.asarray([gnorm.get(k,0.0) for k in keys],float)
    rigidity=1.0-mobility
    # Rank-scale only for combined diagnostics so static and GNM are comparable without fitted coefficients.
    sr=rank01(static); mr=rank01(mobility); rr=1.0-mr
    methods={
      'static_sameSASA':static,
      'gnm_mobility':mobility,
      'gnm_rigidity':rigidity,
      'static_x_gnm_mobility_rank':sr*mr,
      'static_x_gnm_rigidity_rank':sr*rr,
    }
    y=np.asarray([int(ref[k]['seq'] in b['truth']) for k in keys],int)
    for name,s in methods.items():
        rows.append({'condition_id':b['id'],'protein':b['protein'],'pdb':b['pdb'],'method':name,'AUROC':auc(y,s),'AP':ap(y,s),'Recall@5':rec(y,s,5),'Recall@10':rec(y,s,10),'Recall@20':rec(y,s,20)})
    for i,k in enumerate(keys):
        detail.append({'condition_id':b['id'],'protein':b['protein'],'key':k,'res_seq':ref[k]['seq'],'truth':int(y[i]),'static':static[i],'gnm_mobility_norm':mobility[i],'gnm_msf_raw':graw.get(k,np.nan),'static_rank':sr[i]})

res=pd.DataFrame(rows); res.to_csv(OUT/'gnm_method_comparison.csv',index=False)
pd.DataFrame(detail).to_csv(OUT/'gnm_residue_scores.csv',index=False)
summary={'design':{'gnm_cutoff_A':CUTOFF_A,'sasa_points':SASA_POINTS,'combination_policy':'rank products with no fitted coefficients'},'median_metrics_by_method':res.groupby('method')[['AUROC','AP','Recall@5','Recall@10','Recall@20']].median().to_dict(orient='index'),'mean_metrics_by_method':res.groupby('method')[['AUROC','AP','Recall@5','Recall@10','Recall@20']].mean().to_dict(orient='index')}
(OUT/'gnm_summary.json').write_text(json.dumps(summary,indent=2))
print(res.to_string(index=False)); print(json.dumps(summary,indent=2))
