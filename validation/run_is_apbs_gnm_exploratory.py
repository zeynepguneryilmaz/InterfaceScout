from __future__ import annotations
import json, sys, tempfile, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from Bio.PDB import PDBParser

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import backend.main as ism

BASE=Path('validation/packages/all15/results_ready_batch')
PDBDIR=BASE/'pdb_cache'; PDBDIR.mkdir(parents=True,exist_ok=True)
OUT=Path('validation/packages/all15/results_is_apbs_gnm_exploratory'); OUT.mkdir(parents=True,exist_ok=True)
CUTOFF_A=10.0
SASA_POINTS=120
IONIC_MM=150.0
TEMP_K=298.0

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

def ensure_pdb(pid):
    p=PDBDIR/f'{pid}.pdb'
    if not p.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',p)
    return p

def rkey(chain,res): return f"{chain.id}:{int(res.id[1])}:{str(res.id[2]).strip()}"

def gnm_norm(pid):
    p=ensure_pdb(pid); s=PDBParser(QUIET=True).get_structure(pid,str(p)); model=next(s.get_models())
    keys=[]; xyz=[]
    for ch in model:
        for r in ch:
            if 'CA' in r:
                keys.append(rkey(ch,r)); xyz.append(np.asarray(r['CA'].coord,float))
    X=np.vstack(xyz); n=len(X); K=np.zeros((n,n),float)
    for i in range(n-1):
        d=np.linalg.norm(X[i+1:]-X[i],axis=1); js=np.where(d<=CUTOFF_A)[0]+i+1
        for j in js: K[i,j]=K[j,i]=-1.; K[i,i]+=1.; K[j,j]+=1.
    vals,vecs=np.linalg.eigh(K); pos=vals>1e-8
    cov=(vecs[:,pos]*(1./vals[pos]))@vecs[:,pos].T; msf=np.diag(cov)
    lo,hi=float(msf.min()),float(msf.max()); norm=(msf-lo)/(hi-lo) if hi>lo else np.zeros_like(msf)
    return dict(zip(keys,norm))

def get_features(pid,chem,ph):
    p=ensure_pdb(pid); old=ism.SASA_POINTS; ism.SASA_POINTS=SASA_POINTS
    try:
        _,res,atoms,_=ism.build_surface_residues(p,ph)
        env=ism.EnvParams(pH=ph,ionic=IONIC_MM,temp=TEMP_K)
        with tempfile.TemporaryDirectory(prefix='is_apbs_') as td:
            status=ism.attach_apbs_auxiliary(p,res,atoms,env,Path(td))
    finally: ism.SASA_POINTS=old
    surf=[r for r in res if r['surface_exposed']]; D=ism.build_distances(surf); mp=ism.chemistry_map(surf,D,chem,ph)
    sm={x['center_key']:x['multiscale_persistence']/100.0 for x in mp['patch_centers']}
    expected=ism.CHEMISTRIES[chem].get('expected_phi_sign')
    out={}
    for r in surf:
        phi=r.get('phi'); compat=None
        if expected=='positive' and phi is not None: compat=max(float(phi),0.0)
        elif expected=='negative' and phi is not None: compat=max(-float(phi),0.0)
        out[r['key']]={'seq':int(r['res_seq']),'static':sm.get(r['key'],0.0),'phi':phi,'apbs_compat':compat}
    return out,status,expected

def pct(x):
    x=np.asarray(x,float)
    if len(x)<=1:return np.ones_like(x)
    return (rankdata(x,method='average')-1.0)/(len(x)-1.0)

def auc(y,s):
    p=s[y==1]; n=s[y==0]
    if not len(p) or not len(n):return np.nan
    return float(np.mean([(a>b)+.5*(a==b) for a in p for b in n]))

def ap(y,s):
    if not y.sum():return np.nan
    o=np.argsort(-s,kind='stable'); yy=y[o]; pr=np.cumsum(yy)/(np.arange(len(yy))+1)
    return float((pr*yy).sum()/yy.sum())

def rec(y,s,k):
    if not y.sum():return np.nan
    o=np.argsort(-s,kind='stable')[:min(k,len(s))]; return float(y[o].sum()/y.sum())

gcache={}; rows=[]; detail=[]
for b in BENCH:
    if b['pdb'] not in gcache:gcache[b['pdb']]=gnm_norm(b['pdb'])
    feat,status,expected=get_features(b['pdb'],b['chem'],b['pH']); keys=list(feat)
    st=np.asarray([feat[k]['static'] for k in keys]); gn=np.asarray([gcache[b['pdb']].get(k,0.) for k in keys])
    sr=pct(st); gr=pct(gn)
    applicable=expected in ('positive','negative')
    if applicable:
        apraw=np.asarray([feat[k]['apbs_compat'] if feat[k]['apbs_compat'] is not None else 0. for k in keys]); ar=pct(apraw)
    else:
        apraw=np.ones(len(keys)); ar=np.ones(len(keys))
    scores={'IS_static':st,'GNM_only':gn,'APBS_only':apraw if applicable else np.zeros(len(keys)),
            'IS_x_GNM':sr*gr,'IS_x_APBS':sr*ar,'IS_x_GNM_x_APBS':sr*gr*ar}
    y=np.asarray([int(feat[k]['seq'] in b['truth']) for k in keys])
    for name,s in scores.items():
        if name=='APBS_only' and not applicable: continue
        rows.append({'condition_id':b['id'],'protein':b['protein'],'pdb':b['pdb'],'chemistry':b['chem'],'method':name,'apbs_applicable':applicable,'apbs_status':status,'AUROC':auc(y,s),'AP':ap(y,s),'Recall@5':rec(y,s,5),'Recall@10':rec(y,s,10),'Recall@20':rec(y,s,20)})
    for i,k in enumerate(keys):
        detail.append({'condition_id':b['id'],'protein':b['protein'],'key':k,'res_seq':feat[k]['seq'],'truth':int(y[i]),'IS_static':st[i],'GNM_norm':gn[i],'phi':feat[k]['phi'],'APBS_compat':feat[k]['apbs_compat'],'apbs_applicable':applicable,'apbs_status':status})

res=pd.DataFrame(rows); res.to_csv(OUT/'is_apbs_gnm_method_comparison.csv',index=False)
pd.DataFrame(detail).to_csv(OUT/'is_apbs_gnm_residue_scores.csv',index=False)
primary_methods=['IS_static','IS_x_GNM','IS_x_APBS','IS_x_GNM_x_APBS']
summary={'design':{'gnm_cutoff_A':CUTOFF_A,'sasa_points':SASA_POINTS,'ionic_strength_mM':IONIC_MM,'temperature_K':TEMP_K,'APBS_rule':'anionic uses max(+phi,0); cationic uses max(-phi,0); non-electrostatic chemistries use neutral factor 1','combination':'average-tie percentile-rank products; no fitted coefficients'},'all_conditions_median':res[res.method.isin(primary_methods)].groupby('method')[['AUROC','AP','Recall@5','Recall@10','Recall@20']].median().to_dict(orient='index'),'all_conditions_mean':res[res.method.isin(primary_methods)].groupby('method')[['AUROC','AP','Recall@5','Recall@10','Recall@20']].mean().to_dict(orient='index'),'charged_conditions_median':res[(res.apbs_applicable)&res.method.isin(primary_methods)].groupby('method')[['AUROC','AP','Recall@5','Recall@10','Recall@20']].median().to_dict(orient='index'),'apbs_statuses':res[['condition_id','apbs_status']].drop_duplicates().set_index('condition_id')['apbs_status'].to_dict()}
(OUT/'is_apbs_gnm_summary.json').write_text(json.dumps(summary,indent=2))
print(res.to_string(index=False)); print(json.dumps(summary,indent=2))
