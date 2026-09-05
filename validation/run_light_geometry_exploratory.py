from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from Bio.PDB import PDBParser

ROOT=Path('validation/packages/all15')
SRC=ROOT/'results_is_apbs_gnm_exploratory'/'is_apbs_gnm_residue_scores.csv'
PDBDIR=ROOT/'results_ready_batch'/'pdb_cache'
OUT=ROOT/'results_light_geometry_exploratory'; OUT.mkdir(parents=True,exist_ok=True)

DETAIL=pd.read_csv(SRC)
CONDITIONS=DETAIL['condition_id'].drop_duplicates().tolist()
PDB_BY_COND=DETAIL.groupby('condition_id')['protein'].first().to_dict()
PDB_IDS={'Fibronectin III8-10':'1FNF','Cytochrome c':'3NWV','Lysozyme':'1LYZ','Alpha-chymotrypsin':'4CHA','Myoglobin':'1MBN'}
TRUTH=DETAIL[['condition_id','res_seq','truth']]

# lightweight, coordinate-only descriptors
# 1) inverse local CA packing within 10 A (fewer neighbours => more protruding)
# 2) radial prominence relative to protein centroid
# 3) side-chain outwardness: cosine between centroid->CA and CA->sidechain-centroid

def pct(x):
    x=np.asarray(x,float)
    if len(x)<=1:return np.ones_like(x)
    return (rankdata(x,method='average')-1)/(len(x)-1)

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
    o=np.argsort(-s,kind='stable')[:min(k,len(s))]; return float(y[o].sum()/y.sum())

def geom_for_pdb(pid):
    p=PDBDIR/f'{pid}.pdb'; s=PDBParser(QUIET=True).get_structure(pid,str(p)); model=next(s.get_models())
    rows=[]; all_ca=[]
    for ch in model:
        for r in ch:
            if 'CA' in r:
                ca=np.asarray(r['CA'].coord,float); all_ca.append(ca)
                sc=[]
                for a in r.get_atoms():
                    nm=a.get_name().strip().upper()
                    if nm not in {'N','CA','C','O','OXT'} and getattr(a,'element','')!='H': sc.append(np.asarray(a.coord,float))
                sc_cent=np.mean(np.vstack(sc),axis=0) if sc else ca.copy()
                rows.append([f"{ch.id}:{int(r.id[1])}:{str(r.id[2]).strip()}",ca,sc_cent])
    centroid=np.mean(np.vstack(all_ca),axis=0)
    X=np.vstack([r[1] for r in rows]); dist=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2)
    out={}
    for i,(key,ca,sc) in enumerate(rows):
        n10=max(int(np.sum(dist[i]<=10.0))-1,0)
        invpack=1.0/(1.0+n10)
        radial=float(np.linalg.norm(ca-centroid))
        v1=ca-centroid; v2=sc-ca
        nv1=np.linalg.norm(v1); nv2=np.linalg.norm(v2)
        outward=float(np.dot(v1,v2)/(nv1*nv2)) if nv1>0 and nv2>1e-8 else 0.0
        out[key]=(invpack,radial,outward,n10)
    return out

GCACHE={}
rows=[]; drows=[]
for cond in CONDITIONS:
    sub=DETAIL[DETAIL.condition_id==cond].copy()
    protein=sub['protein'].iloc[0]; pid=PDB_IDS[protein]
    if pid not in GCACHE: GCACHE[pid]=geom_for_pdb(pid)
    g=GCACHE[pid]
    keys=sub['key'].astype(str).tolist(); y=sub['truth'].to_numpy(int)
    isr=pct(sub['IS_static'].to_numpy(float)); gr=pct(sub['GNM_norm'].to_numpy(float))
    # conditional APBS: if descriptor unavailable/non-applicable, neutral factor 1
    if bool(sub['apbs_applicable'].iloc[0]):
        ar=pct(np.nan_to_num(sub['APBS_compat'].to_numpy(float),nan=0.0))
    else:
        ar=np.ones(len(sub))
    base=isr*gr*ar
    inv=np.asarray([g.get(k,(0,0,0,0))[0] for k in keys]); rad=np.asarray([g.get(k,(0,0,0,0))[1] for k in keys]); outw=np.asarray([g.get(k,(0,0,0,0))[2] for k in keys])
    ir,rr,orank=pct(inv),pct(rad),pct(outw)
    methods={
        'IS_GNM_APBS_base':base,
        'geometry_invpack_only':inv,
        'geometry_radial_only':rad,
        'geometry_outward_only':outw,
        'base_x_invpack':base*ir,
        'base_x_radial':base*rr,
        'base_x_outward':base*orank,
        'base_x_invpack_x_outward':base*ir*orank,
    }
    for m,s in methods.items():
        rows.append({'condition_id':cond,'protein':protein,'pdb':pid,'method':m,'AUROC':auc(y,s),'AP':ap(y,s),'Recall@5':rec(y,s,5),'Recall@10':rec(y,s,10),'Recall@20':rec(y,s,20)})
    for i,k in enumerate(keys):
        drows.append({'condition_id':cond,'key':k,'res_seq':int(sub.iloc[i]['res_seq']),'truth':int(y[i]),'invpack10':inv[i],'radial_prominence':rad[i],'sidechain_outwardness':outw[i],'n_neighbors10':g.get(k,(0,0,0,0))[3]})

res=pd.DataFrame(rows); res.to_csv(OUT/'light_geometry_method_comparison.csv',index=False)
pd.DataFrame(drows).to_csv(OUT/'light_geometry_residue_scores.csv',index=False)
metrics=['AUROC','AP','Recall@5','Recall@10','Recall@20']
summary={'median_metrics_by_method':res.groupby('method')[metrics].median().to_dict(orient='index'),'mean_metrics_by_method':res.groupby('method')[metrics].mean().to_dict(orient='index'),'notes':{'invpack10':'1/(1+number of CA neighbors within 10 A)','radial_prominence':'CA distance from protein CA centroid','sidechain_outwardness':'cosine between centroid->CA and CA->sidechain centroid','combination':'percentile-rank products; no fitted coefficients; no APBS recomputation'}}
(OUT/'light_geometry_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
