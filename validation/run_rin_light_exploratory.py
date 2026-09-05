from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from Bio.PDB import PDBParser

ROOT=Path('validation/packages/all15')
BASE=ROOT/'results_is_apbs_gnm_exploratory'/'is_apbs_gnm_residue_scores.csv'
GEOM=ROOT/'results_light_geometry_exploratory'/'light_geometry_residue_scores.csv'
PDBDIR=ROOT/'results_ready_batch'/'pdb_cache'
OUT=ROOT/'results_rin_light_exploratory'; OUT.mkdir(parents=True,exist_ok=True)

D=pd.read_csv(BASE)
G=pd.read_csv(GEOM)
D=D.merge(G[['condition_id','key','radial_prominence']],on=['condition_id','key'],how='left')
PDB_IDS={'Fibronectin III8-10':'1FNF','Cytochrome c':'3NWV','Lysozyme':'1LYZ','Alpha-chymotrypsin':'4CHA','Myoglobin':'1MBN'}
CUTOFF=8.0

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

def rin(pid):
    p=PDBDIR/f'{pid}.pdb'; s=PDBParser(QUIET=True).get_structure(pid,str(p)); model=next(s.get_models())
    keys=[]; xyz=[]
    for ch in model:
        for r in ch:
            if 'CA' in r:
                keys.append(f"{ch.id}:{int(r.id[1])}:{str(r.id[2]).strip()}"); xyz.append(np.asarray(r['CA'].coord,float))
    X=np.vstack(xyz); n=len(X); A=np.zeros((n,n),dtype=int)
    dist=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2)
    A[(dist<=CUTOFF)&(dist>0)]=1
    degree=A.sum(1).astype(float)
    # Floyd-Warshall on unweighted graph
    INF=1e9; SP=np.where(A>0,1.0,INF); np.fill_diagonal(SP,0.0)
    for k in range(n): SP=np.minimum(SP,SP[:,k,None]+SP[None,k,:])
    closeness=np.zeros(n)
    for i in range(n):
        finite=(SP[i]<INF)&(np.arange(n)!=i)
        closeness[i]=finite.sum()/SP[i,finite].sum() if finite.any() and SP[i,finite].sum()>0 else 0.0
    # Brandes betweenness, unweighted
    bet=np.zeros(n)
    from collections import deque
    for ss in range(n):
        S=[]; P=[[] for _ in range(n)]; sigma=np.zeros(n); sigma[ss]=1.; dd=-np.ones(n,int); dd[ss]=0; Q=deque([ss])
        while Q:
            v=Q.popleft(); S.append(v)
            for w in np.where(A[v]>0)[0]:
                if dd[w]<0: Q.append(w); dd[w]=dd[v]+1
                if dd[w]==dd[v]+1: sigma[w]+=sigma[v]; P[w].append(v)
        delta=np.zeros(n)
        while S:
            w=S.pop()
            if sigma[w]>0:
                for v in P[w]: delta[v]+=(sigma[v]/sigma[w])*(1+delta[w])
            if w!=ss: bet[w]+=delta[w]
    bet*=0.5
    if n>2: bet/=((n-1)*(n-2)/2)
    return {k:(degree[i],closeness[i],bet[i]) for i,k in enumerate(keys)}

cache={}; rows=[]; detail=[]
for cond,sub in D.groupby('condition_id',sort=False):
    protein=sub['protein'].iloc[0]; pid=PDB_IDS[protein]
    if pid not in cache: cache[pid]=rin(pid)
    r=cache[pid]; keys=sub['key'].astype(str).tolist(); y=sub['truth'].to_numpy(int)
    isr=pct(sub['IS_static'].to_numpy(float)); gr=pct(sub['GNM_norm'].to_numpy(float)); rr=pct(sub['radial_prominence'].fillna(0).to_numpy(float))
    if bool(sub['apbs_applicable'].iloc[0]): ar=pct(np.nan_to_num(sub['APBS_compat'].to_numpy(float),nan=0.0))
    else: ar=np.ones(len(sub))
    base=isr*gr*ar*rr
    deg=np.asarray([r.get(k,(0,0,0))[0] for k in keys]); clo=np.asarray([r.get(k,(0,0,0))[1] for k in keys]); bet=np.asarray([r.get(k,(0,0,0))[2] for k in keys])
    # adsorption candidate hypothesis: peripheral/low-centrality residues favored
    per_deg=1-pct(deg); per_clo=1-pct(clo); per_bet=1-pct(bet)
    methods={'IS_base':base,'RIN_degree_peripheral_only':per_deg,'RIN_closeness_peripheral_only':per_clo,'RIN_betweenness_peripheral_only':per_bet,'IS_x_degree_peripheral':base*per_deg,'IS_x_closeness_peripheral':base*per_clo,'IS_x_betweenness_peripheral':base*per_bet}
    for m,s in methods.items(): rows.append({'condition_id':cond,'protein':protein,'pdb':pid,'method':m,'AUROC':auc(y,s),'AP':ap(y,s),'Recall@5':rec(y,s,5),'Recall@10':rec(y,s,10),'Recall@20':rec(y,s,20)})
    for i,k in enumerate(keys): detail.append({'condition_id':cond,'key':k,'truth':int(y[i]),'degree':deg[i],'closeness':clo[i],'betweenness':bet[i]})
res=pd.DataFrame(rows); res.to_csv(OUT/'rin_method_comparison.csv',index=False); pd.DataFrame(detail).to_csv(OUT/'rin_residue_scores.csv',index=False)
metrics=['AUROC','AP','Recall@5','Recall@10','Recall@20']
summary={'median':res.groupby('method')[metrics].median().to_dict(orient='index'),'mean':res.groupby('method')[metrics].mean().to_dict(orient='index'),'design':{'graph':'C-alpha unweighted contact graph','cutoff_A':CUTOFF,'centralities':['degree','closeness','betweenness'],'combination':'percentile-rank product with current IS base; no fitted coefficients'}}
(OUT/'rin_summary.json').write_text(json.dumps(summary,indent=2)); print(json.dumps(summary,indent=2))
