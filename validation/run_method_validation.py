from __future__ import annotations

import json, math, tempfile, urllib.request
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import backend.main as ism

OUT=Path('validation/results'); OUT.mkdir(parents=True, exist_ok=True)
PANEL=['1CRN','1UBQ','2CI2','1LYZ','5PTI','1AKE','1TIM','1HRC','1MBN','4F5S']
POINTS=[50,100,200,500,1000]
THRESHOLDS=[0.03,0.05,0.075,0.10]
RADII=[4.0,5.0,6.0,8.0,10.0,12.0,15.0]
PH=7.4

class FirstChain(Select):
    def __init__(self,cid): self.cid=cid
    def accept_chain(self,c): return 1 if str(c.id)==self.cid else 0
    def accept_residue(self,r): return 1 if is_aa(r,standard=True) else 0

def download_first_chain(pid:str)->Path:
    raw=OUT/f'{pid}_raw.pdb'; clean=OUT/f'{pid}.pdb'
    if not raw.exists(): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb',raw)
    s=PDBParser(QUIET=True).get_structure(pid,str(raw)); model=next(s.get_models())
    cid=next(c.id for c in model if any(is_aa(r,standard=True) for r in c))
    io=PDBIO(); io.set_structure(s); io.save(str(clean),FirstChain(str(cid)))
    return clean

def build(pid,points,probe=1.4):
    oldp,oldprobe=ism.SASA_POINTS,ism.SASA_PROBE_A
    ism.SASA_POINTS=int(points); ism.SASA_PROBE_A=float(probe)
    try:
        _,res,_,_=ism.build_surface_residues(download_first_chain(pid),PH)
        return res
    finally:
        ism.SASA_POINTS,ism.SASA_PROBE_A=oldp,oldprobe

def dist(surface):
    if not surface:return np.empty((0,0))
    x=np.asarray([[r['x'],r['y'],r['z']] for r in surface],float)
    d=x[:,None,:]-x[None,:,:]
    return np.sqrt(np.sum(d*d,axis=2))

def score(res,threshold,chem,r1=5.0,r2=8.0,operator='min'):
    surf=[r for r in res if r['scrsa_raw']>=threshold]
    d=dist(surf); defs=ism.CHEMISTRIES[chem]['favorable']; n=len(surf)
    local=np.zeros(n)
    for i,r in enumerate(surf):
        if r['res_name'] in defs:
            state=defs[r['res_name']][2]
            local[i]=r['scrsa']*ism.state_availability(r['res_name'],state,PH)
    def dn(R):
        if n==0:return np.zeros(0)
        a=(d<=R).astype(float)@local
        m=a.max() if len(a) else 0
        return a/m if m>0 else np.zeros_like(a)
    a,b=dn(r1),dn(r2)
    if operator=='min': m=np.minimum(a,b)
    elif operator=='mean':m=(a+b)/2
    elif operator=='geo':m=np.sqrt(a*b)
    elif operator=='r1':m=a
    elif operator=='r2':m=b
    else:raise ValueError(operator)
    return {r['key']:float(m[i]) for i,r in enumerate(surf)}, {r['key']:(r['x'],r['y'],r['z']) for r in surf}, (a,b)

def vecmap(a,b):
    keys=sorted(set(a)|set(b)); return np.array([a.get(k,0.) for k in keys]),np.array([b.get(k,0.) for k in keys])
def rho(a,b):
    x,y=vecmap(a,b)
    if len(x)<3 or np.all(x==x[0]) or np.all(y==y[0]):return np.nan
    return float(spearmanr(x,y).statistic)
def topkeys(a,k=10):return [x[0] for x in sorted(a.items(),key=lambda z:(-z[1],z[0]))[:k] if xval(z:=z)>0]
def xval(z): return z[1]
def jacc(a,b,k=10):
    A=set(topkeys(a,k));B=set(topkeys(b,k)); return len(A&B)/len(A|B) if A|B else np.nan
def topcenter(a): return max(a,key=a.get) if a and max(a.values())>0 else None
def center_disp(a,ca,b,cb):
    x,y=topcenter(a),topcenter(b)
    if x is None or y is None:return np.nan
    p=np.asarray(ca[x]);q=np.asarray(cb[y]);return float(np.linalg.norm(p-q))

# Download once
for p in PANEL: download_first_chain(p)

# V1/V2 SASA convergence + final ranking stability
conv=[]; rank=[]
for pid in PANEL:
    by={n:build(pid,n) for n in POINTS}; ref={r['key']:r for r in by[1000]}
    for n in POINTS[:-1]:
        cur={r['key']:r for r in by[n]}; ks=sorted(set(cur)&set(ref))
        dif=np.array([cur[k]['scrsa_raw']-ref[k]['scrsa_raw'] for k in ks])
        surf_dis=np.mean([(cur[k]['scrsa_raw']>=.05)!=(ref[k]['scrsa_raw']>=.05) for k in ks])
        sr=spearmanr([cur[k]['scrsa_raw'] for k in ks],[ref[k]['scrsa_raw'] for k in ks]).statistic
        conv.append([pid,n,float(np.mean(np.abs(dif))),float(np.max(np.abs(dif))),float(sr),float(surf_dis)])
    for chem in ism.CHEMISTRIES:
        s200,c200,_=score(by[200],.05,chem); s1000,c1000,_=score(by[1000],.05,chem)
        rank.append([pid,chem,rho(s200,s1000),jacc(s200,s1000),center_disp(s200,c200,s1000,c1000)])
pd.DataFrame(conv,columns=['pdb','points','mae_scrsa_vs1000','max_abs_scrsa_vs1000','spearman_scrsa','surface_disagreement']).to_csv(OUT/'sasa_convergence.csv',index=False)
pd.DataFrame(rank,columns=['pdb','chemistry','spearman_M_200_vs1000','top10_jaccard_200_vs1000','top_center_displacement_A']).to_csv(OUT/'ranking_convergence.csv',index=False)

# V3 threshold robustness
thr=[]
for pid in PANEL:
    res=build(pid,200); nall=len(res)
    for t in THRESHOLDS:
        retained=sum(r['scrsa_raw']>=t for r in res)/nall
        for chem in ism.CHEMISTRIES:
            st,ct,_=score(res,t,chem); sb,cb,_=score(res,.05,chem)
            thr.append([pid,chem,t,retained,rho(st,sb),jacc(st,sb),center_disp(st,ct,sb,cb)])
pd.DataFrame(thr,columns=['pdb','chemistry','threshold','surface_fraction','spearman_vs_0.05','top10_jaccard_vs_0.05','top_center_displacement_A']).to_csv(OUT/'threshold_robustness.csv',index=False)

# V4 radius-pair audit under threshold perturbation (adsorption-independent)
rad=[]
for pid in PANEL:
    res=build(pid,200)
    for r1,r2 in combinations(RADII,2):
        base,cb,_=score(res,.05,'hydrophobic',r1,r2)
        # aggregate across all chemistries and threshold perturbations
        vals=[]
        for chem in ism.CHEMISTRIES:
            b,cbase,_=score(res,.05,chem,r1,r2)
            for t in [.03,.075]:
                q,cq,_=score(res,t,chem,r1,r2)
                vals.append((rho(q,b),jacc(q,b),center_disp(q,cq,b,cbase)))
        arr=np.array(vals,float)
        rad.append([pid,r1,r2,np.nanmedian(arr[:,0]),np.nanmedian(arr[:,1]),np.nanmedian(arr[:,2])])
pd.DataFrame(rad,columns=['pdb','r1','r2','median_spearman_threshold_perturb','median_top10_jaccard_threshold_perturb','median_top_center_disp_A']).to_csv(OUT/'radius_pair_audit.csv',index=False)

# V5 aggregation operator behavioral audit
abl=[]
for pid in PANEL:
    res=build(pid,200)
    for chem in ism.CHEMISTRIES:
        base,_,ab=score(res,.05,chem,5,8,'min')
        a,b=ab
        imbalance=float(np.mean((np.minimum(a,b)/(np.maximum(a,b)+1e-12))<0.5)) if len(a) else np.nan
        for op in ['min','mean','geo','r1','r2']:
            q,_,_=score(res,.05,chem,5,8,op)
            abl.append([pid,chem,op,rho(q,base),jacc(q,base),imbalance])
pd.DataFrame(abl,columns=['pdb','chemistry','operator','spearman_vs_min','top10_jaccard_vs_min','fraction_centers_scale_imbalanced']).to_csv(OUT/'aggregation_ablation.csv',index=False)

# V6 pH state analytical tests
phrows=[]
for rn,pka in ism.PKA.items():
    for state in ['protonated','deprotonated']:
        for delta in [-1,0,1]:
            obs=ism.state_availability(rn,state,pka+delta)
            exp=(1/(1+10**delta)) if state=='protonated' else (1/(1+10**(-delta)))
            phrows.append([rn,state,delta,obs,exp,abs(obs-exp)])
pd.DataFrame(phrows,columns=['residue','state','pH_minus_pKa','observed','expected','abs_error']).to_csv(OUT/'ph_state_unit_tests.csv',index=False)

# Summary
C=pd.read_csv(OUT/'sasa_convergence.csv');R=pd.read_csv(OUT/'ranking_convergence.csv');T=pd.read_csv(OUT/'threshold_robustness.csv');A=pd.read_csv(OUT/'radius_pair_audit.csv');H=pd.read_csv(OUT/'ph_state_unit_tests.csv')
summary={
 'panel':PANEL,
 'sasa_200_vs_1000':{
   'median_mae_scrsa':float(C[C.points==200].mae_scrsa_vs1000.median()),
   'max_mae_across_proteins':float(C[C.points==200].mae_scrsa_vs1000.max()),
   'median_spearman_scrsa':float(C[C.points==200].spearman_scrsa.median()),
   'median_surface_disagreement':float(C[C.points==200].surface_disagreement.median()),
   'median_M_spearman_all_chemistries':float(R.spearman_M_200_vs1000.median()),
   'median_top10_jaccard':float(R.top10_jaccard_200_vs1000.median()),
   'median_top_center_displacement_A':float(R.top_center_displacement_A.median())},
 'thresholds':T.groupby('threshold').agg({'surface_fraction':'median','spearman_vs_0.05':'median','top10_jaccard_vs_0.05':'median','top_center_displacement_A':'median'}).reset_index().to_dict('records'),
 'radius_pairs':A.groupby(['r1','r2']).agg({'median_spearman_threshold_perturb':'median','median_top10_jaccard_threshold_perturb':'median','median_top_center_disp_A':'median'}).reset_index().sort_values(['median_top_center_disp_A','median_top10_jaccard_threshold_perturb'],ascending=[True,False]).head(10).to_dict('records'),
 'ph_max_abs_error':float(H.abs_error.max())
}
(OUT/'method_validation_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
