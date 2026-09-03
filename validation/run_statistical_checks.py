from __future__ import annotations
import csv, json, random, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from backend import main as IS
from validation.run_external_benchmark import fetch_pdb_text, run_is, score_vectors, exposure_matched_permutation
IS.PDB2PQR=None; IS.APBS=None; IS.MKDSSP=None
OUT=ROOT/'validation_results'; OUT.mkdir(exist_ok=True)
RNG=np.random.default_rng(20260904)

def labels(res,nums):
    s=set(nums); return np.array([1 if int(r['res_seq']) in s else 0 for r in res['surface_residues']],int)

def bootstrap_ci(y,s,n=5000):
    pos=np.where(y==1)[0]; neg=np.where(y==0)[0]
    if not len(pos) or not len(neg): return {}
    auc=[]; ap=[]
    for _ in range(n):
        ip=RNG.choice(pos,len(pos),replace=True); ineg=RNG.choice(neg,len(neg),replace=True); idx=np.r_[ip,ineg]
        yy=y[idx]; ss=s[idx]
        auc.append(roc_auc_score(yy,ss)); ap.append(average_precision_score(yy,ss))
    return {'AUROC_lo':float(np.quantile(auc,.025)),'AUROC_hi':float(np.quantile(auc,.975)),
            'AP_lo':float(np.quantile(ap,.025)),'AP_hi':float(np.quantile(ap,.975))}

def bh(ps):
    p=np.asarray(ps,float); n=len(p); order=np.argsort(p); q=np.empty(n); prev=1.0
    for rank_idx in range(n-1,-1,-1):
        idx=order[rank_idx]; rank=rank_idx+1; val=min(prev,p[idx]*n/rank); q[idx]=val; prev=val
    return q.tolist()

def main():
    pdb={x:fetch_pdb_text(x) for x in ['1MBN','2HHB','2PTN','5H7A']}; rows=[]
    # Exact SpA primary maps.
    specs=[('SpA_Au111','hydrophobic',[221,220,218,33,34]),('SpA_O_rich','anionic',[33,34,35,36,37]),('SpA_Si_rich','cationic',[221,220,219])]
    for ch in ['C','B']:
        res=run_is(pdb['5H7A'],ch,7.0,298.0,20.0)
        for case,chem,nums in specs:
            _,scrsa,_,_,persist,_,_=score_vectors(res,chem); y=labels(res,nums)
            au=float(roc_auc_score(y,persist)); ap=float(average_precision_score(y,persist)); pm=exposure_matched_permutation(y,persist,scrsa,10000,.05)
            rows.append({'group':'SpA_exact','case':f'{case}_chain_{ch}','chemistry':chem,'AUROC':au,'AP':ap,
                         'p_exposure_matched':pm['empirical_p_ge'],**bootstrap_ci(y,persist)})
    # Tavanti region labels; pH 7.4 is an InterfaceScout evaluation condition, not claimed as a reported CG pH.
    tav=[('1MBN','1MBN','A',[43,45]+list(range(96,100))+list(range(146,154))),
         ('2PTN','2PTN','A',[94]+list(range(125,136))+[166,167]+list(range(231,245))),
         ('2HHB_A','2HHB','A',list(range(12,26))+list(range(61,82))),
         ('2HHB_C','2HHB','C',list(range(12,26))+list(range(61,82))),
         ('2HHB_B','2HHB','B',list(range(51,54))),('2HHB_D','2HHB','D',list(range(45,54)))]
    for case,pid,ch,nums in tav:
        res=run_is(pdb[pid],ch,7.4,310.0,150.0)
        for chem in ['anionic','hydrophobic','hbond_acceptor']:
            _,scrsa,_,_,persist,_,_=score_vectors(res,chem); y=labels(res,nums)
            au=float(roc_auc_score(y,persist)); ap=float(average_precision_score(y,persist)); pm=exposure_matched_permutation(y,persist,scrsa,10000,.05)
            rows.append({'group':'Tavanti_region','case':case,'chemistry':chem,'AUROC':au,'AP':ap,
                         'p_exposure_matched':pm['empirical_p_ge'],**bootstrap_ci(y,persist)})
    qs=bh([r['p_exposure_matched'] for r in rows])
    for r,q in zip(rows,qs): r['q_BH_all_tests']=q
    fields=sorted(set().union(*(r.keys() for r in rows)))
    with (OUT/'statistical_checks.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    (OUT/'statistical_checks.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__': main()
