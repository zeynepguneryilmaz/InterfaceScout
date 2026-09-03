from __future__ import annotations

import csv, json, math, random, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from backend import main as IS
from validation.run_external_benchmark import fetch_pdb_text, run_is, score_vectors, recompute_persistence

# Primary score is independent of optional APBS/PDB2PQR/DSSP; disable them for a fast sensitivity audit.
IS.PDB2PQR=None; IS.APBS=None; IS.MKDSSP=None
OUT=ROOT/'validation_results'; OUT.mkdir(exist_ok=True)
RNG=random.Random(20260904)

def auc_ap(y,s):
    if len(set(y.tolist()))<2: return None,None
    return float(roc_auc_score(y,s)), float(average_precision_score(y,s))

def labels(res, nums):
    ns=set(nums)
    return np.array([1 if int(r['res_seq']) in ns else 0 for r in res['surface_residues']],int)

def topk(y,s,k):
    order=np.argsort(-s,kind='mergesort')[:k]
    tp=int(np.sum(y[order])); pred=len(order); pos=int(np.sum(y))
    p=tp/pred if pred else 0.0; r=tp/pos if pos else 0.0
    f=2*p*r/(p+r) if p+r else 0.0
    return {'TP':tp,'precision':p,'recall':r,'F1':f}

def persistence_pair(local,coords,pair):
    return recompute_persistence(local,coords,pair)

def spatial_from_scores(res,s,anchors,k=10):
    surf=res['surface_residues']; coords=np.array([[r['x'],r['y'],r['z']] for r in surf],float)
    nums=np.array([int(r['res_seq']) for r in surf])
    inds=np.argsort(-s,kind='mergesort')[:k]; c=coords[inds]
    ds=[]
    for a in anchors:
        hit=np.where(nums==a)[0]
        if len(hit) and len(c): ds.append(float(np.min(np.linalg.norm(c-coords[hit[0]],axis=1))))
    return {'R5':float(np.mean(np.array(ds)<=5)) if ds else None,
            'R8':float(np.mean(np.array(ds)<=8)) if ds else None,
            'median_nearest_A':float(np.median(ds)) if ds else None}

def main():
    rows=[]; radius=[]; region=[]; threshold=[]
    pdb={x:fetch_pdb_text(x) for x in ['1MBN','2HHB','2PTN','5H7A']}
    # SpA exact anchors; both crystallographic copies C/B retained as a copy-sensitivity audit.
    specs=[('SpA_Au111','hydrophobic',[221,220,218,33,34]),
           ('SpA_O_rich_silica','anionic',[33,34,35,36,37]),
           ('SpA_Si_rich_silica','cationic',[221,220,219])]
    for ch in ['C','B']:
        res=run_is(pdb['5H7A'],ch,7.0,298.0,20.0)
        for case,chem,anc in specs:
            keys,s0,s1,s2,s3,ca,sc=score_vectors(res,chem); y=labels(res,anc)
            for name,s in [('scRSA',s0),('chemistry_exposure',s1),('state_score',s2),('persistence_5_8',s3)]:
                au,ap=auc_ap(y,s); k=len(anc); t=topk(y,s,k); t10=topk(y,s,10)
                rows.append({'case':f'{case}_chain_{ch}','chemistry':chem,'score':name,'AUROC':au,'AP':ap,
                             'n_surface':len(y),'n_positive':int(y.sum()),
                             'P_at_m':t['precision'],'R_at_m':t['recall'],'F1_at_m':t['F1'],
                             'P_at_10':t10['precision'],'R_at_10':t10['recall'],'F1_at_10':t10['F1']})
            # all radius-pair choices from the frozen developmental audit, using the same local score and C-alpha geometry
            pairs=[(5,8),(5,10),(5,12),(5,15),(8,10),(8,12),(8,15),(10,12),(10,15),(12,15)]
            for pair in pairs:
                sp=persistence_pair(s2,ca,pair); au,ap=auc_ap(y,sp); geo=spatial_from_scores(res,sp,anc,10)
                radius.append({'case':f'{case}_chain_{ch}','chemistry':chem,'r1':pair[0],'r2':pair[1],
                               'AUROC':au,'AP':ap,**geo})
    # Tavanti reported binding regions: evaluate region membership as a region-level label, not exact-contact labels.
    tav=[('1MBN_citrate_AuNP','1MBN','A',[43,45]+list(range(96,100))+list(range(146,154))),
         ('2PTN_citrate_AuNP','2PTN','A',[94]+list(range(125,136))+[166,167]+list(range(231,245))),
         ('2HHB_A_citrate_AuNP','2HHB','A',list(range(12,26))+list(range(61,82))),
         ('2HHB_C_citrate_AuNP','2HHB','C',list(range(12,26))+list(range(61,82))),
         ('2HHB_B_citrate_AuNP','2HHB','B',list(range(51,54))),
         ('2HHB_D_citrate_AuNP','2HHB','D',list(range(45,54)))]
    for case,pid,ch,reg in tav:
        res=run_is(pdb[pid],ch,7.4,310.0,150.0)
        for chem in ['anionic','hydrophobic','hbond_acceptor']:
            _,s0,s1,s2,s3,ca,_=score_vectors(res,chem); y=labels(res,reg)
            for name,s in [('scRSA',s0),('chemistry_exposure',s1),('state_score',s2),('persistence_5_8',s3)]:
                au,ap=auc_ap(y,s)
                region.append({'case':case,'chemistry':chem,'score':name,'AUROC':au,'AP':ap,
                               'n_surface':len(y),'n_region_surface':int(y.sum())})
            for pair in [(5,8),(8,10)]:
                sp=persistence_pair(s2,ca,pair); au,ap=auc_ap(y,sp)
                radius.append({'case':case,'chemistry':chem,'r1':pair[0],'r2':pair[1],
                               'AUROC':au,'AP':ap,'R5':None,'R8':None,'median_nearest_A':None})
    # External scRSA threshold sensitivity on representative exact and region-level systems.
    # This changes only the surface inclusion threshold, never chemistry assignments or radii.
    reps=[('SpA_O_rich_C','5H7A','C',7.0,298.0,20.0,'anionic',[33,34,35,36,37],'exact'),
          ('1MBN_citrate','1MBN','A',7.4,310.0,150.0,'anionic',[43,45]+list(range(96,100))+list(range(146,154)),'region'),
          ('2PTN_citrate','2PTN','A',7.4,310.0,150.0,'anionic',[94]+list(range(125,136))+[166,167]+list(range(231,245)),'region')]
    original=IS.SC_RSA_THRESHOLD
    try:
        for th in [0.02,0.03,0.05,0.075,0.10,0.15]:
            IS.SC_RSA_THRESHOLD=th
            for case,pid,ch,ph,temp,ionic,chem,lab,kind in reps:
                res=run_is(pdb[pid],ch,ph,temp,ionic); _,_,_,_,s,_,_=score_vectors(res,chem); y=labels(res,lab)
                au,ap=auc_ap(y,s)
                threshold.append({'case':case,'label_kind':kind,'threshold':th,'chemistry':chem,
                                  'n_surface':len(y),'n_positive_retained':int(y.sum()),'AUROC':au,'AP':ap})
    finally:
        IS.SC_RSA_THRESHOLD=original
    def write(name,data):
        if not data:return
        fields=sorted(set().union(*(d.keys() for d in data)))
        with (OUT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(data)
    write('extended_exact_metrics.csv',rows);write('radius_performance_sensitivity.csv',radius)
    write('tavanti_region_auc_ap.csv',region);write('external_threshold_sensitivity.csv',threshold)
    (OUT/'extended_checks.json').write_text(json.dumps({'exact':rows,'radius':radius,'region':region,'threshold':threshold},indent=2))

if __name__=='__main__': main()
