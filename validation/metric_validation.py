#!/usr/bin/env python3
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'validation'))
from metrics import auroc, average_precision, precision_at_k, recall_at_k, spatial_recovery, nearest_distances, exposure_matched_permutation, permutation_pvalue

checks=[]
def check(name, condition, detail=None):
    checks.append({'name':name,'pass':bool(condition),'detail':detail})

# Independent sklearn cross-checks across fixed random datasets, including ties.
rng=np.random.default_rng(20260904)
max_auc_err=0.0; max_ap_err=0.0
for rep in range(100):
    n=50+rng.integers(0,150)
    y=(rng.random(n)<rng.uniform(0.03,0.35)).astype(int)
    if y.sum()==0: y[0]=1
    if y.sum()==n: y[0]=0
    s=np.round(rng.normal(size=n),2)  # deliberate ties
    a1=auroc(y,s); a2=float(roc_auc_score(y,s))
    p1=average_precision(y,s); p2=float(average_precision_score(y,s))
    max_auc_err=max(max_auc_err,abs(a1-a2)); max_ap_err=max(max_ap_err,abs(p1-p2))
check('AUROC agrees with sklearn on 100 tied datasets',max_auc_err<1e-12,{'max_abs_error':max_auc_err})
check('Average precision agrees with sklearn on 100 tied datasets',max_ap_err<1e-12,{'max_abs_error':max_ap_err})

# Analytic edge cases.
y=np.array([1,1,0,0]); perfect=np.array([4,3,2,1],float); reverse=-perfect; tied=np.ones(4)
check('Perfect AUROC = 1',auroc(y,perfect)==1.0)
check('Reverse AUROC = 0',auroc(y,reverse)==0.0)
check('All-tied AUROC = 0.5',auroc(y,tied)==0.5)
check('Perfect AP = 1',average_precision(y,perfect)==1.0)
check('P@2 perfect = 1',precision_at_k(y,perfect,2)==1.0)
check('R@2 perfect = 1',recall_at_k(y,perfect,2)==1.0)

# Monotonic-score invariance for ranking metrics.
s=rng.normal(size=200); y=(rng.random(200)<0.08).astype(int); y[0]=1; y[1]=0
check('AUROC invariant to positive affine transform',abs(auroc(y,s)-auroc(y,7*s+13))<1e-15)
check('AP invariant to positive affine transform',abs(average_precision(y,s)-average_precision(y,7*s+13))<1e-15)

# Spatial metrics with exact known geometry.
true=np.array([[0,0,0],[5,0,0],[20,0,0]],float); pred=np.array([[0,0,0],[10,0,0]],float)
check('Spatial recovery R=0 exact',abs(spatial_recovery(true,pred,0)-1/3)<1e-15)
check('Spatial recovery R=5 known',abs(spatial_recovery(true,pred,5)-2/3)<1e-15)
check('Spatial recovery monotonic in radius',spatial_recovery(true,pred,8)<=spatial_recovery(true,pred,10))
nd=nearest_distances(true,pred)
check('Nearest distances known',np.allclose(nd,[0,5,10]))

# Exposure-matched permutation preserves positive counts in each stratum in expectation by construction;
# test null calibration on scores independent of labels. Across repeated seeds p-values should not be systematically tiny.
n=240; strata=np.repeat(np.arange(4),60); scores=rng.normal(size=n); y=np.zeros(n,int)
for g in range(4): y[rng.choice(np.where(strata==g)[0],size=5,replace=False)]=1
obs=auroc(y,scores); null=exposure_matched_permutation(scores,y,strata,auroc,n_perm=2000,seed=14)
p=permutation_pvalue(obs,null,'greater')
check('Exposure-matched permutation produces finite 2000-null distribution',len(null)==2000 and np.all(np.isfinite(null)),{'p':p,'null_mean':float(np.mean(null))})
check('Permutation AUROC null centered near 0.5',abs(float(np.mean(null))-0.5)<0.04,{'null_mean':float(np.mean(null))})

# Sparse-positive scenario: explicitly demonstrate why AP/prevalence accompanies AUROC.
y=np.zeros(1000,int); y[:10]=1; s=np.zeros(1000); s[:5]=2; s[5:10]=-2; s[10:]=0
auc=auroc(y,s); ap=average_precision(y,s); prev=float(y.mean())
check('AP enrichment is finite for sparse labels',math.isfinite(ap/prev),{'AUROC':auc,'AP':ap,'prevalence':prev,'AP_over_prevalence':ap/prev})

ok=all(c['pass'] for c in checks)
report={'status':'PASS' if ok else 'FAIL','classification':'metric implementation validation','reference':'scikit-learn + analytic cases','checks':checks}
out=ROOT/'validation'/'metric_validation_report.json'; out.write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
if not ok: raise SystemExit(1)
