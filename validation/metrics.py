"""Validation metrics used by InterfaceScout external benchmarks.

These functions live under validation/ and are not part of the predictive model.
They are deterministic and operate on predeclared ground-truth labels/scores.
"""
from __future__ import annotations
import numpy as np


def _arrays(y, score):
    y=np.asarray(y,dtype=int); s=np.asarray(score,dtype=float)
    if y.ndim!=1 or s.ndim!=1 or len(y)!=len(s): raise ValueError('y and score must be same-length 1D arrays')
    if not set(np.unique(y)).issubset({0,1}): raise ValueError('labels must be binary')
    if not np.all(np.isfinite(s)): raise ValueError('scores must be finite')
    return y,s


def auroc(y, score):
    """Tie-aware Mann-Whitney AUROC."""
    y,s=_arrays(y,score); n1=int(y.sum()); n0=len(y)-n1
    if n1==0 or n0==0: return float('nan')
    order=np.argsort(s,kind='mergesort'); ranks=np.empty(len(s),float)
    i=0
    while i<len(s):
        j=i+1
        while j<len(s) and s[order[j]]==s[order[i]]: j+=1
        ranks[order[i:j]]=(i+1+j)/2.0
        i=j
    u=float(ranks[y==1].sum()-n1*(n1+1)/2.0)
    return u/(n1*n0)


def average_precision(y, score):
    """Average precision using threshold-grouped precision-recall increments."""
    y,s=_arrays(y,score); npos=int(y.sum())
    if npos==0: return float('nan')
    order=np.argsort(-s,kind='mergesort'); ys=y[order]; ss=s[order]
    tp=0; fp=0; ap=0.0; prev_recall=0.0; i=0
    while i<len(y):
        j=i
        while j<len(y) and ss[j]==ss[i]:
            if ys[j]: tp+=1
            else: fp+=1
            j+=1
        recall=tp/npos; precision=tp/(tp+fp)
        ap+=(recall-prev_recall)*precision; prev_recall=recall; i=j
    return float(ap)


def precision_at_k(y, score, k):
    y,s=_arrays(y,score); k=max(1,min(int(k),len(y)))
    idx=np.argsort(-s,kind='mergesort')[:k]
    return float(y[idx].mean())


def recall_at_k(y, score, k):
    y,s=_arrays(y,score); npos=int(y.sum())
    if npos==0: return float('nan')
    k=max(1,min(int(k),len(y))); idx=np.argsort(-s,kind='mergesort')[:k]
    return float(y[idx].sum()/npos)


def spatial_recovery(true_xyz, pred_xyz, radius):
    """Fraction of true contact coordinates within radius of any predicted center."""
    t=np.asarray(true_xyz,float); p=np.asarray(pred_xyz,float)
    if len(t)==0: return float('nan')
    if len(p)==0: return 0.0
    d=np.sqrt(((t[:,None,:]-p[None,:,:])**2).sum(axis=2))
    return float((d.min(axis=1)<=float(radius)).mean())


def nearest_distances(true_xyz, pred_xyz):
    t=np.asarray(true_xyz,float); p=np.asarray(pred_xyz,float)
    if len(t)==0: return np.asarray([],float)
    if len(p)==0: return np.full(len(t),np.inf)
    return np.sqrt(((t[:,None,:]-p[None,:,:])**2).sum(axis=2)).min(axis=1)


def permutation_pvalue(observed, null_values, alternative='greater'):
    null=np.asarray(null_values,float)
    if alternative=='greater': count=int(np.sum(null>=observed))
    elif alternative=='less': count=int(np.sum(null<=observed))
    else: raise ValueError('alternative must be greater or less')
    return float((count+1)/(len(null)+1))


def exposure_matched_permutation(scores, positives, strata, metric_fn, n_perm=10000, seed=20260904):
    """Permute positive labels within exposure strata, preserving positives/stratum."""
    scores=np.asarray(scores,float); y=np.asarray(positives,int); strata=np.asarray(strata)
    if len(scores)!=len(y) or len(y)!=len(strata): raise ValueError('length mismatch')
    rng=np.random.default_rng(seed); out=np.empty(int(n_perm),float)
    groups=[np.where(strata==g)[0] for g in np.unique(strata)]
    counts=[int(y[idx].sum()) for idx in groups]
    for b in range(int(n_perm)):
        yp=np.zeros_like(y)
        for idx,n in zip(groups,counts):
            if n: yp[rng.choice(idx,size=n,replace=False)]=1
        out[b]=metric_fn(yp,scores)
    return out
