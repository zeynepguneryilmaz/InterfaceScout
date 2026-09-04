"""Development-only diagnostic for v5.3 established nonpolar physics.

IMPORTANT: 1AKI/PE and 5H7A/Au are development systems after prior inspection.
These results must not be presented as held-out validation.  The purpose is to
check whether literature-derived continuous hydrophobicity and hydrophobic-dipole
layers repair the identified model failure without fitting weights to labels.
"""
from __future__ import annotations

import json, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend.v52_app import AnalyzeRequest, EnvParams
from backend.v53_app import analyze

OUT = ROOT / "validation" / "v53_established_physics_diagnostic.json"


def auc(y, s):
    y=np.asarray(y,dtype=int); s=np.asarray(s,dtype=float)
    pos=np.where(y==1)[0]; neg=np.where(y==0)[0]
    if len(pos)==0 or len(neg)==0: return None
    wins=0.0
    for i in pos:
        wins += np.sum(s[i] > s[neg]) + 0.5*np.sum(s[i] == s[neg])
    return float(wins/(len(pos)*len(neg)))


def ap(y,s):
    y=np.asarray(y,dtype=int); s=np.asarray(s,dtype=float)
    npos=int(y.sum())
    if npos==0: return None
    order=np.argsort(-s, kind='mergesort')
    yy=y[order]; tp=0; vals=[]
    for k,v in enumerate(yy,1):
        if v:
            tp += 1; vals.append(tp/k)
    return float(np.mean(vals))


def metrics(rows, anchors, score_key):
    anchor=set(anchors)
    y=[1 if int(r['res_seq']) in anchor else 0 for r in rows]
    s=[float(r.get(score_key,0.0) or 0.0) for r in rows]
    order=np.argsort(-np.asarray(s), kind='mergesort')
    k=min(10,len(rows)); top=[rows[i] for i in order[:k]]
    return {
        'n':len(rows),'n_anchor_surface':int(sum(y)),
        'auroc':auc(y,s),'ap':ap(y,s),
        'recall_at_10':float(sum(1 for r in top if int(r['res_seq']) in anchor)/max(sum(y),1)),
        'top10':[[r['res_name'],int(r['res_seq']),round(float(r.get(score_key,0.0) or 0.0),3)] for r in top],
    }


def spatial(rows, anchors, top_rows):
    coords={int(r['res_seq']):np.asarray([r['x'],r['y'],r['z']],float) for r in rows}
    tc=[np.asarray([r['x'],r['y'],r['z']],float) for r in top_rows]
    vals={}
    dists=[]
    for a in anchors:
        if a not in coords or not tc: continue
        d=min(float(np.linalg.norm(coords[a]-x)) for x in tc); dists.append(d)
    for R in [5.0,8.0,10.0]:
        vals[f'recall_{int(R)}A']=float(sum(d<=R for d in dists)/len(dists)) if dists else None
    vals['median_nearest_A']=float(np.median(dists)) if dists else None
    vals['n_anchor_surface']=len(dists)
    return vals


def run_case(pdb, chain, pH, ionic, anchors, label):
    req=AnalyzeRequest(pdb_id=pdb, chain=chain, env=EnvParams(pH=pH,ionic=ionic,temp=298.0), structure_context='deposited_structure', protrusion=True)
    r=analyze(req)
    surface=r['surface_residues']
    old=r['chemistries']['hydrophobic']['residues']
    old_by={x['key']:x for x in old}
    baseline=[]
    for x in surface:
        z=dict(x); z['old_M']=float(old_by.get(x['key'],{}).get('persistence',0.0) or 0.0); baseline.append(z)
    new=r['nonpolar_physics']['primary_local_field']['residues']
    gated=r['nonpolar_physics']['orientation_gated_top_patches']
    new_top=r['nonpolar_physics']['primary_local_field']['top_patches']
    ww=r['nonpolar_physics']['independent_interface_scale_sensitivity']['residues']
    return {
        'label':label,'pdb':pdb,'chain':chain,'anchors':anchors,
        'old_binary_hydrophobic_M':metrics(baseline,anchors,'old_M'),
        'scrsa_baseline':metrics(surface,anchors,'scrsa'),
        'new_eisenberg_field':metrics(new,anchors,'persistence'),
        'ww_interface_sensitivity':metrics(ww,anchors,'persistence'),
        'new_eisenberg_spatial_top10':spatial(surface,anchors,new_top[:10]),
        'orientation_gated_spatial_top10':spatial(surface,anchors,gated[:10]),
        'dipole':r['nonpolar_physics']['hydrophobic_dipole'],
    }


def main():
    cases=[]
    # Development system: successful PE landing/contact region from Wei et al. 2012.
    cases.append(run_case('1AKI','A',7.0,150.0,[67,68,69,70,71,81],'lysozyme_PE_development'))
    # Development systems: neutral hydrophobic Au(111) anchors from Farouq et al.
    for chain in ['B','C']:
        cases.append(run_case('5H7A',chain,7.0,20.0,[33,34,218,220,221],f'SpA_Au111_{chain}_development'))
    report={
        'status':'ok','classification':'development diagnostic; NOT held-out validation',
        'model_changes':'published hydrophobicity scales + hydrophobic dipole only; no label-fitted weights',
        'cases':cases,
    }
    OUT.write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
