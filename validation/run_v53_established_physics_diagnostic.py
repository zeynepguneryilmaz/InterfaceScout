"""Development-only diagnostic for v5.3 established nonpolar physics.

IMPORTANT: 1AKI/PE and 5H7A/Au are development systems after prior inspection.
These results must not be presented as held-out validation. The purpose is to
check literature-derived nonpolar physics without fitting weights to labels.
"""
from __future__ import annotations

import json, sys
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
    order=np.argsort(-s, kind='mergesort'); yy=y[order]; tp=0; vals=[]
    for k,v in enumerate(yy,1):
        if v: tp+=1; vals.append(tp/k)
    return float(np.mean(vals))


def metrics(rows, anchors, score_key):
    anchor=set(anchors); y=[1 if int(r['res_seq']) in anchor else 0 for r in rows]
    s=[float(r.get(score_key,0.0) or 0.0) for r in rows]
    order=np.argsort(-np.asarray(s), kind='mergesort'); k=min(10,len(rows)); top=[rows[i] for i in order[:k]]
    return {
        'n':len(rows),'n_anchor_surface':int(sum(y)),'auroc':auc(y,s),'ap':ap(y,s),
        'recall_at_10':float(sum(1 for r in top if int(r['res_seq']) in anchor)/max(sum(y),1)),
        'top10':[[r['res_name'],int(r['res_seq']),round(float(r.get(score_key,0.0) or 0.0),3)] for r in top],
    }


def canonical_patch_metrics(surface, patch_centers, anchors):
    by={(r['chain'],int(r['res_seq']),r.get('icode','')):r for r in patch_centers}
    rows=[]
    for s in surface:
        z=dict(s); p=by.get((s['chain'],int(s['res_seq']),s.get('icode','')),{})
        z['canonical_M']=float(p.get('multiscale_persistence',0.0) or 0.0); rows.append(z)
    return metrics(rows,anchors,'canonical_M')


def spatial(rows, anchors, top_rows):
    coords={int(r['res_seq']):np.asarray([r['x'],r['y'],r['z']],float) for r in rows}
    tc=[np.asarray([r['x'],r['y'],r['z']],float) for r in top_rows]
    dists=[]
    for a in anchors:
        if a in coords and tc: dists.append(min(float(np.linalg.norm(coords[a]-x)) for x in tc))
    out={}
    for R in [5.0,8.0,10.0]: out[f'recall_{int(R)}A']=float(sum(d<=R for d in dists)/len(dists)) if dists else None
    out['median_nearest_A']=float(np.median(dists)) if dists else None; out['n_anchor_surface']=len(dists)
    return out


def orientation_contact_metrics(scan, anchors):
    aset=set(anchors); tops=scan.get('top_orientations',[])
    def one(o):
        rr=o.get('contact_residues',[]); ids={int(r['res_seq']) for r in rr}
        hits=sorted(aset & ids)
        return {'score':o.get('solvation_contact_score'),'n_contact_residues':len(ids),'hits':hits,
                'recall':len(hits)/len(aset) if aset else None,'precision':len(hits)/len(ids) if ids else None}
    top1=one(tops[0]) if tops else None
    best=None
    for o in tops:
        m=one(o)
        if best is None or (m['recall'],m['precision'] or 0)>(best['recall'],best['precision'] or 0): best=m
    union=set()
    for o in tops[:10]: union.update(int(r['res_seq']) for r in o.get('contact_residues',[]))
    return {'top1':top1,'best_among_top20_descriptive':best,
            'union_top10_recall':len(aset&union)/len(aset) if aset else None,
            'method':scan.get('method'),'best_score':scan.get('best_score'),'median_score':scan.get('median_score')}


def run_case(pdb, chain, pH, ionic, anchors, label):
    req=AnalyzeRequest(pdb_id=pdb, chain=chain, env=EnvParams(pH=pH,ionic=ionic,temp=298.0), structure_context='deposited_structure', protrusion=True)
    r=analyze(req); surface=r['surface_residues']
    new=r['nonpolar_physics']['primary_local_field']['residues']; gated=r['nonpolar_physics']['orientation_gated_top_patches']
    new_top=r['nonpolar_physics']['primary_local_field']['top_patches']; ww=r['nonpolar_physics']['independent_interface_scale_sensitivity']['residues']
    atomscan=r['nonpolar_physics']['atom_level_sasa_orientation']
    return {
        'label':label,'pdb':pdb,'chain':chain,'anchors':anchors,
        'canonical_binary_hydrophobic_M':canonical_patch_metrics(surface,r['chemistries']['hydrophobic']['patch_centers'],anchors),
        'scrsa_baseline':metrics(surface,anchors,'scrsa'),
        'new_eisenberg_field':metrics(new,anchors,'persistence'),
        'ww_interface_sensitivity':metrics(ww,anchors,'persistence'),
        'new_eisenberg_spatial_top10':spatial(surface,anchors,new_top[:10]),
        'orientation_gated_spatial_top10':spatial(surface,anchors,gated[:10]),
        'atom_level_sasa_orientation':orientation_contact_metrics(atomscan,anchors),
        'dipole_summary':{k:r['nonpolar_physics']['hydrophobic_dipole'][k] for k in ['definition','vector','magnitude','hemisphere_rule']},
    }


def main():
    cases=[run_case('1AKI','A',7.0,150.0,[67,68,69,70,71,81],'lysozyme_PE_development')]
    for chain in ['B','C']:
        cases.append(run_case('5H7A',chain,7.0,20.0,[33,34,218,220,221],f'SpA_Au111_{chain}_development'))
    report={'status':'ok','classification':'development diagnostic; NOT held-out validation',
            'model_changes':'published hydrophobicity scales + tertiary hydrophobic vector + Harrison-style SASA solvation orientation component; no label-fitted weights',
            'cases':cases}
    OUT.write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2))

if __name__=='__main__': main()
