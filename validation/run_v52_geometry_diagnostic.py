from __future__ import annotations

import csv
import io
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from backend import main as IS  # noqa: E402

OUT=ROOT/'validation_results'; OUT.mkdir(exist_ok=True)
PAIR=(5.0,8.0)


def fetch(pid):
    with urllib.request.urlopen(f'https://files.rcsb.org/download/{pid}.pdb',timeout=60) as r:
        return r.read().decode('utf-8')


def run(txt,chain,pH,temp=298,ionic=150):
    return IS.analyze(IS.AnalyzeRequest(pdb_text=txt,chain=chain,env=IS.EnvParams(pH=pH,temp=temp,ionic=ionic)))


def local_vector(res,chem):
    members={x['key']:float(x['local_score']) for x in res['chemistries'][chem]['residues']}
    surf=res['surface_residues']
    return np.asarray([members.get(x['key'],0.0) for x in surf],float)


def persistence(res,chem,mode):
    surf=res['surface_residues']; local=local_vector(res,chem)
    if mode=='C_alpha': coords=np.asarray([[x['x'],x['y'],x['z']] for x in surf],float)
    elif mode=='sidechain_centroid': coords=np.asarray([[x['sc_x'],x['sc_y'],x['sc_z']] for x in surf],float)
    else: raise ValueError(mode)
    d=np.linalg.norm(coords[:,None,:]-coords[None,:,:],axis=2)
    norms=[]
    for R in PAIR:
        z=(d<=R).astype(float)@local
        mx=float(z.max()) if z.size else 0.0
        norms.append(z/mx if mx>0 else np.zeros_like(z))
    return 100*np.minimum(norms[0],norms[1]),coords


def exact_metrics(res,chem,anchors,mode):
    p,coords=persistence(res,chem,mode); surf=res['surface_residues']
    y=np.asarray([1 if int(x['res_seq']) in set(anchors) else 0 for x in surf],int)
    eligible=[int(x['res_seq']) for x in surf if int(x['res_seq']) in set(anchors)]
    auc=float(roc_auc_score(y,p)) if len(set(y.tolist()))==2 else None
    ap=float(average_precision_score(y,p)) if y.sum() else None
    top=np.argsort(-p)[:min(10,len(p))]
    coord_by_num={int(x['res_seq']):coords[i] for i,x in enumerate(surf)}
    ds=[]
    for a in anchors:
        if a not in coord_by_num: continue
        ds.append(min(float(np.linalg.norm(coord_by_num[a]-coords[j])) for j in top))
    return {'AUROC':auc,'AP':ap,'R5':float(np.mean(np.asarray(ds)<=5)) if ds else None,
            'R8':float(np.mean(np.asarray(ds)<=8)) if ds else None,'median_nearest_A':float(np.median(ds)) if ds else None,
            'n_anchor_surface':len(eligible),'top10_keys':[surf[j]['key'] for j in top]}


def region_metrics(res,chem,nums,mode):
    p,_=persistence(res,chem,mode); surf=res['surface_residues']; S=set(nums)
    y=np.asarray([1 if int(x['res_seq']) in S else 0 for x in surf],int)
    return {'AUROC':float(roc_auc_score(y,p)) if len(set(y.tolist()))==2 else None,
            'AP':float(average_precision_score(y,p)) if y.sum() else None,'n_region_surface':int(y.sum())}


def main():
    pdb={x:fetch(x) for x in ['1MBN','2PTN','2HHB','5H7A']}
    rows=[]
    # SpA exact-anchor diagnostics. Chains B/C both contain the reported numbering.
    spa=[('Au111','hydrophobic',[221,220,218,33,34]),('O_rich_silica','anionic',[33,34,35,36,37]),('Si_rich_silica','cationic',[221,220,219])]
    for chain in ['B','C']:
        res=run(pdb['5H7A'],chain,7.0,298,20)
        for label,chem,anchors in spa:
            for mode in ['C_alpha','sidechain_centroid']:
                rows.append({'family':'SpA','case':label,'chain':chain,'chemistry':chem,'geometry':mode,**exact_metrics(res,chem,anchors,mode)})
    # Tavanti regions are diagnostic only; no strict residue-level ground truth implied.
    tav=[('1MBN','A',[43,45]+list(range(96,100))+list(range(146,154))),
         ('2PTN','A',[94]+list(range(125,136))+[166,167]+list(range(231,245))),
         ('2HHB','A',list(range(12,26))+list(range(61,82))),
         ('2HHB','C',list(range(12,26))+list(range(61,82))),
         ('2HHB','B',list(range(51,54))),('2HHB','D',list(range(45,54)))]
    for pid,chain,nums in tav:
        res=run(pdb[pid],chain,7.4,310,150)
        for chem in ['anionic','hbond_acceptor','hydrophobic']:
            for mode in ['C_alpha','sidechain_centroid']:
                rows.append({'family':'Tavanti_region','case':pid,'chain':chain,'chemistry':chem,'geometry':mode,**region_metrics(res,chem,nums,mode)})
    fields=sorted(set().union(*(r.keys() for r in rows)))
    with (OUT/'v52_geometry_diagnostic.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    # paired summaries: SC minus CA for matching cases
    pairs=[]
    by={}
    for r in rows:
        key=(r['family'],r['case'],r['chain'],r['chemistry']); by.setdefault(key,{})[r['geometry']]=r
    for key,d in by.items():
        if 'C_alpha' not in d or 'sidechain_centroid' not in d: continue
        a=d['C_alpha']; b=d['sidechain_centroid']
        pairs.append({'family':key[0],'case':key[1],'chain':key[2],'chemistry':key[3],
                      'delta_AUROC_SC_minus_CA':None if a.get('AUROC') is None or b.get('AUROC') is None else b['AUROC']-a['AUROC'],
                      'delta_AP_SC_minus_CA':None if a.get('AP') is None or b.get('AP') is None else b['AP']-a['AP'],
                      'delta_R8_SC_minus_CA':None if a.get('R8') is None or b.get('R8') is None else b['R8']-a['R8'],
                      'delta_median_distance_SC_minus_CA':None if a.get('median_nearest_A') is None or b.get('median_nearest_A') is None else b['median_nearest_A']-a['median_nearest_A']})
    with (OUT/'v52_geometry_paired_changes.csv').open('w',newline='') as f:
        fields=sorted(set().union(*(r.keys() for r in pairs))); w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(pairs)
    report={'candidate_geometry':'sidechain_centroid','candidate_radii_A':[5.0,8.0],
            'diagnostic_only':True,'note':'These systems were previously inspected and are not eligible as final held-out validation.',
            'paired_median_delta_AUROC':float(np.median([x['delta_AUROC_SC_minus_CA'] for x in pairs if x['delta_AUROC_SC_minus_CA'] is not None])),
            'paired_median_delta_AP':float(np.median([x['delta_AP_SC_minus_CA'] for x in pairs if x['delta_AP_SC_minus_CA'] is not None])),
            'SpA_rows':[r for r in rows if r['family']=='SpA']}
    (OUT/'v52_geometry_diagnostic_report.json').write_text(json.dumps(report,indent=2))

if __name__=='__main__': main()
