from __future__ import annotations

import json, sys, tempfile, urllib.request
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from backend import main as IS
from backend.steering import fibonacci_sphere, debye_length_A
from backend.steering_atomistic import parse_pqr_atoms, compute_pqr_steering

OUT=ROOT/'validation_results'; OUT.mkdir(exist_ok=True)

def fetch(pid):
    with urllib.request.urlopen(f'https://files.rcsb.org/download/{pid}.pdb',timeout=60) as r: return r.read().decode()

def make_pqr(pid,chain,pH):
    req=IS.AnalyzeRequest(pdb_text=fetch(pid),chain=chain,env=IS.EnvParams(pH=pH,ionic=20,temp=298))
    td=tempfile.TemporaryDirectory(); root=Path(td.name); pdb,_=IS.prepare_input_pdb(req,root); pqr=root/'protein.pqr'; IS.run_pdb2pqr(pdb,pqr,pH); txt=pqr.read_text(); td.cleanup(); return txt

def fp_set(st,side,shell):
    return {(r['chain'],int(r['res_seq'])) for r in st['surfaces'][side]['surface_facing_footprint'][f'within_{shell}A']}

def recovery(anchors,fp):
    h=[x for x in anchors if x in fp]; return {'recall':len(h)/len(anchors),'n_hit':len(h),'n':len(anchors),'hits':h}

def orientation_null(pqr_text,ionic,temp,sign,anchors,shell):
    atoms=parse_pqr_atoms(pqr_text); lam=debye_length_A(ionic,temp); normals=fibonacci_sphere(4096)
    xyz=np.vstack([a['coord'] for a in atoms]); q=np.asarray([a['charge_e'] for a in atoms]); rad=np.asarray([a['radius_A'] for a in atoms]);
    proj=xyz@normals.T; plane=np.min(proj-rad[:,None],axis=0); z=np.maximum(proj-plane[None,:],0); E=sign*np.sum(q[:,None]*np.exp(-z/lam),axis=0); bi=int(np.argmin(E))
    heavy=[i for i,a in enumerate(atoms) if not str(a['atom_name']).upper().startswith('H')]
    anchor=set(anchors); hits=[]; fps=[]
    for k in range(len(normals)):
        ids=set();
        for i in heavy:
            if proj[i,k]-rad[i]-plane[k] <= shell:
                ids.add((str(atoms[i]['chain']),int(atoms[i]['res_seq'])))
        hits.append(len(anchor&ids)); fps.append(len(ids))
    obs=hits[bi]; p=(1+sum(h>=obs for h in hits))/(len(hits)+1)
    return {'obs_hits':obs,'n':len(anchor),'recall':obs/len(anchor),'random_mean_hits':float(np.mean(hits)),'empirical_p_hit_ge':p,'footprint_size':fps[bi]}

def main():
    out={'model':'PDB2PQR atom-charge steering diagnostic','cases':[]}
    for chain in ['B','C']:
        pqr=make_pqr('5H7A',chain,7.0); st=compute_pqr_steering(pqr,20,298,4096,(2.,5.,8.))
        pos=[(chain,i) for i in [219,220,221]]; neg=[(chain,i) for i in [33,34,35,36,37]]
        c={'case':f'5H7A_{chain}','net_pqr_charge_e':st['net_pqr_charge_e'],'debye_A':st['debye_length_A'],'positive':{},'negative':{}}
        for sh in [2,5,8]:
            c['positive'][f'R{sh}']={'recovery':recovery(pos,fp_set(st,'positive',sh)),'null':orientation_null(pqr,20,298,+1,pos,sh)}
            c['negative'][f'R{sh}']={'recovery':recovery(neg,fp_set(st,'negative',sh)),'null':orientation_null(pqr,20,298,-1,neg,sh)}
        out['cases'].append(c)
    (OUT/'atomistic_steering_diagnostic.json').write_text(json.dumps(out,indent=2))
    lines=['# Atom-charge steering diagnostic','']
    for c in out['cases']:
        lines.append(f"## {c['case']} (net PQR charge {c['net_pqr_charge_e']})")
        for side in ['positive','negative']:
            for k,v in c[side].items():
                n=v['null']; lines.append(f"- {side} {k}: recall={n['recall']:.3f}, random mean hits={n['random_mean_hits']:.3f}, p={n['empirical_p_hit_ge']:.4f}, footprint n={n['footprint_size']}")
        lines.append('')
    (OUT/'ATOMISTIC_STEERING_DIAGNOSTIC.md').write_text('\n'.join(lines))
if __name__=='__main__': main()
