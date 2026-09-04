from __future__ import annotations

import json
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend import main as IS
from backend.steering import debye_length_A, fibonacci_sphere

OUT = ROOT / "validation_results"
OUT.mkdir(exist_ok=True)


def fetch(pid):
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=60) as r:
        return r.read().decode("utf-8")


def rows(pid, chain, pH):
    req=IS.AnalyzeRequest(pdb_text=fetch(pid),chain=chain,env=IS.EnvParams(pH=pH,ionic=20,temp=298))
    with tempfile.TemporaryDirectory() as td:
        pdb,_=IS.prepare_input_pdb(req,Path(td)); _,allr,_,_=IS.build_surface_residues(pdb,pH)
    return allr,[r for r in allr if r["surface_exposed"]]


def coords(rs): return np.asarray([[r["x"],r["y"],r["z"]] for r in rs],float)
def charges(rs): return np.asarray([r["charge_descriptor"] for r in rs],float)


def evaluate(pid,chain,pH,ionic,temp,surface_sign,anchors,depth):
    allr,surf=rows(pid,chain,pH)
    ac=coords(allr); sc=coords(surf); q=charges(allr); mask=np.abs(q)>1e-8; ac=ac[mask];q=q[mask]
    normals=fibonacci_sphere(4096); lam=debye_length_A(ionic,temp)
    sp=sc@normals.T; plane=sp.min(axis=0)
    z=np.maximum(ac@normals.T-plane[None,:],0); energy=surface_sign*np.sum(q[:,None]*np.exp(-z/lam),axis=0)
    best=int(np.argmin(energy))
    anchor_set=set(anchors); surf_ids=[(r["chain"],int(r["res_seq"])) for r in surf]
    recalls=[]; hits=[]; footprint_sizes=[]
    for k in range(len(normals)):
        dep=np.maximum(sp[:,k]-plane[k],0); ids={surf_ids[i] for i in np.where(dep<=depth)[0]}
        h=len(anchor_set&ids); hits.append(h); recalls.append(h/len(anchor_set)); footprint_sizes.append(len(ids))
    obs=recalls[best]; obs_hits=hits[best]
    p=(1+sum(h>=obs_hits for h in hits))/(len(hits)+1)
    return {"case":f"{pid}_{chain}_{'pos' if surface_sign>0 else 'neg'}_R{depth}","depth_A":depth,
            "observed_recall":obs,"observed_hits":obs_hits,"n_anchors":len(anchor_set),
            "observed_footprint_size":footprint_sizes[best],"random_orientation_mean_recall":float(np.mean(recalls)),
            "random_orientation_median_recall":float(np.median(recalls)),"empirical_p_hit_ge":p,
            "best_energy_percentile":float(np.mean(energy>=energy[best])),"debye_A":lam}


def main():
    out=[]
    for chain in ["B","C"]:
        pos=[(chain,219),(chain,220),(chain,221)]; neg=[(chain,i) for i in [33,34,35,36,37]]
        for d in [5,8,10]:
            out.append(evaluate("5H7A",chain,7.0,20.0,298.0,+1,pos,d))
            out.append(evaluate("5H7A",chain,7.0,20.0,298.0,-1,neg,d))
    (OUT/"steering_orientation_null.json").write_text(json.dumps(out,indent=2))
    lines=["# Steering orientation-null diagnostic",""]
    for r in out:
        lines.append(f"- {r['case']}: recall {r['observed_recall']:.3f}, random mean {r['random_orientation_mean_recall']:.3f}, p(hit>=obs)={r['empirical_p_hit_ge']:.4f}, footprint n={r['observed_footprint_size']}")
    (OUT/"STEERING_ORIENTATION_NULL.md").write_text("\n".join(lines))

if __name__=='__main__': main()
