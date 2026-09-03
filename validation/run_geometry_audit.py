from __future__ import annotations

import csv
import io
import itertools
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend import main as IS  # noqa: E402

OUT = ROOT / "validation_results"
OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20260904)

# Adsorption-independent structural-development panel. None of the published
# protein-surface anchor/contact benchmark systems is included.
PANEL = [
    "1UBQ","1CRN","1BTA","1VII","1PGB","1R69","1CSP","2PTL","1TEN","1AKE",
    "1LYZ","1HRC","1MBO","1EMA","1FKJ","1D3Z","1TIT","1WIT","2CI2","2HPR",
]
RADIUS_VALUES = [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0]
RADIUS_PAIRS = list(itertools.combinations(RADIUS_VALUES, 2))
CHEMISTRIES = list(IS.CHEMISTRIES.keys())
PH_VALUES = [5.0, 7.4, 9.0]
NOISE_SIGMAS = [0.25, 0.50]
NOISE_REPS = 3


def fetch_pdb(pid):
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=60) as r:
        return r.read().decode("utf-8")


def residue_table(pdb_text):
    s = PDBParser(QUIET=True).get_structure("x", io.StringIO(pdb_text))
    ShrakeRupley(probe_radius=IS.SASA_PROBE_A, n_points=IS.SASA_POINTS).compute(s, level="A")
    model = next(s.get_models())
    chosen = next((c for c in model if any(is_aa(r, standard=True) and r.has_id("CA") for r in c)), None)
    if chosen is None: return []
    rows=[]
    for res in chosen:
        if not is_aa(res, standard=True) or not res.has_id("CA"): continue
        rn=res.get_resname().strip().upper()
        if rn not in IS.SIDECHAIN_REF_ASA: continue
        if rn == "GLY": sc_atoms=[res["CA"]]
        else:
            backbone={"N","CA","C","O","OXT"}
            sc_atoms=[a for a in res.get_atoms() if a.get_name().strip().upper() not in backbone and a.element != "H"] or [res["CA"]]
        sasa=np.asarray([max(float(getattr(a,"sasa",0.0)),0.0) for a in sc_atoms])
        scrsa_raw=float(sasa.sum())/IS.SIDECHAIN_REF_ASA[rn]
        if scrsa_raw < IS.SC_RSA_THRESHOLD: continue
        coords=np.vstack([np.asarray(a.coord,float) for a in sc_atoms])
        sc=coords.mean(axis=0)
        exp=np.average(coords,axis=0,weights=sasa) if sasa.sum()>1e-12 else sc.copy()
        rows.append({"rn":rn,"scrsa":min(max(scrsa_raw,0.0),1.0),"ca":np.asarray(res["CA"].coord,float),"sc":sc,"exp":exp})
    return rows


def dmat(coords):
    x=coords[:,None,:]-coords[None,:,:]
    return np.sqrt(np.sum(x*x,axis=2))


def radius_masks(dm):
    return {R:(dm<=R) for R in RADIUS_VALUES}


def normalized_radius_maps(local,masks):
    out={}
    for R,m in masks.items():
        x=m@local
        mx=float(np.max(x)) if x.size else 0.0
        out[R]=x/mx if mx>0 else np.zeros_like(x)
    return out


def pair_map(norms,pair):
    return 100.0*np.minimum(norms[pair[0]],norms[pair[1]])


def local_scores(rows,chem,pH):
    defs=IS.CHEMISTRIES[chem]["favorable"]
    v=np.zeros(len(rows),float)
    for i,r in enumerate(rows):
        if r["rn"] in defs:
            _,_,state=defs[r["rn"]]
            v[i]=r["scrsa"]*IS.state_availability(r["rn"],state,pH)
    return v


def topk(x,k=10):
    return set(np.argsort(-x)[:min(k,len(x))].tolist()) if len(x) else set()


def jac(a,b):
    u=a|b
    return len(a&b)/len(u) if u else 1.0


def main():
    raw=[]; failures=[]
    for pid in PANEL:
        try:
            rows=residue_table(fetch_pdb(pid))
            if len(rows)<20: raise ValueError(f"too few exposed residues ({len(rows)})")
        except Exception as e:
            failures.append({"pdb":pid,"reason":str(e)}); continue

        reps={"C_alpha":np.vstack([r["ca"] for r in rows]),
              "sidechain_centroid":np.vstack([r["sc"] for r in rows]),
              "SASA_weighted_sidechain":np.vstack([r["exp"] for r in rows])}
        geom={}
        for rep,coords in reps.items():
            base_dm=dmat(coords)
            base_masks=radius_masks(base_dm)
            noise_masks={sig:[radius_masks(dmat(coords+RNG.normal(0,sig,size=coords.shape))) for _ in range(NOISE_REPS)] for sig in NOISE_SIGMAS}
            geom[rep]=(coords,base_dm,base_masks,noise_masks)

        for chem in CHEMISTRIES:
            for pH in PH_VALUES:
                local=local_scores(rows,chem,pH)
                if np.count_nonzero(local)<2: continue
                for rep,(coords,base_dm,base_masks,noise_masks) in geom.items():
                    base_norms=normalized_radius_maps(local,base_masks)
                    base_maps={pair:pair_map(base_norms,pair) for pair in RADIUS_PAIRS}
                    tops={pair:topk(x) for pair,x in base_maps.items()}
                    noise_pair_maps={}
                    for sig in NOISE_SIGMAS:
                        noise_pair_maps[sig]=[]
                        for masks in noise_masks[sig]:
                            nm=normalized_radius_maps(local,masks)
                            noise_pair_maps[sig].append({pair:pair_map(nm,pair) for pair in RADIUS_PAIRS})

                    for pair in RADIUS_PAIRS:
                        base=base_maps[pair]; bt=tops[pair]
                        bi=int(np.argmax(base)); bcoord=coords[bi]
                        vals={}
                        for sig in NOISE_SIGMAS:
                            js=[]; dis=[]
                            for maps in noise_pair_maps[sig]:
                                x=maps[pair]
                                js.append(jac(bt,topk(x)))
                                xi=int(np.argmax(x)); dis.append(float(np.linalg.norm(bcoord-coords[xi])))
                            vals[sig]=(float(np.median(js)),float(np.median(dis)))
                        nearby=[]
                        for alt in RADIUS_PAIRS:
                            if alt==pair: continue
                            if abs(alt[0]-pair[0])+abs(alt[1]-pair[1])<=3.0:
                                nearby.append(jac(bt,tops[alt]))
                        rad=float(np.median(nearby)) if nearby else None
                        idx=np.where((base_dm[bi]<=pair[1])&(local>0))[0]
                        diam=0.0 if len(idx)<2 else float(np.max(base_dm[np.ix_(idx,idx)]))
                        raw.append({"pdb":pid,"n_surface":len(rows),"chemistry":chem,"pH":pH,"representation":rep,
                                    "r1":pair[0],"r2":pair[1],
                                    "noise025_top10_jaccard":vals[0.25][0],"noise025_top1_disp_A":vals[0.25][1],
                                    "noise050_top10_jaccard":vals[0.50][0],"noise050_top1_disp_A":vals[0.50][1],
                                    "radius_neighbor_top10_jaccard":rad,"top_patch_compatible_members":float(len(idx)),"top_patch_diameter_A":diam})

    fields=sorted(set().union(*(r.keys() for r in raw))) if raw else []
    with (OUT/"geometry_audit_raw.csv").open("w",newline="") as f:
        if fields:
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(raw)

    groups={}
    for r in raw: groups.setdefault((r["representation"],r["r1"],r["r2"]),[]).append(r)
    summary=[]
    for key,rs in groups.items():
        def med(k):
            v=[x[k] for x in rs if x[k] is not None]
            return float(np.median(v)) if v else None
        j25,j50,rad=med("noise025_top10_jaccard"),med("noise050_top10_jaccard"),med("radius_neighbor_top10_jaccard")
        d25,d50=med("noise025_top1_disp_A"),med("noise050_top1_disp_A")
        score=.30*(j25 or 0)+.30*(j50 or 0)+.30*(rad or 0)+.05/(1+(d25 if d25 is not None else 999))+.05/(1+(d50 if d50 is not None else 999))
        summary.append({"representation":key[0],"r1":key[1],"r2":key[2],"n_conditions":len(rs),
                        "median_noise025_top10_jaccard":j25,"median_noise050_top10_jaccard":j50,"median_radius_neighbor_top10_jaccard":rad,
                        "median_noise025_top1_disp_A":d25,"median_noise050_top1_disp_A":d50,
                        "median_top_patch_members":med("top_patch_compatible_members"),"median_top_patch_diameter_A":med("top_patch_diameter_A"),
                        "audit_rank_score":float(score)})
    summary.sort(key=lambda x:-x["audit_rank_score"])
    fields=sorted(set().union(*(r.keys() for r in summary))) if summary else []
    with (OUT/"geometry_audit_summary.csv").open("w",newline="") as f:
        if fields:
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(summary)

    best={}
    for x in summary: best.setdefault(x["representation"],x)
    report={"panel_requested":PANEL,"panel_successful":sorted(set(r["pdb"] for r in raw)),"failures":failures,
            "selection_is_adsorption_independent":True,"noise_reps":NOISE_REPS,
            "best_by_representation":best,"overall_best":summary[0] if summary else None,"top15":summary[:15]}
    (OUT/"geometry_audit_report.json").write_text(json.dumps(report,indent=2))
    lines=["# Adsorption-independent InterfaceScout geometry audit","",
           "No adsorption labels, anchor residues, MD contacts, or adsorption outcomes were used.","",
           f"Successful structural panel: {len(report['panel_successful'])}/{len(PANEL)} proteins.","",
           "## Best pair by representation"]
    for rep,x in best.items():
        lines.append(f"- **{rep}**: {x['r1']}/{x['r2']} Å; audit={x['audit_rank_score']:.4f}; J0.25={x['median_noise025_top10_jaccard']:.3f}; J0.50={x['median_noise050_top10_jaccard']:.3f}; radius-J={x['median_radius_neighbor_top10_jaccard']:.3f}")
    if summary:
        x=summary[0]; lines += ["","## Overall winner",f"**{x['representation']} at {x['r1']}/{x['r2']} Å** (audit score {x['audit_rank_score']:.4f})."]
    (OUT/"GEOMETRY_AUDIT.md").write_text("\n".join(lines))

if __name__=="__main__": main()
