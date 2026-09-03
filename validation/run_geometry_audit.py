from __future__ import annotations

import csv
import io
import itertools
import json
import math
import random
import sys
import urllib.request
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend import main as IS  # noqa: E402

OUT = ROOT / "validation_results"
OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20260904)

# Adsorption-independent structural panel. These proteins were chosen only to span
# small/medium/large soluble folds and oligomeric states; no protein-surface contact labels
# are used anywhere in this audit.
PANEL = [
    "1UBQ","1CRN","1BTA","1L2Y","1VII","1PGB","1R69","1CSP","2PTL","1TEN",
    "1AKE","1TIM","1LYZ","1HRC","1MBO","1A3N","1FNF","1EMA","1AON","1GFL",
    "1FKJ","1APS","1D3Z","1E0L","1EAZ","1EJG","1HHP","1KTE","1MJC","1N55",
    "1O6X","1PHT","1QYS","1R2R","1SMD","1TIT","1VCC","1WIT","2CI2","2HPR",
]

RADIUS_VALUES = [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0]
RADIUS_PAIRS = list(itertools.combinations(RADIUS_VALUES, 2))
CHEMISTRIES = list(IS.CHEMISTRIES.keys())
PH_VALUES = [5.0, 7.4, 9.0]


def fetch_pdb(pid: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=60) as r:
        return r.read().decode("utf-8")


def atom_sasa_model(pdb_text: str):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", io.StringIO(pdb_text))
    sr = ShrakeRupley(probe_radius=IS.SASA_PROBE_A, n_points=IS.SASA_POINTS)
    sr.compute(s, level="A")
    return s


def residue_table(pdb_text: str):
    s = atom_sasa_model(pdb_text)
    model = next(s.get_models())
    rows = []
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=True) or not res.has_id("CA"):
                continue
            rn = res.get_resname().strip().upper()
            if rn not in IS.SIDECHAIN_REF_ASA:
                continue
            if rn == "GLY":
                sc_atoms = [res["CA"]]
            else:
                backbone = {"N","CA","C","O","OXT"}
                sc_atoms = [a for a in res.get_atoms() if a.get_name().strip().upper() not in backbone and a.element != "H"]
            if not sc_atoms:
                sc_atoms = [res["CA"]]
            sc_sasa = sum(max(float(getattr(a, "sasa", 0.0)), 0.0) for a in sc_atoms)
            scrsa_raw = sc_sasa / IS.SIDECHAIN_REF_ASA[rn]
            if scrsa_raw < IS.SC_RSA_THRESHOLD:
                continue
            scrsa = min(max(float(scrsa_raw), 0.0), 1.0)
            ca = np.asarray(res["CA"].coord, float)
            coords = np.vstack([np.asarray(a.coord, float) for a in sc_atoms])
            sc_centroid = coords.mean(axis=0)
            weights = np.asarray([max(float(getattr(a, "sasa", 0.0)), 0.0) for a in sc_atoms], float)
            if weights.sum() > 1e-12:
                exposed_centroid = np.average(coords, axis=0, weights=weights)
            else:
                exposed_centroid = sc_centroid.copy()
            rows.append({
                "key": f"{chain.id}:{int(res.id[1])}:{str(res.id[2]).strip()}",
                "rn": rn,
                "scrsa": scrsa,
                "ca": ca,
                "sc": sc_centroid,
                "exp": exposed_centroid,
            })
    return rows


def local_scores(rows, chemistry, pH):
    defs = IS.CHEMISTRIES[chemistry]["favorable"]
    out = np.zeros(len(rows), float)
    for i,r in enumerate(rows):
        if r["rn"] in defs:
            _,_,state = defs[r["rn"]]
            out[i] = r["scrsa"] * IS.state_availability(r["rn"], state, pH)
    return out


def persistence(local, coords, pair):
    if len(local) == 0:
        return np.zeros(0)
    dif = coords[:,None,:] - coords[None,:,:]
    dmat = np.sqrt(np.sum(dif*dif, axis=2))
    ns = []
    for R in pair:
        d = (dmat <= R).astype(float) @ local
        mx = float(np.max(d)) if d.size else 0.0
        ns.append(d/mx if mx > 0 else np.zeros_like(d))
    return 100.0*np.minimum(ns[0], ns[1])


def topk(arr, k=10):
    if len(arr) == 0:
        return set()
    return set(np.argsort(-arr)[:min(k,len(arr))].tolist())


def jacc(a,b):
    u = a|b
    return len(a&b)/len(u) if u else 1.0


def safe_spearman(a,b):
    if len(a)<3 or np.allclose(a,a[0]) or np.allclose(b,b[0]):
        return None
    x = spearmanr(a,b).statistic
    return None if x is None or math.isnan(x) else float(x)


def noise_stability(local, coords, pair, sigma, reps=20):
    base = persistence(local, coords, pair)
    js=[]; rs=[]; disp=[]
    if len(base)==0:
        return None,None,None
    btop=topk(base); bi=int(np.argmax(base)); bcoord=coords[bi]
    for _ in range(reps):
        pert = coords + RNG.normal(0.0, sigma, size=coords.shape)
        x = persistence(local, pert, pair)
        js.append(jacc(btop, topk(x)))
        rr=safe_spearman(base,x)
        if rr is not None: rs.append(rr)
        xi=int(np.argmax(x)); disp.append(float(np.linalg.norm(bcoord-coords[xi])))
    return float(np.median(js)), (float(np.median(rs)) if rs else None), float(np.median(disp))


def radius_neighbor_stability(local, coords, pair):
    # Compare a pair with nearby radius-pair choices; no adsorption labels.
    base = persistence(local, coords, pair); bt=topk(base)
    scores=[]
    for alt in RADIUS_PAIRS:
        if alt==pair: continue
        # only local perturbations in radius space
        dist=abs(alt[0]-pair[0])+abs(alt[1]-pair[1])
        if dist <= 3.0:
            scores.append(jacc(bt, topk(persistence(local,coords,alt))))
    return float(np.median(scores)) if scores else None


def patch_coherence(local, coords, p):
    # Geometric compactness of compatible residues around top center at the larger radius.
    if len(p)==0 or np.max(p)<=0:
        return None,None
    i=int(np.argmax(p)); R=8.0
    d=np.linalg.norm(coords-coords[i],axis=1)
    idx=np.where((d<=R)&(local>0))[0]
    if len(idx)<2:
        return float(len(idx)), 0.0
    pts=coords[idx]
    diam=float(np.max(np.linalg.norm(pts[:,None,:]-pts[None,:,:],axis=2)))
    return float(len(idx)), diam


def main():
    raw=[]; failures=[]
    for pid in PANEL:
        try:
            txt=fetch_pdb(pid); rows=residue_table(txt)
            if len(rows)<20:
                failures.append({"pdb":pid,"reason":f"too few exposed residues ({len(rows)})"}); continue
        except Exception as e:
            failures.append({"pdb":pid,"reason":str(e)}); continue
        coords_by_rep={"C_alpha":np.vstack([r["ca"] for r in rows]),
                       "sidechain_centroid":np.vstack([r["sc"] for r in rows]),
                       "SASA_weighted_sidechain":np.vstack([r["exp"] for r in rows])}
        for chem in CHEMISTRIES:
            for pH in PH_VALUES:
                local=local_scores(rows,chem,pH)
                if np.count_nonzero(local)<2:
                    continue
                for rep,coords in coords_by_rep.items():
                    for pair in RADIUS_PAIRS:
                        p=persistence(local,coords,pair)
                        j025,r025,d025=noise_stability(local,coords,pair,0.25,reps=10)
                        j050,r050,d050=noise_stability(local,coords,pair,0.50,reps=10)
                        rad=radius_neighbor_stability(local,coords,pair)
                        nmem,diam=patch_coherence(local,coords,p)
                        raw.append({"pdb":pid,"n_surface":len(rows),"chemistry":chem,"pH":pH,
                                    "representation":rep,"r1":pair[0],"r2":pair[1],
                                    "noise025_top10_jaccard":j025,"noise025_spearman":r025,"noise025_top1_disp_A":d025,
                                    "noise050_top10_jaccard":j050,"noise050_spearman":r050,"noise050_top1_disp_A":d050,
                                    "radius_neighbor_top10_jaccard":rad,"top_patch_compatible_members":nmem,
                                    "top_patch_diameter_A":diam})
    fields=sorted(set().union(*(r.keys() for r in raw))) if raw else []
    with (OUT/"geometry_audit_raw.csv").open("w",newline="") as f:
        if fields:
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(raw)

    # Aggregate robustly over structures/chemistries/pH. Selection score rewards coordinate-noise
    # stability and local radius robustness, while penalizing large top-1 jumps. It deliberately
    # contains no adsorption/anchor information.
    groups={}
    for r in raw:
        key=(r["representation"],r["r1"],r["r2"])
        groups.setdefault(key,[]).append(r)
    summary=[]
    for key,rs in groups.items():
        def med(k):
            vals=[x[k] for x in rs if x[k] is not None]
            return float(np.median(vals)) if vals else None
        j025=med("noise025_top10_jaccard"); j050=med("noise050_top10_jaccard"); rad=med("radius_neighbor_top10_jaccard")
        d025=med("noise025_top1_disp_A"); d050=med("noise050_top1_disp_A")
        # unitless audit score; used only for developmental ranking, not reported as a physical quantity.
        score=(0.30*(j025 or 0)+0.30*(j050 or 0)+0.30*(rad or 0)
               +0.05*(1/(1+(d025 or 999)))+0.05*(1/(1+(d050 or 999))))
        summary.append({"representation":key[0],"r1":key[1],"r2":key[2],"n_conditions":len(rs),
                        "median_noise025_top10_jaccard":j025,"median_noise050_top10_jaccard":j050,
                        "median_radius_neighbor_top10_jaccard":rad,
                        "median_noise025_top1_disp_A":d025,"median_noise050_top1_disp_A":d050,
                        "median_top_patch_members":med("top_patch_compatible_members"),
                        "median_top_patch_diameter_A":med("top_patch_diameter_A"),
                        "audit_rank_score":float(score)})
    summary.sort(key=lambda x:-x["audit_rank_score"])
    fields=sorted(set().union(*(r.keys() for r in summary))) if summary else []
    with (OUT/"geometry_audit_summary.csv").open("w",newline="") as f:
        if fields:
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(summary)

    # Representation-level winner across radius pairs using best developmental pair per representation.
    reps={}
    for r in summary:
        reps.setdefault(r["representation"],r)
    report={"panel_requested":PANEL,"panel_successful":sorted(set(r["pdb"] for r in raw)),"failures":failures,
            "selection_is_adsorption_independent":True,"best_by_representation":reps,
            "overall_best":summary[0] if summary else None,"top15":summary[:15]}
    (OUT/"geometry_audit_report.json").write_text(json.dumps(report,indent=2))
    md=["# Adsorption-independent InterfaceScout geometry audit","",
        "No protein-surface adsorption labels or anchor residues were used in this audit.","",
        f"Successful structural panel: {len(report['panel_successful'])}/{len(PANEL)} proteins.","",
        "## Best radius pair within each representation"]
    for rep,x in reps.items():
        md.append(f"- **{rep}**: {x['r1']}/{x['r2']} Å; audit score {x['audit_rank_score']:.4f}; noise0.25 J={x['median_noise025_top10_jaccard']:.3f}; noise0.50 J={x['median_noise050_top10_jaccard']:.3f}; radius-neighbor J={x['median_radius_neighbor_top10_jaccard']:.3f}")
    if summary:
        x=summary[0]; md += ["","## Overall audit winner",f"**{x['representation']} at {x['r1']}/{x['r2']} Å** (developmental audit score {x['audit_rank_score']:.4f})."]
    (OUT/"GEOMETRY_AUDIT.md").write_text("\n".join(md))

if __name__=="__main__":
    main()
