from __future__ import annotations

import copy
import csv
import io
import json
import math
import random
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from Bio.PDB import PDBParser
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from backend import main as IS  # noqa: E402

OUT = ROOT / "validation_results"
OUT.mkdir(exist_ok=True)
RNG = random.Random(20260903)


def fetch_pdb_text(pdb_id: str) -> str:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8")


def chain_inventory(pdb_text: str) -> Dict[str, dict]:
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("x", io.StringIO(pdb_text))
    out = {}
    model = next(s.get_models())
    for chain in model:
        ids = []
        names = {}
        for res in chain:
            if not IS.is_aa(res, standard=True) or not res.has_id("CA"):
                continue
            seq = int(res.id[1]); ic = str(res.id[2]).strip()
            ids.append((seq, ic)); names[(seq, ic)] = res.get_resname().strip().upper()
        out[str(chain.id)] = {"n": len(ids), "ids": ids, "names": names,
                              "min_resseq": min([x[0] for x in ids]) if ids else None,
                              "max_resseq": max([x[0] for x in ids]) if ids else None}
    return out


def run_is(pdb_text: str, chain: str | None, pH: float, temp: float = 298.0, ionic: float = 150.0):
    req = IS.AnalyzeRequest(pdb_text=pdb_text, chain=chain,
                            env=IS.EnvParams(pH=pH, temp=temp, ionic=ionic))
    return IS.analyze(req)


def score_vectors(result: dict, chemistry: str):
    surf = result["surface_residues"]
    keys = [r["key"] for r in surf]
    scrsa = np.array([float(r["scrsa"]) for r in surf], dtype=float)
    chem = result["chemistries"][chemistry]
    member = {r["key"]: r for r in chem["residues"]}
    m1 = np.array([float(r["scrsa"]) if r["key"] in member else 0.0 for r in surf])
    m2 = np.array([float(member[r["key"]]["local_score"]) if r["key"] in member else 0.0 for r in surf])
    pmap = {r["center_key"]: float(r["multiscale_persistence"]) for r in chem["patch_centers"]}
    m3 = np.array([pmap.get(r["key"], 0.0) for r in surf])
    coords = np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in surf])
    sccoords = np.array([[float(r["sc_x"]), float(r["sc_y"]), float(r["sc_z"])] for r in surf])
    return keys, scrsa, m1, m2, m3, coords, sccoords


def key_nums(result: dict):
    return {(r["chain"], int(r["res_seq"]), str(r["icode"]).strip()): r["key"] for r in result["surface_residues"]}


def exact_labels(result: dict, anchor_numbers: List[int]):
    pos = set(anchor_numbers)
    y = np.array([1 if int(r["res_seq"]) in pos else 0 for r in result["surface_residues"]], dtype=int)
    eligible = sorted(set(int(r["res_seq"]) for r in result["surface_residues"]) & pos)
    missing = sorted(pos - set(eligible))
    return y, eligible, missing


def safe_auc(y, s):
    return float(roc_auc_score(y, s)) if len(set(y.tolist())) == 2 else None


def safe_ap(y, s):
    return float(average_precision_score(y, s)) if np.sum(y) > 0 else None


def top_centers(result: dict, chemistry: str, k: int = 10):
    return result["chemistries"][chemistry]["patch_centers"][:k]


def spatial_recovery(result: dict, chemistry: str, anchor_numbers: List[int], k: int = 10):
    surf = result["surface_residues"]
    coord_by_num = {}
    for r in surf:
        coord_by_num.setdefault(int(r["res_seq"]), np.array([r["x"], r["y"], r["z"]], float))
    center_keys = [x["center_key"] for x in top_centers(result, chemistry, k)]
    coord_by_key = {r["key"]: np.array([r["x"], r["y"], r["z"]], float) for r in surf}
    centers = [coord_by_key[x] for x in center_keys if x in coord_by_key]
    ds = []
    for a in anchor_numbers:
        if a not in coord_by_num or not centers:
            continue
        ds.append(min(float(np.linalg.norm(coord_by_num[a] - c)) for c in centers))
    return {
        "n_eligible": len(ds),
        "R5": float(np.mean(np.array(ds) <= 5.0)) if ds else None,
        "R8": float(np.mean(np.array(ds) <= 8.0)) if ds else None,
        "median_nearest_A": float(np.median(ds)) if ds else None,
        "distances_A": ds,
        "top_centers": center_keys,
    }


def region_enrichment(result: dict, chemistry: str, regions: Dict[str, List[int]], label: str):
    surf = result["surface_residues"]
    _, scrsa, m1, m2, m3, _, _ = score_vectors(result, chemistry)
    rows = []
    all_nums = np.array([int(r["res_seq"]) for r in surf], dtype=int)
    for region_name, nums in regions.items():
        mask = np.isin(all_nums, np.array(nums, dtype=int))
        for name, s in [("scRSA", scrsa), ("chemistry_exposure", m1), ("state_score", m2), ("persistence", m3)]:
            obs = float(np.median(s[mask])) if np.any(mask) else None
            # region-size-matched random exposed-residue null
            vals = []
            m = int(np.sum(mask))
            if m > 0:
                inds = list(range(len(s)))
                for _ in range(10000):
                    sel = RNG.sample(inds, m)
                    vals.append(float(np.median(s[sel])))
                p = (1 + sum(v >= obs for v in vals)) / (len(vals) + 1)
                pct = float(np.mean(np.asarray(vals) <= obs))
            else:
                p = None; pct = None
            rows.append({"case": label, "region": region_name, "metric": name,
                         "n_region_surface": m, "observed_median": obs,
                         "random_empirical_p_ge": p, "random_percentile": pct})
    return rows


def exposure_matched_permutation(y, score, scrsa, nperm=10000, tol=0.05):
    pos_idx = np.where(y == 1)[0].tolist()
    if not pos_idx:
        return None
    obs = float(np.mean(score[pos_idx]))
    all_idx = list(range(len(y)))
    neg_idx = np.where(y == 0)[0].tolist()
    null = []
    for _ in range(nperm):
        chosen = []
        used = set()
        for i in pos_idx:
            cand = [j for j in neg_idx if j not in used and abs(float(scrsa[j]) - float(scrsa[i])) <= tol]
            if not cand:
                cand = [j for j in neg_idx if j not in used]
                cand.sort(key=lambda j: abs(float(scrsa[j]) - float(scrsa[i])))
                cand = cand[:max(1, min(10, len(cand)))]
            j = RNG.choice(cand); chosen.append(j); used.add(j)
        null.append(float(np.mean(score[chosen])))
    p = (1 + sum(v >= obs for v in null)) / (len(null) + 1)
    return {"observed_mean": obs, "null_mean": float(np.mean(null)), "null_sd": float(np.std(null, ddof=1)),
            "empirical_p_ge": p, "percentile": float(np.mean(np.asarray(null) <= obs))}


def recompute_persistence(local: np.ndarray, coords: np.ndarray, radii=(5.0, 8.0)):
    dif = coords[:, None, :] - coords[None, :, :]
    dmat = np.sqrt(np.sum(dif * dif, axis=2))
    norms = []
    for R in radii:
        d = (dmat <= R).astype(float) @ local
        mx = float(np.max(d)) if len(d) else 0.0
        norms.append(d / mx if mx > 0 else np.zeros_like(d))
    return 100.0 * np.minimum(norms[0], norms[1])


def geometry_sensitivity(result: dict, chemistry: str):
    _, _, _, local, canonical, ca, sc = score_vectors(result, chemistry)
    alt = recompute_persistence(local, sc, (5.0, 8.0))
    rho = spearmanr(canonical, alt).statistic if len(canonical) > 2 else None
    top_ca = set(np.argsort(-canonical)[:10].tolist())
    top_sc = set(np.argsort(-alt)[:10].tolist())
    jac = len(top_ca & top_sc) / len(top_ca | top_sc) if (top_ca | top_sc) else 1.0
    i0 = int(np.argmax(canonical)); i1 = int(np.argmax(alt))
    disp = float(np.linalg.norm(ca[i0] - ca[i1]))
    return {"spearman_Calpha_vs_sidechain_centroid": None if rho is None or math.isnan(rho) else float(rho),
            "top10_jaccard": float(jac), "top1_center_displacement_A": disp}


def mathematical_checks():
    checks = []
    ph_grid = np.linspace(0, 14, 141)
    vals = []
    for rn in IS.PKA:
        for pH in ph_grid:
            vals.append(IS.charged_fraction(rn, float(pH)))
            vals.append(IS.state_availability(rn, "protonated", float(pH)))
            vals.append(IS.state_availability(rn, "deprotonated", float(pH)))
    checks.append({"check": "Henderson-Hasselbalch fractions bounded [0,1]", "pass": bool(min(vals) >= 0 and max(vals) <= 1)})
    z = IS.normalize_to_100(np.zeros(5))
    checks.append({"check": "zero-map normalization finite and zero", "pass": bool(np.all(np.isfinite(z)) and np.all(z == 0))})
    checks.append({"check": "canonical radii frozen at 5/8 A", "pass": tuple(IS.PATCH_RADII_A) == (5.0, 8.0)})
    checks.append({"check": "scRSA threshold frozen at 0.05", "pass": float(IS.SC_RSA_THRESHOLD) == 0.05})
    checks.append({"check": "SASA points frozen at 200", "pass": int(IS.SASA_POINTS) == 200})
    return checks


def ebase_apbs_invariance(result: dict, chemistry: str):
    # Rebuild chemistry map from already computed surface geometry. First mutate only Ebase metadata.
    surf = copy.deepcopy(result["surface_residues"])
    dmat = IS.build_distances(surf)
    base = IS.chemistry_map(surf, dmat, chemistry, result["settings"]["pH"])
    base_vec = {x["center_key"]: x["multiscale_persistence"] for x in base["patch_centers"]}
    original = copy.deepcopy(IS.CHEMISTRIES[chemistry])
    try:
        for channel in ("favorable", "repulsive"):
            for rn, (e, mech, state) in list(IS.CHEMISTRIES[chemistry][channel].items()):
                IS.CHEMISTRIES[chemistry][channel][rn] = (e * (3.17 if e != 0 else 1.0), mech, state)
        mutated = IS.chemistry_map(surf, dmat, chemistry, result["settings"]["pH"])
        mut_vec = {x["center_key"]: x["multiscale_persistence"] for x in mutated["patch_centers"]}
    finally:
        IS.CHEMISTRIES[chemistry] = original
    ebase_same = base_vec == mut_vec
    # APBS: overwrite phi only; canonical map must be invariant.
    for i, r in enumerate(surf):
        r["phi"] = float((i % 17) - 8) * 9.37
    phi_mut = IS.chemistry_map(surf, dmat, chemistry, result["settings"]["pH"])
    phi_vec = {x["center_key"]: x["multiscale_persistence"] for x in phi_mut["patch_centers"]}
    return {"Ebase_rescaling_exact_invariance": bool(ebase_same), "APBS_phi_exact_invariance": bool(base_vec == phi_vec)}


def exact_case(label, result, chemistry, anchors):
    keys, scrsa, m1, m2, m3, _, _ = score_vectors(result, chemistry)
    y, eligible, missing = exact_labels(result, anchors)
    rows = []
    for name, s in [("scRSA", scrsa), ("chemistry_exposure", m1), ("state_score", m2), ("persistence", m3)]:
        pm = exposure_matched_permutation(y, s, scrsa) if np.sum(y) else None
        rows.append({"case": label, "chemistry": chemistry, "score": name,
                     "n_surface": len(y), "n_anchor_total": len(anchors), "n_anchor_surface": int(np.sum(y)),
                     "anchors_surface": eligible, "anchors_below_surface_threshold": missing,
                     "AUROC": safe_auc(y, s), "AP": safe_ap(y, s),
                     "exposure_matched_null": pm})
    sp = spatial_recovery(result, chemistry, anchors, k=10)
    return rows, sp


def main():
    summary = {"model_version": IS.APP_VERSION, "mathematical_checks": mathematical_checks(),
               "pdb_inventory": {}, "exact_anchor_metrics": [], "spatial": [], "region_enrichment": [],
               "geometry_sensitivity": [], "invariance": [], "notes": []}

    pdbs = {pid: fetch_pdb_text(pid) for pid in ["1MBN", "2HHB", "2PTN", "5H7A"]}
    for pid, txt in pdbs.items():
        inv = chain_inventory(txt)
        summary["pdb_inventory"][pid] = {c: {"n": v["n"], "min_resseq": v["min_resseq"], "max_resseq": v["max_resseq"]} for c, v in inv.items()}

    # Tavanti 2019: binding REGIONS, not treated as strict binary anchors because exact persistent residues are graphically encoded.
    tavanti = [
        ("1MBN_citrate_AuNP", "1MBN", "A", {"MB_regions": [43,45] + list(range(96,100)) + list(range(146,154))}),
        ("2PTN_citrate_AuNP", "2PTN", "A", {"TRP_regions": [94] + list(range(125,136)) + [166,167] + list(range(231,245))}),
        ("2HHB_A_citrate_AuNP", "2HHB", "A", {"HB_A_regions": list(range(12,26)) + list(range(61,82))}),
        ("2HHB_C_citrate_AuNP", "2HHB", "C", {"HB_C_regions": list(range(12,26)) + list(range(61,82))}),
        ("2HHB_B_citrate_AuNP", "2HHB", "B", {"HB_B_region": list(range(51,54))}),
        ("2HHB_D_citrate_AuNP", "2HHB", "D", {"HB_D_region": list(range(45,54))}),
    ]
    for label, pid, chain, regions in tavanti:
        inv = summary["pdb_inventory"][pid]
        if chain not in inv:
            summary["notes"].append(f"{label}: requested chain {chain} absent; skipped")
            continue
        res = run_is(pdbs[pid], chain, pH=7.4, temp=310.0, ionic=150.0)
        # Predeclared citrate mechanistic panel; no cross-chemistry score combination.
        for chemistry in ["anionic", "hydrophobic", "hbond_acceptor"]:
            summary["region_enrichment"].extend(region_enrichment(res, chemistry, regions, label))
            summary["geometry_sensitivity"].append({"case": label, "chemistry": chemistry, **geometry_sensitivity(res, chemistry)})
        summary["invariance"].append({"case": label, "chemistry": "anionic", **ebase_apbs_invariance(res, "anionic")})

    # SpA 5H7A: exact Table-3 anchors. Paper does not state which crystallographic copy was extracted.
    # Run every chain that contains all anchor residue numbers and summarize chain-copy sensitivity instead of choosing post hoc.
    spa_specs = [
        ("SpA_Au111", "hydrophobic", [221,220,218,33,34]),
        ("SpA_O_rich_silica", "anionic", [33,34,35,36,37]),
        ("SpA_Si_rich_silica", "cationic", [221,220,219]),
    ]
    spa_inv = chain_inventory(pdbs["5H7A"])
    needed = set([33,34,35,36,37,218,219,220,221])
    candidate_chains = []
    for c, meta in spa_inv.items():
        nums = set(x[0] for x in meta["ids"])
        if needed.issubset(nums):
            candidate_chains.append(c)
    summary["notes"].append(f"5H7A candidate crystallographic chains containing all benchmark residue numbers: {candidate_chains}")
    for chain in candidate_chains:
        res = run_is(pdbs["5H7A"], chain, pH=7.0, temp=298.0, ionic=20.0)
        for label, chemistry, anchors in spa_specs:
            case = f"{label}_5H7A_chain_{chain}"
            rows, sp = exact_case(case, res, chemistry, anchors)
            summary["exact_anchor_metrics"].extend(rows)
            summary["spatial"].append({"case": case, "chemistry": chemistry, **sp})
            summary["geometry_sensitivity"].append({"case": case, "chemistry": chemistry, **geometry_sensitivity(res, chemistry)})
        # Gold-class sensitivity comparator only, never substituted for the predeclared hydrophobic primary map.
        rows, sp = exact_case(f"SpA_Au111_gold_comparator_5H7A_chain_{chain}", res, "gold", [221,220,218,33,34])
        summary["exact_anchor_metrics"].extend(rows)
        summary["spatial"].append({"case": f"SpA_Au111_gold_comparator_5H7A_chain_{chain}", "chemistry": "gold", **sp})
        summary["invariance"].append({"case": f"5H7A_chain_{chain}", "chemistry": "anionic", **ebase_apbs_invariance(res, "anionic")})

    # Determinism and bounds on one independent benchmark structure.
    mb1 = run_is(pdbs["1MBN"], "A", pH=7.4, temp=310.0, ionic=150.0)
    mb2 = run_is(pdbs["1MBN"], "A", pH=7.4, temp=310.0, ionic=150.0)
    for chem in IS.CHEMISTRIES:
        v1 = [(r["key"], r["propensity"], r["multiscale_persistence"]) for r in mb1["chemistries"][chem]["residues"]]
        v2 = [(r["key"], r["propensity"], r["multiscale_persistence"]) for r in mb2["chemistries"][chem]["residues"]]
        bounded = all(0 <= x[1] <= 100 and 0 <= x[2] <= 100 for x in v1)
        summary["mathematical_checks"].append({"check": f"1MBN {chem}: deterministic canonical outputs", "pass": v1 == v2})
        summary["mathematical_checks"].append({"check": f"1MBN {chem}: propensity/persistence bounded [0,100]", "pass": bounded})

    (OUT / "benchmark_results.json").write_text(json.dumps(summary, indent=2))

    # Flat CSVs for manuscript plotting/statistics.
    def flatten(rows):
        out = []
        for r in rows:
            rr = dict(r)
            for k,v in list(rr.items()):
                if isinstance(v, (dict,list)):
                    rr[k] = json.dumps(v, sort_keys=True)
            out.append(rr)
        return out
    for name, rows in [("exact_anchor_metrics.csv", summary["exact_anchor_metrics"]),
                       ("spatial_recovery.csv", summary["spatial"]),
                       ("region_enrichment.csv", summary["region_enrichment"]),
                       ("geometry_sensitivity.csv", summary["geometry_sensitivity"]),
                       ("invariance_checks.csv", summary["invariance"]),
                       ("mathematical_checks.csv", summary["mathematical_checks"])]:
        rows = flatten(rows)
        if not rows:
            continue
        fields = sorted(set().union(*(r.keys() for r in rows)))
        with (OUT / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

    # Compact markdown report.
    lines = ["# InterfaceScout frozen-v5.1 external validation", "",
             f"Model: `{IS.APP_VERSION}`", "",
             "## Mathematical/code checks"]
    for c in summary["mathematical_checks"]:
        lines.append(f"- {'PASS' if c['pass'] else 'FAIL'} — {c['check']}")
    lines += ["", "## PDB inventory", "```json", json.dumps(summary["pdb_inventory"], indent=2), "```",
              "", "## Notes"]
    for n in summary["notes"]: lines.append(f"- {n}")
    (OUT / "REPORT.md").write_text("\n".join(lines))

if __name__ == "__main__":
    main()
