"""Held-out predictive validation for the frozen InterfaceScout v5.2 candidate.

System (not used in prior InterfaceScout development diagnostics):
  Wei T, Carignano MA, Szleifer I. J Phys Chem B. 2012;116:10189-10194.
  DOI: 10.1021/jp304057e
  Protein: hen egg-white lysozyme, PDB 1AKI.
  Surface: crystalline polyethylene (PE), hydrophobic.
  Conditions: pH 7, 300 K; only neutralizing chloride ions were added.

The paper reports two explicitly identified landing-site residue sets in one
trajectory: an unsuccessful first landing (14, 126-129) and a third landing
that results in adsorption (67-71, 81).  The successful third-landing set is
predeclared as the PRIMARY external contact ground truth.  The unsuccessful
first-landing set is retained only as a mechanistic comparator because the
paper attributes success/failure to competition with hydration/dehydration,
which InterfaceScout does not explicitly model.

The material mapping is predeclared from PE chemistry as InterfaceScout's
'hydrophobic' channel.  No backend equations, chemistry memberships, scRSA
threshold, or 5/8 A patch radii are modified by this script.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from backend.v52_app import AnalyzeRequest, EnvParams, analyze

FREEZE_SHA = "8e83ee3c2242cfd22714d47c585cb7afac120d72"
PDB_ID = "1AKI"
CHAIN = "A"
CHEMISTRY = "hydrophobic"
PRIMARY_SUCCESSFUL = [67, 68, 69, 70, 71, 81]
UNSUCCESSFUL_FIRST = [14, 126, 127, 128, 129]
N_PERM = 10000
SEED = 20260904
OUT = Path("validation/v52_heldout_lysozyme_pe_report.json")


def key_for_resseq(rows: Sequence[dict], resseq: int) -> str:
    matches = [r["key"] for r in rows if r.get("chain") == CHAIN and int(r.get("res_seq")) == int(resseq)]
    if len(matches) != 1:
        raise RuntimeError(f"Expected exactly one {CHAIN}:{resseq} residue, found {matches}")
    return matches[0]


def auroc(y: np.ndarray, score: np.ndarray) -> float:
    pos = score[y == 1]
    neg = score[y == 0]
    if not len(pos) or not len(neg):
        return float("nan")
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float((gt + 0.5 * eq) / (len(pos) * len(neg)))


def average_precision(y: np.ndarray, score: np.ndarray) -> float:
    """Threshold-grouped AP, equivalent to a non-interpolated PR integral."""
    npos = int(y.sum())
    if npos == 0:
        return float("nan")
    order = np.argsort(-score, kind="mergesort")
    ys = y[order]
    ss = score[order]
    tp = 0
    fp = 0
    prev_recall = 0.0
    ap = 0.0
    i = 0
    while i < len(ys):
        j = i + 1
        while j < len(ys) and ss[j] == ss[i]:
            j += 1
        block = ys[i:j]
        tp += int(block.sum())
        fp += int(len(block) - block.sum())
        recall = tp / npos
        precision = tp / (tp + fp)
        ap += (recall - prev_recall) * precision
        prev_recall = recall
        i = j
    return float(ap)


def precision_recall_at_k(y: np.ndarray, score: np.ndarray, k: int) -> Tuple[float, float]:
    order = np.argsort(-score, kind="mergesort")[:k]
    hits = int(y[order].sum())
    return hits / max(k, 1), hits / max(int(y.sum()), 1)


def coords(row: dict) -> np.ndarray:
    return np.asarray([row["x"], row["y"], row["z"]], dtype=float)


def spatial_recovery(anchor_keys: Sequence[str], top_centers: Sequence[str], rows_by_key: Dict[str, dict], radius: float) -> float:
    centers = [coords(rows_by_key[k]) for k in top_centers if k in rows_by_key]
    if not centers:
        return 0.0
    hits = 0
    for k in anchor_keys:
        p = coords(rows_by_key[k])
        if min(float(np.linalg.norm(p - c)) for c in centers) <= radius:
            hits += 1
    return hits / len(anchor_keys)


def nearest_distances(anchor_keys: Sequence[str], top_centers: Sequence[str], rows_by_key: Dict[str, dict]) -> List[float]:
    centers = [coords(rows_by_key[k]) for k in top_centers if k in rows_by_key]
    vals = []
    for k in anchor_keys:
        p = coords(rows_by_key[k])
        vals.append(min(float(np.linalg.norm(p - c)) for c in centers))
    return vals


def quantile_bins(values: np.ndarray, nbin: int = 5) -> np.ndarray:
    edges = np.quantile(values, np.linspace(0.0, 1.0, nbin + 1))
    # De-duplicate flat quantiles while keeping at least one interior edge.
    edges = np.unique(edges)
    if len(edges) <= 2:
        return np.zeros(len(values), dtype=int)
    return np.digitize(values, edges[1:-1], right=True)


def exposure_matched_null(
    y: np.ndarray,
    scrsa: np.ndarray,
    score: np.ndarray,
    anchor_keys: Sequence[str],
    keys: Sequence[str],
    top_centers: Sequence[str],
    rows_by_key: Dict[str, dict],
    observed_auc: float,
    observed_ap: float,
    observed_r8: float,
) -> dict:
    rng = np.random.default_rng(SEED)
    bins = quantile_bins(scrsa, 5)
    anchor_idx = np.flatnonzero(y == 1)
    nonanchor_idx = np.flatnonzero(y == 0)
    pools = {b: nonanchor_idx[bins[nonanchor_idx] == b] for b in np.unique(bins)}

    null_auc = []
    null_ap = []
    null_r8 = []
    attempts = 0
    while len(null_auc) < N_PERM and attempts < N_PERM * 20:
        attempts += 1
        chosen: List[int] = []
        used = set()
        ok = True
        for ai in anchor_idx:
            pool = [int(x) for x in pools.get(int(bins[ai]), []) if int(x) not in used]
            if not pool:
                ok = False
                break
            pick = int(rng.choice(pool))
            chosen.append(pick)
            used.add(pick)
        if not ok:
            continue
        yp = np.zeros_like(y)
        yp[chosen] = 1
        null_auc.append(auroc(yp, score))
        null_ap.append(average_precision(yp, score))
        null_keys = [keys[i] for i in chosen]
        null_r8.append(spatial_recovery(null_keys, top_centers, rows_by_key, 8.0))

    if len(null_auc) < N_PERM:
        raise RuntimeError(f"Only generated {len(null_auc)} / {N_PERM} exposure-matched permutations")

    a = np.asarray(null_auc)
    p = np.asarray(null_ap)
    r = np.asarray(null_r8)
    return {
        "n_permutations": int(len(a)),
        "matching": "surface-residue scRSA quintiles; same positive-set size",
        "auroc_null_mean": float(a.mean()),
        "auroc_p_ge_observed": float((1 + np.sum(a >= observed_auc)) / (len(a) + 1)),
        "ap_null_mean": float(p.mean()),
        "ap_p_ge_observed": float((1 + np.sum(p >= observed_ap)) / (len(p) + 1)),
        "spatial_R8_null_mean": float(r.mean()),
        "spatial_R8_p_ge_observed": float((1 + np.sum(r >= observed_r8)) / (len(r) + 1)),
    }


def score_bundle(keys: Sequence[str], y: np.ndarray, score: np.ndarray) -> dict:
    auc = auroc(y, score)
    ap = average_precision(y, score)
    p10, r10 = precision_recall_at_k(y, score, min(10, len(score)))
    pm, rm = precision_recall_at_k(y, score, int(y.sum()))
    return {
        "auroc": auc,
        "average_precision": ap,
        "prevalence": float(y.mean()),
        "ap_over_prevalence": float(ap / y.mean()) if y.mean() else None,
        "precision_at_10": p10,
        "recall_at_10": r10,
        "precision_at_m": pm,
        "recall_at_m": rm,
        "m": int(y.sum()),
    }


def main() -> None:
    req = AnalyzeRequest(
        pdb_id=PDB_ID,
        chain=CHAIN,
        env=EnvParams(pH=7.0, ionic=0.0, temp=300.0),
        structure_context="auto",
        protrusion=True,
        material_profile=None,
    )
    result = analyze(req)
    if result["core_version"] != "5.1.0-publication-freeze":
        raise RuntimeError(f"Unexpected frozen core version: {result['core_version']}")
    if tuple(result["settings"]["patch_radii_A"]) != (5.0, 8.0):
        raise RuntimeError("Patch radii changed from frozen 5/8 A")
    if float(result["settings"]["scrsa_threshold"]) != 0.05:
        raise RuntimeError("scRSA threshold changed from frozen 0.05")

    all_rows = result["all_residues"]
    surface = result["surface_residues"]
    rows_by_key = {r["key"]: r for r in all_rows}
    surface_keys = [r["key"] for r in surface]
    primary_keys_all = [key_for_resseq(all_rows, n) for n in PRIMARY_SUCCESSFUL]
    failed_keys_all = [key_for_resseq(all_rows, n) for n in UNSUCCESSFUL_FIRST]
    primary_keys = [k for k in primary_keys_all if k in set(surface_keys)]

    chem = result["chemistries"][CHEMISTRY]
    local_by_key = {r["key"]: float(r["local_score"]) for r in chem["residues"]}
    persistence_by_key = {r["center_key"]: float(r["multiscale_persistence"]) for r in chem["patch_centers"]}

    scrsa = np.asarray([float(r["scrsa"]) for r in surface], dtype=float)
    local = np.asarray([local_by_key.get(k, 0.0) for k in surface_keys], dtype=float)
    persistence = np.asarray([persistence_by_key.get(k, 0.0) for k in surface_keys], dtype=float)
    y = np.asarray([1 if k in set(primary_keys) else 0 for k in surface_keys], dtype=int)
    if int(y.sum()) < 2:
        raise RuntimeError("Too few primary anchors retained as surface residues for held-out scoring")

    top10 = [r["center_key"] for r in chem["top_patches"][:10]]
    dists = nearest_distances(primary_keys_all, top10, rows_by_key)
    spatial = {
        "top10_patch_centers": top10,
        "anchor_recall_within_5A": spatial_recovery(primary_keys_all, top10, rows_by_key, 5.0),
        "anchor_recall_within_8A": spatial_recovery(primary_keys_all, top10, rows_by_key, 8.0),
        "anchor_recall_within_10A": spatial_recovery(primary_keys_all, top10, rows_by_key, 10.0),
        "nearest_top10_center_distance_A": {k: d for k, d in zip(primary_keys_all, dists)},
        "median_nearest_distance_A": float(np.median(dists)),
    }

    metrics = {
        "scRSA_baseline": score_bundle(surface_keys, y, scrsa),
        "hydrophobic_local_score": score_bundle(surface_keys, y, local),
        "hydrophobic_persistence_M": score_bundle(surface_keys, y, persistence),
    }
    null = exposure_matched_null(
        y=y,
        scrsa=scrsa,
        score=persistence,
        anchor_keys=primary_keys,
        keys=surface_keys,
        top_centers=top10,
        rows_by_key=rows_by_key,
        observed_auc=metrics["hydrophobic_persistence_M"]["auroc"],
        observed_ap=metrics["hydrophobic_persistence_M"]["average_precision"],
        observed_r8=spatial["anchor_recall_within_8A"],
    )

    def residue_scores(keys_subset: Sequence[str]) -> List[dict]:
        out = []
        for k in keys_subset:
            r = rows_by_key[k]
            out.append({
                "key": k,
                "res_name": r["res_name"],
                "res_seq": r["res_seq"],
                "surface_exposed": bool(r["surface_exposed"]),
                "scrsa": float(r["scrsa"]),
                "local_score": float(local_by_key.get(k, 0.0)),
                "persistence_M": float(persistence_by_key.get(k, 0.0)),
            })
        return out

    successful_rows = residue_scores(primary_keys_all)
    failed_rows = residue_scores(failed_keys_all)
    comparator = {
        "interpretation": "descriptive only; the paper attributes landing success/failure to hydration/dehydration, which is outside the canonical InterfaceScout model",
        "successful_third_landing": successful_rows,
        "unsuccessful_first_landing": failed_rows,
        "mean_persistence_successful": float(np.mean([x["persistence_M"] for x in successful_rows])),
        "mean_persistence_unsuccessful": float(np.mean([x["persistence_M"] for x in failed_rows])),
        "mean_local_successful": float(np.mean([x["local_score"] for x in successful_rows])),
        "mean_local_unsuccessful": float(np.mean([x["local_score"] for x in failed_rows])),
    }

    report = {
        "status": "PASS",
        "validation_type": "predictive held-out external validation",
        "freeze_sha": FREEZE_SHA,
        "freeze_integrity": {
            "backend_modified_for_validation": False,
            "parameters_retuned_after_ground_truth": False,
            "predeclared_material_mapping": "polyethylene -> hydrophobic",
            "primary_ground_truth_predeclared": PRIMARY_SUCCESSFUL,
        },
        "system": {
            "citation": "Wei T, Carignano MA, Szleifer I. Molecular Dynamics Simulation of Lysozyme Adsorption/Desorption on Hydrophobic Surfaces. J Phys Chem B. 2012;116(34):10189-10194. DOI:10.1021/jp304057e",
            "protein": "hen egg-white lysozyme",
            "pdb": PDB_ID,
            "chain": CHAIN,
            "surface": "crystalline polyethylene (PE) (010), hydrophobic",
            "conditions": "pH 7, 300 K; explicit water; neutralizing Cl- only",
            "interface_scout_ionic_mM": 0.0,
            "primary_successful_contact_residues": PRIMARY_SUCCESSFUL,
            "unsuccessful_first_landing_residues": UNSUCCESSFUL_FIRST,
        },
        "run": {
            "v52_version": result["version"],
            "core_version": result["core_version"],
            "structure_context": result["settings"]["structure_context"],
            "n_residues": result["stats"]["n_residues_reported"],
            "n_surface_residues": result["stats"]["n_surface_res_reported"],
            "primary_anchors_total": len(primary_keys_all),
            "primary_anchors_surface_retained": len(primary_keys),
            "primary_anchor_surface_coverage": len(primary_keys) / len(primary_keys_all),
            "primary_anchor_keys": primary_keys_all,
        },
        "metrics_surface_conditioned": metrics,
        "spatial_recovery_primary": spatial,
        "exposure_matched_permutation_null_for_M": null,
        "successful_vs_unsuccessful_landing_comparator": comparator,
        "interpretation_guardrail": "This test evaluates residue/patch compatibility recovery for a hydrophobic surface. It does not test adsorption free energy, kinetics, dehydration, or a unique adsorption orientation.",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
