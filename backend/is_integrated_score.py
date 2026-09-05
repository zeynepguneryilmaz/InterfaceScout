"""Integrated InterfaceScout vNext score.

This module folds the validated lightweight descriptors into a single InterfaceScout
ranking while retaining each component for ablation/provenance.

Final residue/patch score components
------------------------------------
1. Legacy InterfaceScout core: chemistry + scRSA + pH state + 5/8 A patch persistence.
2. GNM intrinsic mobility: C-alpha Gaussian Network Model, 10 A contact cutoff.
3. Conditional APBS electrostatic complementarity: used only when the chemistry
   declares an expected protein-potential sign; neutral when APBS is unavailable or
   electrostatics are not designated as a primary mechanism.
4. Radial prominence: C-alpha distance from the protein C-alpha centroid.

No fitted coefficients are used. Components are converted to within-protein percentile
ranks and combined multiplicatively. The product is normalized to 0-100 within each
chemistry map. This file intentionally preserves the component values so all ablations
remain reproducible.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional
import numpy as np

GNM_CUTOFF_A = 10.0


def _percentile_rank(values: np.ndarray) -> np.ndarray:
    """Average-tie percentile ranks in [0,1], dependency-free."""
    x = np.asarray(values, dtype=float)
    n = len(x)
    if n <= 1:
        return np.ones(n, dtype=float)
    order = np.argsort(x, kind="stable")
    ranks = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i + 1
        while j < n and x[order[j]] == x[order[i]]:
            j += 1
        avg_rank = 0.5 * (i + j - 1)
        ranks[order[i:j]] = avg_rank
        i = j
    return ranks / float(n - 1)


def gnm_mobility_by_key(all_residues: List[Dict[str, Any]], cutoff_A: float = GNM_CUTOFF_A) -> Dict[str, float]:
    """Normalized residue MSF from a simple C-alpha Gaussian Network Model."""
    if not all_residues:
        return {}
    keys = [r["key"] for r in all_residues]
    xyz = np.asarray([[r["x"], r["y"], r["z"]] for r in all_residues], dtype=float)
    n = len(xyz)
    kirchhoff = np.zeros((n, n), dtype=float)
    for i in range(n - 1):
        d = np.linalg.norm(xyz[i + 1:] - xyz[i], axis=1)
        js = np.where(d <= cutoff_A)[0] + i + 1
        for j in js:
            kirchhoff[i, j] = kirchhoff[j, i] = -1.0
            kirchhoff[i, i] += 1.0
            kirchhoff[j, j] += 1.0
    vals, vecs = np.linalg.eigh(kirchhoff)
    positive = vals > 1e-8
    if not np.any(positive):
        return {k: 0.0 for k in keys}
    covariance = (vecs[:, positive] * (1.0 / vals[positive])) @ vecs[:, positive].T
    msf = np.maximum(np.diag(covariance), 0.0)
    lo, hi = float(np.min(msf)), float(np.max(msf))
    norm = (msf - lo) / (hi - lo) if hi > lo else np.zeros_like(msf)
    return {k: float(v) for k, v in zip(keys, norm)}


def radial_prominence_by_key(all_residues: List[Dict[str, Any]]) -> Dict[str, float]:
    """C-alpha radial prominence relative to the protein C-alpha centroid."""
    if not all_residues:
        return {}
    xyz = np.asarray([[r["x"], r["y"], r["z"]] for r in all_residues], dtype=float)
    centroid = np.mean(xyz, axis=0)
    d = np.linalg.norm(xyz - centroid, axis=1)
    return {r["key"]: float(v) for r, v in zip(all_residues, d)}


def electrostatic_compatibility(phi: Optional[float], expected_phi_sign: Optional[str]) -> Optional[float]:
    """Favorable magnitude of APBS potential for the selected surface chemistry.

    expected_phi_sign describes the desired protein potential at the interface:
    positive for an anionic material, negative for a cationic material.
    """
    if phi is None or expected_phi_sign not in {"positive", "negative"}:
        return None
    p = float(phi)
    if expected_phi_sign == "positive":
        return max(p, 0.0)
    return max(-p, 0.0)


def integrate_interfacescout_score(
    all_residues: List[Dict[str, Any]],
    surface_residues: List[Dict[str, Any]],
    chemistry_result: Dict[str, Any],
    expected_phi_sign: Optional[str],
) -> Dict[str, Any]:
    """Attach a unified `interfacescout_score` to every patch centre.

    APBS is conditional. If the chemistry is not electrostatic, or if APBS potential is
    unavailable, the APBS rank factor is set to 1.0 so the score degrades gracefully.
    """
    centres = chemistry_result.get("patch_centers", [])
    if not centres:
        chemistry_result["score_model"] = "InterfaceScout integrated score"
        chemistry_result["score_components"] = []
        return chemistry_result

    by_key = {r["key"]: r for r in surface_residues}
    gnm = gnm_mobility_by_key(all_residues)
    radial = radial_prominence_by_key(all_residues)

    core = np.asarray([float(c.get("multiscale_persistence", 0.0)) for c in centres], dtype=float)
    mob = np.asarray([float(gnm.get(c["center_key"], 0.0)) for c in centres], dtype=float)
    rad = np.asarray([float(radial.get(c["center_key"], 0.0)) for c in centres], dtype=float)

    apbs_applicable = expected_phi_sign in {"positive", "negative"}
    apbs_available = False
    apbs_raw = np.zeros(len(centres), dtype=float)
    if apbs_applicable:
        for i, c in enumerate(centres):
            phi = by_key.get(c["center_key"], {}).get("phi")
            ec = electrostatic_compatibility(phi, expected_phi_sign)
            if ec is not None:
                apbs_raw[i] = ec
                apbs_available = True

    core_r = _percentile_rank(core)
    gnm_r = _percentile_rank(mob)
    radial_r = _percentile_rank(rad)
    apbs_r = _percentile_rank(apbs_raw) if (apbs_applicable and apbs_available) else np.ones(len(centres), dtype=float)

    product = core_r * gnm_r * radial_r * apbs_r
    mx = float(np.max(product)) if len(product) else 0.0
    final = 100.0 * product / mx if mx > 0 else np.zeros_like(product)

    for i, c in enumerate(centres):
        phi = by_key.get(c["center_key"], {}).get("phi")
        ec = electrostatic_compatibility(phi, expected_phi_sign)
        c["legacy_core_score"] = round(float(core[i]), 3)
        c["gnm_mobility"] = round(float(mob[i]), 6)
        c["radial_prominence_A"] = round(float(rad[i]), 4)
        c["apbs_phi"] = phi
        c["apbs_compatibility"] = round(float(ec), 6) if ec is not None else None
        c["rank_core"] = round(float(core_r[i]), 6)
        c["rank_gnm"] = round(float(gnm_r[i]), 6)
        c["rank_radial"] = round(float(radial_r[i]), 6)
        c["rank_apbs"] = round(float(apbs_r[i]), 6)
        c["interfacescout_score"] = round(float(final[i]), 3)

    centres.sort(key=lambda c: (-c["interfacescout_score"], -c["legacy_core_score"], c["center_key"]))
    chemistry_result["top_patch"] = centres[0] if centres else None
    chemistry_result["top_patches"] = centres[:10]
    chemistry_result["patch_centers"] = centres
    chemistry_result["score_model"] = "InterfaceScout integrated score"
    chemistry_result["score_formula"] = "normalized percentile-rank product: core × GNM mobility × radial prominence × conditional APBS complementarity"
    chemistry_result["score_components"] = [
        "chemistry/pH/scRSA/5-8A patch persistence (legacy core)",
        "GNM intrinsic mobility (10 A C-alpha network)",
        "radial prominence",
        "conditional APBS electrostatic complementarity",
    ]
    chemistry_result["apbs_in_final_score"] = bool(apbs_applicable and apbs_available)
    chemistry_result["apbs_policy"] = "conditional; neutral factor when non-electrostatic or unavailable"
    return chemistry_result
