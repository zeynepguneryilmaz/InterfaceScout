"""Protein-only InterfaceScout prototype for material-blind ablation experiments.

Scientific question
-------------------
How much information about an experimentally observed protein-material contact
region is already encoded in the native protein itself, before the identity of
the contacting material is supplied?

The predictor therefore accepts protein structure + solution environment only.
Material identity and experimental interface labels are deliberately withheld.

Four predeclared ranking variants are exposed for ablation:
    M0  local surface physicochemistry + geometry
    M1  M0 + residue-interaction-network (RIN) context
    M2  M0 + GNM dynamic coherence
    M3  M0 + RIN + GNM

RIN/GNM do not change patch membership. They only refine the ordering of the
same candidate local surface patches. No adsorption-label-fitted weights are
used. The directional hypotheses tested here are explicitly exploratory:
structurally connected/bridging surface patches and dynamically coherent local
patches may be more interface-ready than otherwise similar exposed patches.
"""

from __future__ import annotations

import importlib
from typing import Any, Dict, Iterable, List, Optional, Sequence

import numpy as np

from .geometry import (
    build_surface_geometry,
    ca_distance,
    same_face,
    patch_orientation_coherence,
    patch_diameter_A,
)
from .gnm import solve_gnm
from .prepare import prepare_pdb_text
from .rin import build_rin, annotate_rin_percentiles, summarize_patch_rin

PATCH_SCALE_A = 8.0  # same frozen local spatial scale used by publication V1
GENERIC_FEATURES = (
    "positive_charge",
    "negative_charge",
    "hydrophobic",
    "aromatic",
    "hbond_donor",
    "hbond_acceptor",
)

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
AROMATIC = {"PHE", "TRP", "TYR", "HIS"}
HBOND_DONOR = {"SER", "THR", "TYR", "ASN", "GLN", "TRP", "LYS", "ARG", "HIS"}
HBOND_ACCEPTOR = {"ASP", "GLU", "ASN", "GLN", "SER", "THR", "TYR", "HIS"}


def _load_v1():
    return importlib.import_module("main")


def _feature_vector(row: dict) -> Dict[str, float]:
    rn = str(row.get("res_name", "")).upper()
    q = float(row.get("charge_descriptor", 0.0) or 0.0)
    return {
        "positive_charge": max(q, 0.0),
        "negative_charge": max(-q, 0.0),
        "hydrophobic": 1.0 if rn in HYDROPHOBIC else 0.0,
        "aromatic": 1.0 if rn in AROMATIC else 0.0,
        "hbond_donor": 1.0 if rn in HBOND_DONOR else 0.0,
        "hbond_acceptor": 1.0 if rn in HBOND_ACCEPTOR else 0.0,
    }


def _mean_profile(keys: Iterable[str], rows: Dict[str, dict]) -> Dict[str, float]:
    valid = [k for k in keys if k in rows]
    if not valid:
        return {f: 0.0 for f in GENERIC_FEATURES}
    vecs = [_feature_vector(rows[k]) for k in valid]
    return {f: float(np.mean([v[f] for v in vecs])) for f in GENERIC_FEATURES}


def _pair_corr(keys: Sequence[str], gnm: dict) -> List[float]:
    index = gnm["index"]
    matrix = gnm["correlation_matrix"]
    vals: List[float] = []
    for i in range(len(keys)):
        if keys[i] not in index:
            continue
        ii = index[keys[i]]
        for j in range(i + 1, len(keys)):
            if keys[j] not in index:
                continue
            jj = index[keys[j]]
            vals.append(float(matrix[ii, jj]))
    return vals


def _pareto_fronts(rows: List[dict], fields: Sequence[str]) -> Dict[int, int]:
    """Return row-index -> Pareto front, maximizing all named fields."""
    if not rows:
        return {}
    if not fields:
        return {i: 1 for i in range(len(rows))}

    def dominates(a: dict, b: dict) -> bool:
        av = [float(a.get(f, 0.0) or 0.0) for f in fields]
        bv = [float(b.get(f, 0.0) or 0.0) for f in fields]
        return all(x >= y for x, y in zip(av, bv)) and any(x > y for x, y in zip(av, bv))

    remaining = set(range(len(rows)))
    result: Dict[int, int] = {}
    front = 1
    while remaining:
        current = []
        for i in remaining:
            if not any(dominates(rows[j], rows[i]) for j in remaining if j != i):
                current.append(i)
        for i in current:
            result[i] = front
        remaining -= set(current)
        front += 1
    return result


def _secondary_fronts(rows: List[dict], base_fronts: Dict[int, int], fields: Sequence[str]) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for bf in sorted(set(base_fronts.values())):
        idx = [i for i, f in base_fronts.items() if f == bf]
        subset = [rows[i] for i in idx]
        subfront = _pareto_fronts(subset, fields)
        for local_i, global_i in enumerate(idx):
            out[global_i] = int(subfront.get(local_i, 1))
    return out


def _local_patch(center: str, geometry: dict) -> List[str]:
    members = [
        k for k in geometry["coords"]
        if ca_distance(center, k, geometry) <= PATCH_SCALE_A and same_face(center, k, geometry)
    ]
    if center not in members:
        members.append(center)
    return sorted(set(members))


def _candidate_patches(*, v1_result: dict, gnm: dict, rin: dict) -> List[dict]:
    geometry = build_surface_geometry(v1_result, gnm)
    surface_rows = {
        str(r["key"]): r
        for r in v1_result.get("surface_residues", [])
        if r.get("key") and str(r["key"]) in geometry["coords"]
    }
    surface_keys = sorted(surface_rows)
    surface_profile = _mean_profile(surface_keys, surface_rows)

    patches: List[dict] = []
    for center in surface_keys:
        members = _local_patch(center, geometry)
        profile = _mean_profile(members, surface_rows)
        enrichment = {f: float(profile[f] - surface_profile[f]) for f in GENERIC_FEATURES}
        chemical_contrast = max([max(v, 0.0) for v in enrichment.values()] + [0.0])

        access = [float(geometry["scrsa"].get(k, 0.0)) for k in members]
        corr = _pair_corr(members, gnm)
        fluctuations = [
            float(gnm["residue_metrics"][k]["normalized_fluctuation"])
            for k in members if k in gnm["residue_metrics"]
        ]
        rin_summary = summarize_patch_rin(members, rin)

        patches.append({
            "center_key": center,
            "center_residue": geometry["meta"].get(center),
            "members": members,
            "member_residues": [geometry["meta"][k] for k in members if k in geometry["meta"]],
            "n_members": len(members),
            "diameter_A": patch_diameter_A(members, geometry),
            "center_accessibility": float(geometry["scrsa"].get(center, 0.0)),
            "mean_accessibility": float(np.mean(access)) if access else 0.0,
            "orientation_coherence": patch_orientation_coherence(members, geometry),
            "chemical_profile": profile,
            "surface_chemical_profile": surface_profile,
            "chemical_enrichment": enrichment,
            "chemical_contrast": float(chemical_contrast),
            "dynamic_coupling_abs": float(np.mean(np.abs(corr))) if corr else 0.0,
            "dynamic_coupling_signed": float(np.mean(corr)) if corr else 0.0,
            "mean_normalized_fluctuation": float(np.mean(fluctuations)) if fluctuations else None,
            "rin": rin_summary,
            "rin_mean_degree_percentile_surface": float(rin_summary.get("mean_degree_percentile_surface") or 0.0),
            "rin_max_betweenness_percentile_surface": float(rin_summary.get("max_betweenness_percentile_surface") or 0.0),
            "rin_internal_edge_density": float(rin_summary.get("internal_edge_density") or 0.0),
            "rin_boundary_edge_fraction": float(rin_summary.get("boundary_edge_fraction") or 0.0),
        })
    return patches


BASE_FIELDS = (
    "center_accessibility",
    "mean_accessibility",
    "orientation_coherence",
    "chemical_contrast",
)
RIN_FIELDS = (
    "rin_mean_degree_percentile_surface",
    "rin_max_betweenness_percentile_surface",
)
GNM_FIELDS = ("dynamic_coupling_abs",)

VARIANT_FIELDS = {
    "M0_surface": (),
    "M1_surface_rin": RIN_FIELDS,
    "M2_surface_gnm": GNM_FIELDS,
    "M3_surface_rin_gnm": RIN_FIELDS + GNM_FIELDS,
}


def _rank_variant(patches: List[dict], variant: str, geometry: dict) -> List[dict]:
    if variant not in VARIANT_FIELDS:
        raise ValueError(f"Unknown ablation variant: {variant}")
    base_front = _pareto_fronts(patches, BASE_FIELDS)
    secondary = _secondary_fronts(patches, base_front, VARIANT_FIELDS[variant])

    indexed = []
    for i, patch in enumerate(patches):
        row = dict(patch)
        row["base_pareto_front"] = int(base_front[i])
        row["secondary_pareto_front"] = int(secondary[i])
        indexed.append(row)

    # Predeclared lexicographic tie-breakers; no fitted numerical weights.
    indexed.sort(key=lambda p: (
        int(p["base_pareto_front"]),
        int(p["secondary_pareto_front"]),
        -float(p["chemical_contrast"]),
        -float(p["mean_accessibility"]),
        -float(p["orientation_coherence"]),
        str(p["center_key"]),
    ))

    # Keep distinct local hypotheses rather than multiple nearly identical
    # centres from the same 8 A surface neighbourhood.
    kept: List[dict] = []
    for patch in indexed:
        center = str(patch["center_key"])
        if any(
            ca_distance(center, str(prev["center_key"]), geometry) <= PATCH_SCALE_A
            and same_face(center, str(prev["center_key"]), geometry)
            for prev in kept
        ):
            continue
        kept.append(patch)

    for rank, patch in enumerate(kept, start=1):
        patch["display_rank"] = rank
        patch["variant"] = variant
    return kept


def analyze_protein_only(
    *,
    pdb_text: str,
    chain: Optional[str] = None,
    pH: float = 7.4,
    ionic_mM: float = 150.0,
    temp_K: float = 298.0,
    gnm_cutoff_A: float = 7.3,
    rin_cutoff_A: float = 4.5,
) -> Dict[str, Any]:
    """Generate material-blind candidate biointerface patches and ablations."""
    if not (0.0 <= float(pH) <= 14.0):
        raise ValueError("pH must be between 0 and 14")
    if float(ionic_mM) < 0:
        raise ValueError("ionic_mM must be non-negative")

    prepared, prep_report = prepare_pdb_text(pdb_text, chain=chain)
    v1 = _load_v1()
    request = v1.AnalyzeRequest(
        pdb_text=prepared,
        chain=None,
        env=v1.EnvParams(pH=float(pH), ionic=float(ionic_mM), temp=float(temp_K)),
    )
    v1_result = v1.analyze(request)
    if hasattr(v1_result, "body"):
        raise RuntimeError("Unexpected HTTP response object returned by V1 analyze()")

    gnm = solve_gnm(prepared, cutoff_A=float(gnm_cutoff_A))
    geometry = build_surface_geometry(v1_result, gnm)
    rin = build_rin(prepared, cutoff_A=float(rin_cutoff_A))
    rin = annotate_rin_percentiles(rin, geometry["coords"].keys())
    candidates = _candidate_patches(v1_result=v1_result, gnm=gnm, rin=rin)

    variants = {
        name: _rank_variant(candidates, name, geometry)
        for name in VARIANT_FIELDS
    }

    return {
        "engine": "InterfaceScout protein-only ablation prototype",
        "version": "0.1.0",
        "scope": {
            "material_identity_used": False,
            "experimental_labels_used_for_prediction": False,
            "prediction_unit": "coarse native/native-like protein surface patch",
            "models_adsorption_induced_unfolding": False,
            "predicts_adsorption_energy": False,
            "predicts_unique_material_specific_orientation": False,
        },
        "input": {
            "chain": prep_report.get("selected_chain", "ALL"),
            "pH": float(pH),
            "ionic_mM": float(ionic_mM),
            "temp_K": float(temp_K),
        },
        "method": {
            "patch_scale_A": PATCH_SCALE_A,
            "base_fields": list(BASE_FIELDS),
            "rin_fields": list(RIN_FIELDS),
            "gnm_fields": list(GNM_FIELDS),
            "ranking": "base Pareto front; optional RIN/GNM Pareto refinement within each base front; predeclared lexicographic tie-breakers; no numerical weights",
            "chemical_contrast": "maximum positive enrichment of a generic protein-side physicochemical feature relative to that protein's exposed surface background",
            "material_blinding": "material identity is deliberately absent from predictor input",
        },
        "structure_preparation": prep_report,
        "n_surface_residues": len(geometry["coords"]),
        "n_candidate_centres": len(candidates),
        "variants": variants,
    }
