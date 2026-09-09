"""InterfaceScout V2: coarse, weight-free protein-material interface prediction."""

from __future__ import annotations

import importlib
import os
import urllib.request
from typing import Any, Dict, Optional

import numpy as np

from .chemistry_freeze import apply_publication_chemistry
from .coarse_patch import build_coarse_patches, PATCH_SCALE_A
from .gnm import extract_ca_nodes, solve_gnm
from .model_settings import (
    MODEL_VERSION,
    SASA_POINTS,
    SASA_PROBE_A,
    SC_RSA_THRESHOLD,
    MULTISCALE_RADII_A,
    COARSE_PATCH_RADIUS_A,
    PARAMETER_SELECTION,
)
from .prepare import prepare_pdb_text
from .rin import build_rin, annotate_rin_percentiles, summarize_patch_rin
from .surface_modes import get_surface_mode


def _load_v1():
    module = importlib.import_module("main")
    apply_publication_chemistry(module)
    return module


def _obtain_pdb_text(pdb_id: Optional[str], pdb_text: Optional[str]) -> str:
    if pdb_text:
        return pdb_text
    pid = (pdb_id or "").strip().upper()
    if not pid:
        raise ValueError("Provide pdb_id or pdb_text")
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


def _geometry_only_gnm(prepared: str, cutoff_A: float) -> dict:
    """Return only the GNM fields required by coarse-patch geometry."""
    nodes = extract_ca_nodes(prepared)
    n = len(nodes)
    return {
        "cutoff_A": float(cutoff_A),
        "n_nodes": n,
        "nodes": nodes,
        "index": {x["key"]: i for i, x in enumerate(nodes)},
        "correlation_matrix": np.zeros((n, n), dtype=float),
    }


def analyze_interface_v2(
    *,
    surface: str,
    pH: float = 7.4,
    ionic_mM: float = 150.0,
    temp_K: float = 298.0,
    pdb_id: Optional[str] = None,
    pdb_text: Optional[str] = None,
    chain: Optional[str] = None,
    gnm_cutoff_A: float = 7.3,
) -> Dict[str, Any]:
    """Predict coarse plausible protein-material interface regions."""
    if not (0.0 <= float(pH) <= 14.0):
        raise ValueError("pH must be between 0 and 14")
    if float(ionic_mM) < 0:
        raise ValueError("ionic_mM must be non-negative")
    if float(temp_K) <= 0:
        raise ValueError("temp_K must be positive")

    mode = get_surface_mode(surface)
    raw = _obtain_pdb_text(pdb_id, pdb_text)
    prepared, prep_report = prepare_pdb_text(raw, chain=chain)

    v1 = _load_v1()
    request = v1.AnalyzeRequest(
        pdb_text=prepared,
        chain=None,
        env=v1.EnvParams(pH=float(pH), ionic=float(ionic_mM), temp=float(temp_K)),
    )
    v1_result = v1.analyze(request)
    if hasattr(v1_result, "body"):
        raise RuntimeError("Unexpected HTTP response object returned by V1 analyze()")

    canonical_only = os.environ.get("INTERFACESCOUT_VALIDATION_CANONICAL_ONLY") == "1"
    gnm = _geometry_only_gnm(prepared, float(gnm_cutoff_A)) if canonical_only else solve_gnm(prepared, cutoff_A=float(gnm_cutoff_A))
    patches = build_coarse_patches(v1_result=v1_result, chemistry=mode.chemistry, gnm=gnm)

    surface_keys = [str(r["key"]) for r in v1_result.get("surface_residues", []) if r.get("key")]
    if canonical_only:
        rin = {"cutoff_A": None, "n_nodes": 0, "n_edges": 0}
        for patch in patches:
            patch["rin_context"] = {"status": "disabled_in_canonical_only_validation"}
    else:
        rin = annotate_rin_percentiles(build_rin(prepared), surface_keys)
        for patch in patches:
            patch["rin_context"] = summarize_patch_rin(patch.get("members", []), rin)

    pareto_primary = [p for p in patches if int(p.get("pareto_front", 999)) == 1]

    return {
        "engine": "InterfaceScout V2",
        "version": MODEL_VERSION,
        "scope": {
            "prediction_unit": "coarse protein surface region / interface patch",
            "predicts_absolute_adsorption_free_energy": False,
            "predicts_adsorption_amount": False,
            "predicts_unique_orientation": False,
            "models_adsorption_induced_unfolding": False,
            "benchmark_fitted_weights": False,
            "residue_precision_claim": False,
            "rin_changes_patch_prediction": False,
            "gnm_changes_patch_prediction": False,
        },
        "input": {
            "pdb_id": (pdb_id or "").strip().upper() or None,
            "chain": prep_report.get("selected_chain", "ALL"),
            "surface": mode.key,
            "surface_label": mode.label,
            "primary_chemistry": mode.chemistry,
            "pH": float(pH),
            "ionic_mM": float(ionic_mM),
            "temp_K": float(temp_K),
        },
        "structure_preparation": prep_report,
        "method": {
            "core_question": "Where on the native folded protein is a plausible material-contact region under the defined surface chemistry and environment?",
            "chemistry_source": "InterfaceScout compatibility channel with publication chemistry corrections",
            "accessibility_source": "side-chain relative solvent accessibility",
            "sasa_probe_A": SASA_PROBE_A,
            "sasa_points_per_atom": SASA_POINTS,
            "scrsa_threshold": SC_RSA_THRESHOLD,
            "multiscale_radii_A": list(MULTISCALE_RADII_A),
            "multiscale_radius_basis": "selected on an independent adsorption-label-free development panel before external experimental validation",
            "parameter_selection": PARAMETER_SELECTION,
            "patch_radius_A": PATCH_SCALE_A,
            "patch_radius_basis": "separate 8 A structural neighbourhood rule; not equated to the selected outer multiscale aggregation radius",
            "patch_definition": "non-transitive local surface neighbourhood around a chemistry-patch maximum",
            "orientation": "coarse outward C-alpha face consistency",
            "ranking": "Pareto fronts across chemistry support, accessibility, patch coherence and orientation coherence; no weighted sum",
            "dynamics": "unweighted C-alpha GNM; downstream descriptive context only",
            "gnm_cutoff_A": float(gnm_cutoff_A),
            "rin": "heavy-atom-contact residue interaction network; downstream structural-network context only",
            "rin_heavy_atom_cutoff_A": rin["cutoff_A"],
        },
        "surface_mode": {
            "key": mode.key,
            "label": mode.label,
            "chemistry": mode.chemistry,
            "description": mode.description,
        },
        "n_patches": len(patches),
        "n_pareto_primary_patches": len(pareto_primary),
        "primary_patches": pareto_primary,
        "patches": patches,
        "network": {
            "n_rin_nodes": rin["n_nodes"],
            "n_rin_edges": rin["n_edges"],
            "rin_cutoff_A": rin["cutoff_A"],
            "interpretation": "RIN centrality/context is descriptive and is not used to improve or tune interface localization.",
        },
        "diagnostics": {
            "n_surface_residues": len(surface_keys),
            "surface_residue_keys": surface_keys,
            "canonical_only_validation": canonical_only,
        },
        "method_notes": [
            "Experimental interface labels are not inputs to patch construction or ranking.",
            "SASA sampling density and the 6/9 A multiscale pair were selected before external experimental validation on a separate label-free development panel.",
            "The 8 A coarse-patch radius is a separate structural-neighbourhood rule and was not selected by the 6/9 A multiscale analysis.",
            "Patch membership is intentionally coarse; individual residues are not claimed as precise adsorption contacts.",
            "Patch growth is non-transitive to prevent surface percolation into unrealistically large regions.",
            "The V2 publication chemistry excludes protonated Lys/Arg from protein-side acceptor eligibility for H-bond-donor surfaces.",
            "Numerical literature Ebase values inherited from V1 remain metadata only and never change V2 ranking.",
            "GNM is excluded from patch ranking and is retained only as native-state dynamic context.",
            "RIN is excluded from patch prediction and is used only to characterize the structural-network location of a predicted patch.",
            "Canonical-only validation may skip GNM/RIN numerical descriptors because they do not affect patch membership or ranking.",
            "Multiple Pareto-optimal patches are allowed because protein adsorption may have alternative plausible encounter interfaces.",
        ],
    }
