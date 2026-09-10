"""InterfaceScout coarse protein-material interface prediction engine.

The public analysis path prepares the protein and computes solvent exposure once,
then derives all chemistry maps from that shared state. A single run therefore
supports instant switching among chemistry maps in the user interface.
"""
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
    PARAMETER_SELECTION,
)
from .prepare import prepare_pdb_text
from .rin import build_rin, annotate_rin_percentiles, summarize_patch_rin
from .surface_modes import SURFACE_MODES, get_surface_mode

MAP_ORDER = [
    "anionic",
    "cationic",
    "hydrophobic",
    "pi_carbon",
    "hbond_donor",
    "hbond_acceptor",
    "oxide",
    "hydroxyapatite",
    "metal_coord",
    "gold",
    "phosphate",
]


def _load_core():
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


def _prepare_shared_context(
    *,
    pH: float,
    ionic_mM: float,
    temp_K: float,
    pdb_id: Optional[str],
    pdb_text: Optional[str],
    chain: Optional[str],
    gnm_cutoff_A: float,
) -> dict:
    if not (0.0 <= float(pH) <= 14.0):
        raise ValueError("pH must be between 0 and 14")
    if float(ionic_mM) < 0:
        raise ValueError("ionic_mM must be non-negative")
    if float(temp_K) <= 0:
        raise ValueError("temp_K must be positive")

    raw = _obtain_pdb_text(pdb_id, pdb_text)
    prepared, prep_report = prepare_pdb_text(raw, chain=chain)

    core = _load_core()
    request = core.AnalyzeRequest(
        pdb_text=prepared,
        chain=None,
        env=core.EnvParams(pH=float(pH), ionic=float(ionic_mM), temp=float(temp_K)),
    )
    core_result = core.analyze(request)
    if hasattr(core_result, "body"):
        raise RuntimeError("Unexpected HTTP response object returned by the chemistry core")

    canonical_only = os.environ.get("INTERFACESCOUT_VALIDATION_CANONICAL_ONLY") == "1"
    gnm = _geometry_only_gnm(prepared, float(gnm_cutoff_A)) if canonical_only else solve_gnm(prepared, cutoff_A=float(gnm_cutoff_A))

    surface_keys = [str(r["key"]) for r in core_result.get("surface_residues", []) if r.get("key")]
    if canonical_only:
        rin = {"cutoff_A": None, "n_nodes": 0, "n_edges": 0}
    else:
        rin = annotate_rin_percentiles(build_rin(prepared), surface_keys)

    return {
        "raw": raw,
        "prepared": prepared,
        "prep_report": prep_report,
        "core_result": core_result,
        "canonical_only": canonical_only,
        "gnm": gnm,
        "rin": rin,
        "surface_keys": surface_keys,
    }


def _build_map(context: dict, chemistry: str) -> dict:
    core_result = context["core_result"]
    patches = build_coarse_patches(v1_result=core_result, chemistry=chemistry, gnm=context["gnm"])

    if context["canonical_only"]:
        for patch in patches:
            patch["rin_context"] = {"status": "disabled_in_canonical_only_validation"}
    else:
        for patch in patches:
            patch["rin_context"] = summarize_patch_rin(patch.get("members", []), context["rin"])

    chemistry_meta = core_result.get("chemistries", {}).get(chemistry, {})
    primary = [p for p in patches if int(p.get("pareto_front", 999)) == 1]
    presets = [
        {
            "key": mode.key,
            "label": mode.label,
            "description": mode.description,
        }
        for mode in SURFACE_MODES.values()
        if mode.chemistry == chemistry
    ]

    return {
        "key": chemistry,
        "label": chemistry_meta.get("label", chemistry.replace("_", " ").title()),
        "surface_group": chemistry_meta.get("surface_group", ""),
        "description": chemistry_meta.get("description", ""),
        "material_presets": presets,
        "n_patches": len(patches),
        "n_pareto_primary_patches": len(primary),
        "primary_patches": primary,
        "patches": patches,
    }


def _shared_response(context: dict, *, pdb_id: Optional[str], pH: float, ionic_mM: float, temp_K: float, gnm_cutoff_A: float) -> dict:
    return {
        "engine": "InterfaceScout",
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
            "chain": context["prep_report"].get("selected_chain", "ALL"),
            "pH": float(pH),
            "ionic_mM": float(ionic_mM),
            "temp_K": float(temp_K),
        },
        "structure_preparation": context["prep_report"],
        "method": {
            "core_question": "Where on the native folded protein are plausible material-contact regions under the defined chemistry and environment?",
            "chemistry_source": "InterfaceScout publication-frozen compatibility channels",
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
            "rin_heavy_atom_cutoff_A": context["rin"]["cutoff_A"],
            "single_pass_all_maps": True,
        },
        "network": {
            "n_rin_nodes": context["rin"]["n_nodes"],
            "n_rin_edges": context["rin"]["n_edges"],
            "rin_cutoff_A": context["rin"]["cutoff_A"],
            "interpretation": "RIN centrality/context is descriptive and is not used to improve or tune interface localization.",
        },
        "diagnostics": {
            "n_surface_residues": len(context["surface_keys"]),
            "surface_residue_keys": context["surface_keys"],
            "canonical_only_validation": context["canonical_only"],
        },
        "method_notes": [
            "All chemistry maps in one analysis share the same prepared structure, SASA/scRSA calculation, pH state calculation, GNM and RIN context.",
            "Experimental interface labels are not inputs to patch construction or ranking.",
            "SASA sampling density and the 6/9 A multiscale pair were selected before external experimental validation on a separate label-free development panel.",
            "The 8 A coarse-patch radius is a separate structural-neighbourhood rule and was not selected by the 6/9 A multiscale analysis.",
            "Patch membership is intentionally coarse; individual residues are not claimed as precise adsorption contacts.",
            "Patch growth is non-transitive to prevent surface percolation into unrealistically large regions.",
            "GNM and RIN are excluded from patch prediction/ranking and retained only as native-state context.",
            "Multiple Pareto-optimal patches are allowed because protein adsorption may have alternative plausible encounter interfaces.",
        ],
    }


def analyze_all_maps(
    *,
    pH: float = 7.4,
    ionic_mM: float = 150.0,
    temp_K: float = 298.0,
    pdb_id: Optional[str] = None,
    pdb_text: Optional[str] = None,
    chain: Optional[str] = None,
    initial_surface: Optional[str] = None,
    gnm_cutoff_A: float = 7.3,
) -> Dict[str, Any]:
    """Compute every publication chemistry map from one shared protein analysis."""
    context = _prepare_shared_context(
        pH=pH,
        ionic_mM=ionic_mM,
        temp_K=temp_K,
        pdb_id=pdb_id,
        pdb_text=pdb_text,
        chain=chain,
        gnm_cutoff_A=gnm_cutoff_A,
    )

    available = context["core_result"].get("chemistries", {})
    order = [key for key in MAP_ORDER if key in available]
    order.extend(key for key in available if key not in order)
    maps = {chemistry: _build_map(context, chemistry) for chemistry in order}

    initial_map = order[0] if order else None
    if initial_surface:
        mode = get_surface_mode(initial_surface)
        if mode.chemistry in maps:
            initial_map = mode.chemistry

    out = _shared_response(
        context,
        pdb_id=pdb_id,
        pH=pH,
        ionic_mM=ionic_mM,
        temp_K=temp_K,
        gnm_cutoff_A=gnm_cutoff_A,
    )
    out.update({
        "n_maps": len(order),
        "map_order": order,
        "initial_map": initial_map,
        "maps": maps,
        "surface_presets": {
            key: {
                "label": mode.label,
                "chemistry": mode.chemistry,
                "description": mode.description,
            }
            for key, mode in sorted(SURFACE_MODES.items())
        },
    })
    return out


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
    """Compatibility wrapper returning one material-selected map for validation scripts."""
    mode = get_surface_mode(surface)
    context = _prepare_shared_context(
        pH=pH,
        ionic_mM=ionic_mM,
        temp_K=temp_K,
        pdb_id=pdb_id,
        pdb_text=pdb_text,
        chain=chain,
        gnm_cutoff_A=gnm_cutoff_A,
    )
    map_payload = _build_map(context, mode.chemistry)
    out = _shared_response(
        context,
        pdb_id=pdb_id,
        pH=pH,
        ionic_mM=ionic_mM,
        temp_K=temp_K,
        gnm_cutoff_A=gnm_cutoff_A,
    )
    out.update({
        "input": {
            **out["input"],
            "surface": mode.key,
            "surface_label": mode.label,
            "primary_chemistry": mode.chemistry,
        },
        "surface_mode": {
            "key": mode.key,
            "label": mode.label,
            "chemistry": mode.chemistry,
            "description": mode.description,
        },
        "n_patches": map_payload["n_patches"],
        "n_pareto_primary_patches": map_payload["n_pareto_primary_patches"],
        "primary_patches": map_payload["primary_patches"],
        "patches": map_payload["patches"],
    })
    return out
