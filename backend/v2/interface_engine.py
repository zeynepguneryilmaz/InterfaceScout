"""InterfaceScout protein-side surface-chemistry mapping engine.

A single analysis prepares the protein and computes solvent exposure once, then
returns every canonical surface-chemistry map. Material names, GNM dynamics and
residue-interaction-network descriptors are not part of the public model.
"""
from __future__ import annotations

import importlib
import urllib.request
from typing import Any, Dict, Optional

from .chemistry_freeze import apply_publication_chemistry
from .coarse_patch import build_coarse_patches, PATCH_SCALE_A
from .geometry import extract_ca_nodes
from .model_settings import (
    MODEL_VERSION,
    SASA_POINTS,
    SASA_PROBE_A,
    SC_RSA_THRESHOLD,
    MULTISCALE_RADII_A,
    PARAMETER_SELECTION,
)
from .prepare import prepare_pdb_text
from .surface_modes import get_surface_mode

PUBLIC_MAP_SPECS = [
    ("anionic", "anionic", "Anionic surface", "Negatively charged surface chemistry, including exposed anionic or carboxylate-rich groups."),
    ("cationic", "cationic", "Cationic surface", "Positively charged surface chemistry, including protonated amine-rich groups."),
    ("hydrophobic", "hydrophobic", "Hydrophobic surface", "Nonpolar surface chemistry supporting hydrophobic contacts."),
    ("pi_carbon", "pi_carbon", "π / aromatic surface", "Aromatic or graphitic-like surface chemistry supporting π-associated and cation–π contacts."),
    ("hbond_donor", "hbond_donor", "H-bond donor surface", "Surface groups capable of donating hydrogen bonds to exposed protein side chains."),
    ("hbond_acceptor", "hbond_acceptor", "H-bond acceptor surface", "Surface groups capable of accepting hydrogen bonds from exposed protein side chains."),
    ("oxide", "oxide", "Metal-oxide surface", "Oxide-like surface chemistry represented by the frozen carboxylate-compatible oxide channel."),
    ("calcium_phosphate", "hydroxyapatite", "Calcium/phosphate charged sites", "Calcium- and phosphate-rich charged surface sites with complementary protein-side interactions."),
    ("metal_coord", "metal_coord", "Transition-metal coordination", "Accessible transition-metal sites supporting coordination by exposed protein side chains."),
    ("soft_metal_sulfur", "gold", "Soft-metal sulfur affinity", "Soft-metal-like surface affinity dominated by accessible sulfur-containing side chains."),
    ("phosphate", "phosphate", "Phosphate-rich surface", "Phosphate-rich surface chemistry supporting electrostatic and hydrogen-bond interactions."),
]

PUBLIC_TO_INTERNAL = {public: internal for public, internal, _label, _description in PUBLIC_MAP_SPECS}


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


def _prepare_shared_context(
    *,
    pH: float,
    ionic_mM: float,
    temp_K: float,
    pdb_id: Optional[str],
    pdb_text: Optional[str],
    chain: Optional[str],
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

    ca_nodes = extract_ca_nodes(prepared)
    surface_keys = [str(r["key"]) for r in core_result.get("surface_residues", []) if r.get("key")]
    return {
        "prepared": prepared,
        "prep_report": prep_report,
        "core_result": core_result,
        "ca_nodes": ca_nodes,
        "surface_keys": surface_keys,
    }


def _residue_rows(channel: dict) -> list[dict]:
    rows = []
    for row in channel.get("residues", []):
        key = row.get("key")
        if not key:
            continue
        rows.append({
            "key": str(key),
            "chain": row.get("chain"),
            "res_seq": row.get("res_seq"),
            "icode": row.get("icode", ""),
            "res_name": row.get("res_name"),
            "scrsa": row.get("scrsa", row.get("scrsa_raw")),
            "local_score": row.get("local_score"),
            "propensity": row.get("propensity"),
            "multiscale_persistence": row.get("multiscale_persistence"),
        })
    return rows


def _build_map(context: dict, public_key: str, internal_key: str, label: str, description: str) -> dict:
    core_result = context["core_result"]
    patches = build_coarse_patches(
        v1_result=core_result,
        chemistry=internal_key,
        ca_nodes=context["ca_nodes"],
    )
    channel = core_result.get("chemistries", {}).get(internal_key, {})
    primary = [p for p in patches if int(p.get("pareto_front", 999)) == 1]
    return {
        "key": public_key,
        "label": label,
        "description": description,
        "internal_chemistry_key": internal_key,
        "n_patches": len(patches),
        "n_pareto_primary_patches": len(primary),
        "primary_patches": primary,
        "patches": patches,
        "residues": _residue_rows(channel),
    }


def _shared_response(context: dict, *, pdb_id: Optional[str], pH: float, ionic_mM: float, temp_K: float) -> dict:
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
            "accessibility_source": "side-chain relative solvent accessibility",
            "sasa_probe_A": SASA_PROBE_A,
            "sasa_points_per_atom": SASA_POINTS,
            "scrsa_threshold": SC_RSA_THRESHOLD,
            "multiscale_radii_A": list(MULTISCALE_RADII_A),
            "parameter_selection": PARAMETER_SELECTION,
            "patch_radius_A": PATCH_SCALE_A,
            "patch_definition": "non-transitive local surface neighbourhood around a chemistry-patch maximum",
            "orientation": "coarse outward C-alpha face consistency",
            "ranking": "Pareto fronts across chemistry support, accessibility, patch coherence and orientation coherence; no weighted sum",
            "single_pass_all_maps": True,
        },
        "diagnostics": {
            "n_surface_residues": len(context["surface_keys"]),
            "surface_residue_keys": context["surface_keys"],
        },
        "method_notes": [
            "All single-channel chemistry maps use the same prepared structure and environmental conditions.",
            "Experimental interface labels are not inputs to patch construction or ranking.",
            "Multi-chemistry composite views are derived overlays of canonical single-channel maps; they do not introduce fitted material weights or constitute a separately validated combined score.",
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
) -> Dict[str, Any]:
    """Compute all public chemistry maps in one run."""
    context = _prepare_shared_context(
        pH=pH,
        ionic_mM=ionic_mM,
        temp_K=temp_K,
        pdb_id=pdb_id,
        pdb_text=pdb_text,
        chain=chain,
    )

    available = context["core_result"].get("chemistries", {})
    specs = [spec for spec in PUBLIC_MAP_SPECS if spec[1] in available]
    maps = {public: _build_map(context, public, internal, label, desc) for public, internal, label, desc in specs}
    order = [spec[0] for spec in specs]

    out = _shared_response(context, pdb_id=pdb_id, pH=pH, ionic_mM=ionic_mM, temp_K=temp_K)
    out.update({
        "n_maps": len(order),
        "map_order": order,
        "initial_map": order[0] if order else None,
        "maps": maps,
        "composite_policy": {
            "type": "derived_overlay",
            "validated_as_combined_model": False,
            "default_residue_rule": "maximum normalized support across selected canonical chemistry maps",
            "double_counting": False,
            "weights": "none unless externally known surface composition is supplied in a future explicit model extension",
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
    """Internal compatibility wrapper retained for frozen validation scripts.

    gnm_cutoff_A remains accepted only so archived validation calls do not break;
    it has no effect on the current model.
    """
    _ = gnm_cutoff_A
    mode = get_surface_mode(surface)
    context = _prepare_shared_context(
        pH=pH,
        ionic_mM=ionic_mM,
        temp_K=temp_K,
        pdb_id=pdb_id,
        pdb_text=pdb_text,
        chain=chain,
    )
    internal = mode.chemistry
    spec = next((x for x in PUBLIC_MAP_SPECS if x[1] == internal), (internal, internal, internal.replace("_", " ").title(), ""))
    map_payload = _build_map(context, *spec)
    out = _shared_response(context, pdb_id=pdb_id, pH=pH, ionic_mM=ionic_mM, temp_K=temp_K)
    out.update({
        "input": {**out["input"], "surface": mode.key, "primary_chemistry": internal},
        "n_patches": map_payload["n_patches"],
        "n_pareto_primary_patches": map_payload["n_pareto_primary_patches"],
        "primary_patches": map_payload["primary_patches"],
        "patches": map_payload["patches"],
    })
    return out
