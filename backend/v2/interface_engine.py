"""InterfaceScout protein-side surface-chemistry mapping engine.

A single analysis prepares the protein and computes solvent exposure once, then
returns every canonical surface-chemistry map plus optional structural-context
maps. Material names are not part of the public model. Structural-context maps
never enter patch construction or ranking.
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
from .surface_modes import get_surface_mode

# Public labels describe generalized surface chemistry rather than named materials.
# The second item is the frozen internal chemistry key retained for reproducibility.
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


def _geometry_only_gnm(prepared: str, cutoff_A: float) -> dict:
    nodes = extract_ca_nodes(prepared)
    n = len(nodes)
    return {
        "cutoff_A": float(cutoff_A),
        "n_nodes": n,
        "nodes": nodes,
        "index": {x["key"]: i for i, x in enumerate(nodes)},
        "correlation_matrix": np.zeros((n, n), dtype=float),
        "residue_metrics": {},
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
        rin = {"cutoff_A": None, "n_nodes": 0, "n_edges": 0, "residue_metrics": {}}
    else:
        rin = annotate_rin_percentiles(build_rin(prepared), surface_keys)

    return {
        "prepared": prepared,
        "prep_report": prep_report,
        "core_result": core_result,
        "canonical_only": canonical_only,
        "gnm": gnm,
        "rin": rin,
        "surface_keys": surface_keys,
    }


def _residue_rows(channel: dict) -> list[dict]:
    """Expose the residue-level chemistry fields needed for inspection/overlay."""
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
    patches = build_coarse_patches(v1_result=core_result, chemistry=internal_key, gnm=context["gnm"])

    if context["canonical_only"]:
        for patch in patches:
            patch["rin_context"] = {"status": "disabled_in_canonical_only_validation"}
    else:
        for patch in patches:
            patch["rin_context"] = summarize_patch_rin(patch.get("members", []), context["rin"])

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


def _range(rows: list[dict]) -> tuple[float | None, float | None]:
    vals = [float(r["value"]) for r in rows if r.get("value") is not None and np.isfinite(float(r["value"]))]
    return (min(vals), max(vals)) if vals else (None, None)


def _build_property_maps(context: dict) -> dict:
    if context["canonical_only"]:
        return {
            "available": False,
            "effect_on_prediction": False,
            "reason": "Structural-context calculations are disabled in canonical-only validation mode.",
            "order": [],
            "maps": {},
        }

    surface = set(context["surface_keys"])
    gnm_metrics = context["gnm"].get("residue_metrics", {})
    rin_metrics = context["rin"].get("residue_metrics", {})

    gnm_rows = [{
        "key": key,
        "chain": m["chain"],
        "res_seq": int(m["res_seq"]),
        "icode": m.get("icode", ""),
        "res_name": m["res_name"],
        "value": float(m["normalized_fluctuation"]),
        "surface": key in surface,
    } for key, m in gnm_metrics.items()]

    def rin_rows(field: str, percentile_field: str) -> list[dict]:
        return [{
            "key": key,
            "chain": m["chain"],
            "res_seq": int(m["res_seq"]),
            "icode": m.get("icode", ""),
            "res_name": m["res_name"],
            "value": float(m[field]),
            "surface": key in surface,
            "surface_percentile": m.get(percentile_field),
        } for key, m in rin_metrics.items()]

    definitions = {
        "gnm_fluctuation": {
            "label": "GNM fluctuation",
            "description": "Normalized native-state Cα fluctuation. Values >1 indicate above-average mobility within the analyzed protein.",
            "value_label": "normalized fluctuation",
            "rows": gnm_rows,
        },
        "rin_degree": {
            "label": "RIN degree",
            "description": "Normalized residue degree in the 4.5 Å heavy-atom contact network; higher values indicate more direct structural contacts.",
            "value_label": "normalized degree",
            "rows": rin_rows("degree_normalized", "degree_normalized_percentile_surface"),
        },
        "rin_betweenness": {
            "label": "RIN betweenness",
            "description": "Normalized betweenness centrality; higher values indicate residues lying on more shortest network paths.",
            "value_label": "betweenness",
            "rows": rin_rows("betweenness", "betweenness_percentile_surface"),
        },
        "rin_closeness": {
            "label": "RIN closeness",
            "description": "Closeness centrality; higher values indicate shorter network distance to the rest of the protein.",
            "value_label": "closeness",
            "rows": rin_rows("closeness", "closeness_percentile_surface"),
        },
    }
    for item in definitions.values():
        lo, hi = _range(item["rows"])
        item["min_value"] = lo
        item["max_value"] = hi
        item["effect_on_prediction"] = False
        item["role"] = "descriptive structural context only"

    return {
        "available": True,
        "effect_on_prediction": False,
        "order": ["gnm_fluctuation", "rin_degree", "rin_betweenness", "rin_closeness"],
        "maps": definitions,
        "gnm": {"cutoff_A": context["gnm"]["cutoff_A"]},
        "rin": {"cutoff_A": context["rin"]["cutoff_A"]},
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
            "gnm_cutoff_A": float(gnm_cutoff_A),
            "rin_heavy_atom_cutoff_A": context["rin"]["cutoff_A"],
            "single_pass_all_maps": True,
        },
        "diagnostics": {
            "n_surface_residues": len(context["surface_keys"]),
            "surface_residue_keys": context["surface_keys"],
            "canonical_only_validation": context["canonical_only"],
        },
        "method_notes": [
            "All single-channel chemistry maps use the same prepared structure and environmental conditions.",
            "Experimental interface labels are not inputs to patch construction or ranking.",
            "GNM and RIN are excluded from patch prediction/ranking and retained only as structural context.",
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
    gnm_cutoff_A: float = 7.3,
) -> Dict[str, Any]:
    """Compute all public chemistry maps and structural-context maps in one run."""
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
    specs = [spec for spec in PUBLIC_MAP_SPECS if spec[1] in available]
    maps = {public: _build_map(context, public, internal, label, desc) for public, internal, label, desc in specs}
    order = [spec[0] for spec in specs]

    out = _shared_response(context, pdb_id=pdb_id, pH=pH, ionic_mM=ionic_mM, temp_K=temp_K, gnm_cutoff_A=gnm_cutoff_A)
    out.update({
        "n_maps": len(order),
        "map_order": order,
        "initial_map": order[0] if order else None,
        "maps": maps,
        "structural_context": _build_property_maps(context),
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
    """Internal compatibility wrapper retained for frozen validation scripts."""
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
    internal = mode.chemistry
    spec = next((x for x in PUBLIC_MAP_SPECS if x[1] == internal), (internal, internal, internal.replace("_", " ").title(), ""))
    map_payload = _build_map(context, *spec)
    out = _shared_response(context, pdb_id=pdb_id, pH=pH, ionic_mM=ionic_mM, temp_K=temp_K, gnm_cutoff_A=gnm_cutoff_A)
    out.update({
        "input": {**out["input"], "surface": mode.key, "primary_chemistry": internal},
        "n_patches": map_payload["n_patches"],
        "n_pareto_primary_patches": map_payload["n_pareto_primary_patches"],
        "primary_patches": map_payload["primary_patches"],
        "patches": map_payload["patches"],
    })
    return out
