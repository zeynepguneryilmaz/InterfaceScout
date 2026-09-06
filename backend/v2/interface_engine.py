"""InterfaceScout V2: coarse, weight-free protein-material interface prediction."""

from __future__ import annotations

import importlib
import urllib.request
from typing import Any, Dict, Optional

from .coarse_patch import build_coarse_patches, PATCH_SCALE_A
from .gnm import solve_gnm
from .prepare import prepare_pdb_text
from .surface_modes import get_surface_mode


def _load_v1():
    return importlib.import_module("main")


def _obtain_pdb_text(pdb_id: Optional[str], pdb_text: Optional[str]) -> str:
    if pdb_text:
        return pdb_text
    pid = (pdb_id or "").strip().upper()
    if not pid:
        raise ValueError("Provide pdb_id or pdb_text")
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


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

    # The input is already reduced to the requested chain subset. Passing
    # chain=None prevents the frozen V1 parser from interpreting a multi-chain
    # selector such as "C,E" as one literal chain ID.
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
    patches = build_coarse_patches(v1_result=v1_result, chemistry=mode.chemistry, gnm=gnm)
    pareto_primary = [p for p in patches if int(p.get("pareto_front", 999)) == 1]

    return {
        "engine": "InterfaceScout V2",
        "version": "2.2.1-coarse-prototype",
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
            "chemistry_source": "frozen InterfaceScout V1 compatibility channel",
            "accessibility_source": "V1 side-chain relative solvent accessibility",
            "patch_radius_A": PATCH_SCALE_A,
            "patch_radius_basis": "frozen V1 8 A patch scale; not adsorption-label fitted",
            "patch_definition": "non-transitive local surface neighbourhood around a V1 chemistry-patch maximum",
            "orientation": "coarse outward C-alpha face consistency; descriptive",
            "dynamics": "unweighted C-alpha GNM; descriptive",
            "gnm_cutoff_A": float(gnm_cutoff_A),
            "ranking": "Pareto fronts across chemistry support, accessibility, patch coherence, dynamic coupling, and orientation coherence; no weighted sum",
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
        "diagnostics": {
            "n_surface_residues": len(v1_result.get("surface_residues", [])),
            "surface_residue_keys": [str(r["key"]) for r in v1_result.get("surface_residues", []) if r.get("key")],
        },
        "method_notes": [
            "Experimental interface labels are not inputs to patch construction or ranking.",
            "Patch membership is intentionally coarse; individual residues are not claimed as precise adsorption contacts.",
            "Patch growth is non-transitive to prevent surface percolation into unrealistically large regions.",
            "GNM and orientation provide physical context without benchmark-tuned membership thresholds.",
            "Multiple Pareto-optimal patches are allowed because protein adsorption may have alternative plausible interfaces.",
        ],
    }
