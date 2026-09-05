"""InterfaceScout V2-alpha analysis engine."""

from __future__ import annotations

import importlib
import urllib.request
from typing import Any, Dict, Optional

from .patches import assemble_candidate_patches, nonredundant_top_patches
from .prepare import prepare_pdb_text
from .scoring import rank_patches, summarize_competition
from .surface_profiles import get_surface_profile


def _load_v1():
    # backend/main.py remains the publication-frozen V1 implementation.
    return importlib.import_module("main")


def _obtain_pdb_text(pdb_id: Optional[str], pdb_text: Optional[str]) -> str:
    if pdb_text:
        return pdb_text
    pid = (pdb_id or "").strip().upper()
    if not pid:
        raise ValueError("Provide pdb_id or pdb_text")
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


def analyze_v2(
    *,
    surface: str,
    pH: float = 7.4,
    ionic_mM: float = 150.0,
    temp_K: float = 298.0,
    pdb_id: Optional[str] = None,
    pdb_text: Optional[str] = None,
    chain: Optional[str] = None,
    top_n: int = 3,
) -> Dict[str, Any]:
    """Run V2-alpha patch-level material-aware analysis.

    V2-alpha deliberately reuses the frozen V1 residue chemistry maps as input
    descriptors.  It does not alter V1 and it does not fit any coefficients.
    """
    if not (0.0 <= float(pH) <= 14.0):
        raise ValueError("pH must be between 0 and 14")
    if float(ionic_mM) < 0:
        raise ValueError("ionic_mM must be non-negative")
    if float(temp_K) <= 0:
        raise ValueError("temp_K must be positive")
    if int(top_n) < 1 or int(top_n) > 20:
        raise ValueError("top_n must be between 1 and 20")

    profile = get_surface_profile(surface)
    weights = profile.normalized_weights()
    raw_text = _obtain_pdb_text(pdb_id, pdb_text)
    prepared_text, prep_report = prepare_pdb_text(raw_text, chain=chain)

    v1 = _load_v1()
    request = v1.AnalyzeRequest(
        pdb_text=prepared_text,
        chain=chain,
        env=v1.EnvParams(pH=float(pH), ionic=float(ionic_mM), temp=float(temp_K)),
    )
    v1_result = v1.analyze(request)
    if hasattr(v1_result, "body"):
        raise RuntimeError("Unexpected HTTP response object returned by V1 analyze()")
    chemistries = v1_result.get("chemistries", {})

    candidates = assemble_candidate_patches(chemistries, weights)
    ranked = rank_patches(candidates)
    top = nonredundant_top_patches(ranked, top_n=int(top_n))
    # Re-label retained patch ranks after redundancy filtering while preserving raw rank.
    for final_rank, row in enumerate(top, start=1):
        row["raw_rank"] = row.get("rank")
        row["rank"] = final_rank

    return {
        "engine": "InterfaceScout V2-alpha",
        "version": "2.0.0-alpha.1",
        "scope": {
            "prediction_unit": "protein surface patch",
            "orientation_search": False,
            "adsorption_free_energy": False,
            "desolvation_model": False,
            "fitted_to_external_benchmark": False,
        },
        "input": {
            "pdb_id": (pdb_id or "").strip().upper() or None,
            "chain": (chain or "").strip() or "ALL",
            "pH": float(pH),
            "ionic_mM": float(ionic_mM),
            "temp_K": float(temp_K),
            "surface_profile": profile.key,
        },
        "structure_preparation": prep_report,
        "surface_profile": {
            "key": profile.key,
            "label": profile.label,
            "description": profile.description,
            "normalized_weights": weights,
        },
        "n_candidate_patch_centers": len(candidates),
        "top_patches": top,
        "competition": summarize_competition(ranked),
        "method_notes": [
            "V1 publication-freeze chemistry maps are reused without modification.",
            "Patch score is the normalized surface-profile-weighted mean of V1 5/8 Å multiscale persistence.",
            "Near-duplicate Top patches are suppressed by compatible-member-set overlap.",
            "V2-alpha is a patch-localization prototype, not an adsorption free-energy or rigid-body docking model.",
        ],
    }
