"""InterfaceScout V2 dynamic engine: frozen V1 anchors + native-state GNM context."""

from __future__ import annotations

import importlib
import urllib.request
from typing import Any, Dict, Optional

from .dynamic_patch import analyze_anchor_dynamics
from .gnm import solve_gnm
from .prepare import prepare_pdb_text


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


def analyze_v2_dynamic(
    *,
    chemistry: str,
    pH: float = 7.4,
    ionic_mM: float = 150.0,
    temp_K: float = 298.0,
    pdb_id: Optional[str] = None,
    pdb_text: Optional[str] = None,
    chain: Optional[str] = None,
    top_n_anchors: int = 10,
    gnm_cutoff_A: float = 7.3,
    extension_radius_A: float = 12.0,
    min_abs_corr: float = 0.35,
    max_extensions: int = 12,
) -> Dict[str, Any]:
    """Run frozen V1 chemistry scoring, then annotate top anchors with GNM dynamics.

    GNM does not alter the V1 ranking. This separation is deliberate so that the
    added dynamics can be validated independently before any composite model is
    considered.
    """
    raw = _obtain_pdb_text(pdb_id, pdb_text)
    prepared, prep_report = prepare_pdb_text(raw, chain=chain)

    v1 = _load_v1()
    request = v1.AnalyzeRequest(
        pdb_text=prepared,
        chain=chain,
        env=v1.EnvParams(pH=float(pH), ionic=float(ionic_mM), temp=float(temp_K)),
    )
    v1_result = v1.analyze(request)
    chemistries = v1_result.get("chemistries", {})
    if chemistry not in chemistries:
        raise ValueError(f"Unknown V1 chemistry {chemistry!r}. Available: {', '.join(sorted(chemistries))}")

    rows = sorted(
        chemistries[chemistry].get("patch_centers", []),
        key=lambda r: float(r.get("multiscale_persistence", 0.0)),
        reverse=True,
    )[: int(top_n_anchors)]

    gnm = solve_gnm(prepared, cutoff_A=float(gnm_cutoff_A))
    dynamic = analyze_anchor_dynamics(
        anchor_rows=rows,
        gnm=gnm,
        v1_result=v1_result,
        radius_A=float(extension_radius_A),
        min_abs_corr=float(min_abs_corr),
        max_extensions=int(max_extensions),
    )

    return {
        "engine": "InterfaceScout V2-dynamic-prototype",
        "version": "2.1.0-prototype",
        "scope": {
            "v1_anchor_ranking_modified": False,
            "gnm_predicts_adsorption_energy": False,
            "gnm_creates_new_anchors": False,
            "dynamic_layer_role": "native-state context and anchor-neighborhood expansion",
        },
        "input": {
            "pdb_id": (pdb_id or "").strip().upper() or None,
            "chain": (chain or "").strip() or "ALL",
            "chemistry": chemistry,
            "pH": float(pH),
            "ionic_mM": float(ionic_mM),
            "temp_K": float(temp_K),
        },
        "structure_preparation": prep_report,
        "gnm": {
            "cutoff_A": gnm["cutoff_A"],
            "n_nodes": gnm["n_nodes"],
            "n_zero_modes": gnm["n_zero_modes"],
            "n_nonzero_modes": gnm["n_nonzero_modes"],
        },
        "parameters": {
            "top_n_anchors": int(top_n_anchors),
            "extension_radius_A": float(extension_radius_A),
            "min_abs_corr": float(min_abs_corr),
            "max_extensions": int(max_extensions),
        },
        "v1_top_anchors": rows,
        "dynamic_anchor_analysis": dynamic,
        "method_notes": [
            "V1 chemistry ranking is kept frozen and unchanged.",
            "GNM uses a standard unweighted C-alpha Kirchhoff network.",
            "Normalized fluctuation >1 denotes above-average native-state mobility.",
            "Dynamic extensions must be V1-solvent-exposed, spatially near the anchor, and sufficiently GNM-correlated.",
            "No benchmark-fitted coefficient or V1+GNM weighted score is used in this prototype.",
        ],
    }
