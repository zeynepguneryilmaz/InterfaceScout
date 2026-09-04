"""InterfaceScout v5.3 development wrapper.

Runs the v5.2 structural-context candidate unchanged, then adds independent
established-physics descriptors for nonpolar interfaces. The frozen canonical
v5.1 chemistry/state and C-alpha 5/8 A scores remain available and unchanged.
"""
from __future__ import annotations

import shutil
import tempfile
from pathlib import Path
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

try:
    from . import main as core
    from .v52_app import AnalyzeRequest, analyze as analyze_v52, _prepare_context
    from .physics_refinement import enrich_nonpolar_physics
    from .nonpolar_sasa_orientation import scan as scan_nonpolar_sasa
except ImportError:
    import main as core
    from v52_app import AnalyzeRequest, analyze as analyze_v52, _prepare_context
    from physics_refinement import enrich_nonpolar_physics
    from nonpolar_sasa_orientation import scan as scan_nonpolar_sasa

APP_VERSION = "5.3.0-established-physics-development"


def _atom_level_orientation(req: AnalyzeRequest):
    """Rebuild the selected structural context to retain atom-level SASA."""
    workdir = Path(tempfile.mkdtemp(prefix="interfacescout_v53_atom_"))
    try:
        pdb, _, _, _ = _prepare_context(req, workdir)
        struct, _, _, _ = core.build_surface_residues(pdb, req.env.pH)
        return scan_nonpolar_sasa(struct)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def analyze(req: AnalyzeRequest):
    result = dict(analyze_v52(req))
    parent_version = result.get("version")
    result["version"] = APP_VERSION
    result["parent_structural_version"] = parent_version
    result["nonpolar_physics"] = enrich_nonpolar_physics(
        result.get("surface_residues", []), result.get("all_residues", [])
    )
    result["nonpolar_physics"]["atom_level_sasa_orientation"] = _atom_level_orientation(req)
    result.setdefault("applicability", {}).setdefault("included_in_this_run", []).extend([
        "continuous Eisenberg hydrophobic surface field",
        "Wimley-White interfacial-scale sensitivity descriptor",
        "exposure-weighted tertiary hydrophobic vector / preferred-hemisphere descriptor",
        "Harrison-style atom-level SASA nonpolar orientation descriptor",
    ])
    result["applicability"].setdefault("not_included_or_interpretation_limits", []).extend([
        "the atom-level SASA orientation layer implements the published solvation component but not the CHARMM van der Waals term or full Metropolis trajectory",
        "explicit molecular water and adsorption-induced large conformational rearrangement remain outside the lightweight model",
    ])
    return result


app = FastAPI(title="InterfaceScout v5.3 Development", version=APP_VERSION)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return {
        "status": "ok",
        "version": APP_VERSION,
        "canonical_v51_changed": False,
        "new_nonpolar_model": "continuous hydrophobic field + tertiary hydrophobic vector + atom-level SASA orientation; no fitted weighted sum",
    }


@app.post("/analyze_surface")
def analyze_surface(req: AnalyzeRequest):
    try:
        return analyze(req)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(500, str(exc))
