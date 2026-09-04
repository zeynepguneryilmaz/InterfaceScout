"""InterfaceScout 2.0 development application.

InterfaceScout 2.0 retains the frozen 1.0 canonical compatibility core for
backward-comparable residue/patch maps and adds structural-context and
surface-type-specific physics layers under separate, explicitly qualified
outputs.
"""
from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse

try:
    from . import main as core
    from .structural_context import AnalyzeRequest, analyze_structural, prepare_context, available_material_profiles
    from .physics_refinement import enrich_nonpolar_physics
    from .nonpolar_sasa_orientation import scan as scan_nonpolar_sasa
except ImportError:
    import main as core
    from structural_context import AnalyzeRequest, analyze_structural, prepare_context, available_material_profiles
    from physics_refinement import enrich_nonpolar_physics
    from nonpolar_sasa_orientation import scan as scan_nonpolar_sasa

APP_VERSION = "2.0.0-dev"
CORE_RELEASE_VERSION = "1.0.0"


def _atom_level_nonpolar_orientation(req: AnalyzeRequest):
    workdir = Path(tempfile.mkdtemp(prefix="interfacescout_v2_nonpolar_"))
    try:
        pdb, _, _, _ = prepare_context(req, workdir)
        struct, _, _, _ = core.build_surface_residues(pdb, req.env.pH)
        return scan_nonpolar_sasa(struct)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def analyze(req: AnalyzeRequest):
    result = dict(analyze_structural(req))
    result["version"] = APP_VERSION
    result["core_version"] = CORE_RELEASE_VERSION
    result["model"] = "InterfaceScout 2.0 development"

    nonpolar = enrich_nonpolar_physics(
        result.get("surface_residues", []), result.get("all_residues", [])
    )
    nonpolar["atom_level_sasa_orientation"] = _atom_level_nonpolar_orientation(req)
    result["nonpolar_physics"] = nonpolar

    result.setdefault("applicability", {}).setdefault("included_in_this_run", []).extend([
        "continuous Eisenberg hydrophobic surface field",
        "Wimley-White interfacial preference sensitivity descriptor",
        "exposure-weighted tertiary hydrophobic vector",
        "atom-level SASA nonpolar orientation descriptor",
    ])
    result["applicability"].setdefault("not_included_or_interpretation_limits", []).extend([
        "the current atom-level nonpolar orientation layer contains the published SASA-solvation component but not yet the full van der Waals term",
        "explicit molecular water and adsorption-induced large conformational rearrangement remain outside the lightweight model",
    ])
    return result


app = FastAPI(title="InterfaceScout", version=APP_VERSION)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return {
        "status": "ok",
        "version": APP_VERSION,
        "release_status": "development",
        "frozen_reference": CORE_RELEASE_VERSION,
        "canonical_v1_score_changed": False,
        "structural_context": "biological assembly aware",
        "nonpolar_physics": "continuous hydrophobic fields + hydrophobic vector + SASA orientation component",
    }


@app.get("/material_profiles")
def material_profiles():
    return {"profiles": available_material_profiles(), "combination_rule": "none; channels remain separate"}


@app.post("/analyze_surface")
def analyze_surface(req: AnalyzeRequest):
    try:
        return analyze(req)
    except HTTPException:
        raise
    except Exception as exc:
        core.log.exception("InterfaceScout 2.0 analysis failed")
        raise HTTPException(500, str(exc))


@app.get("/model_spec")
def model_spec():
    return {
        "version": APP_VERSION,
        "status": "development",
        "frozen_reference": CORE_RELEASE_VERSION,
        "canonical_v1_core_unchanged": True,
        "structure_context_modes": ["auto", "biological_assembly_1", "deposited_structure", "selected_chain_legacy"],
        "nonpolar_layer": {
            "continuous_scales": ["Eisenberg", "Wimley-White interface"],
            "orientation_descriptors": ["tertiary hydrophobic vector", "atom-level SASA orientation component"],
            "vdw_term_complete": False,
        },
    }


@app.get("/")
def root():
    path = core.PROJECT_DIR / "frontend" / "index.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding="utf-8"))
    return JSONResponse({"name": "InterfaceScout", "version": APP_VERSION, "frontend": "not found"}, status_code=404)


@app.get("/favicon.png")
def favicon():
    path = core.PROJECT_DIR / "frontend" / "favicon.png"
    if not path.exists():
        raise HTTPException(404, "favicon.png not found")
    return FileResponse(path)


@app.get("/logo.png")
def logo():
    path = core.PROJECT_DIR / "frontend" / "logo.png"
    if not path.exists():
        raise HTTPException(404, "logo.png not found")
    return FileResponse(path)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
