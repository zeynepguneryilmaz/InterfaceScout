"""InterfaceScout 2.0 development application.

InterfaceScout 2.0 is protein-centered. The program derives a target interface
profile from the supplied protein structure and reports which protein patches
could engage each generalized surface property. It does not use a named-material
library and it does not recommend materials by name.

The frozen InterfaceScout 1.0 residue/patch equations remain the reference
compatibility core. Material-specific physics modules are retained only as
validation/research modules and are not run by the default application path.
"""
from __future__ import annotations

from typing import List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse

try:
    from . import main as core
    from .structural_context import AnalyzeRequest, analyze_structural
    from .functional_sites import normalize_protected_residues, pdb_site_annotations
    from .target_profile import build_target_interface_profile
    from .ui import inject_ui
except ImportError:
    import main as core
    from structural_context import AnalyzeRequest, analyze_structural
    from functional_sites import normalize_protected_residues, pdb_site_annotations
    from target_profile import build_target_interface_profile
    from ui import inject_ui

APP_VERSION = "2.0.0-dev"
CORE_RELEASE_VERSION = "1.0.0"


class ProteinAnalyzeRequest(AnalyzeRequest):
    """Default v2 request: protein structure plus optional protected residues."""
    protected_residue_keys: Optional[List[str]] = None


def analyze(req: ProteinAnalyzeRequest):
    result = dict(analyze_structural(req))
    result["version"] = APP_VERSION
    result["core_version"] = CORE_RELEASE_VERSION
    result["model"] = "InterfaceScout 2.0 protein-centered development"

    protected = normalize_protected_residues(req.protected_residue_keys)
    sites = pdb_site_annotations(req.pdb_id, req.pdb_text)
    target = build_target_interface_profile(
        chemistries=result.get("chemistries", {}),
        surface_residues=result.get("surface_residues", []),
        all_residues=result.get("all_residues", []),
        site_annotations=sites,
        protected_residue_keys=protected,
    )
    result["protein_derived_target_interface_profile"] = target

    # Named-material profiles are intentionally suppressed from the primary v2
    # output. Any legacy field that survives in the structural regression layer
    # is neutralized here.
    result["material_profile"] = None
    result["available_material_profiles"] = []
    if "settings" in result:
        result["settings"]["material_profile"] = None
        result["settings"]["protected_residue_keys"] = protected

    app = result.setdefault("applicability", {})
    included = app.setdefault("included_in_this_run", [])
    included.extend([
        "protein-derived target interface profile",
        "patch-resolved generalized surface-property requirements",
        "PDB SITE proximity/overlap annotation when present",
        "user-specified protected functional-residue overlap when provided",
    ])
    limits = app.setdefault("not_included_or_interpretation_limits", [])
    limits.extend([
        "no named-material library or named-material recommendation is used",
        "PDB SITE records are generic structural annotations and are not automatically interpreted as catalytic active sites",
        "material-specific adsorption physics is not part of the default protein-only prediction path",
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
        "primary_output": "protein-derived target interface profile",
        "material_library": False,
        "named_material_recommendation": False,
    }


@app.post("/analyze_surface")
def analyze_surface(req: ProteinAnalyzeRequest):
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
        "primary_direction": "protein -> target interface properties -> patch-level interpretation",
        "material_library": False,
        "named_material_recommendation": False,
        "target_profile_channels": list(core.CHEMISTRIES.keys()),
        "functional_site_annotation": {
            "PDB_SITE": "generic annotation only; not automatically catalytic",
            "user_protected_residues": "supported; used only to annotate patch overlap/proximity",
            "changes_compatibility_score": False,
        },
        "material_specific_physics": "research/validation modules only; excluded from default protein-only analysis",
    }


@app.get("/")
def root():
    path = core.PROJECT_DIR / "frontend" / "index.html"
    if path.exists():
        return HTMLResponse(inject_ui(path.read_text(encoding="utf-8")))
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
