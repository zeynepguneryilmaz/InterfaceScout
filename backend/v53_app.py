"""InterfaceScout v5.3 development wrapper.

Runs the v5.2 structural-context candidate unchanged, then adds independent
established-physics descriptors for nonpolar interfaces.  The frozen canonical
v5.1 chemistry/state and C-alpha 5/8 A scores remain available and unchanged.
"""
from __future__ import annotations

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

try:
    from .v52_app import AnalyzeRequest, analyze as analyze_v52
    from .physics_refinement import enrich_nonpolar_physics
except ImportError:
    from v52_app import AnalyzeRequest, analyze as analyze_v52
    from physics_refinement import enrich_nonpolar_physics

APP_VERSION = "5.3.0-established-physics-development"


def analyze(req: AnalyzeRequest):
    result = analyze_v52(req)
    result = dict(result)
    result["version"] = APP_VERSION
    result["parent_structural_version"] = result.get("version")
    result["nonpolar_physics"] = enrich_nonpolar_physics(result.get("surface_residues", []))
    result.setdefault("applicability", {}).setdefault("included_in_this_run", []).extend([
        "continuous Eisenberg hydrophobic surface field",
        "Wimley-White interfacial-scale sensitivity descriptor",
        "3-D hydrophobic dipole / preferred-hemisphere descriptor",
    ])
    result["applicability"].setdefault("not_included_or_interpretation_limits", []).append(
        "explicit molecular water and adsorption-induced large conformational rearrangement remain outside the lightweight model"
    )
    return result


app = FastAPI(title="InterfaceScout v5.3 Development", version=APP_VERSION)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return {
        "status": "ok",
        "version": APP_VERSION,
        "canonical_v51_changed": False,
        "new_nonpolar_model": "continuous hydrophobic field + hydrophobic dipole; no fitted weighted sum",
    }


@app.post("/analyze_surface")
def analyze_surface(req: AnalyzeRequest):
    try:
        return analyze(req)
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(500, str(exc))
