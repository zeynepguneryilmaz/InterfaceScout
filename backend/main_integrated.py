"""InterfaceScout vNext integrated backend entrypoint.

Run with:
    uvicorn backend.main_integrated:app

The public score is named `interfacescout_score` and includes the legacy core,
GNM mobility, conditional APBS electrostatic complementarity, and radial prominence.
The legacy core is retained only as an ablation/provenance component.
"""
from __future__ import annotations

from pathlib import Path
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse

from backend import main as core
from backend.is_integrated_score import integrate_interfacescout_score, GNM_CUTOFF_A

APP_VERSION = "5.2.0-integrated-candidate"


def analyze(req: core.AnalyzeRequest):
    result = core.analyze(req)
    all_residues = result["all_residues"]
    surface_residues = result["surface_residues"]
    for key, chemistry_result in result["chemistries"].items():
        expected = core.CHEMISTRIES[key].get("expected_phi_sign")
        integrate_interfacescout_score(all_residues, surface_residues, chemistry_result, expected)

    result["version"] = APP_VERSION
    result["model"] = "InterfaceScout integrated interface-propensity score"
    result["score_name"] = "InterfaceScout score"
    result["score_field"] = "interfacescout_score"
    result["score_components"] = [
        "surface exposure (side-chain scRSA)",
        "pH-dependent residue state",
        "material-specific chemistry compatibility",
        "5/8 A multiscale patch persistence",
        "GNM intrinsic mobility",
        "conditional APBS electrostatic complementarity",
        "radial prominence",
    ]
    result["settings"]["gnm_cutoff_A"] = GNM_CUTOFF_A
    result["settings"]["apbs_score_policy"] = "conditional by chemistry; neutral factor if not applicable/unavailable"
    result["settings"]["geometry_descriptor"] = "C-alpha radial prominence from protein C-alpha centroid"
    return result


app = FastAPI(title="InterfaceScout Backend", version=APP_VERSION)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return {
        "status": "ok",
        "version": APP_VERSION,
        "pdb2pqr": bool(core.PDB2PQR),
        "apbs": bool(core.APBS),
        "dssp": bool(core.MKDSSP),
        "primary_score": "InterfaceScout score",
        "components": "legacy core × GNM mobility × radial prominence × conditional APBS complementarity",
    }


@app.post("/analyze_surface")
def analyze_surface(req: core.AnalyzeRequest):
    try:
        return analyze(req)
    except HTTPException:
        raise
    except Exception as exc:
        core.log.exception("Integrated analysis failed")
        raise HTTPException(500, str(exc))


@app.get("/model_spec")
def model_spec():
    return JSONResponse({
        "version": APP_VERSION,
        "primary_score_name": "InterfaceScout score",
        "primary_score_field": "interfacescout_score",
        "components": {
            "core": "chemistry × scRSA × pH state, aggregated as 5/8 A multiscale patch persistence",
            "dynamics": f"GNM C-alpha mobility; cutoff {GNM_CUTOFF_A:g} A",
            "electrostatics": "APBS residue potential complementarity, conditional on chemistry class",
            "geometry": "C-alpha radial prominence relative to protein centroid",
        },
        "combination": "within-protein percentile-rank product, normalized to 0-100; no fitted coefficients",
        "legacy_core": "retained for ablation/provenance only",
        "scope": "relative protein-side interface propensity; not adsorption free energy or unique orientation",
    })


BACKEND_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BACKEND_DIR.parent
FRONTEND_INDEX = PROJECT_DIR / "frontend" / "index.html"
SAME_DIR_INDEX = BACKEND_DIR / "index.html"


@app.get("/")
def root():
    if FRONTEND_INDEX.exists():
        return FileResponse(FRONTEND_INDEX)
    if SAME_DIR_INDEX.exists():
        return FileResponse(SAME_DIR_INDEX)
    return JSONResponse({"name":"InterfaceScout","version":APP_VERSION,"frontend":"not found"}, status_code=404)


@app.get("/favicon.png")
def favicon():
    path = PROJECT_DIR / "frontend" / "favicon.png"
    if not path.exists():
        raise HTTPException(status_code=404, detail="favicon.png not found")
    return FileResponse(path)


@app.get("/logo.png")
def logo():
    path = PROJECT_DIR / "frontend" / "logo.png"
    if not path.exists():
        raise HTTPException(status_code=404, detail="logo.png not found")
    return FileResponse(path)
