"""InterfaceScout web application."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from interface_engine import PUBLIC_MAP_SPECS, analyze_all_maps
from model_settings import MODEL_VERSION
from research_workflows import list_research_goals, compare_results, combine_maps
from research_workflows import list_research_goals, compare_results, combine_maps

BACKEND = Path(__file__).resolve().parent
ROOT = BACKEND.parent
FRONTEND = ROOT / "frontend"

PUBLIC_VERSION = "1.0-publication"
app = FastAPI(title="InterfaceScout", version=PUBLIC_VERSION)


class AnalyzeRequest(BaseModel):
    pdb_id: Optional[str] = None
    pdb_text: Optional[str] = None
    chain: Optional[str] = None
    pH: float = Field(7.4, ge=0.0, le=14.0)
    ionic_mM: float = Field(150.0, ge=0.0)
    temp_K: float = Field(298.0, gt=0.0)

class CompareRequest(BaseModel):
    structure_a: AnalyzeRequest
    structure_b: AnalyzeRequest


class CombineRequest(BaseModel):
    analysis: dict
    selected_maps: list[str] = Field(min_length=2, max_length=11)
    max_regions: int = Field(25, ge=1, le=100)


@app.get("/research-goals")
def research_goals():
    return list_research_goals()


@app.post("/compare")
def compare(req: CompareRequest):
    try:
        a = analyze_all_maps(**req.structure_a.model_dump())
        b = analyze_all_maps(**req.structure_b.model_dump())
        return compare_results(a, b)
    except (KeyError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/combine")
def combine(req: CombineRequest):
    try:
        return combine_maps(req.analysis, req.selected_maps, req.max_regions)
    except (KeyError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/health")
def health():
    return {
        "status": "ok",
        "engine": "InterfaceScout",
        "version": PUBLIC_VERSION,
        "model_version": MODEL_VERSION,
        "single_pass_all_maps": True,
        "material_presets": False,
        "benchmark_fitted_weights": False,
    }


@app.get("/chemistries")
def chemistries():
    return [
        {"key": public_key, "label": label, "description": description}
        for public_key, _internal_key, label, description in PUBLIC_MAP_SPECS
    ]


@app.post("/analyze")
def analyze(req: AnalyzeRequest):
    try:
        out = analyze_all_maps(
            pH=req.pH,
            ionic_mM=req.ionic_mM,
            temp_K=req.temp_K,
            pdb_id=req.pdb_id,
            pdb_text=req.pdb_text,
            chain=req.chain,
        )
        out["public_version"] = PUBLIC_VERSION
        return out
    except (KeyError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"InterfaceScout analysis failed: {exc}") from exc


if FRONTEND.exists():
    app.mount("/", StaticFiles(directory=str(FRONTEND), html=True), name="frontend")
