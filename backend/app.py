"""InterfaceScout web application."""
from __future__ import annotations

from pathlib import Path
from typing import Optional
from io import BytesIO

from fastapi import FastAPI, HTTPException
from fastapi.responses import StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from interface_engine import PUBLIC_MAP_SPECS, analyze_all_maps
from model_settings import MODEL_VERSION
from research_workflows import list_research_goals, compare_results, combine_maps
from ensemble import analyze_structure_set, build_consensus
from export_bundle import build_export_zip

BACKEND = Path(__file__).resolve().parent
ROOT = BACKEND.parent
FRONTEND = ROOT / "frontend"

PUBLIC_VERSION = "1.1.0-rc1"
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
    max_regions: Optional[int] = Field(25, ge=1, le=10000)



class StructureInput(BaseModel):
    name: str = Field(..., min_length=1, max_length=85)
    pdb_id: Optional[str] = None
    pdb_text: Optional[str] = None
    chain: Optional[str] = None


class StructureSetRequest(BaseModel):
    structures: list[StructureInput] = Field(..., min_length=1, max_length=20)
    pH: float = Field(7.4, ge=0.0, le=14.0)
    ionic_mM: float = Field(150.0, ge=0.0)
    temp_K: float = Field(298.0, gt=0.0)
    combined_maps: Optional[list[str]] = None


class DeriveConsensusRequest(BaseModel):
    structures: list[dict] = Field(..., min_length=2, max_length=20)
    selected_maps: list[str] = Field(default_factory=list, max_length=11)


@app.post("/derive-consensus")
def derive_consensus(req: DeriveConsensusRequest):
    """Recompute descriptive consensus from existing analyses (no re-analysis)."""
    try:
        selected = list(dict.fromkeys(req.selected_maps))
        if len(selected) == 1:
            raise ValueError("Select either zero or at least two chemistry maps")
        runs = []
        for item in req.structures:
            analysis = item["analysis"]
            runs.append({
                "name": item["name"],
                "analysis": analysis,
                "combined_maps": combine_maps(analysis, selected, max_regions=None)
                if selected else None,
            })
        return build_consensus(runs)
    except (KeyError, ValueError, TypeError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


class ExportRequest(BaseModel):
    structures: list[dict] = Field(..., min_length=1, max_length=20)
    consensus: Optional[dict] = None
    comparison: list[dict] = Field(default_factory=list)
    consensus_screenshots: dict[str, str] = Field(default_factory=dict)
    comparison_skipped: list[dict] = Field(default_factory=list)
    consensus_unavailable_reason: Optional[str] = None


@app.post("/structure-set")
def structure_set(req: StructureSetRequest):
    """One structure or multiple PDB structures; all existing map rules unchanged."""
    try:
        return analyze_structure_set(
            [s.model_dump() for s in req.structures],
            pH=req.pH, ionic_mM=req.ionic_mM, temp_K=req.temp_K,
            combined_maps=req.combined_maps,
        )
    except (KeyError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/export-package")
def export_package(req: ExportRequest):
    """Export full data plus actual client-captured 3D viewer PNGs."""
    try:
        data = build_export_zip(req.model_dump())
    except (KeyError, ValueError, TypeError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return StreamingResponse(
        BytesIO(data), media_type="application/zip",
        headers={"Content-Disposition": 'attachment; filename="InterfaceScout_complete_results.zip"'},
    )


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
