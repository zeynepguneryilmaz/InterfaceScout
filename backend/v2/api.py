"""Canonical InterfaceScout web application."""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

_BACKEND = Path(__file__).resolve().parents[1]
_ROOT = _BACKEND.parent
_FRONTEND = _ROOT / "frontend"
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

from .interface_engine import analyze_all_maps
from .model_settings import MODEL_VERSION
from .surface_modes import SURFACE_MODES

PUBLIC_VERSION = "1.0-publication"
app = FastAPI(title="InterfaceScout", version=PUBLIC_VERSION)


class AnalyzeRequest(BaseModel):
    pdb_id: Optional[str] = None
    pdb_text: Optional[str] = None
    chain: Optional[str] = None
    pH: float = Field(7.4, ge=0.0, le=14.0)
    ionic_mM: float = Field(150.0, ge=0.0)
    temp_K: float = Field(298.0, gt=0.0)
    initial_surface: Optional[str] = None


@app.get("/health")
def health():
    return {
        "status": "ok",
        "engine": "InterfaceScout",
        "version": PUBLIC_VERSION,
        "model_version": MODEL_VERSION,
        "single_pass_all_maps": True,
        "benchmark_fitted_weights": False,
    }


@app.get("/surfaces")
def surfaces():
    return {
        key: {
            "label": mode.label,
            "primary_chemistry": mode.chemistry,
            "description": mode.description,
        }
        for key, mode in sorted(SURFACE_MODES.items())
    }


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
            initial_surface=req.initial_surface,
        )
        out["public_version"] = PUBLIC_VERSION
        return out
    except (KeyError, ValueError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"InterfaceScout analysis failed: {exc}") from exc


if _FRONTEND.exists():
    app.mount("/", StaticFiles(directory=str(_FRONTEND), html=True), name="frontend")
