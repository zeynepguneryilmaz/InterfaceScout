"""Standalone FastAPI surface for InterfaceScout V2 coarse biointerface prediction.

Run from repository root with:
    uvicorn backend.v2.api:app --reload

V1 remains served by backend/main.py and is not modified by this module.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

_BACKEND = Path(__file__).resolve().parents[1]
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

from . import V2_VERSION
from .interface_engine import analyze_interface_v2
from .surface_modes import SURFACE_MODES

app = FastAPI(title="InterfaceScout V2", version=V2_VERSION)


class V2AnalyzeRequest(BaseModel):
    surface: str
    pdb_id: Optional[str] = None
    pdb_text: Optional[str] = None
    chain: Optional[str] = None
    pH: float = Field(7.4, ge=0.0, le=14.0)
    ionic_mM: float = Field(150.0, ge=0.0)
    temp_K: float = Field(298.0, gt=0.0)


@app.get("/v2/health")
def health():
    return {
        "status": "ok",
        "engine": "InterfaceScout V2 coarse interface predictor",
        "version": V2_VERSION,
        "v1_modified": False,
        "benchmark_fitted_weights": False,
    }


@app.get("/v2/surfaces")
def surfaces():
    return {
        key: {
            "label": mode.label,
            "primary_chemistry": mode.chemistry,
            "description": mode.description,
        }
        for key, mode in sorted(SURFACE_MODES.items())
    }


@app.post("/v2/analyze")
def analyze(req: V2AnalyzeRequest):
    try:
        return analyze_interface_v2(
            surface=req.surface,
            pH=req.pH,
            ionic_mM=req.ionic_mM,
            temp_K=req.temp_K,
            pdb_id=req.pdb_id,
            pdb_text=req.pdb_text,
            chain=req.chain,
        )
    except KeyError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"V2 analysis failed: {exc}") from exc
