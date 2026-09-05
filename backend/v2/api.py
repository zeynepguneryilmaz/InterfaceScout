"""Standalone FastAPI surface for InterfaceScout V2-alpha.

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

# Allow engine.py to import the frozen backend/main.py as module ``main``.
_BACKEND = Path(__file__).resolve().parents[1]
if str(_BACKEND) not in sys.path:
    sys.path.insert(0, str(_BACKEND))

from .engine import analyze_v2
from .surface_profiles import SURFACE_PROFILES

app = FastAPI(title="InterfaceScout V2-alpha", version="2.0.0-alpha.1")


class V2AnalyzeRequest(BaseModel):
    surface: str
    pdb_id: Optional[str] = None
    pdb_text: Optional[str] = None
    chain: Optional[str] = None
    pH: float = Field(7.4, ge=0.0, le=14.0)
    ionic_mM: float = Field(150.0, ge=0.0)
    temp_K: float = Field(298.0, gt=0.0)
    top_n: int = Field(3, ge=1, le=20)


@app.get("/v2/health")
def health():
    return {
        "status": "ok",
        "engine": "InterfaceScout V2-alpha",
        "version": "2.0.0-alpha.1",
        "v1_modified": False,
    }


@app.get("/v2/surfaces")
def surfaces():
    return {
        key: {
            "label": profile.label,
            "description": profile.description,
            "weights": profile.normalized_weights(),
        }
        for key, profile in sorted(SURFACE_PROFILES.items())
    }


@app.post("/v2/analyze")
def analyze(req: V2AnalyzeRequest):
    try:
        return analyze_v2(
            surface=req.surface,
            pH=req.pH,
            ionic_mM=req.ionic_mM,
            temp_K=req.temp_K,
            pdb_id=req.pdb_id,
            pdb_text=req.pdb_text,
            chain=req.chain,
            top_n=req.top_n,
        )
    except KeyError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"V2-alpha analysis failed: {exc}") from exc
