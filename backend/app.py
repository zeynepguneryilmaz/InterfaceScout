"""InterfaceScout web application."""
from __future__ import annotations

from pathlib import Path
from typing import Optional
from io import BytesIO
import os
import time
import uuid
import tempfile
import threading
import zipfile

from fastapi import FastAPI, HTTPException
from fastapi.responses import StreamingResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from interface_engine import PUBLIC_MAP_SPECS, analyze_all_maps
from model_settings import MODEL_VERSION, PUBLIC_VERSION
from research_workflows import list_research_goals, compare_results, combine_maps
from ensemble import analyze_structure_set, build_consensus
from export_bundle import build_export_zip
from optional_electrostatics import ElectrostaticSettings, ElectrostaticsError, compute as compute_optional_electrostatics, status as electrostatics_status

BACKEND = Path(__file__).resolve().parent
ROOT = BACKEND.parent
FRONTEND = ROOT / "frontend"

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


# Prepared exports are stored under a private temporary directory on the
# LOCAL application server, not in a browser Blob URL. The browser can request
# a normal GET/download directly, even after long-running screenshot work.
_EXPORT_DIRECTORY = Path(tempfile.mkdtemp(prefix="interfacescout_exports_"))
_EXPORT_TTL_SECONDS = 60 * 60
_EXPORT_LIMIT = 20
_EXPORT_LOCK = threading.Lock()
_EXPORT_FILES: dict[str, tuple[Path, float]] = {}


def _clean_expired_exports(now: float) -> None:
    # Called only while holding _EXPORT_LOCK.
    for token, (path, created) in list(_EXPORT_FILES.items()):
        if now - created > _EXPORT_TTL_SECONDS or not path.exists():
            path.unlink(missing_ok=True)
            _EXPORT_FILES.pop(token, None)


@app.post("/export-package/prepare")
def prepare_export_package(req: ExportRequest):
    """Generate the full archive on disk and issue a temporary GET download."""
    try:
        archive = build_export_zip(req.model_dump())
    except (KeyError, ValueError, TypeError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    token = uuid.uuid4().hex
    path = _EXPORT_DIRECTORY / (token + ".zip")
    # Write atomically: never advertise a partially written ZIP.
    draft = _EXPORT_DIRECTORY / (token + ".tmp")
    try:
        draft.write_bytes(archive)
        os.replace(draft, path)
    finally:
        draft.unlink(missing_ok=True)

    with _EXPORT_LOCK:
        now = time.time()
        _clean_expired_exports(now)
        while len(_EXPORT_FILES) >= _EXPORT_LIMIT:
            oldest = min(_EXPORT_FILES, key=lambda t: _EXPORT_FILES[t][1])
            _EXPORT_FILES.pop(oldest)[0].unlink(missing_ok=True)
        _EXPORT_FILES[token] = (path, now)
    return {
        "download_url": "/export-package/download/" + token,
        "filename": "InterfaceScout_Complete_Results.zip",
        "bytes": path.stat().st_size,
        "expires_in_seconds": _EXPORT_TTL_SECONDS,
        "message": "Ready for direct browser download",
    }


@app.get("/export-package/download/{token}")
def download_prepared_export(token: str):
    """Serve a prepared ZIP with native browser attachment semantics."""
    with _EXPORT_LOCK:
        _clean_expired_exports(time.time())
        entry = _EXPORT_FILES.get(token)
    if entry is None or not entry[0].is_file():
        raise HTTPException(status_code=404, detail="This export has expired. Prepare a new download.")
    return FileResponse(
        path=entry[0],
        media_type="application/zip",
        filename="InterfaceScout_Complete_Results.zip",
        content_disposition_type="attachment",
        headers={"Cache-Control": "no-store", "X-Content-Type-Options": "nosniff"},
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



class ElectrostaticsRequest(BaseModel):
    pdb_text: str = Field(...,min_length=20,max_length=16_000_000)
    chain: Optional[str] = None
    pH: float = Field(7.4,ge=0,le=14)
    ionic_mM: float = Field(150.0,ge=0,le=1000)
    temperature_K: float = Field(298.15,ge=273.15,le=373.15)
    electrolyte: str = Field("NaCl",pattern="^(NaCl|KCl)$")
    forcefield: str = Field("PARSE",pattern="^(PARSE|AMBER|CHARMM)$")
    protein_dielectric: float = Field(2.0,ge=1,le=30)
    solvent_dielectric: float = Field(78.5,ge=20,le=110)
    apbs_enabled: bool = True
    propka_enabled: bool = True

@app.get("/electrostatics/status")
def optional_electrostatics_status():
    return electrostatics_status()

@app.post("/electrostatics/prepare")
def optional_electrostatics_prepare(req: ElectrostaticsRequest):
    """Calculate once, retain evidence ZIP briefly, return viewer-compatible 11 maps."""
    try:
        values=req.model_dump()
        pdb_text=values.pop("pdb_text")
        chain=values.pop("chain")
        settings=ElectrostaticSettings(**values)
        result, archive=compute_optional_electrostatics(
            pdb_text=pdb_text,chain=chain,settings=settings)
    except ElectrostaticsError as exc:
        raise HTTPException(status_code=422,detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400,detail=str(exc)) from exc
    physical=result.pop("physics_aware_view",None)
    canonical=result.pop("canonical_view",None)
    token=uuid.uuid4().hex
    path=_EXPORT_DIRECTORY/(token+".zip")
    draft=_EXPORT_DIRECTORY/(token+".tmp")
    try:
        draft.write_bytes(archive)
        os.replace(draft,path)
    finally:
        draft.unlink(missing_ok=True)
    with _EXPORT_LOCK:
        now=time.time()
        _clean_expired_exports(now)
        while len(_EXPORT_FILES)>=_EXPORT_LIMIT:
            oldest=min(_EXPORT_FILES,key=lambda k:_EXPORT_FILES[k][1])
            _EXPORT_FILES.pop(oldest)[0].unlink(missing_ok=True)
        _EXPORT_FILES[token]=(path,now)
    return {
        "status":"complete",
        "result":result,
        "analysis":physical,
        "canonical_analysis":canonical,
        "download_url":"/electrostatics/download/"+token,
        "potential_url":"/electrostatics/potential/"+token if physical is not None else None,
        "filename":"InterfaceScout_PhysicsAware_Electrostatics.zip",
        "expires_in_seconds":_EXPORT_TTL_SECONDS,
    }

@app.get("/electrostatics/potential/{token}")
def electrostatics_potential_dx(token:str):
    """Private temporary OpenDX visualization source for the current local run."""
    with _EXPORT_LOCK:
        _clean_expired_exports(time.time())
        entry=_EXPORT_FILES.get(token)
    if not entry or not entry[0].is_file():
        raise HTTPException(status_code=404,detail="Potential file expired; recalculate")
    with zipfile.ZipFile(entry[0]) as archive:
        if "data/potential.dx" not in archive.namelist():
            raise HTTPException(status_code=404,detail="APBS potential unavailable")
        raw=archive.read("data/potential.dx")
    if len(raw)>100_000_000:
        raise HTTPException(status_code=413,detail="Potential grid too large for browser")
    return StreamingResponse(BytesIO(raw),media_type="text/plain",headers={
        "Cache-Control":"no-store","X-Content-Type-Options":"nosniff"})

@app.get("/electrostatics/download/{token}")
def download_electrostatics(token:str):
    with _EXPORT_LOCK:
        _clean_expired_exports(time.time())
        entry=_EXPORT_FILES.get(token)
    if not entry or not entry[0].is_file():
        raise HTTPException(status_code=404,detail="Electrostatics result expired; recalculate")
    return FileResponse(entry[0],media_type="application/zip",
            filename="InterfaceScout_PhysicsAware_Electrostatics.zip",
            content_disposition_type="attachment",
            headers={"Cache-Control":"no-store","X-Content-Type-Options":"nosniff"})


@app.post("/electrostatics/export")
def optional_electrostatics_export(req: ElectrostaticsRequest):
    """Isolated optional physically parameterized calculation and evidence ZIP."""
    try:
        values=req.model_dump()
        pdb_text=values.pop("pdb_text")
        chain=values.pop("chain")
        settings=ElectrostaticSettings(**values)
        _, archive=compute_optional_electrostatics(pdb_text=pdb_text,chain=chain,settings=settings)
    except ElectrostaticsError as exc:
        raise HTTPException(status_code=422,detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400,detail=str(exc)) from exc
    return StreamingResponse(BytesIO(archive),media_type="application/zip",headers={
        "Content-Disposition":'attachment; filename="InterfaceScout_Electrostatics_Optional.zip"',
        "X-Content-Type-Options":"nosniff"})

if FRONTEND.exists():
    app.mount("/", StaticFiles(directory=str(FRONTEND), html=True), name="frontend")
