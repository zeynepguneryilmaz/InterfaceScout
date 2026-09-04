"""Runtime entry point for the InterfaceScout v5.2 structural candidate."""

from __future__ import annotations

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse

try:
    from . import v52_app as model
    from .v52_ui import inject_v52_ui
except ImportError:
    import v52_app as model
    from v52_ui import inject_v52_ui

app = FastAPI(title="InterfaceScout Backend", version=model.APP_VERSION)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return model.health()


@app.get("/material_profiles")
def material_profiles():
    return model.material_profiles()


@app.post("/analyze_surface")
def analyze_surface(req: model.AnalyzeRequest):
    return model.analyze_surface(req)


@app.get("/model_spec")
def model_spec():
    return model.model_spec()


@app.get("/")
def root():
    path = model.core.PROJECT_DIR / "frontend" / "index.html"
    if not path.exists():
        return JSONResponse({"name": "InterfaceScout", "version": model.APP_VERSION, "frontend": "not found"}, status_code=404)
    return HTMLResponse(inject_v52_ui(path.read_text(encoding="utf-8")))


@app.get("/favicon.png")
def favicon():
    path = model.core.PROJECT_DIR / "frontend" / "favicon.png"
    if not path.exists():
        raise HTTPException(404, "favicon.png not found")
    return FileResponse(path)


@app.get("/logo.png")
def logo():
    path = model.core.PROJECT_DIR / "frontend" / "logo.png"
    if not path.exists():
        raise HTTPException(404, "logo.png not found")
    return FileResponse(path)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
