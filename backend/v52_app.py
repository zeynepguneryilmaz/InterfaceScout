"""InterfaceScout v5.2 candidate application.

This module wraps the frozen v5.1 computational core without changing its
canonical chemistry/state or C-alpha 5/8 A scoring equations.  It adds three
lightweight, separately reported capabilities:

1) structural-context handling: RCSB biological assembly 1 when available,
   otherwise the deposited structure; a selected chain is reported *within*
   that structural context rather than being isolated before SASA calculation;
2) Pintar-style CX protrusion descriptors (10 A sphere, 20.1 A^3 mean heavy-
   atom volume), reported as geometry descriptors only;
3) material mechanism profiles that expose several predeclared chemistry maps
   side-by-side without weighted score combination.

The frozen v5.1 backend remains in backend/main.py for regression comparisons.
"""

from __future__ import annotations

import gzip
import math
import re
import shutil
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB import NeighborSearch, PDBIO, PDBParser
from Bio.PDB.Polypeptide import is_aa
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from pydantic import BaseModel, Field

try:  # package import
    from . import main as core
    from .options import MATERIAL_PROFILES, available_material_profiles, material_profile
except ImportError:  # direct execution from backend/
    import main as core
    from options import MATERIAL_PROFILES, available_material_profiles, material_profile

APP_VERSION = "5.2.0-structural-candidate"
CX_RADIUS_A = 10.0
CX_MEAN_ATOM_VOLUME_A3 = 20.1


class EnvParams(BaseModel):
    pH: float = Field(7.4, ge=0.0, le=14.0)
    ionic: float = Field(150.0, ge=0.0, description="mM")
    temp: float = Field(298.0, gt=0.0, description="K")


class AnalyzeRequest(BaseModel):
    pdb_id: Optional[str] = None
    pdb_text: Optional[str] = None
    chain: Optional[str] = None
    env: EnvParams = EnvParams()
    # auto = biological assembly 1 for identifiable RCSB entries, otherwise
    # deposited/full uploaded structure.  selected_chain_legacy reproduces the
    # historical pre-v5.2 isolated-chain behaviour for regression testing.
    structure_context: str = "auto"
    protrusion: bool = True
    material_profile: Optional[str] = None


class FirstModelStandardAA(core.Select):
    """Write only standard amino acids from the first model."""

    def __init__(self) -> None:
        super().__init__()
        self.first_model_id: Any = None

    def accept_model(self, model):
        if self.first_model_id is None:
            self.first_model_id = model.id
        return 1 if model.id == self.first_model_id else 0

    def accept_residue(self, residue):
        return 1 if is_aa(residue, standard=True) else 0


class FirstModelChainAA(FirstModelStandardAA):
    def __init__(self, chain_id: str) -> None:
        super().__init__()
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return 1 if str(chain.id) == self.chain_id else 0


def _detect_pdb_id(text: str) -> Optional[str]:
    """Best-effort RCSB ID extraction from a PDB text file."""
    for line in text.splitlines()[:80]:
        if line.startswith("HEADER") and len(line) >= 66:
            candidate = line[62:66].strip().upper()
            if re.fullmatch(r"[0-9][A-Z0-9]{3}", candidate):
                return candidate
    return None


def _download_text(url: str, target: Path) -> None:
    urllib.request.urlretrieve(url, target)


def _download_biological_assembly_1(pdb_id: str, target: Path) -> None:
    """Download legacy PDB-format biological assembly 1 from RCSB."""
    gz = target.with_suffix(target.suffix + ".gz")
    _download_text(f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb1.gz", gz)
    with gzip.open(gz, "rb") as src, target.open("wb") as dst:
        shutil.copyfileobj(src, dst)


def _raw_input(req: AnalyzeRequest, workdir: Path) -> Tuple[Path, Optional[str]]:
    raw = workdir / "raw.pdb"
    if req.pdb_text:
        raw.write_text(req.pdb_text)
        return raw, (req.pdb_id or _detect_pdb_id(req.pdb_text))
    if req.pdb_id:
        pdb_id = req.pdb_id.strip().upper()
        if not re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb_id):
            raise HTTPException(400, f"Invalid PDB ID: {req.pdb_id}")
        try:
            _download_text(f"https://files.rcsb.org/download/{pdb_id}.pdb", raw)
        except Exception as exc:
            raise HTTPException(400, f"Could not download PDB {pdb_id}: {exc}")
        return raw, pdb_id
    raise HTTPException(400, "Provide pdb_id or pdb_text")


def _prepare_context(req: AnalyzeRequest, workdir: Path) -> Tuple[Path, str, str, Optional[str]]:
    """Return context PDB, report chain, context label and detected PDB ID."""
    raw, pdb_id = _raw_input(req, workdir)
    mode = (req.structure_context or "auto").strip().lower()
    allowed = {"auto", "biological_assembly_1", "deposited_structure", "selected_chain_legacy"}
    if mode not in allowed:
        raise HTTPException(400, f"Unknown structure_context '{req.structure_context}'. Allowed: {sorted(allowed)}")

    source = raw
    context_label = "deposited_structure"
    if mode in {"auto", "biological_assembly_1"} and pdb_id:
        assembly = workdir / "assembly1.pdb"
        try:
            _download_biological_assembly_1(pdb_id, assembly)
            source = assembly
            context_label = "biological_assembly_1"
        except Exception as exc:
            if mode == "biological_assembly_1":
                raise HTTPException(400, f"Biological assembly 1 could not be retrieved for {pdb_id}: {exc}")
            core.log.info("Assembly 1 unavailable for %s; using deposited structure: %s", pdb_id, exc)

    parser = PDBParser(QUIET=True)
    try:
        struct = parser.get_structure("context", str(source))
        first_model = next(struct.get_models())
    except Exception as exc:
        raise HTTPException(400, f"Could not parse structural context: {exc}")

    available = [str(c.id) for c in first_model]
    requested_chain = (req.chain or "").strip()
    if requested_chain and requested_chain not in available:
        raise HTTPException(400, f"Chain '{requested_chain}' not found in {context_label}. Available chains: {', '.join(available)}")

    out = workdir / "analysis_context.pdb"
    io = PDBIO()
    io.set_structure(struct)
    if mode == "selected_chain_legacy" and requested_chain:
        io.save(str(out), FirstModelChainAA(requested_chain))
        context_label = "selected_chain_legacy"
    else:
        io.save(str(out), FirstModelStandardAA())

    return out, requested_chain or "ALL", context_label, pdb_id


def _heavy_atoms_first_model(struct) -> List[Any]:
    model = next(struct.get_models())
    atoms: List[Any] = []
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=True):
                continue
            for atom in res.get_atoms():
                if str(getattr(atom, "element", "")).upper() == "H":
                    continue
                atoms.append(atom)
    return atoms


def compute_cx_descriptors(struct, residues: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Attach Pintar CX descriptors to residue dictionaries.

    Original approximation:
      V_int = N_atom * 20.1 A^3
      V_ext = 4/3*pi*R^3 - V_int
      CX    = V_ext / V_int
    with R=10 A.  The central atom is included in N_atom.  Residue values are
    descriptive summaries; CX is never multiplied into L, P, D or M.
    """
    atoms = _heavy_atoms_first_model(struct)
    if not atoms:
        return {"status": "unavailable", "n_atoms": 0}

    ns = NeighborSearch(atoms)
    sphere_v = 4.0 / 3.0 * math.pi * (CX_RADIUS_A ** 3)
    atom_cx: Dict[int, float] = {}
    for atom in atoms:
        n = len(ns.search(atom.coord, CX_RADIUS_A, level="A"))
        vint = max(float(n) * CX_MEAN_ATOM_VOLUME_A3, 1e-12)
        vext = max(sphere_v - vint, 0.0)
        atom_cx[id(atom)] = vext / vint

    by_key = {r["key"]: r for r in residues}
    for model in struct:
        for chain in model:
            for res in chain:
                if not is_aa(res, standard=True):
                    continue
                key = f"{chain.id}:{int(res.id[1])}:{str(res.id[2]).strip()}"
                row = by_key.get(key)
                if row is None:
                    continue
                vals: List[float] = []
                sc_vals: List[float] = []
                ca_val: Optional[float] = None
                backbone = {"N", "CA", "C", "O", "OXT"}
                for atom in res.get_atoms():
                    if str(getattr(atom, "element", "")).upper() == "H":
                        continue
                    v = atom_cx.get(id(atom))
                    if v is None:
                        continue
                    vals.append(v)
                    name = atom.get_name().strip().upper()
                    if name == "CA":
                        ca_val = v
                    if name not in backbone:
                        sc_vals.append(v)
                row["cx_residue_mean"] = round(float(np.mean(vals)), 5) if vals else None
                row["cx_sidechain_mean"] = round(float(np.mean(sc_vals)), 5) if sc_vals else ca_val
                row["cx_max"] = round(float(np.max(vals)), 5) if vals else None
                row["cx_ca"] = round(float(ca_val), 5) if ca_val is not None else None
        break

    values = [v for v in atom_cx.values() if math.isfinite(v)]
    return {
        "status": "ok",
        "definition": "Pintar CX = external/internal volume in a 10 A sphere",
        "radius_A": CX_RADIUS_A,
        "mean_heavy_atom_volume_A3": CX_MEAN_ATOM_VOLUME_A3,
        "n_atoms": len(atoms),
        "atom_cx_min": round(float(min(values)), 5),
        "atom_cx_median": round(float(np.median(values)), 5),
        "atom_cx_max": round(float(max(values)), 5),
        "used_in_primary_score": False,
    }


def _filter_chemistry_to_chain(group: Dict[str, Any], chain: str) -> Dict[str, Any]:
    if chain == "ALL":
        return group
    out = dict(group)
    out["residues"] = [r for r in group.get("residues", []) if r.get("chain") == chain]
    out["repulsive_residues"] = [r for r in group.get("repulsive_residues", []) if r.get("chain") == chain]
    pcs = [r for r in group.get("patch_centers", []) if r.get("chain") == chain]
    out["patch_centers"] = pcs
    out["top_patches"] = pcs[:10]
    out["top_patch"] = pcs[0] if pcs else None
    out["n_favorable_residues"] = len(out["residues"])
    out["n_repulsive_residues"] = len(out["repulsive_residues"])
    return out


def _filter_feature_to_chain(group: Dict[str, Any], chain: str) -> Dict[str, Any]:
    if chain == "ALL":
        return group
    out = dict(group)
    rows = [r for r in group.get("residues", []) if r.get("chain") == chain]
    out["residues"] = rows
    out["n"] = len(rows)
    return out


def _applicability(context_label: str, protrusion: bool, selected_profile: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    included = [
        "side-chain solvent accessibility (scRSA)",
        "reference-pKa state availability at the requested bulk pH",
        "C-alpha 5/8 A multiscale persistence",
    ]
    limits: List[str] = []
    if context_label == "biological_assembly_1":
        included.append("biological-assembly structural context")
    elif context_label == "selected_chain_legacy":
        limits.append("oligomeric shielding and inter-chain neighborhoods are not represented")
    else:
        included.append("full deposited/uploaded structural context")
        limits.append("for uploaded structures, biological-assembly correctness depends on the supplied coordinates")
    if protrusion:
        included.append("Pintar CX protrusion descriptors (auxiliary; not scored)")
    else:
        limits.append("accessible residues are not distinguished by geometric protrusion")
    if selected_profile:
        included.append("material mechanism profile with unweighted, separately reported chemistry channels")
    else:
        limits.append("each chemistry map represents an idealized single interface chemistry unless a material profile is selected")

    limits.extend([
        "structure-specific pKa shifts are not included in the canonical score",
        "whole-protein electrostatic steering is not included in the canonical score",
        "explicit interfacial hydration/desolvation is not modeled",
        "adsorption-induced large conformational change/unfolding is not modeled",
        "multi-protein crowding and mature corona organization are not modeled",
        "absolute adsorption free energy and adsorption capacity are not predicted",
    ])
    return {"included_in_this_run": included, "not_included_or_interpretation_limits": limits}


def analyze(req: AnalyzeRequest) -> Dict[str, Any]:
    workdir = Path(tempfile.mkdtemp(prefix="interfacescout_v52_"))
    try:
        pdb, selected_chain, context_label, detected_id = _prepare_context(req, workdir)
        struct, all_residues, atoms_by_key, _ = core.build_surface_residues(pdb, req.env.pH)
        if not all_residues:
            raise HTTPException(400, "No standard amino-acid residues were found")

        # CX is computed after the canonical structural parse and is auxiliary.
        cx_meta = compute_cx_descriptors(struct, all_residues) if req.protrusion else {"status": "disabled", "used_in_primary_score": False}

        core_env = core.EnvParams(pH=req.env.pH, ionic=req.env.ionic, temp=req.env.temp)
        electrostatics = core.attach_apbs_auxiliary(pdb, all_residues, atoms_by_key, core_env, workdir)
        surface_all = [r for r in all_residues if r["surface_exposed"]]
        distances = core.build_distances(surface_all)
        if len(surface_all):
            for i, r in enumerate(surface_all):
                r["n_neighbors_8A"] = int(np.sum((distances[i] <= core.CONTEXT_RADIUS_A) & (distances[i] > 0.0)))

        # Canonical chemistry calculations are performed in the full structural
        # context.  If a chain is requested, only the reported centres/residues
        # are filtered after calculation; cross-chain shielding/neighborhoods
        # remain part of the structural context.
        chem_all = {k: core.chemistry_map(surface_all, distances, k, req.env.pH) for k in core.CHEMISTRIES}
        features_all = core.feature_map(surface_all, req.env.pH)
        chem = {k: _filter_chemistry_to_chain(v, selected_chain) for k, v in chem_all.items()}
        features = {k: _filter_feature_to_chain(v, selected_chain) for k, v in features_all.items()}
        all_report = all_residues if selected_chain == "ALL" else [r for r in all_residues if r["chain"] == selected_chain]
        surface_report = surface_all if selected_chain == "ALL" else [r for r in surface_all if r["chain"] == selected_chain]

        profile = material_profile(req.material_profile)
        if req.material_profile and profile is None:
            raise HTTPException(400, f"Unknown material_profile '{req.material_profile}'.")
        profile_result = None
        if profile:
            profile_result = {
                "key": req.material_profile,
                **profile,
                "channel_results": {ch: chem[ch] for ch in profile["channels"] if ch in chem},
                "combination_rule": "none; chemistry channels are reported separately",
            }

        n_atoms_context = sum(1 for m in struct for c in m for res in c if is_aa(res, standard=True) for _ in res.get_atoms())
        result = {
            "status": "ok",
            "version": APP_VERSION,
            "core_version": core.APP_VERSION,
            "model": "InterfaceScout v5.2 structural-context candidate; frozen v5.1 scoring core",
            "scope": {
                "predicts": "protein-side residue/patch compatibility hypotheses for generalized interface chemistries",
                "does_not_predict": [
                    "adsorption capacity", "absolute adsorption free energy", "unique adsorption orientation",
                    "adsorption-induced conformational change", "explicit interfacial hydration", "multi-protein corona organization",
                ],
            },
            "settings": {
                "chain": selected_chain,
                "structure_context": context_label,
                "detected_pdb_id": detected_id,
                "pH": req.env.pH, "ionic_mM": req.env.ionic, "temperature_K": req.env.temp,
                "sasa_probe_A": core.SASA_PROBE_A, "sasa_points_per_atom": core.SASA_POINTS,
                "scrsa_threshold": core.SC_RSA_THRESHOLD,
                "patch_radii_A": list(core.PATCH_RADII_A),
                "patch_pair_selection": core.PATCH_PAIR_AUDIT,
                "context_radius_A_auxiliary": core.CONTEXT_RADIUS_A,
                "protrusion_enabled": bool(req.protrusion),
                "material_profile": req.material_profile,
            },
            "stats": {
                "n_atoms_context": n_atoms_context,
                "n_residues_context": len(all_residues),
                "n_surface_res_context": len(surface_all),
                "n_residues_reported": len(all_report),
                "n_surface_res_reported": len(surface_report),
                # aliases retained for the existing frontend
                "n_atoms": n_atoms_context,
                "n_residues": len(all_report),
                "n_surface_res": len(surface_report),
                "electrostatics": electrostatics,
                "pdb2pqr": bool(core.PDB2PQR), "apbs": bool(core.APBS), "dssp": bool(core.MKDSSP),
            },
            "protrusion": cx_meta,
            "material_profile": profile_result,
            "available_material_profiles": available_material_profiles(),
            "applicability": _applicability(context_label, req.protrusion, profile),
            "chemistry_list": list(core.CHEMISTRIES.keys()),
            "chemistries": chem,
            "feature_list": list(core.FEATURE_RESIDUES.keys()),
            "features": features,
            "all_residues": all_report,
            "surface_residues": surface_report,
            "reference_sidechain_asa": core.SIDECHAIN_REF_ASA,
        }
        return result
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


app = FastAPI(title="InterfaceScout Backend", version=APP_VERSION)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return {
        "status": "ok",
        "version": APP_VERSION,
        "core_version": core.APP_VERSION,
        "canonical_model": "frozen v5.1: scRSA x membership x state availability; C-alpha 5/8 A persistence",
        "structural_context": "auto biological assembly 1 when an RCSB ID is detectable",
        "cx_protrusion": "auxiliary only",
        "pdb2pqr": bool(core.PDB2PQR), "apbs": bool(core.APBS), "dssp": bool(core.MKDSSP),
    }


@app.get("/material_profiles")
def material_profiles():
    return {"profiles": available_material_profiles(), "combination_rule": "none; channels remain separate"}


@app.post("/analyze_surface")
def analyze_surface(req: AnalyzeRequest):
    try:
        return analyze(req)
    except HTTPException:
        raise
    except Exception as exc:
        core.log.exception("v5.2 analysis failed")
        raise HTTPException(500, str(exc))


@app.get("/model_spec")
def model_spec():
    spec = core.model_spec().body.decode("utf-8") if hasattr(core.model_spec(), "body") else None
    return {
        "version": APP_VERSION,
        "core_version": core.APP_VERSION,
        "canonical_core_unchanged": True,
        "cx": {"radius_A": CX_RADIUS_A, "mean_atom_volume_A3": CX_MEAN_ATOM_VOLUME_A3, "used_in_primary_score": False},
        "structure_context_modes": ["auto", "biological_assembly_1", "deposited_structure", "selected_chain_legacy"],
        "material_profiles": list(MATERIAL_PROFILES.keys()),
        "core_spec_json": spec,
    }


@app.get("/")
def root():
    path = core.PROJECT_DIR / "frontend" / "index.html"
    if path.exists():
        return HTMLResponse(path.read_text(encoding="utf-8"))
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
