"""Structural-context layer for InterfaceScout 2.0 development.

This layer preserves the frozen InterfaceScout 1.0 chemistry/state and C-alpha
5/8 A scoring equations while adding protein-side structural context:

- RCSB biological assembly 1 when available,
- Pintar-style CX protrusion descriptors (auxiliary only),
- optional user-specified function-critical residue annotations.

No named-material library is used here.  The chemistry maps are generalized
protein-side interaction channels that feed the protein-derived target interface
profile in ``target_profile.py``.
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
from fastapi import HTTPException
from pydantic import BaseModel, Field

try:
    from . import main as core
except ImportError:
    import main as core

STRUCTURAL_LAYER_VERSION = "2.0.0-dev"
CORE_RELEASE_VERSION = "1.0.0"
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
    structure_context: str = "auto"
    protrusion: bool = True
    protected_residue_keys: Optional[List[str]] = None


class FirstModelStandardAA(core.Select):
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
    for line in text.splitlines()[:80]:
        if line.startswith("HEADER") and len(line) >= 66:
            candidate = line[62:66].strip().upper()
            if re.fullmatch(r"[0-9][A-Z0-9]{3}", candidate):
                return candidate
    return None


def _download(url: str, target: Path) -> None:
    urllib.request.urlretrieve(url, target)


def _download_biological_assembly_1(pdb_id: str, target: Path) -> None:
    gz = target.with_suffix(target.suffix + ".gz")
    _download(f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb1.gz", gz)
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
            _download(f"https://files.rcsb.org/download/{pdb_id}.pdb", raw)
        except Exception as exc:
            raise HTTPException(400, f"Could not download PDB {pdb_id}: {exc}")
        return raw, pdb_id
    raise HTTPException(400, "Provide pdb_id or pdb_text")


def prepare_context(req: AnalyzeRequest, workdir: Path) -> Tuple[Path, str, str, Optional[str]]:
    """Prepare structural coordinates and return path, report chain, context, PDB ID."""
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
                if str(getattr(atom, "element", "")).upper() != "H":
                    atoms.append(atom)
    return atoms


def compute_cx_descriptors(struct, residues: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Attach Pintar-style assembly-context CX descriptors without changing scores."""
    atoms = _heavy_atoms_first_model(struct)
    if not atoms:
        return {"status": "unavailable", "n_atoms": 0, "used_in_primary_score": False}

    ns = NeighborSearch(atoms)
    sphere_v = 4.0 / 3.0 * math.pi * CX_RADIUS_A**3
    atom_cx: Dict[int, float] = {}
    for atom in atoms:
        n = len(ns.search(atom.coord, CX_RADIUS_A, level="A"))
        vint = max(float(n) * CX_MEAN_ATOM_VOLUME_A3, 1e-12)
        atom_cx[id(atom)] = max(sphere_v - vint, 0.0) / vint

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
        "definition": "Pintar-style CX = external/internal volume in a 10 A sphere",
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


def applicability_notes(context_label: str, protrusion: bool) -> Dict[str, Any]:
    included = [
        "side-chain solvent accessibility (scRSA)",
        "reference-pKa state availability at the requested bulk pH",
        "C-alpha 5/8 A multiscale persistence",
        "protein-derived generalized interaction channels",
    ]
    limits: List[str] = []
    if context_label == "biological_assembly_1":
        included.append("biological-assembly structural context")
    elif context_label == "selected_chain_legacy":
        limits.append("oligomeric shielding and inter-chain neighborhoods are not represented")
    else:
        included.append("full deposited/uploaded structural context")
        limits.append("for uploaded structures, biological-assembly correctness depends on supplied coordinates")

    if protrusion:
        included.append("Pintar-style assembly-context CX protrusion descriptors (auxiliary; not scored)")
    else:
        limits.append("accessible residues are not distinguished by geometric protrusion")

    limits.extend([
        "no named-material library is used",
        "structure-specific pKa shifts are not included in the canonical score",
        "whole-protein electrostatic steering is not included in the canonical score",
        "explicit interfacial hydration/desolvation is not modeled by the canonical core",
        "adsorption-induced large conformational change/unfolding is not modeled",
        "multi-protein crowding and mature corona organization are not modeled",
        "absolute adsorption free energy and adsorption capacity are not predicted",
    ])
    return {"included_in_this_run": included, "not_included_or_interpretation_limits": limits}


def analyze_structural(req: AnalyzeRequest) -> Dict[str, Any]:
    workdir = Path(tempfile.mkdtemp(prefix="interfacescout_v2_"))
    try:
        pdb, selected_chain, context_label, detected_id = prepare_context(req, workdir)
        struct, all_residues, atoms_by_key, _ = core.build_surface_residues(pdb, req.env.pH)
        if not all_residues:
            raise HTTPException(400, "No standard amino-acid residues were found")

        cx_meta = compute_cx_descriptors(struct, all_residues) if req.protrusion else {"status": "disabled", "used_in_primary_score": False}
        core_env = core.EnvParams(pH=req.env.pH, ionic=req.env.ionic, temp=req.env.temp)
        electrostatics = core.attach_apbs_auxiliary(pdb, all_residues, atoms_by_key, core_env, workdir)

        surface_all = [r for r in all_residues if r["surface_exposed"]]
        distances = core.build_distances(surface_all)
        if len(surface_all):
            for i, r in enumerate(surface_all):
                r["n_neighbors_8A"] = int(np.sum((distances[i] <= core.CONTEXT_RADIUS_A) & (distances[i] > 0.0)))

        chem_all = {k: core.chemistry_map(surface_all, distances, k, req.env.pH) for k in core.CHEMISTRIES}
        features_all = core.feature_map(surface_all, req.env.pH)
        chem = {k: _filter_chemistry_to_chain(v, selected_chain) for k, v in chem_all.items()}
        features = {k: _filter_feature_to_chain(v, selected_chain) for k, v in features_all.items()}
        all_report = all_residues if selected_chain == "ALL" else [r for r in all_residues if r["chain"] == selected_chain]
        surface_report = surface_all if selected_chain == "ALL" else [r for r in surface_all if r["chain"] == selected_chain]

        n_atoms_context = sum(1 for m in struct for c in m for res in c if is_aa(res, standard=True) for _ in res.get_atoms())
        return {
            "status": "ok",
            "version": STRUCTURAL_LAYER_VERSION,
            "core_version": CORE_RELEASE_VERSION,
            "model": "InterfaceScout 2.0 material-agnostic structural layer over the InterfaceScout 1.0 canonical scoring core",
            "scope": {
                "predicts": "protein-side residue/patch compatibility hypotheses for generalized interface properties",
                "does_not_predict": [
                    "named material identity", "adsorption capacity", "absolute adsorption free energy", "unique adsorption orientation",
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
            },
            "stats": {
                "n_atoms_context": n_atoms_context,
                "n_residues_context": len(all_residues),
                "n_surface_res_context": len(surface_all),
                "n_residues_reported": len(all_report),
                "n_surface_res_reported": len(surface_report),
                "n_atoms": n_atoms_context,
                "n_residues": len(all_report),
                "n_surface_res": len(surface_report),
                "electrostatics": electrostatics,
                "pdb2pqr": bool(core.PDB2PQR), "apbs": bool(core.APBS), "dssp": bool(core.MKDSSP),
            },
            "protrusion": cx_meta,
            "applicability": applicability_notes(context_label, req.protrusion),
            "chemistry_list": list(core.CHEMISTRIES.keys()),
            "chemistries": chem,
            "feature_list": list(core.FEATURE_RESIDUES.keys()),
            "features": features,
            "all_residues": all_report,
            "surface_residues": surface_report,
            "reference_sidechain_asa": core.SIDECHAIN_REF_ASA,
        }
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
