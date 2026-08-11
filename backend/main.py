"""
InterfaceScout — backend
==================================================
Canonical model implemented to match the manuscript/SI formulation:

    L_i,c = I_i,c * scRSA_i * f_state,i,c(pH)
    P_i,c = 100 * L_i,c / max(L_c,fav)

Patch enrichment is evaluated independently at 5 and 8 Å:

    D_i,c(R) = sum_j[d_ij <= R] L_j,c
    M_i,c = 100 * min(D5_norm, D8_norm)

Important scope boundaries
--------------------------
* Literature Ebase values are mechanistic metadata only; they are NOT score weights.
* APBS electrostatic potential is an auxiliary descriptor only; it is NOT a score weight.
* The historical 8 Å context statistic is descriptive only; it is NOT a score weight.
* No adsorption capacity, unique orientation, conformational change, or adsorption free
  energy is predicted.

Stack: FastAPI + Biopython + optional pdb2pqr/APBS + optional DSSP.
"""

from __future__ import annotations

import json
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field
from Bio.PDB import PDBParser, PDBIO, Select, DSSP
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
log = logging.getLogger("interfacescout")

APP_VERSION = "5.1.0-publication-freeze"
SC_RSA_THRESHOLD = 0.05
SASA_PROBE_A = 1.40
SASA_POINTS = 200
PATCH_RADII_A = (5.0, 8.0)
CONTEXT_RADIUS_A = 8.0  # auxiliary statistic only
PKA_SENSITIVITY_WINDOW = 1.0

# Frozen spatial-pair selection. The 5/8 Å pair was selected from a
# developmental, adsorption-independent radius-pair audit as a compromise
# between local chemical resolution and top-patch positional stability.
# The 8/10 Å pair showed the highest broad top-10 spatial coherence but
# substantially poorer top-patch-centre stability; it is therefore retained
# as an audit comparator rather than the canonical pair.
PATCH_PAIR_AUDIT = {
    "tested_radii_A": [5.0, 8.0, 10.0, 12.0, 15.0],
    "selected_pair_A": [5.0, 8.0],
    "selection_basis": "developmental adsorption-independent spatial robustness audit",
    "broad_coherence_comparator_A": [8.0, 10.0],
}

# NACCESS-style extended Ala-X-Ala side-chain reference areas used in the SI (Å²)
SIDECHAIN_REF_ASA: Dict[str, float] = {
    "ALA": 69.23, "ARG": 200.35, "ASN": 106.25, "ASP": 102.06,
    "CYS": 96.69, "GLN": 140.58, "GLU": 134.61, "GLY": 32.28,
    "HIS": 147.00, "ILE": 137.91, "LEU": 140.76, "LYS": 162.50,
    "MET": 156.08, "PHE": 163.90, "PRO": 119.65, "SER": 78.16,
    "THR": 101.67, "TRP": 210.89, "TYR": 176.61, "VAL": 114.14,
}

# Side-chain pKa reference values used by the manuscript/SI.
PKA: Dict[str, float] = {
    "ASP": 3.9, "GLU": 4.1, "HIS": 6.5, "CYS": 8.3,
    "TYR": 10.1, "LYS": 10.5, "ARG": 12.5,
}

# Generic charged-state descriptor sign: +1 basic/protonated; -1 acidic/deprotonated.
ION_SIGN: Dict[str, int] = {
    "ARG": +1, "LYS": +1, "HIS": +1,
    "ASP": -1, "GLU": -1, "CYS": -1, "TYR": -1,
}

# Literature reference energies are retained strictly as metadata.
# Each entry: residue -> (Ebase kcal/mol, mechanism, required_state)
# required_state is one of: None, "protonated", "deprotonated", "auxiliary_only".
CHEMISTRIES: Dict[str, Dict[str, Any]] = {
    "cationic": {
        "label": "Cationic surface",
        "surface_group": "Amine / cationic",
        "description": "Protein-side compatibility with positively charged/amine-rich interfaces.",
        "favorable": {
            "ASP": (-4.0, "electrostatic + H-bond", "deprotonated"),
            "GLU": (-4.0, "electrostatic + H-bond", "deprotonated"),
            "TYR": (-1.0, "H-bond", None),
            "SER": (-0.8, "H-bond", None), "THR": (-0.8, "H-bond", None),
            "ASN": (-0.8, "H-bond", None), "GLN": (-0.8, "H-bond", None),
        },
        "repulsive": {},
        "expected_phi_sign": "negative",
    },
    "anionic": {
        "label": "Anionic surface",
        "surface_group": "Carboxyl / anionic",
        "description": "Protein-side compatibility with negatively charged/carboxyl-rich interfaces.",
        "favorable": {
            "LYS": (-4.0, "electrostatic + H-bond", "protonated"),
            "ARG": (-4.0, "electrostatic + H-bond", "protonated"),
            "HIS": (-2.0, "electrostatic (if protonated)", "protonated"),
            "SER": (-0.8, "H-bond donor", None), "THR": (-0.8, "H-bond donor", None),
        },
        "repulsive": {
            "ASP": (+2.0, "like-charge repulsion", "deprotonated"),
            "GLU": (+2.0, "like-charge repulsion", "deprotonated"),
        },
        "expected_phi_sign": "positive",
    },
    "hbond_donor": {
        "label": "H-bond donor surface",
        "surface_group": "Hydroxyl",
        "description": "Compatibility with surfaces capable of donating hydrogen bonds.",
        "favorable": {
            "TYR": (-1.0, "H-bond", None), "ASP": (-1.0, "H-bond acceptor", None),
            "GLU": (-1.0, "H-bond acceptor", None),
            "SER": (-0.8, "H-bond", None), "THR": (-0.8, "H-bond", None),
            "ASN": (-0.8, "H-bond", None), "GLN": (-0.8, "H-bond", None),
            "LYS": (-0.8, "H-bond", None), "ARG": (-0.8, "H-bond", None),
        },
        "repulsive": {}, "expected_phi_sign": None,
    },
    "hbond_acceptor": {
        "label": "H-bond acceptor surface",
        "surface_group": "Carbonyl",
        "description": "Compatibility with surfaces capable of accepting hydrogen bonds.",
        "favorable": {
            "LYS": (-1.2, "H-bond", None), "ARG": (-1.2, "H-bond", None),
            "SER": (-1.0, "H-bond donor", None), "THR": (-1.0, "H-bond donor", None),
            "TRP": (-1.0, "H-bond (indole NH)", None),
        },
        "repulsive": {}, "expected_phi_sign": None,
    },
    "pi_carbon": {
        "label": "π / carbon-like surface",
        "surface_group": "Aromatic / graphitic π",
        "description": "Graphitic/aromatic compatibility through π-associated, cation–π, and related contacts.",
        "favorable": {
            "PHE": (-4.0, "π–π stacking", None), "TRP": (-5.0, "π–π stacking", None),
            "TYR": (-4.0, "π–π stacking", None),
            "ARG": (-3.0, "cation–π (guanidinium)", "protonated"),
            "LYS": (-2.0, "cation–π", "protonated"),
            # Histidine protonated fraction is exported as auxiliary descriptor only;
            # mixed π/cation–π mechanisms are not given an arbitrary primary weight.
            "HIS": (-2.5, "π–π / cation–π", "auxiliary_only"),
        },
        "repulsive": {}, "expected_phi_sign": None,
    },
    "hydrophobic": {
        "label": "Hydrophobic surface",
        "surface_group": "Methyl / nonpolar",
        "description": "Compatibility with nonpolar interfaces through hydrophobic/CH–π contacts.",
        "favorable": {
            "LEU": (-2.2, "hydrophobic contact", None), "ILE": (-2.2, "hydrophobic contact", None),
            "VAL": (-1.8, "hydrophobic contact", None), "MET": (-1.8, "hydrophobic contact", None),
            "PHE": (-2.0, "hydrophobic + CH–π", None), "TRP": (-2.2, "hydrophobic + CH–π", None),
            "ALA": (-1.0, "hydrophobic contact", None), "PRO": (-1.0, "hydrophobic contact", None),
        },
        "repulsive": {}, "expected_phi_sign": None,
    },
    "oxide": {
        "label": "Oxide surface",
        "surface_group": "Metal oxide (Ti–O/Zr–O-like)",
        "description": "Compatibility with oxide interfaces through carboxylate coordination and H-bonding.",
        "favorable": {
            "ASP": (-8.0, "carboxylate–metal coordination", "deprotonated"),
            "GLU": (-8.0, "carboxylate–metal coordination", "deprotonated"),
            "SER": (-3.0, "hydroxyl–oxide H-bond", None), "THR": (-3.0, "hydroxyl–oxide H-bond", None),
            "TYR": (-3.0, "phenol–oxide H-bond", None),
        },
        "repulsive": {}, "expected_phi_sign": None,
    },
    "hydroxyapatite": {
        "label": "Hydroxyapatite / Ca²⁺ surface",
        "surface_group": "Ca²⁺ / calcium phosphate",
        "description": "Compatibility with calcium-rich interfaces through carboxylate coordination and complementary contacts.",
        "favorable": {
            "ASP": (-10.0, "carboxylate–Ca coordination", "deprotonated"),
            "GLU": (-10.0, "carboxylate–Ca coordination", "deprotonated"),
            "SER": (-6.0, "hydroxyl/phospho-Ser binding", None),
            "TRP": (-6.0, "cation–π", None), "TYR": (-5.0, "cation–π", None),
            "PHE": (-4.0, "cation–π", None), "THR": (-4.0, "H-bond", None),
        },
        "repulsive": {}, "expected_phi_sign": None,
    },
    "metal_coord": {
        "label": "Transition-metal coordination",
        "surface_group": "Transition-metal sites",
        "description": "Compatibility with accessible transition-metal sites through coordinating side chains.",
        "favorable": {
            "HIS": (-12.0, "imidazole–metal coordination", None),
            "CYS": (-10.0, "thiolate–metal coordination", None),
            "ASP": (-6.0, "carboxylate–metal coordination", "deprotonated"),
            "GLU": (-6.0, "carboxylate–metal coordination", "deprotonated"),
            "MET": (-4.0, "thioether–metal (soft)", None),
        },
        "repulsive": {}, "expected_phi_sign": None,
    },
    "gold": {
        "label": "Gold affinity",
        "surface_group": "Au / soft-metal",
        "description": "Compatibility with gold/soft-metal interfaces, dominated by sulfur-containing residues.",
        "favorable": {
            "CYS": (-45.0, "Au–S semi-covalent", None),
            "MET": (-5.0, "weak S coordination", None),
        },
        "repulsive": {}, "expected_phi_sign": None,
    },
    "phosphate": {
        "label": "Phosphate surface",
        "surface_group": "PO₄",
        "description": "Compatibility with phosphate-rich interfaces through electrostatic and H-bond contacts.",
        "favorable": {
            "ARG": (-7.0, "bidentate phosphate–guanidinium", "protonated"),
            "LYS": (-6.0, "phosphate–ammonium electrostatic", "protonated"),
            "HIS": (-3.0, "electrostatic", "protonated"),
            "SER": (-1.0, "H-bond", None),
        },
        "repulsive": {}, "expected_phi_sign": "positive",
    },
}

FEATURE_RESIDUES: Dict[str, set] = {
    "charge_pos": {"ARG", "LYS", "HIS"},
    "charge_neg": {"ASP", "GLU", "CYS", "TYR"},
    "hbond_donor": {"SER", "THR", "TYR", "ASN", "GLN", "TRP", "LYS", "ARG", "HIS"},
    "hbond_acceptor": {"ASP", "GLU", "ASN", "GLN", "SER", "THR", "TYR", "HIS"},
    "hydrophobic": {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"},
    "aromatic": {"PHE", "TRP", "TYR", "HIS"},
    "metal_binding": {"HIS", "CYS", "ASP", "GLU", "MET"},
    "thiol": {"CYS"},
    "carboxyl": {"ASP", "GLU"},
    "amine": {"LYS", "ARG", "HIS"},
}

FEATURE_INFO = {
    "charge_pos": ("Positive charge", "Protonated basic side chains"),
    "charge_neg": ("Negative charge", "Deprotonated acidic side chains"),
    "hbond_donor": ("H-bond donor", "Side chains capable of donating H-bonds"),
    "hbond_acceptor": ("H-bond acceptor", "Side chains capable of accepting H-bonds"),
    "hydrophobic": ("Hydrophobic", "Nonpolar/aliphatic/aromatic side chains"),
    "aromatic": ("Aromatic", "Aromatic/imidazole ring systems"),
    "metal_binding": ("Metal-binding", "Common transition-metal coordinating side chains"),
    "thiol": ("Thiol", "Cysteine sulfur chemistry"),
    "carboxyl": ("Carboxyl", "Asp/Glu carboxylate chemistry"),
    "amine": ("Amine/basic", "Basic nitrogen-containing side chains"),
}


def _find_binary(name: str) -> Optional[str]:
    env_key = f"{name.upper().replace('-', '_')}_PATH"
    p = os.environ.get(env_key)
    if p and Path(p).exists():
        return p
    p = shutil.which(name)
    if p:
        return p
    # apbs_binary package support
    if name == "apbs":
        try:
            from apbs_binary import APBS_BIN_PATH
            if Path(str(APBS_BIN_PATH)).exists():
                return str(APBS_BIN_PATH)
        except Exception:
            pass
    return None

PDB2PQR = _find_binary("pdb2pqr")
APBS = _find_binary("apbs")
MKDSSP = _find_binary("mkdssp") or _find_binary("dssp")


class EnvParams(BaseModel):
    pH: float = Field(7.4, ge=0.0, le=14.0)
    ionic: float = Field(150.0, ge=0.0, description="mM")
    temp: float = Field(298.0, gt=0.0, description="K")


class AnalyzeRequest(BaseModel):
    pdb_id: Optional[str] = None
    pdb_text: Optional[str] = None
    chain: Optional[str] = None  # one chain ID, e.g. A; blank/None = all chains
    env: EnvParams = EnvParams()


def charged_fraction(res_name: str, pH: float) -> float:
    """Generic fraction in the charged side-chain state."""
    rn = res_name.upper()
    if rn not in PKA or rn not in ION_SIGN:
        return 0.0
    pka = PKA[rn]
    if ION_SIGN[rn] > 0:
        return 1.0 / (1.0 + 10.0 ** (pH - pka))
    return 1.0 / (1.0 + 10.0 ** (pka - pH))


def state_availability(res_name: str, state: Optional[str], pH: float) -> float:
    if state in (None, "auxiliary_only"):
        return 1.0
    rn = res_name.upper()
    pka = PKA.get(rn)
    if pka is None:
        return 1.0
    if state == "protonated":
        return 1.0 / (1.0 + 10.0 ** (pH - pka))
    if state == "deprotonated":
        return 1.0 / (1.0 + 10.0 ** (pka - pH))
    return 1.0


def residue_charge_descriptor(res_name: str, pH: float) -> float:
    rn = res_name.upper()
    if rn not in ION_SIGN:
        return 0.0
    f = charged_fraction(rn, pH)
    return float(ION_SIGN[rn] * f)


class ChainSelect(Select):
    def __init__(self, chain_id: str):
        super().__init__()
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return 1 if str(chain.id) == self.chain_id else 0

    def accept_residue(self, residue):
        return 1 if is_aa(residue, standard=True) else 0


class StandardAASelect(Select):
    def accept_residue(self, residue):
        return 1 if is_aa(residue, standard=True) else 0


def prepare_input_pdb(req: AnalyzeRequest, workdir: Path) -> Tuple[Path, str]:
    raw = workdir / "raw.pdb"
    if req.pdb_text:
        raw.write_text(req.pdb_text)
    elif req.pdb_id:
        pdb_id = req.pdb_id.strip().upper()
        if not pdb_id:
            raise HTTPException(400, "Empty PDB ID")
        try:
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", raw)
        except Exception as exc:
            raise HTTPException(400, f"Could not download PDB {pdb_id}: {exc}")
    else:
        raise HTTPException(400, "Provide pdb_id or pdb_text")

    parser = PDBParser(QUIET=True)
    try:
        struct = parser.get_structure("input", str(raw))
    except Exception as exc:
        raise HTTPException(400, f"Could not parse PDB: {exc}")

    available = [str(c.id) for c in next(struct.get_models())]
    chain = (req.chain or "").strip()
    if chain and chain not in available:
        raise HTTPException(400, f"Chain '{chain}' not found. Available chains: {', '.join(available)}")

    out = workdir / "analysis.pdb"
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(out), ChainSelect(chain) if chain else StandardAASelect())
    return out, chain or "ALL"


def compute_dssp(struct, pdb: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not MKDSSP:
        return out
    try:
        dssp = DSSP(struct[0], str(pdb), dssp=MKDSSP)
        for key in dssp.keys():
            chain = str(key[0])
            resid = key[1]
            resseq = int(resid[1])
            icode = str(resid[2]).strip()
            ss = dssp[key][2]
            ss_simple = "H" if ss in {"H", "G", "I"} else ("E" if ss in {"E", "B"} else "C")
            out[f"{chain}:{resseq}:{icode}"] = ss_simple
    except Exception as exc:
        log.warning("DSSP unavailable for this structure: %s", exc)
    return out


def sidechain_atoms(residue) -> List[Any]:
    rn = residue.get_resname().strip().upper()
    if rn == "GLY":
        return [residue["CA"]] if residue.has_id("CA") else []
    backbone = {"N", "CA", "C", "O", "OXT"}
    return [a for a in residue.get_atoms() if a.get_name().strip().upper() not in backbone and a.element != "H"]


def build_surface_residues(pdb: Path, pH: float) -> Tuple[Any, List[Dict[str, Any]], Dict[str, List[Any]], Dict[str, str]]:
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("protein", str(pdb))
    # Atom-level SASA is required to isolate side-chain accessibility.
    sr = ShrakeRupley(probe_radius=SASA_PROBE_A, n_points=SASA_POINTS)
    sr.compute(struct, level="A")
    dssp = compute_dssp(struct, pdb)

    residues: List[Dict[str, Any]] = []
    atoms_by_key: Dict[str, List[Any]] = {}
    for model in struct:
        for chain in model:
            for res in chain:
                if not is_aa(res, standard=True) or not res.has_id("CA"):
                    continue
                rn = res.get_resname().strip().upper()
                if rn not in SIDECHAIN_REF_ASA:
                    continue
                resseq = int(res.id[1])
                icode = str(res.id[2]).strip()
                cid = str(chain.id)
                key = f"{cid}:{resseq}:{icode}"
                sc_atoms = sidechain_atoms(res)
                sc_sasa = float(sum(float(getattr(a, "sasa", 0.0)) for a in sc_atoms))
                total_sasa = float(sum(float(getattr(a, "sasa", 0.0)) for a in res.get_atoms()))
                scrsa_raw = sc_sasa / SIDECHAIN_REF_ASA[rn] if SIDECHAIN_REF_ASA[rn] > 0 else 0.0
                scrsa_clip = min(max(scrsa_raw, 0.0), 1.0)
                ca = res["CA"]
                sc_coords = [np.asarray(a.coord, dtype=float) for a in sc_atoms]
                if sc_coords:
                    centroid = np.mean(np.vstack(sc_coords), axis=0)
                else:
                    centroid = np.asarray(ca.coord, dtype=float)
                residues.append({
                    "key": key,
                    "res_name": rn,
                    "res_seq": resseq,
                    "icode": icode,
                    "chain": cid,
                    "x": float(ca.coord[0]), "y": float(ca.coord[1]), "z": float(ca.coord[2]),
                    "sc_x": float(centroid[0]), "sc_y": float(centroid[1]), "sc_z": float(centroid[2]),
                    "total_sasa": round(total_sasa, 3),
                    "sidechain_sasa": round(sc_sasa, 3),
                    "scrsa_raw": round(scrsa_raw, 5),
                    "scrsa": round(scrsa_clip, 5),
                    "surface_exposed": bool(scrsa_raw >= SC_RSA_THRESHOLD),
                    "charge_fraction": round(charged_fraction(rn, pH), 5),
                    "charge_descriptor": round(residue_charge_descriptor(rn, pH), 5),
                    "pka": PKA.get(rn),
                    "ionization_sensitive": bool(rn in PKA and abs(pH - PKA[rn]) <= PKA_SENSITIVITY_WINDOW),
                    "ss": dssp.get(key, "C"),
                    "bfactor": round(float(ca.bfactor), 3),
                    "phi": None,
                    "n_neighbors_8A": 0,
                })
                atoms_by_key[key] = sc_atoms
        break  # publication workflow uses first MODEL only

    return struct, residues, atoms_by_key, dssp


def run_pdb2pqr(pdb: Path, pqr: Path, pH: float) -> None:
    if not PDB2PQR:
        raise RuntimeError("pdb2pqr not found")
    cmd = [PDB2PQR, "--ff=PARSE", f"--with-ph={pH:.2f}", "--titration-state-method=propka",
           "--drop-water", "--keep-chain", str(pdb), str(pqr)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if r.returncode != 0 or not pqr.exists():
        raise RuntimeError(f"pdb2pqr failed: {r.stderr[-500:]}")


def build_apbs_input(pqr: Path, dx_stem: Path, ionic_M: float, temp: float) -> Path:
    xs: List[float] = []; ys: List[float] = []; zs: List[float] = []
    for line in pqr.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
        except Exception:
            continue
    if not xs:
        raise RuntimeError("PQR contains no atoms")

    pad = 15.0
    spacing = 0.5
    def dim(lo: float, hi: float) -> Tuple[int, float]:
        n = max(33, int(math.ceil((hi - lo + 2.0 * pad) / spacing)))
        if n % 2 == 0:
            n += 1
        return n, (lo + hi) / 2.0

    nx, cx = dim(min(xs), max(xs)); ny, cy = dim(min(ys), max(ys)); nz, cz = dim(min(zs), max(zs))
    gx, gy, gz = nx * spacing, ny * spacing, nz * spacing
    inp = pqr.parent / "apbs.in"
    inp.write_text(f"""read
    mol pqr {pqr}
end
elec
    mg-manual
    mol 1
    dime {nx} {ny} {nz}
    glen {gx:.2f} {gy:.2f} {gz:.2f}
    gcent {cx:.3f} {cy:.3f} {cz:.3f}
    lpbe
    bcfl sdh
    ion charge +1 conc {ionic_M:.5f} radius 2.0
    ion charge -1 conc {ionic_M:.5f} radius 2.0
    pdie 2.0
    sdie 78.54
    srfm smol
    sdens 10.0
    chgm spl2
    srad 1.4
    swin 0.3
    temp {temp:.2f}
    calcenergy no
    calcforce no
    write pot dx {dx_stem}
end
quit
""")
    return inp


def _apbs_env() -> Dict[str, str]:
    env = dict(os.environ)
    try:
        import apbs_binary
        lib_dir = str(getattr(apbs_binary, "LIB_DIR", ""))
        if lib_dir:
            key = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
            env[key] = lib_dir + os.pathsep + env.get(key, "")
    except Exception:
        pass
    if APBS:
        env["PATH"] = str(Path(APBS).parent) + os.pathsep + env.get("PATH", "")
    return env


def run_apbs(inp: Path, workdir: Path) -> Path:
    if not APBS:
        raise RuntimeError("APBS not found")
    r = subprocess.run([APBS, str(inp)], capture_output=True, text=True, cwd=str(workdir),
                       timeout=600, env=_apbs_env())
    dxs = list(workdir.glob("*.dx"))
    if not dxs:
        raise RuntimeError(f"APBS produced no .dx (exit={r.returncode}): {r.stderr[-400:]}")
    return dxs[0]


def parse_dx(dx: Path) -> Dict[str, Any]:
    lines = dx.read_text(errors="replace").splitlines()
    origin = None; deltas: List[List[float]] = []; counts = None; data: List[float] = []
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith("object 1 class gridpositions counts"):
            p = ln.split(); counts = [int(p[-3]), int(p[-2]), int(p[-1])]
        elif ln.startswith("origin"):
            origin = list(map(float, ln.split()[1:4]))
        elif ln.startswith("delta") and len(deltas) < 3:
            deltas.append(list(map(float, ln.split()[1:4])))
        elif ln.startswith("object 3"):
            i += 1
            while i < len(lines):
                row = lines[i].strip()
                if row.startswith(("attribute", "object")):
                    break
                if row:
                    try:
                        data.extend(float(x) for x in row.split())
                    except ValueError:
                        pass
                i += 1
            continue
        i += 1
    if origin is None or counts is None or len(deltas) < 3:
        raise RuntimeError("DX parse failed")
    return {"origin": origin, "delta": deltas, "counts": counts, "data": np.asarray(data, dtype=np.float64)}


def phi_at(grid: Dict[str, Any], x: float, y: float, z: float) -> Optional[float]:
    o = grid["origin"]; d = grid["delta"]; nx, ny, nz = grid["counts"]; dat = grid["data"]
    dx, dy, dz = d[0][0], d[1][1], d[2][2]
    if dx == 0 or dy == 0 or dz == 0:
        return None
    fx, fy, fz = (x - o[0]) / dx, (y - o[1]) / dy, (z - o[2]) / dz
    ix, iy, iz = int(math.floor(fx)), int(math.floor(fy)), int(math.floor(fz))
    if ix < 0 or iy < 0 or iz < 0 or ix >= nx - 1 or iy >= ny - 1 or iz >= nz - 1:
        return None
    tx, ty, tz = fx - ix, fy - iy, fz - iz
    def idx(a: int, b: int, c: int) -> int:
        return a * ny * nz + b * nz + c
    vals = [
        dat[idx(ix,iy,iz)], dat[idx(ix+1,iy,iz)], dat[idx(ix,iy+1,iz)], dat[idx(ix,iy,iz+1)],
        dat[idx(ix+1,iy+1,iz)], dat[idx(ix+1,iy,iz+1)], dat[idx(ix,iy+1,iz+1)], dat[idx(ix+1,iy+1,iz+1)],
    ]
    ws = [
        (1-tx)*(1-ty)*(1-tz), tx*(1-ty)*(1-tz), (1-tx)*ty*(1-tz), (1-tx)*(1-ty)*tz,
        tx*ty*(1-tz), tx*(1-ty)*tz, (1-tx)*ty*tz, tx*ty*tz,
    ]
    return float(sum(v*w for v, w in zip(vals, ws)))


def attach_apbs_auxiliary(pdb: Path, residues: List[Dict[str, Any]], atoms_by_key: Dict[str, List[Any]],
                          env: EnvParams, workdir: Path) -> str:
    if not PDB2PQR:
        return "unavailable"
    try:
        pqr = workdir / "protein.pqr"
        run_pdb2pqr(pdb, pqr, env.pH)
        if not APBS:
            return "pdb2pqr_only"
        inp = build_apbs_input(pqr, workdir / "potential", env.ionic / 1000.0, env.temp)
        dx = run_apbs(inp, workdir)
        grid = parse_dx(dx)
        for r in residues:
            atoms = atoms_by_key.get(r["key"], [])
            weighted = 0.0; wsum = 0.0
            for atom in atoms:
                w = max(float(getattr(atom, "sasa", 0.0)), 0.0)
                if w <= 0:
                    continue
                ph = phi_at(grid, float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2]))
                if ph is None:
                    continue
                weighted += w * ph; wsum += w
            if wsum > 0:
                r["phi"] = round(weighted / wsum, 5)
            else:
                # glycine or fully buried side chain: use side-chain centroid as fallback descriptor
                ph = phi_at(grid, r["sc_x"], r["sc_y"], r["sc_z"])
                r["phi"] = round(ph, 5) if ph is not None else None
        return "APBS_LinearPB_auxiliary"
    except Exception as exc:
        log.warning("APBS auxiliary descriptor unavailable: %s", exc)
        return "pdb2pqr_only"


def build_distances(surface: List[Dict[str, Any]]) -> np.ndarray:
    coords = np.asarray([[r["x"], r["y"], r["z"]] for r in surface], dtype=float)
    if len(coords) == 0:
        return np.empty((0, 0), dtype=float)
    dif = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.sum(dif * dif, axis=2))


def normalize_to_100(values: np.ndarray) -> np.ndarray:
    mx = float(np.max(values)) if values.size else 0.0
    if mx <= 0:
        return np.zeros_like(values, dtype=float)
    return values / mx * 100.0


def electrostatic_relation(phi: Optional[float], expected: Optional[str]) -> Optional[str]:
    if phi is None or expected is None:
        return None
    if abs(phi) < 0.25:
        return "near-neutral"
    if expected == "negative":
        return "complementary" if phi < 0 else "opposing"
    if expected == "positive":
        return "complementary" if phi > 0 else "opposing"
    return None


def chemistry_map(surface: List[Dict[str, Any]], distances: np.ndarray, key: str, pH: float) -> Dict[str, Any]:
    meta = CHEMISTRIES[key]
    fav_defs = meta["favorable"]
    rep_defs = meta["repulsive"]
    n = len(surface)
    local = np.zeros(n, dtype=float)
    rep_local = np.zeros(n, dtype=float)

    # Canonical residue-level scores.
    for i, r in enumerate(surface):
        rn = r["res_name"]
        if rn in fav_defs:
            _, _, state = fav_defs[rn]
            fstate = state_availability(rn, state, pH)
            local[i] = r["scrsa"] * fstate
        if rn in rep_defs:
            _, _, state = rep_defs[rn]
            fstate = state_availability(rn, state, pH)
            rep_local[i] = r["scrsa"] * fstate

    prop = normalize_to_100(local)
    rep_prop = normalize_to_100(rep_local)

    # 5 Å and 8 Å patch-density maps over all surface-exposed residues.
    dens: Dict[float, np.ndarray] = {}
    dens_norm: Dict[float, np.ndarray] = {}
    for R in PATCH_RADII_A:
        if n:
            mask = distances <= R
            d = mask.astype(float) @ local
        else:
            d = np.zeros(0, dtype=float)
        dens[R] = d
        mx = float(np.max(d)) if d.size else 0.0
        dens_norm[R] = (d / mx) if mx > 0 else np.zeros_like(d)

    if n:
        persistence = 100.0 * np.minimum(dens_norm[5.0], dens_norm[8.0])
        geom = 100.0 * np.sqrt(dens_norm[5.0] * dens_norm[8.0])
    else:
        persistence = np.zeros(0); geom = np.zeros(0)

    # Top patch centres are surface residues, not necessarily members of the chemistry class.
    patch_centers = []
    for i, r in enumerate(surface):
        if persistence[i] <= 0:
            continue
        members8 = [surface[j]["key"] for j in range(n) if distances[i, j] <= 8.0 and local[j] > 0]
        patch_centers.append({
            "center_key": r["key"], "res_name": r["res_name"], "res_seq": r["res_seq"],
            "icode": r["icode"], "chain": r["chain"],
            "density_5A_raw": round(float(dens[5.0][i]), 6),
            "density_8A_raw": round(float(dens[8.0][i]), 6),
            "density_5A_norm": round(float(dens_norm[5.0][i] * 100.0), 3),
            "density_8A_norm": round(float(dens_norm[8.0][i] * 100.0), 3),
            "multiscale_persistence": round(float(persistence[i]), 3),
            "multiscale_geomean": round(float(geom[i]), 3),
            "compatible_members_8A": members8,
        })
    patch_centers.sort(key=lambda x: (-x["multiscale_persistence"], -x["multiscale_geomean"], x["center_key"]))

    members = []
    repulsive = []
    for i, r in enumerate(surface):
        rn = r["res_name"]
        if rn in fav_defs:
            ebase, mechanism, state = fav_defs[rn]
            fstate = state_availability(rn, state, pH)
            aux_state = charged_fraction(rn, pH) if state == "auxiliary_only" else None
            members.append({
                "key": r["key"], "res_name": rn, "res_seq": r["res_seq"], "icode": r["icode"], "chain": r["chain"],
                "sidechain_sasa": r["sidechain_sasa"], "scrsa_raw": r["scrsa_raw"], "scrsa": r["scrsa"],
                "state_requirement": state or "none", "state_availability": round(float(fstate), 5),
                "auxiliary_charged_fraction": round(float(aux_state), 5) if aux_state is not None else None,
                "local_score": round(float(local[i]), 6),
                "propensity": round(float(prop[i]), 3),
                "patch_density_5A": round(float(dens_norm[5.0][i] * 100.0), 3),
                "patch_density_8A": round(float(dens_norm[8.0][i] * 100.0), 3),
                "multiscale_persistence": round(float(persistence[i]), 3),
                "multiscale_geomean": round(float(geom[i]), 3),
                "phi": r["phi"],
                "electrostatic_relation": electrostatic_relation(r["phi"], meta.get("expected_phi_sign")),
                "ebase_metadata_kcal_mol": ebase,
                "mechanism": mechanism,
                "ionization_sensitive": r["ionization_sensitive"],
                "ss": r["ss"], "bfactor": r["bfactor"], "n_neighbors_8A": r["n_neighbors_8A"],
            })
        if rn in rep_defs:
            ebase, mechanism, state = rep_defs[rn]
            fstate = state_availability(rn, state, pH)
            repulsive.append({
                "key": r["key"], "res_name": rn, "res_seq": r["res_seq"], "icode": r["icode"], "chain": r["chain"],
                "scrsa": r["scrsa"], "state_requirement": state or "none",
                "state_availability": round(float(fstate), 5),
                "repulsion_score": round(float(rep_local[i]), 6),
                "repulsion_propensity": round(float(rep_prop[i]), 3),
                "phi": r["phi"], "ebase_metadata_kcal_mol": ebase, "mechanism": mechanism,
            })

    members.sort(key=lambda x: (-x["propensity"], -x["multiscale_persistence"], x["key"]))
    repulsive.sort(key=lambda x: (-x["repulsion_propensity"], x["key"]))

    return {
        "kind": "compatibility",
        "label": meta["label"],
        "surface_group": meta["surface_group"],
        "description": meta["description"],
        "n_favorable_residues": len(members),
        "n_repulsive_residues": len(repulsive),
        "top_patch": patch_centers[0] if patch_centers else None,
        "top_patches": patch_centers[:10],
        "residues": members,
        "repulsive_residues": repulsive,
        "notes": {
            "primary_score": "L = scRSA × binary membership × required-state availability",
            "normalization": "within-map maximum = 100",
            "patch": "canonical multiscale persistence = 100 × min(normalized 5 Å density, normalized 8 Å density)",
            "ebase": "metadata only; excluded from ranking",
            "apbs": "auxiliary descriptor only; excluded from ranking",
        },
    }


def feature_map(surface: List[Dict[str, Any]], pH: float) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for feat, memberset in FEATURE_RESIDUES.items():
        rows = []
        for r in surface:
            rn = r["res_name"]
            if rn not in memberset:
                continue
            # Charge features are condition-aware.
            if feat == "charge_pos" and residue_charge_descriptor(rn, pH) <= 0.05:
                continue
            if feat == "charge_neg" and residue_charge_descriptor(rn, pH) >= -0.05:
                continue
            rows.append({
                "key": r["key"], "res_name": rn, "res_seq": r["res_seq"], "icode": r["icode"], "chain": r["chain"],
                "sidechain_sasa": r["sidechain_sasa"], "scrsa_raw": r["scrsa_raw"], "scrsa": r["scrsa"],
                "charge_descriptor": r["charge_descriptor"], "phi": r["phi"], "ss": r["ss"],
            })
        rows.sort(key=lambda x: (-x["scrsa"], x["key"]))
        label, mechanism = FEATURE_INFO[feat]
        out[feat] = {"kind": "feature", "label": label, "mechanism": mechanism, "n": len(rows), "residues": rows}
    return out


def analyze(req: AnalyzeRequest) -> Dict[str, Any]:
    workdir = Path(tempfile.mkdtemp(prefix="interfacescout_"))
    try:
        pdb, selected_chain = prepare_input_pdb(req, workdir)
        struct, all_residues, atoms_by_key, _ = build_surface_residues(pdb, req.env.pH)
        if not all_residues:
            raise HTTPException(400, "No standard amino-acid residues were found in the selected structure/chain")

        electrostatics = attach_apbs_auxiliary(pdb, all_residues, atoms_by_key, req.env, workdir)
        surface = [r for r in all_residues if r["surface_exposed"]]
        distances = build_distances(surface)

        # Historical 8 Å local context retained only as an exported descriptor.
        if len(surface):
            for i, r in enumerate(surface):
                r["n_neighbors_8A"] = int(np.sum((distances[i] <= CONTEXT_RADIUS_A) & (distances[i] > 0.0)))

        chem = {k: chemistry_map(surface, distances, k, req.env.pH) for k in CHEMISTRIES}
        features = feature_map(surface, req.env.pH)

        n_atoms = sum(1 for m in struct for c in m for res in c if is_aa(res, standard=True) for _ in res.get_atoms())
        n_res = len(all_residues)
        result = {
            "status": "ok",
            "version": APP_VERSION,
            "model": "InterfaceScout publication-freeze compatibility mapping",
            "scope": {
                "predicts": "protein-side residue/patch compatibility hypotheses for generalized interface chemistries",
                "does_not_predict": ["adsorption capacity", "absolute adsorption free energy", "unique adsorption orientation",
                                     "adsorption-induced conformational change", "material-side transport/porosity/hydration"],
            },
            "settings": {
                "chain": selected_chain,
                "pH": req.env.pH, "ionic_mM": req.env.ionic, "temperature_K": req.env.temp,
                "sasa_probe_A": SASA_PROBE_A, "sasa_points_per_atom": SASA_POINTS,
                "scrsa_threshold": SC_RSA_THRESHOLD,
                "patch_radii_A": list(PATCH_RADII_A),
                "patch_pair_selection": PATCH_PAIR_AUDIT,
                "context_radius_A_auxiliary": CONTEXT_RADIUS_A,
            },
            "stats": {
                "n_atoms": n_atoms, "n_residues": n_res, "n_surface_res": len(surface),
                "electrostatics": electrostatics,
                "pdb2pqr": bool(PDB2PQR), "apbs": bool(APBS), "dssp": bool(MKDSSP),
            },
            "chemistry_list": list(CHEMISTRIES.keys()),
            "chemistries": chem,
            "feature_list": list(FEATURE_RESIDUES.keys()),
            "features": features,
            "surface_residues": surface,
            "reference_sidechain_asa": SIDECHAIN_REF_ASA,
        }
        return result
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


app = FastAPI(title="InterfaceScout Backend", version=APP_VERSION)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])


@app.get("/health")
def health():
    return {
        "status": "ok", "version": APP_VERSION,
        "pdb2pqr": bool(PDB2PQR), "apbs": bool(APBS), "dssp": bool(MKDSSP),
        "canonical_model": "scRSA × membership × state availability; 5/8 Å persistence",
    }


@app.post("/analyze_surface")
def analyze_surface(req: AnalyzeRequest):
    try:
        return analyze(req)
    except HTTPException:
        raise
    except Exception as exc:
        log.exception("Analysis failed")
        raise HTTPException(500, str(exc))


@app.get("/model_spec")
def model_spec():
    """Machine-readable publication-freeze model specification."""
    spec = {
        "version": APP_VERSION,
        "sasa": {"algorithm": "Shrake-Rupley", "probe_A": SASA_PROBE_A, "points_per_atom": SASA_POINTS},
        "exposure": {"metric": "side-chain relative solvent accessibility", "threshold": SC_RSA_THRESHOLD,
                     "glycine": "Cα proxy"},
        "residue_score": "L_i,c = I_i,c * scRSA_i * f_state,i,c(pH)",
        "propensity": "P_i,c = 100 * L_i,c / max(L_c,fav)",
        "patch": {"radii_A": list(PATCH_RADII_A),
                  "pair_selection": PATCH_PAIR_AUDIT,
                  "persistence": "100 * min(D5_norm, D8_norm)",
                  "geometric_mean": "secondary diagnostic"},
        "excluded_from_primary_score": ["Ebase magnitude", "APBS potential", "8 Å context statistic"],
        "chemistry_classes": list(CHEMISTRIES.keys()),
    }
    return JSONResponse(spec)


@app.get("/")
def root():
    # Serve the final frontend automatically when index.html is placed next to this backend.
    index = Path(__file__).with_name("index.html")
    if index.exists():
        return FileResponse(index)
    return {"name": "InterfaceScout", "version": APP_VERSION, "frontend": "Place index.html beside this file."}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
