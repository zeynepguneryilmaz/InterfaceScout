"""Canonical InterfaceScout chemistry/SASA core.

This module contains only the shared protein preparation-independent numerical
kernel used by the public InterfaceScout application and publication
reproducibility scripts. It implements the current 200-point Shrake-Rupley
side-chain accessibility calculation and 6/9 A multiscale chemistry mapping.
No GNM, RIN, APBS, DSSP, material presets, or legacy 5/8 A model are used.
"""
from __future__ import annotations

from io import StringIO
import urllib.request
from typing import Any, Dict, List, Optional

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
from pydantic import BaseModel, Field

APP_VERSION = "1.0-publication"
SC_RSA_THRESHOLD = 0.05
SASA_PROBE_A = 1.40
SASA_POINTS = 200
PATCH_RADII_A = (6.0, 9.0)
COARSE_PATCH_RADIUS_A = 8.0
PKA_SENSITIVITY_WINDOW = 1.0

SIDECHAIN_REF_ASA: Dict[str, float] = {
    "ALA": 69.23, "ARG": 200.35, "ASN": 106.25, "ASP": 102.06,
    "CYS": 96.69, "GLN": 140.58, "GLU": 134.61, "GLY": 32.28,
    "HIS": 147.00, "ILE": 137.91, "LEU": 140.76, "LYS": 162.50,
    "MET": 156.08, "PHE": 163.90, "PRO": 119.65, "SER": 78.16,
    "THR": 101.67, "TRP": 210.89, "TYR": 176.61, "VAL": 114.14,
}

PKA: Dict[str, float] = {
    "ASP": 3.9, "GLU": 4.1, "HIS": 6.5, "CYS": 8.3,
    "TYR": 10.1, "LYS": 10.5, "ARG": 12.5,
}
ION_SIGN: Dict[str, int] = {
    "ARG": +1, "LYS": +1, "HIS": +1,
    "ASP": -1, "GLU": -1, "CYS": -1, "TYR": -1,
}

# Numerical energy magnitudes are intentionally zero. Mechanism and required
# ionization state are retained as transparent chemistry metadata only.
CHEMISTRIES: Dict[str, Dict[str, Any]] = {
    "cationic": {
        "label": "Cationic surface",
        "surface_group": "cationic",
        "description": "Compatibility with positively charged surface groups.",
        "favorable": {
            "ASP": (0.0, "electrostatic + H-bond", "deprotonated"),
            "GLU": (0.0, "electrostatic + H-bond", "deprotonated"),
            "TYR": (0.0, "H-bond", None), "SER": (0.0, "H-bond", None),
            "THR": (0.0, "H-bond", None), "ASN": (0.0, "H-bond", None),
            "GLN": (0.0, "H-bond", None),
        },
        "repulsive": {},
    },
    "anionic": {
        "label": "Anionic surface",
        "surface_group": "anionic",
        "description": "Compatibility with negatively charged surface groups.",
        "favorable": {
            "LYS": (0.0, "electrostatic + H-bond", "protonated"),
            "ARG": (0.0, "electrostatic + H-bond", "protonated"),
            "HIS": (0.0, "electrostatic", "protonated"),
            "SER": (0.0, "H-bond donor", None), "THR": (0.0, "H-bond donor", None),
        },
        "repulsive": {
            "ASP": (0.0, "like-charge repulsion", "deprotonated"),
            "GLU": (0.0, "like-charge repulsion", "deprotonated"),
        },
    },
    "hbond_donor": {
        "label": "H-bond donor surface",
        "surface_group": "H-bond donor",
        "description": "Compatibility with surfaces capable of donating hydrogen bonds.",
        "favorable": {
            "TYR": (0.0, "H-bond acceptor", None),
            "ASP": (0.0, "H-bond acceptor", None), "GLU": (0.0, "H-bond acceptor", None),
            "SER": (0.0, "H-bond", None), "THR": (0.0, "H-bond", None),
            "ASN": (0.0, "H-bond", None), "GLN": (0.0, "H-bond", None),
        },
        "repulsive": {},
    },
    "hbond_acceptor": {
        "label": "H-bond acceptor surface",
        "surface_group": "H-bond acceptor",
        "description": "Compatibility with surfaces capable of accepting hydrogen bonds.",
        "favorable": {
            "LYS": (0.0, "H-bond donor", "protonated"),
            "ARG": (0.0, "H-bond donor", "protonated"),
            "SER": (0.0, "H-bond donor", None), "THR": (0.0, "H-bond donor", None),
            "TRP": (0.0, "indole-NH donor", None),
        },
        "repulsive": {},
    },
    "pi_carbon": {
        "label": "π / aromatic surface",
        "surface_group": "aromatic / graphitic",
        "description": "Compatibility with aromatic or graphitic-like surface chemistry.",
        "favorable": {
            "PHE": (0.0, "π-associated", None), "TRP": (0.0, "π-associated", None),
            "TYR": (0.0, "π-associated", None),
            "ARG": (0.0, "cation-π", "protonated"), "LYS": (0.0, "cation-π", "protonated"),
            "HIS": (0.0, "π-associated / cation-π", "auxiliary_only"),
        },
        "repulsive": {},
    },
    "hydrophobic": {
        "label": "Hydrophobic surface",
        "surface_group": "nonpolar",
        "description": "Compatibility with nonpolar surface chemistry.",
        "favorable": {r: (0.0, "hydrophobic contact", None) for r in ("ALA","VAL","LEU","ILE","MET","PHE","TRP","PRO")},
        "repulsive": {},
    },
    "oxide": {
        "label": "Metal-oxide surface",
        "surface_group": "oxide",
        "description": "Primary compatibility with exposed oxide metal sites through deprotonated carboxylates.",
        "favorable": {
            "ASP": (0.0, "carboxylate-metal coordination", "deprotonated"),
            "GLU": (0.0, "carboxylate-metal coordination", "deprotonated"),
        },
        "repulsive": {},
    },
    "hydroxyapatite": {
        "label": "Calcium/phosphate charged sites",
        "surface_group": "calcium/phosphate charged sites",
        "description": "Charged-site compatibility with calcium/phosphate-rich interfaces.",
        "favorable": {
            "ASP": (0.0, "carboxylate-Ca complementarity", "deprotonated"),
            "GLU": (0.0, "carboxylate-Ca complementarity", "deprotonated"),
            "ARG": (0.0, "charged-site complementarity", "protonated"),
            "LYS": (0.0, "charged-site complementarity", "protonated"),
            "HIS": (0.0, "charged-site complementarity", "protonated"),
        },
        "repulsive": {},
    },
    "metal_coord": {
        "label": "Transition-metal coordination",
        "surface_group": "transition-metal coordination",
        "description": "Compatibility with accessible transition-metal coordination sites.",
        "favorable": {
            "HIS": (0.0, "imidazole-metal coordination", None),
            "CYS": (0.0, "sulfur-metal coordination", None),
            "ASP": (0.0, "carboxylate-metal coordination", "deprotonated"),
            "GLU": (0.0, "carboxylate-metal coordination", "deprotonated"),
            "MET": (0.0, "thioether-metal coordination", None),
        },
        "repulsive": {},
    },
    "gold": {
        "label": "Soft-metal sulfur affinity",
        "surface_group": "soft-metal sulfur affinity",
        "description": "Soft-metal-like surface affinity dominated by exposed sulfur-containing side chains.",
        "favorable": {
            "CYS": (0.0, "soft-metal sulfur affinity", None),
            "MET": (0.0, "weak sulfur coordination", None),
        },
        "repulsive": {},
    },
    "phosphate": {
        "label": "Phosphate-rich surface",
        "surface_group": "phosphate-rich",
        "description": "Compatibility with phosphate-rich surface chemistry.",
        "favorable": {
            "ARG": (0.0, "phosphate-guanidinium", "protonated"),
            "LYS": (0.0, "phosphate-ammonium", "protonated"),
            "HIS": (0.0, "electrostatic", "protonated"),
            "SER": (0.0, "H-bond", None),
        },
        "repulsive": {},
    },
}


class EnvParams(BaseModel):
    pH: float = Field(7.4, ge=0.0, le=14.0)
    ionic: float = Field(150.0, ge=0.0)
    temp: float = Field(298.0, gt=0.0)


class AnalyzeRequest(BaseModel):
    pdb_id: Optional[str] = None
    pdb_text: Optional[str] = None
    chain: Optional[str] = None
    env: EnvParams = EnvParams()


def charged_fraction(res_name: str, pH: float) -> float:
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
    pka = PKA.get(res_name.upper())
    if pka is None:
        return 1.0
    if state == "protonated":
        return 1.0 / (1.0 + 10.0 ** (pH - pka))
    if state == "deprotonated":
        return 1.0 / (1.0 + 10.0 ** (pka - pH))
    return 1.0


def residue_charge_descriptor(res_name: str, pH: float) -> float:
    rn = res_name.upper()
    return float(ION_SIGN.get(rn, 0) * charged_fraction(rn, pH))


def normalize_to_100(values: np.ndarray) -> np.ndarray:
    maximum = float(np.max(values)) if values.size else 0.0
    return values / maximum * 100.0 if maximum > 0 else np.zeros_like(values, dtype=float)


def _pdb_text(req: AnalyzeRequest) -> str:
    if req.pdb_text:
        return req.pdb_text
    pid = (req.pdb_id or "").strip().upper()
    if not pid:
        raise ValueError("Provide pdb_id or pdb_text")
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=30) as response:
        return response.read().decode("utf-8", errors="replace")


def _sidechain_atoms(residue) -> List[Any]:
    if residue.get_resname().strip().upper() == "GLY":
        return [residue["CA"]] if "CA" in residue else []
    backbone = {"N", "CA", "C", "O", "OXT"}
    return [a for a in residue.get_atoms() if a.get_name().strip().upper() not in backbone and (getattr(a, "element", "") or "").upper() != "H"]


def _surface_residues(pdb_text: str, pH: float) -> List[Dict[str, Any]]:
    structure = PDBParser(QUIET=True).get_structure("interfacescout", StringIO(pdb_text))
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("No structural model found")
    ShrakeRupley(probe_radius=SASA_PROBE_A, n_points=SASA_POINTS).compute(structure, level="A")
    rows: List[Dict[str, Any]] = []
    for chain in model:
        for residue in chain:
            if not is_aa(residue, standard=True) or "CA" not in residue:
                continue
            rn = residue.get_resname().strip().upper()
            if rn not in SIDECHAIN_REF_ASA:
                continue
            sc_atoms = _sidechain_atoms(residue)
            sc_sasa = float(sum(float(getattr(a, "sasa", 0.0)) for a in sc_atoms))
            scrsa_raw = sc_sasa / SIDECHAIN_REF_ASA[rn]
            ca = residue["CA"]
            seq = int(residue.id[1]); icode = str(residue.id[2]).strip(); cid = str(chain.id)
            rows.append({
                "key": f"{cid}:{seq}:{icode}", "res_name": rn, "res_seq": seq, "icode": icode, "chain": cid,
                "x": float(ca.coord[0]), "y": float(ca.coord[1]), "z": float(ca.coord[2]),
                "sidechain_sasa": round(sc_sasa, 3), "scrsa_raw": round(scrsa_raw, 5),
                "scrsa": round(min(max(scrsa_raw, 0.0), 1.0), 5),
                "surface_exposed": bool(scrsa_raw >= SC_RSA_THRESHOLD),
                "charge_fraction": round(charged_fraction(rn, pH), 5),
                "charge_descriptor": round(residue_charge_descriptor(rn, pH), 5),
                "pka": PKA.get(rn),
                "ionization_sensitive": bool(rn in PKA and abs(pH - PKA[rn]) <= PKA_SENSITIVITY_WINDOW),
                "bfactor": round(float(ca.bfactor), 3), "phi": None, "ss": None, "n_neighbors_8A": 0,
            })
    return rows


def build_distances(surface: List[Dict[str, Any]]) -> np.ndarray:
    coords = np.asarray([[r["x"], r["y"], r["z"]] for r in surface], dtype=float)
    if len(coords) == 0:
        return np.empty((0, 0), dtype=float)
    delta = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def chemistry_map(surface: List[Dict[str, Any]], distances: np.ndarray, key: str, pH: float) -> Dict[str, Any]:
    meta = CHEMISTRIES[key]
    fav_defs = meta["favorable"]
    rep_defs = meta.get("repulsive", {})
    n = len(surface)
    local = np.zeros(n, dtype=float)
    rep_local = np.zeros(n, dtype=float)
    for i, row in enumerate(surface):
        rn = row["res_name"]
        if rn in fav_defs:
            local[i] = row["scrsa"] * state_availability(rn, fav_defs[rn][2], pH)
        if rn in rep_defs:
            rep_local[i] = row["scrsa"] * state_availability(rn, rep_defs[rn][2], pH)

    prop = normalize_to_100(local)
    rep_prop = normalize_to_100(rep_local)
    inner, outer = PATCH_RADII_A
    dens: Dict[float, np.ndarray] = {}
    dens_norm: Dict[float, np.ndarray] = {}
    for radius in PATCH_RADII_A:
        values = (distances <= radius).astype(float) @ local if n else np.zeros(0, dtype=float)
        dens[radius] = values
        maximum = float(np.max(values)) if values.size else 0.0
        dens_norm[radius] = values / maximum if maximum > 0 else np.zeros_like(values)
    persistence = 100.0 * np.minimum(dens_norm[inner], dens_norm[outer]) if n else np.zeros(0)
    geomean = 100.0 * np.sqrt(dens_norm[inner] * dens_norm[outer]) if n else np.zeros(0)

    patch_centers = []
    for i, row in enumerate(surface):
        if persistence[i] <= 0:
            continue
        patch_centers.append({
            "center_key": row["key"], "res_name": row["res_name"], "res_seq": row["res_seq"],
            "icode": row["icode"], "chain": row["chain"],
            f"density_{int(inner)}A_raw": round(float(dens[inner][i]), 6),
            f"density_{int(outer)}A_raw": round(float(dens[outer][i]), 6),
            f"density_{int(inner)}A_norm": round(float(dens_norm[inner][i] * 100.0), 3),
            f"density_{int(outer)}A_norm": round(float(dens_norm[outer][i] * 100.0), 3),
            "multiscale_persistence": round(float(persistence[i]), 3),
            "multiscale_geomean": round(float(geomean[i]), 3),
            "compatible_members_8A": [surface[j]["key"] for j in range(n) if distances[i, j] <= COARSE_PATCH_RADIUS_A and local[j] > 0],
        })
    patch_centers.sort(key=lambda x: (-x["multiscale_persistence"], -x["multiscale_geomean"], x["center_key"]))

    members = []
    for i, row in enumerate(surface):
        rn = row["res_name"]
        if rn not in fav_defs:
            continue
        ebase, mechanism, state = fav_defs[rn]
        members.append({
            "key": row["key"], "res_name": rn, "res_seq": row["res_seq"], "icode": row["icode"], "chain": row["chain"],
            "sidechain_sasa": row["sidechain_sasa"], "scrsa_raw": row["scrsa_raw"], "scrsa": row["scrsa"],
            "state_requirement": state or "none", "state_availability": round(state_availability(rn, state, pH), 5),
            "local_score": round(float(local[i]), 6), "propensity": round(float(prop[i]), 3),
            f"patch_density_{int(inner)}A": round(float(dens_norm[inner][i] * 100.0), 3),
            f"patch_density_{int(outer)}A": round(float(dens_norm[outer][i] * 100.0), 3),
            "multiscale_persistence": round(float(persistence[i]), 3), "multiscale_geomean": round(float(geomean[i]), 3),
            "ebase_metadata_kcal_mol": ebase, "mechanism": mechanism,
            "ionization_sensitive": row["ionization_sensitive"], "bfactor": row["bfactor"],
        })
    members.sort(key=lambda x: (-x["propensity"], -x["multiscale_persistence"], x["key"]))

    repulsive = []
    for i, row in enumerate(surface):
        rn = row["res_name"]
        if rn in rep_defs:
            repulsive.append({
                "key": row["key"], "res_name": rn, "res_seq": row["res_seq"], "chain": row["chain"],
                "repulsion_propensity": round(float(rep_prop[i]), 3),
            })

    return {
        "kind": "compatibility", "label": meta["label"], "surface_group": meta["surface_group"],
        "description": meta["description"], "n_favorable_residues": len(members),
        "n_repulsive_residues": len(repulsive), "top_patch": patch_centers[0] if patch_centers else None,
        "top_patches": patch_centers[:10], "patch_centers": patch_centers,
        "residues": members, "repulsive_residues": repulsive,
        "notes": {
            "primary_score": "L = scRSA x binary membership x required-state availability",
            "normalization": "within-map maximum = 100",
            "patch": f"multiscale persistence = 100 x min(normalized {int(inner)} A density, normalized {int(outer)} A density)",
        },
    }


def analyze(req: AnalyzeRequest) -> Dict[str, Any]:
    rows = _surface_residues(_pdb_text(req), req.env.pH)
    if not rows:
        raise ValueError("No standard amino-acid residues found")
    surface = [r for r in rows if r["surface_exposed"]]
    distances = build_distances(surface)
    if len(surface):
        for i, row in enumerate(surface):
            row["n_neighbors_8A"] = int(np.sum((distances[i] <= COARSE_PATCH_RADIUS_A) & (distances[i] > 0.0)))
    chemistries = {key: chemistry_map(surface, distances, key, req.env.pH) for key in CHEMISTRIES}
    return {
        "status": "ok", "version": APP_VERSION,
        "settings": {
            "pH": req.env.pH, "ionic_mM": req.env.ionic, "temperature_K": req.env.temp,
            "sasa_probe_A": SASA_PROBE_A, "sasa_points_per_atom": SASA_POINTS,
            "scrsa_threshold": SC_RSA_THRESHOLD, "patch_radii_A": list(PATCH_RADII_A),
            "coarse_patch_radius_A": COARSE_PATCH_RADIUS_A,
        },
        "chemistry_list": list(CHEMISTRIES), "chemistries": chemistries,
        "all_residues": rows, "surface_residues": surface,
        "reference_sidechain_asa": SIDECHAIN_REF_ASA,
    }
