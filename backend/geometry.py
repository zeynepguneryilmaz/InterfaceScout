"""Geometry utilities for coarse InterfaceScout interface patches.

The goal is not atomistic docking. We only ask whether residues belong to the
same exposed protein face and can plausibly participate in one coarse contact
region. No adsorption benchmark labels are used here.
"""

from __future__ import annotations

from io import StringIO
from typing import Dict, Iterable, List

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n <= 1e-12:
        return np.zeros(3, dtype=float)
    return np.asarray(v, dtype=float) / n


def _residue_key(chain_id: str, residue) -> str:
    seq = int(residue.id[1])
    icode = str(residue.id[2]).strip()
    return f"{chain_id}:{seq}:{icode}"


def extract_ca_nodes(pdb_text: str) -> List[dict]:
    """Extract standard-amino-acid C-alpha coordinates from the first PDB model."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("interfacescout_geometry", StringIO(pdb_text))
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("No model found in prepared structure")

    nodes: List[dict] = []
    for chain in model:
        for residue in chain:
            if not is_aa(residue, standard=True) or "CA" not in residue:
                continue
            nodes.append(
                {
                    "key": _residue_key(str(chain.id), residue),
                    "chain": str(chain.id),
                    "res_seq": int(residue.id[1]),
                    "icode": str(residue.id[2]).strip(),
                    "res_name": str(residue.resname).strip(),
                    "coord": np.asarray(residue["CA"].coord, dtype=float),
                }
            )
    if len(nodes) < 1:
        raise ValueError("No standard-amino-acid C-alpha coordinates found")
    return nodes


def build_surface_geometry(v1_result: dict, ca_nodes: List[dict]) -> dict:
    """Return C-alpha coordinates and coarse outward directions for surface residues.

    The outward direction is the vector from the protein C-alpha centroid to the
    residue C-alpha position. It is a coarse face descriptor rather than a true
    molecular-surface normal.
    """
    node_by_key = {str(n["key"]): n for n in ca_nodes}
    protein_centroid = np.mean(np.vstack([n["coord"] for n in ca_nodes]), axis=0)

    surface_rows = {
        str(r["key"]): r
        for r in v1_result.get("surface_residues", [])
        if r.get("key") and str(r["key"]) in node_by_key
    }

    coords: Dict[str, np.ndarray] = {}
    normals: Dict[str, np.ndarray] = {}
    scrsa: Dict[str, float] = {}
    meta: Dict[str, dict] = {}

    for key, row in surface_rows.items():
        node = node_by_key[key]
        coord = np.asarray(node["coord"], dtype=float)
        coords[key] = coord
        normals[key] = _unit(coord - protein_centroid)
        scrsa[key] = float(row.get("scrsa", row.get("scrsa_raw", 0.0)) or 0.0)
        meta[key] = {
            "key": key,
            "chain": node["chain"],
            "res_seq": int(node["res_seq"]),
            "icode": node["icode"],
            "res_name": node["res_name"],
        }

    return {
        "protein_centroid": protein_centroid,
        "coords": coords,
        "normals": normals,
        "scrsa": scrsa,
        "meta": meta,
    }


def ca_distance(key_a: str, key_b: str, geometry: dict) -> float:
    return float(np.linalg.norm(geometry["coords"][key_a] - geometry["coords"][key_b]))


def same_face(key_a: str, key_b: str, geometry: dict) -> bool:
    """Return whether two coarse outward directions lie in the same hemisphere."""
    return float(np.dot(geometry["normals"][key_a], geometry["normals"][key_b])) > 0.0


def patch_orientation_coherence(keys: Iterable[str], geometry: dict) -> float:
    """Resultant length of unit outward directions; range 0..1."""
    valid = [k for k in keys if k in geometry["normals"]]
    if not valid:
        return 0.0
    vec = np.mean(np.vstack([geometry["normals"][k] for k in valid]), axis=0)
    return float(np.linalg.norm(vec))


def patch_diameter_A(keys: Iterable[str], geometry: dict) -> float:
    valid = [k for k in keys if k in geometry["coords"]]
    if len(valid) < 2:
        return 0.0
    pts = np.vstack([geometry["coords"][k] for k in valid])
    delta = pts[:, None, :] - pts[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", delta, delta)
    return float(np.sqrt(np.max(d2)))
