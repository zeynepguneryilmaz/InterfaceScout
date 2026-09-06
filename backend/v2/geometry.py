"""Geometry utilities for coarse InterfaceScout V2 interface patches.

The goal is not atomistic docking.  We only ask whether residues belong to the
same exposed protein face and can plausibly participate in one coarse contact
region.  No adsorption benchmark labels are used here.
"""

from __future__ import annotations

from typing import Dict, Iterable

import numpy as np


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n <= 1e-12:
        return np.zeros(3, dtype=float)
    return np.asarray(v, dtype=float) / n


def build_surface_geometry(v1_result: dict, gnm: dict) -> dict:
    """Return C-alpha coordinates and coarse outward directions for surface residues.

    The outward direction is a deliberately coarse face descriptor: the vector
    from the protein C-alpha centroid to the residue C-alpha position.  For the
    folded/globular proteins in the primary validation scope this is stable,
    transparent, and parameter-free.  It is not claimed to be a molecular
    surface normal for highly concave or intrinsically disordered structures.
    """
    nodes = gnm["nodes"]
    node_by_key = {str(n["key"]): n for n in nodes}
    protein_centroid = np.mean(np.vstack([n["coord"] for n in nodes]), axis=0)

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
    """Co-facing rule with no fitted angular coefficient.

    Positive dot product means the two coarse outward directions lie in the same
    hemisphere (<90 degrees apart).  This is used as a permissive geometric gate,
    not as a score.
    """
    return float(np.dot(geometry["normals"][key_a], geometry["normals"][key_b])) > 0.0


def patch_orientation_coherence(keys: Iterable[str], geometry: dict) -> float:
    """Resultant length of unit outward directions; descriptive range 0..1."""
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
