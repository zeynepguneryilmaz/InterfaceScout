"""Established-physics refinements for InterfaceScout v5.3 development.

This module does not alter the frozen v5.1 chemistry/state equations. It adds a
separate nonpolar-interface model built from published hydrophobicity scales and
an exposure-weighted 3-D hydrophobic vector.

Primary nonpolar descriptor
----------------------------
Eisenberg normalized consensus hydrophobicity (Eisenberg et al., J Mol Biol
179, 125-142, 1984) is treated as a signed continuous residue property instead
of binary hydrophobic membership. For exposed residue j,

    h_j = H_j * scRSA_j

and local fields are

    Q_i(R) = sum_{j:d_ij<=R} h_j , R in {5, 8 A}.

Each Q(R) is min-max normalized over the analysed surface and the two-scale
persistence descriptor is the minimum of the two normalized fields. Negative
hydrophilicity therefore penalizes a region; residues with negative individual
hydrophobicity can still be passenger contacts inside a favorable patch.

Hydrophobic vector / approach hemisphere
----------------------------------------
The global exposure-weighted hydrophobic vector follows the solvent-accessibility
form described for tertiary protein structures by Silverman (Proteins 53,
880-888, 2003; see also the published implementation description):

    H = (1/n) * sum_i h_i * p_i * u_i

where h_i is residue hydrophobicity, p_i is solvent exposure (scRSA here), and
u_i is the unit vector from the centroid of residue centroids to residue i.
Side-chain centroids are used when available, with C-alpha fallback already
provided by the structural core. The positive half-space along H defines the
hydrophobic-facing hemisphere. This is reported as an orientation descriptor,
not an adsorption free energy or unique adsorbed pose.

Independent sensitivity descriptor
----------------------------------
The experimentally determined Wimley-White interfacial hydrophobicity scale
(Nat Struct Biol 3, 842-848, 1996) is evaluated with the same patch geometry.
For convenient interpretation the published transfer-free-energy sign is
reversed so that larger values mean greater interfacial preference. It is
reported independently and is not averaged with the Eisenberg descriptor.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
import math
import numpy as np

EISENBERG = {
    "ALA": 0.62, "ARG": -2.53, "ASN": -0.78, "ASP": -0.90,
    "CYS": 0.29, "GLN": -0.85, "GLU": -0.74, "GLY": 0.48,
    "HIS": -0.40, "ILE": 1.38, "LEU": 1.06, "LYS": -1.50,
    "MET": 0.64, "PHE": 1.19, "PRO": 0.12, "SER": -0.18,
    "THR": -0.05, "TRP": 0.81, "TYR": 0.26, "VAL": 1.08,
}

# Negative of the commonly tabulated Wimley-White water->interface transfer
# free energies, so positive means stronger interfacial preference.
WIMLEY_WHITE_INTERFACE = {
    "ALA": -0.17, "ARG": -0.81, "ASN": -0.42, "ASP": -1.23,
    "CYS": 0.24, "GLN": -0.58, "GLU": -2.02, "GLY": -0.01,
    "HIS": -0.17, "ILE": 0.31, "LEU": 0.56, "LYS": -0.99,
    "MET": 0.23, "PHE": 1.13, "PRO": -0.45, "SER": -0.13,
    "THR": -0.14, "TRP": 1.85, "TYR": 0.94, "VAL": -0.07,
}

RADII_A = (5.0, 8.0)


def _ca_coords(rows: List[Dict[str, Any]]) -> np.ndarray:
    return np.asarray([[float(r["x"]), float(r["y"]), float(r["z"])] for r in rows], dtype=float)


def _sc_coords(rows: List[Dict[str, Any]]) -> np.ndarray:
    return np.asarray([[float(r.get("sc_x", r["x"])), float(r.get("sc_y", r["y"])), float(r.get("sc_z", r["z"]))] for r in rows], dtype=float)


def _distances(rows: List[Dict[str, Any]]) -> np.ndarray:
    c = _ca_coords(rows)
    if len(c) == 0:
        return np.zeros((0, 0), dtype=float)
    d = c[:, None, :] - c[None, :, :]
    return np.sqrt(np.sum(d * d, axis=2))


def _minmax(v: np.ndarray) -> np.ndarray:
    if v.size == 0:
        return v
    lo = float(np.min(v)); hi = float(np.max(v))
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo + 1e-15:
        return np.zeros_like(v)
    return 100.0 * (v - lo) / (hi - lo)


def _scale_patch(rows: List[Dict[str, Any]], scale: Dict[str, float], label: str) -> Dict[str, Any]:
    if not rows:
        return {"label": label, "residues": [], "top_patches": [], "radii_A": list(RADII_A)}
    d = _distances(rows)
    local = np.asarray([float(scale.get(r["res_name"], 0.0)) * float(r.get("scrsa", 0.0)) for r in rows], dtype=float)
    fields = {}; norms = {}
    for R in RADII_A:
        q = np.asarray([float(np.sum(local[d[i] <= R])) for i in range(len(rows))], dtype=float)
        fields[R] = q; norms[R] = _minmax(q)
    persistence = np.minimum(norms[RADII_A[0]], norms[RADII_A[1]])
    out_rows = []
    for i, r in enumerate(rows):
        z = dict(r)
        z.update({
            "scale_value": round(float(scale.get(r["res_name"], 0.0)), 5),
            "exposure_weighted_value": round(float(local[i]), 5),
            "field_5A": round(float(fields[5.0][i]), 5),
            "field_8A": round(float(fields[8.0][i]), 5),
            "field_5A_norm": round(float(norms[5.0][i]), 5),
            "field_8A_norm": round(float(norms[8.0][i]), 5),
            "persistence": round(float(persistence[i]), 5),
            "driver_sign": "favorable" if local[i] > 0 else ("unfavorable" if local[i] < 0 else "neutral"),
        })
        out_rows.append(z)
    ranked = sorted(out_rows, key=lambda x: (-x["persistence"], -x["field_8A"], x["chain"], x["res_seq"], x.get("icode", "")))
    return {
        "label": label, "radii_A": list(RADII_A), "signed_continuous": True,
        "used_in_frozen_v51_score": False, "residues": out_rows,
        "top_patches": ranked[:10], "top_patch": ranked[0] if ranked else None,
    }


def hydrophobic_vector(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Exposure-weighted tertiary hydrophobic vector using Silverman Eq.-12 form."""
    if not rows:
        return {"status": "unavailable"}
    c = _sc_coords(rows)
    center = np.mean(c, axis=0)
    rel = c - center[None, :]
    unit_radial = np.zeros_like(rel)
    norms = np.linalg.norm(rel, axis=1)
    good = norms > 1e-12
    unit_radial[good] = rel[good] / norms[good, None]
    weights = np.asarray([float(EISENBERG.get(r["res_name"], 0.0)) * float(r.get("scrsa", 0.0)) for r in rows], dtype=float)
    mu = np.sum(weights[:, None] * unit_radial, axis=0) / max(len(rows), 1)
    mag = float(np.linalg.norm(mu))
    unit = mu / mag if mag > 1e-12 else np.zeros(3, dtype=float)

    face_rows = []
    for i, r in enumerate(rows):
        align = float(np.dot(unit_radial[i], unit)) if mag > 1e-12 and good[i] else 0.0
        face_rows.append({
            "key": r["key"], "chain": r["chain"], "res_seq": r["res_seq"], "icode": r.get("icode", ""),
            "res_name": r["res_name"], "alignment_cosine": round(align, 5),
            "hydrophobic_hemisphere": bool(align >= 0.0),
        })
    return {
        "status": "ok",
        "definition": "H=(1/n) sum h_i*scRSA_i*u_i; exposure-weighted tertiary hydrophobic vector",
        "coordinate_reference": "side-chain centroid with C-alpha fallback; centroid of residue centroids as origin",
        "vector": [round(float(x), 6) for x in mu],
        "unit_vector": [round(float(x), 6) for x in unit],
        "magnitude": round(mag, 6),
        "hemisphere_rule": "alignment_cosine >= 0; no fitted angular cutoff",
        "used_as_adsorption_free_energy": False,
        "residues": face_rows,
    }


def enrich_nonpolar_physics(surface_rows: List[Dict[str, Any]], all_rows: Optional[List[Dict[str, Any]]] = None) -> Dict[str, Any]:
    eis = _scale_patch(surface_rows, EISENBERG, "Eisenberg normalized-consensus hydrophobic field")
    ww = _scale_patch(surface_rows, WIMLEY_WHITE_INTERFACE, "Wimley-White interfacial preference field")
    dip = hydrophobic_vector(all_rows if all_rows is not None else surface_rows)
    align = {r["key"]: r for r in dip.get("residues", [])}
    gated = []
    for r in eis.get("residues", []):
        a = align.get(r["key"], {})
        z = dict(r); z["alignment_cosine"] = a.get("alignment_cosine"); z["hydrophobic_hemisphere"] = a.get("hydrophobic_hemisphere")
        if a.get("hydrophobic_hemisphere"):
            gated.append(z)
    gated.sort(key=lambda x: (-x["persistence"], -x["field_8A"], x["chain"], x["res_seq"], x.get("icode", "")))

    return {
        "model": "established nonpolar-interface descriptors; no fitted cross-term weights",
        "primary_local_field": eis,
        "independent_interface_scale_sensitivity": ww,
        "hydrophobic_dipole": dip,
        "orientation_gated_top_patches": gated[:10],
        "interpretation": {
            "driver_vs_passenger": "individual hydrophilic/neutral residues may be passenger contacts inside a favorable local hydrophobic field",
            "hydration": "published transfer/hydrophobicity scales provide implicit solvation information; explicit interfacial water is still not simulated",
            "orientation": "the tertiary hydrophobic vector supplies a preferred hemisphere, not a unique adsorbed pose",
        },
    }
