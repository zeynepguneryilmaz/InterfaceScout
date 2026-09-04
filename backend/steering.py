"""Screened electrostatic steering for charged planar interfaces.

This module is deliberately separate from InterfaceScout's canonical local
compatibility score. It does not alter residue propensities or 5/8 Å patch
persistence. Instead it estimates which protein face is electrostatically
favoured to approach a homogeneous charged plane before short-range anchoring.

Model
-----
For a monovalent electrolyte in the linearized Poisson-Boltzmann (Debye-Hückel)
limit, the potential above a planar interface decays as

    psi(z) = psi0 * exp(-kappa z)

For residue-level partial charges q_i (in elementary-charge units), the
orientation-dependent interaction energy is therefore proportional to

    U_tilde(n; s) = s * sum_i q_i exp[-z_i(n)/lambda_D]

where s = +1 for a positively charged plane, s = -1 for a negatively charged
plane, lambda_D = kappa^-1, and z_i(n) is the C-alpha depth of residue i from a
plane tangent to the exposed C-alpha envelope along orientation n.

U_tilde is dimensionless and equals U / (e |psi0|). Because |psi0| is not
specified, only orientation ranking is interpreted; no absolute adsorption
energy is reported.

Important scope boundaries
--------------------------
* Homogeneous planar charged interface.
* Linearized screening profile; no charge regulation or nonlinear PB.
* Bulk-pH residue charge descriptors, not constant-pH titration near surface.
* C-alpha coordinates define rigid residue topology for steering geometry.
* Steering predicts an electrostatic surface-facing preference, not a unique
  adsorbed state and not final short-range anchoring.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# CODATA exact/standard constants in SI where applicable.
_E_CHARGE_C = 1.602176634e-19
_KB_J_K = 1.380649e-23
_NA_MOL = 6.02214076e23
_EPS0_F_M = 8.8541878128e-12
_EPS_R_WATER = 78.54


def debye_length_A(ionic_mM: float, temperature_K: float, eps_r: float = _EPS_R_WATER) -> Optional[float]:
    """Return Debye length in Å for a symmetric monovalent electrolyte.

    Uses
        kappa^2 = 2 e^2 N_A (1000 I_M) / (eps_r eps0 k_B T)
    with I_M in mol/L. Returns None at zero ionic strength because the
    linearized screened-plane form used by this module is then undefined.
    """
    ionic_mM = float(ionic_mM)
    temperature_K = float(temperature_K)
    if ionic_mM <= 0.0 or temperature_K <= 0.0 or eps_r <= 0.0:
        return None
    ionic_M = ionic_mM / 1000.0
    kappa_m = math.sqrt(
        2.0 * (_E_CHARGE_C ** 2) * _NA_MOL * 1000.0 * ionic_M
        / (eps_r * _EPS0_F_M * _KB_J_K * temperature_K)
    )
    if not math.isfinite(kappa_m) or kappa_m <= 0.0:
        return None
    return float((1.0 / kappa_m) * 1.0e10)


def fibonacci_sphere(n: int = 2048) -> np.ndarray:
    """Deterministic approximately uniform unit-vector sampling of the sphere."""
    n = max(int(n), 32)
    i = np.arange(n, dtype=float)
    z = 1.0 - 2.0 * (i + 0.5) / n
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    phi = golden_angle * i
    radial = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    x = radial * np.cos(phi)
    y = radial * np.sin(phi)
    out = np.stack([x, y, z], axis=1)
    out /= np.linalg.norm(out, axis=1, keepdims=True)
    return out


def _coords(rows: Sequence[Dict[str, Any]]) -> np.ndarray:
    return np.asarray([[float(r["x"]), float(r["y"]), float(r["z"])] for r in rows], dtype=float)


def _charges(rows: Sequence[Dict[str, Any]]) -> np.ndarray:
    return np.asarray([float(r.get("charge_descriptor", 0.0) or 0.0) for r in rows], dtype=float)


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n <= 0.0:
        return np.zeros(3, dtype=float)
    return v / n


def _footprint(surface: Sequence[Dict[str, Any]], surface_coords: np.ndarray, normal: np.ndarray,
               depths_A: Iterable[float]) -> Dict[str, Any]:
    proj = surface_coords @ normal
    plane = float(np.min(proj))
    depth = np.maximum(proj - plane, 0.0)
    out: Dict[str, Any] = {}
    for d in depths_A:
        key = f"within_{int(round(float(d)))}A"
        rows = []
        for i, r in enumerate(surface):
            if depth[i] <= float(d) + 1e-12:
                rows.append({
                    "key": r["key"],
                    "res_name": r["res_name"],
                    "res_seq": r["res_seq"],
                    "icode": r["icode"],
                    "chain": r["chain"],
                    "depth_A": round(float(depth[i]), 4),
                    "scrsa": r.get("scrsa"),
                    "charge_descriptor": r.get("charge_descriptor"),
                })
        rows.sort(key=lambda x: (x["depth_A"], x["chain"], x["res_seq"], x["icode"]))
        out[key] = rows
    return out


def _sampled_orientation_rows(order: np.ndarray, normals: np.ndarray, energy: np.ndarray,
                              limit: int = 10) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for rank, idx in enumerate(order[: max(1, int(limit))], start=1):
        n = normals[int(idx)]
        rows.append({
            "rank": rank,
            "normal_plane_to_protein": [round(float(x), 6) for x in n],
            "protein_facing_direction_to_plane": [round(float(x), 6) for x in (-n)],
            "reduced_energy": round(float(energy[int(idx)]), 8),
        })
    return rows


def _one_surface_sign(
    sign: int,
    label: str,
    linked_local_map: str,
    all_residues: Sequence[Dict[str, Any]],
    surface_residues: Sequence[Dict[str, Any]],
    charge_coords: np.ndarray,
    surface_coords: np.ndarray,
    q: np.ndarray,
    normals: np.ndarray,
    lambda_A: float,
    footprint_depths_A: Tuple[float, ...],
) -> Dict[str, Any]:
    # Plane is tangent to the exposed C-alpha envelope for each sampled orientation.
    surface_proj = surface_coords @ normals.T  # n_surface x n_orient
    plane_proj = np.min(surface_proj, axis=0)
    z = charge_coords @ normals.T - plane_proj[None, :]
    z = np.maximum(z, 0.0)
    screened = np.exp(-z / lambda_A)
    base = np.sum(q[:, None] * screened, axis=0)
    energy = float(sign) * base
    order = np.argsort(energy)
    best_idx = int(order[0])
    best_n = normals[best_idx]

    return {
        "surface_charge": label,
        "surface_sign": int(sign),
        "linked_local_compatibility_map": linked_local_map,
        "best_normal_plane_to_protein": [round(float(x), 6) for x in best_n],
        "best_protein_facing_direction_to_plane": [round(float(x), 6) for x in (-best_n)],
        "best_reduced_energy": round(float(energy[best_idx]), 8),
        "worst_reduced_energy": round(float(np.max(energy)), 8),
        "reduced_energy_span": round(float(np.max(energy) - np.min(energy)), 8),
        "sampled_orientations": int(len(normals)),
        "top_sampled_orientations": _sampled_orientation_rows(order, normals, energy, limit=10),
        "surface_facing_footprint": _footprint(
            surface_residues, surface_coords, best_n, footprint_depths_A
        ),
    }


def compute_electrostatic_steering(
    all_residues: Sequence[Dict[str, Any]],
    surface_residues: Sequence[Dict[str, Any]],
    ionic_mM: float,
    temperature_K: float,
    n_orientations: int = 2048,
    footprint_depths_A: Tuple[float, ...] = (5.0, 8.0),
) -> Dict[str, Any]:
    """Compute charged-plane electrostatic steering as a separate auxiliary layer."""
    lambda_A = debye_length_A(ionic_mM, temperature_K)
    if lambda_A is None:
        return {
            "available": False,
            "reason": "ionic strength must be > 0 for the screened-plane Debye-Huckel model",
            "model": "linearized Poisson-Boltzmann charged-plane orientation ranking",
        }
    if not all_residues or not surface_residues:
        return {
            "available": False,
            "reason": "protein contains no analyzable residue/surface coordinates",
            "model": "linearized Poisson-Boltzmann charged-plane orientation ranking",
        }

    charge_coords_all = _coords(all_residues)
    q_all = _charges(all_residues)
    charged_mask = np.abs(q_all) > 1.0e-8
    if not np.any(charged_mask):
        return {
            "available": False,
            "reason": "no ionizable residue charge remained at the selected pH",
            "model": "linearized Poisson-Boltzmann charged-plane orientation ranking",
            "debye_length_A": round(float(lambda_A), 5),
            "net_charge_descriptor_e": 0.0,
        }

    charge_coords = charge_coords_all[charged_mask]
    q = q_all[charged_mask]
    surface_coords = _coords(surface_residues)
    normals = fibonacci_sphere(n_orientations)

    ca_centroid = np.mean(charge_coords_all, axis=0)
    charge_moment = np.sum(q_all[:, None] * (charge_coords_all - ca_centroid[None, :]), axis=0)

    positive = _one_surface_sign(
        +1, "positive", "cationic", all_residues, surface_residues,
        charge_coords, surface_coords, q, normals, lambda_A, footprint_depths_A,
    )
    negative = _one_surface_sign(
        -1, "negative", "anionic", all_residues, surface_residues,
        charge_coords, surface_coords, q, normals, lambda_A, footprint_depths_A,
    )

    return {
        "available": True,
        "kind": "electrostatic_steering",
        "model": "linearized Poisson-Boltzmann charged-plane orientation ranking",
        "canonical_compatibility_scores_modified": False,
        "debye_length_A": round(float(lambda_A), 5),
        "ionic_mM": float(ionic_mM),
        "temperature_K": float(temperature_K),
        "relative_dielectric": _EPS_R_WATER,
        "net_charge_descriptor_e": round(float(np.sum(q_all)), 6),
        "n_charged_residues": int(np.sum(charged_mask)),
        "charge_first_moment_eA": [round(float(x), 6) for x in charge_moment],
        "charge_first_moment_magnitude_eA": round(float(np.linalg.norm(charge_moment)), 6),
        "geometry": {
            "charge_positions": "C-alpha residue coordinates",
            "plane_tangency": "minimum projection of exposed C-alpha envelope for each sampled normal",
            "footprint": "surface residues within the requested C-alpha depth from the tangent plane",
            "orientation_sampling": "deterministic Fibonacci sphere",
            "n_orientations": int(len(normals)),
        },
        "energy_definition": {
            "reduced_energy": "U_tilde = U / (e |psi0|) = s * sum_i q_i exp(-z_i/lambda_D)",
            "interpretation": "lower reduced energy = more favorable screened electrostatic approach",
            "absolute_energy": "not reported because surface-potential magnitude is unspecified",
        },
        "positive_surface": positive,
        "negative_surface": negative,
        "scope": {
            "predicts": "whole-protein electrostatic surface-facing preference for a homogeneous charged plane",
            "does_not_predict": [
                "absolute adsorption free energy",
                "charge regulation",
                "nonlinear Poisson-Boltzmann effects",
                "heterogeneous surface charge domains",
                "final short-range anchoring",
                "adsorption-induced conformational change",
            ],
        },
    }
