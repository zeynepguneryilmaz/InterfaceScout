"""Explicit-graphene nonpolar orientation physics for InterfaceScout 2.0.

This module is a validation prototype intended to reproduce, more faithfully,
the physical ingredients used by Harrison et al. (Biointerphases 2017,
DOI 10.1116/1.4971381) while keeping InterfaceScout deterministic and
lightweight enough for screening.

For each rigid-body protein orientation and protein-surface separation:

    DeltaE = E_LJ(protein, explicit graphene)
             + sum_i sigma_i [SASA_i(complex) - SASA_i(isolated)]

where the surface-tension signs/magnitudes follow Harrison et al. and SASA is
recomputed for the actual protein+graphene complex with a 1.4 A probe using
FreeSASA's Shrake-Rupley implementation. The graphene sheet is an explicit
finite honeycomb lattice large enough that the protein-contact region is far
from sheet edges. Protein CHARMM36 LJ parameters are obtained from the same
OpenMM decoder used by ``nonpolar_energy``.

This is not yet wired into the production API. It exists so that the explicit
physics can be validated against literature before replacing the earlier
continuum approximation.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    import freesasa
except Exception:
    freesasa = None

try:
    from .nonpolar_energy import (
        CONTACT_A,
        DEFAULT_N_ORIENTATIONS,
        GRAPHENE_CC_A,
        HARRISON_SIGMA_KJ_MOL_A2,
        SASA_PROBE_A,
        SURFACE_EPSILON_KJ_MOL,
        SURFACE_RMIN2_A,
        SURFACE_SIGMA_A,
        _biopython_sasa_map,
        _openmm_atoms,
        fibonacci_sphere,
    )
except ImportError:
    from nonpolar_energy import (
        CONTACT_A,
        DEFAULT_N_ORIENTATIONS,
        GRAPHENE_CC_A,
        HARRISON_SIGMA_KJ_MOL_A2,
        SASA_PROBE_A,
        SURFACE_EPSILON_KJ_MOL,
        SURFACE_RMIN2_A,
        SURFACE_SIGMA_A,
        _biopython_sasa_map,
        _openmm_atoms,
        fibonacci_sphere,
    )

EXPLICIT_SEPARATIONS_A = (2.8, 3.2, 3.6, 4.0, 4.4)
LJ_CUTOFF_A = 12.0
SASA_POINTS = 100
GRAPHENE_MARGIN_A = LJ_CUTOFF_A + 2.0 * (SURFACE_RMIN2_A + SASA_PROBE_A) + 4.0


def _orthonormal_basis(normal: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = np.asarray(normal, dtype=float)
    n = n / np.linalg.norm(n)
    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(n, ref))) > 0.90:
        ref = np.array([0.0, 1.0, 0.0])
    u = np.cross(ref, n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    v /= np.linalg.norm(v)
    return u, v


def graphene_patch(half_width_A: float, half_height_A: Optional[float] = None) -> np.ndarray:
    """Generate a centered finite graphene honeycomb sheet in the xy plane."""
    hy = float(half_height_A if half_height_A is not None else half_width_A)
    hx = float(half_width_A)
    d = GRAPHENE_CC_A
    a1 = np.array([math.sqrt(3.0) * d, 0.0])
    a2 = np.array([math.sqrt(3.0) * d / 2.0, 1.5 * d])
    basis = (np.array([0.0, 0.0]), np.array([math.sqrt(3.0) * d / 2.0, 0.5 * d]))
    n1 = int(math.ceil((2.0 * hx + 4.0 * d) / np.linalg.norm(a1))) + 6
    n2 = int(math.ceil((2.0 * hy + 4.0 * d) / np.linalg.norm(a2))) + 6
    pts: List[Tuple[float, float, float]] = []
    for i in range(-n1, n1 + 1):
        for j in range(-n2, n2 + 1):
            origin = i * a1 + j * a2
            for b in basis:
                q = origin + b
                if -hx <= q[0] <= hx and -hy <= q[1] <= hy:
                    pts.append((float(q[0]), float(q[1]), 0.0))
    arr = np.asarray(pts, dtype=float)
    if arr.size == 0:
        raise RuntimeError("graphene patch generation returned no atoms")
    # Center the finite patch without changing bond geometry.
    arr[:, 0] -= 0.5 * (float(arr[:, 0].min()) + float(arr[:, 0].max()))
    arr[:, 1] -= 0.5 * (float(arr[:, 1].min()) + float(arr[:, 1].max()))
    return arr


def _freesasa_areas(coords_A: np.ndarray, radii_A: np.ndarray, n_points: int = SASA_POINTS) -> np.ndarray:
    if freesasa is None:
        raise RuntimeError("FreeSASA is not installed")
    coords = np.asarray(coords_A, dtype=float)
    radii = np.asarray(radii_A, dtype=float)
    if coords.shape != (len(radii), 3):
        raise ValueError("coordinate/radius dimensions are inconsistent")
    params = freesasa.Parameters({
        "algorithm": freesasa.ShrakeRupley,
        "probe-radius": float(SASA_PROBE_A),
        "n-points": int(n_points),
    })
    result = freesasa.calcCoord(coords.reshape(-1).tolist(), radii.tolist(), params)
    return np.asarray([float(result.atomArea(i)) for i in range(len(radii))], dtype=float)


def _explicit_lj(
    protein_xyz_A: np.ndarray,
    graphene_xyz_A: np.ndarray,
    protein_sigma_A: np.ndarray,
    protein_epsilon_kj: np.ndarray,
    cutoff_A: float = LJ_CUTOFF_A,
) -> float:
    """Direct protein-carbon 12-6 LJ sum with Lorentz-Berthelot mixing."""
    total = 0.0
    cutoff2 = float(cutoff_A) ** 2
    surf_sigma = float(SURFACE_SIGMA_A)
    surf_eps = float(SURFACE_EPSILON_KJ_MOL)
    for i in range(len(protein_xyz_A)):
        d = graphene_xyz_A - protein_xyz_A[i]
        r2 = np.einsum("ij,ij->i", d, d)
        keep = (r2 > 1e-12) & (r2 <= cutoff2)
        if not np.any(keep):
            continue
        r = np.sqrt(r2[keep])
        sig = 0.5 * (float(protein_sigma_A[i]) + surf_sigma)
        eps = math.sqrt(float(protein_epsilon_kj[i]) * surf_eps)
        sr6 = (sig / r) ** 6
        total += float(np.sum(4.0 * eps * (sr6 * sr6 - sr6)))
    return total


def _pose(xyz: np.ndarray, normal: np.ndarray, separation_A: float) -> np.ndarray:
    u, v = _orthonormal_basis(normal)
    x = xyz @ u
    y = xyz @ v
    z0 = xyz @ normal
    z = z0 - float(z0.min()) + float(separation_A)
    # Center the lateral footprint over the finite graphene patch.
    x = x - 0.5 * (float(x.min()) + float(x.max()))
    y = y - 0.5 * (float(y.min()) + float(y.max()))
    return np.column_stack([x, y, z])


def scan_explicit(
    pdb_path: Path,
    struct,
    pH: float = 7.4,
    n_orientations: int = 256,
    separations_A: Optional[Sequence[float]] = None,
    sasa_points: int = SASA_POINTS,
) -> Dict[str, Any]:
    """Explicit-graphene rigid-body orientation scan using actual complex SASA."""
    sep_grid = tuple(float(x) for x in (separations_A or EXPLICIT_SEPARATIONS_A))
    try:
        atoms, source, ljdiag = _openmm_atoms(Path(pdb_path), pH, _biopython_sasa_map(struct))
        xyz = np.vstack([a["coord"] for a in atoms]).astype(float)
        sigma = np.asarray([a["sigma_A"] for a in atoms], dtype=float)
        epsilon = np.asarray([a["epsilon_kj_mol"] for a in atoms], dtype=float)
        radii = np.asarray([a["rmin2_A"] for a in atoms], dtype=float)
        signs = np.asarray([a["solvation_sign"] for a in atoms], dtype=float)
        protein_iso = _freesasa_areas(xyz, radii, sasa_points)
    except Exception as exc:
        return {"status": "unavailable", "method": "explicit graphene + complex SASA", "reason": str(exc)}

    # A single patch is sized from the maximum protein lateral span. Rotation
    # cannot make the footprint wider than the protein diameter.
    centered = xyz - xyz.mean(axis=0)
    diameter = 2.0 * float(np.max(np.linalg.norm(centered, axis=1)))
    half = 0.5 * diameter + GRAPHENE_MARGIN_A
    graphene = graphene_patch(half)
    graphene_radii = np.full(len(graphene), float(SURFACE_RMIN2_A), dtype=float)
    graphene_signs = np.ones(len(graphene), dtype=float)
    try:
        graphene_iso = _freesasa_areas(graphene, graphene_radii, sasa_points)
    except Exception as exc:
        return {"status": "unavailable", "method": "explicit graphene + complex SASA", "reason": str(exc)}

    normals = fibonacci_sphere(int(n_orientations))
    bestE = np.full(len(normals), np.inf, dtype=float)
    bestZ = np.zeros(len(normals), dtype=float)
    bestV = np.zeros(len(normals), dtype=float)
    bestS = np.zeros(len(normals), dtype=float)
    bestPB = np.zeros(len(normals), dtype=float)
    bestGB = np.zeros(len(normals), dtype=float)

    all_radii = np.concatenate([radii, graphene_radii])
    all_signs = np.concatenate([signs, graphene_signs])
    iso = np.concatenate([protein_iso, graphene_iso])

    for k, n in enumerate(normals):
        for sep in sep_grid:
            pxyz = _pose(xyz, n, sep)
            vdw = _explicit_lj(pxyz, graphene, sigma, epsilon)
            complex_xyz = np.vstack([pxyz, graphene])
            areas = _freesasa_areas(complex_xyz, all_radii, sasa_points)
            delta = areas - iso
            solv = float(HARRISON_SIGMA_KJ_MOL_A2 * np.dot(all_signs, delta))
            total = vdw + solv
            if total < bestE[k]:
                bestE[k] = total
                bestZ[k] = sep
                bestV[k] = vdw
                bestS[k] = solv
                bestPB[k] = float(np.sum(protein_iso - areas[:len(atoms)]))
                bestGB[k] = float(np.sum(graphene_iso - areas[len(atoms):]))

    order = np.argsort(bestE, kind="mergesort")
    top: List[Dict[str, Any]] = []
    for rank, idx in enumerate(order[:20], 1):
        pxyz = _pose(xyz, normals[idx], bestZ[idx])
        residues: Dict[str, Dict[str, Any]] = {}
        for j, a in enumerate(atoms):
            # Harrison contact convention: any protein atom within 6 A of the
            # graphene plane. Because the explicit sheet lies at z=0, this is
            # equivalent to the minimum protein atom height criterion here.
            if pxyz[j, 2] > CONTACT_A:
                continue
            key = f"{a['chain']}:{a['res_seq']}:{a['icode']}"
            row = residues.setdefault(key, {
                "key": key, "chain": a["chain"], "res_seq": a["res_seq"],
                "icode": a["icode"], "res_name": a["res_name"],
                "min_atom_height_A": float(pxyz[j, 2]), "contact_atom_count": 0,
            })
            row["min_atom_height_A"] = min(row["min_atom_height_A"], float(pxyz[j, 2]))
            row["contact_atom_count"] += 1
        rlist = sorted(residues.values(), key=lambda r: (r["min_atom_height_A"], r["chain"], r["res_seq"]))
        for r in rlist:
            r["min_atom_height_A"] = round(float(r["min_atom_height_A"]), 4)
        top.append({
            "rank": rank,
            "orientation_index": int(idx),
            "normal": [round(float(x), 6) for x in normals[idx]],
            "minimum_separation_A": round(float(bestZ[idx]), 4),
            "total_energy_change_kj_mol": round(float(bestE[idx]), 4),
            "vdw_energy_kj_mol": round(float(bestV[idx]), 4),
            "solvation_energy_change_kj_mol": round(float(bestS[idx]), 4),
            "protein_buried_sasa_A2": round(float(bestPB[idx]), 4),
            "graphene_buried_sasa_A2": round(float(bestGB[idx]), 4),
            "contact_residues": rlist,
        })

    return {
        "status": "ok",
        "method": "explicit finite graphene + CHARMM36 LJ + actual complex Shrake-Rupley SASA",
        "energy_definition": "DeltaE = direct protein-graphene LJ + sum sigma_i*(SASA_complex-SASA_isolated)",
        "parameter_source": source,
        "lj_parameter_diagnostics": ljdiag,
        "n_orientations": int(n_orientations),
        "separation_grid_A": list(sep_grid),
        "contact_distance_A": CONTACT_A,
        "sasa": {"engine": "FreeSASA Shrake-Rupley", "probe_A": SASA_PROBE_A, "n_points": int(sasa_points)},
        "graphene": {
            "n_atoms": int(len(graphene)),
            "half_width_A": float(half),
            "C_C_A": GRAPHENE_CC_A,
            "finite_patch": True,
            "edge_margin_A": GRAPHENE_MARGIN_A,
        },
        "surface_tension_kj_mol_A2": HARRISON_SIGMA_KJ_MOL_A2,
        "no_fitted_weights": True,
        "absolute_adsorption_free_energy": False,
        "best_energy_change_kj_mol": round(float(bestE[order[0]]), 4),
        "median_best_energy_change_kj_mol": round(float(np.median(bestE)), 4),
        "top_orientations": top,
    }
