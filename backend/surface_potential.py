"""Numerical sampling of APBS potential on an approximate solvent-accessible shell.

Each PQR heavy-atom vdW sphere is expanded by the probe and a small outward
offset. Fibonacci sphere points occluded by any *other* expanded sphere are
discarded. APBS potentials are trilinearly interpolated at retained points,
then averaged for each residue using the sampling area weights.

This is a discrete probe-center SAS shell, NOT an analytic solvent-excluded
surface (SES), a binding affinity, or a material-specific electrostatic model.
"""
from __future__ import annotations

from collections import defaultdict
import math
import re
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

PROBE_A = 1.4
OUTWARD_OFFSET_A = 0.25
SPHERE_DIRECTIONS = 42


def fibonacci_directions(n: int = SPHERE_DIRECTIONS) -> np.ndarray:
    if n < 6:
        raise ValueError("Too few directions for a reproducible sphere sampling")
    i = np.arange(n, dtype=np.float64)
    z = 1.0 - 2.0 * (i + .5) / n
    theta = i * (math.pi * (3.0 - math.sqrt(5.0)))
    xy = np.sqrt(1.0 - z * z)
    return np.stack((np.cos(theta) * xy, np.sin(theta) * xy, z), axis=-1)


def _is_heavy(atom: Mapping[str, Any]) -> bool:
    name = str(atom["name"]).upper()
    return float(atom["radius"]) >= .7 and re.match(r"^[0-9]*H", name) is None


def _shell_points(atoms: Sequence[Mapping[str, Any]], valid_keys: set[str],
                  *, n_directions: int = SPHERE_DIRECTIONS):
    heavy = [atom for atom in atoms if _is_heavy(atom)]
    if len(heavy) < 10:
        raise ValueError("Insufficient PQR heavy atoms for surface construction")
    center = np.asarray([a["xyz"] for a in heavy], dtype=np.float64)
    radii = np.asarray(
        [float(a["radius"]) + PROBE_A + OUTWARD_OFFSET_A for a in heavy], dtype=np.float64
    )
    if not np.isfinite(center).all() or not np.isfinite(radii).all():
        raise ValueError("Nonfinite PQR center/radius")
    keys = [f'{a["chain"]}:{a["seq"]}:{a["icode"]}' for a in heavy]
    dirs = fibonacci_directions(n_directions)
    tree = cKDTree(center)
    max_radius = float(np.max(radii))
    coords: list[np.ndarray] = []
    member_keys: list[str] = []
    areas: list[float] = []
    for i, (xyz, r, key) in enumerate(zip(center, radii, keys)):
        if key not in valid_keys:
            continue
        candidates = xyz + r * dirs
        neighbors = tree.query_ball_point(candidates, r=max_radius + 1.e-5)
        weight = 4.0 * math.pi * r * r / n_directions
        for point, nearby in zip(candidates, neighbors):
            others = [j for j in nearby if j != i]
            if others:
                dv = center[others] - point
                d2 = np.einsum("ij,ij->i", dv, dv)
                if np.any(d2 < np.square(radii[others] - 1.e-6)):
                    continue
            coords.append(point)
            member_keys.append(key)
            areas.append(weight)
    if not coords:
        raise ValueError("No solvent-accessible sampling points found")
    return np.asarray(coords, np.float64), member_keys, np.asarray(areas, np.float64), len(heavy)


def sample_residue_sas_potential(
    atoms: Sequence[Mapping[str, Any]],
    grid: np.ndarray,
    origin: np.ndarray,
    spacing: np.ndarray,
    surface_residue_keys: set[str],
    *,
    n_directions: int = SPHERE_DIRECTIONS,
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    """Return per-residue area-weighted potential and auditable method metadata.

    The source atom/PQR coordinate set is unmodified. An absent surface point
    yields an explicit null field in the caller, never a synthetic zero.
    """
    coords, keys, weights, n_heavy = _shell_points(
        atoms, surface_residue_keys, n_directions=n_directions
    )
    positions = (coords - np.asarray(origin, dtype=np.float64)) / np.asarray(spacing, dtype=np.float64)
    grid_shape = np.asarray(grid.shape, dtype=np.int64)
    inside = np.all((positions >= 0) & (positions < grid_shape - 1), axis=1)
    if not np.all(inside):
        raise ValueError(
            f"APBS potential grid excludes {int((~inside).sum())}/{len(inside)} "
            "solvent-facing sampling points"
        )
    vals = map_coordinates(grid, positions.T, order=1, mode="constant",
                           cval=np.nan, prefilter=False)
    if not np.isfinite(vals).all():
        raise ValueError("Nonfinite interpolated solvent-facing APBS potential")

    totals: dict[str, list[float]] = defaultdict(lambda: [0., 0., 0.])
    for key, area, potential in zip(keys, weights, vals):
        entry = totals[key]
        entry[0] += float(area)
        entry[1] += float(area) * float(potential)
        entry[2] += 1
    residues = {
        key: {
            "surface_potential_kT_per_e": round(s[1] / s[0], 7),
            "surface_sampling_area_A2": round(s[0], 4),
            "surface_potential_samples": int(s[2]),
        }
        for key, s in totals.items() if s[0] > 0
    }
    meta = {
        "method": "area-weighted APBS potential on discrete solvent-accessible probe-center shell",
        "surface_type": "approximate SAS shell, not molecular SES",
        "probe_A": PROBE_A,
        "outward_offset_A": OUTWARD_OFFSET_A,
        "directions_per_heavy_atom": n_directions,
        "n_heavy_atoms": n_heavy,
        "n_surface_samples": len(keys),
        "n_exposed_residues_requested": len(surface_residue_keys),
        "n_exposed_residues_with_surface_potential": len(residues),
        "n_exposed_residues_without_shell_samples": len(surface_residue_keys) - len(residues),
        "estimated_total_sampling_area_A2": round(float(weights.sum()), 3),
        "unit": "kT/e",
        "interpretation": "electrostatic descriptor only; not adsorption free energy",
    }
    return residues, meta
