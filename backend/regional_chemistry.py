"""Descriptive regional chemistry profiles derived from frozen InterfaceScout maps.

No adsorption labels, optimized weights, fit coefficients, or new ranking are used.
All descriptors reuse canonical single-map residues and the existing 8 Å patches.
"""
from __future__ import annotations

from itertools import combinations
from math import dist, isfinite
from statistics import median
from typing import Any, Mapping, Sequence


METHOD_VERSION = "regional-chemistry-descriptive-v1"


def _local_support(analysis: Mapping[str, Any], name: str) -> dict[str, float]:
    """Canonical local_score = scRSA * compatible-class * state availability."""
    scores: dict[str, float] = {}
    for row in analysis["maps"][name].get("residues", []):
        value = float(row.get("local_score") or 0.0)
        if not isfinite(value) or value < 0.0 or value > 1.000001:
            raise ValueError("Invalid canonical local_score for map " + name)
        key = str(row["key"])
        if key in scores:
            raise ValueError("Duplicate residue key in chemistry map: " + key)
        scores[key] = value
    return scores


def _distance_coordinates(analysis: Mapping[str, Any]) -> dict[str, tuple[float, float, float]]:
    raw = analysis.get("diagnostics", {}).get("surface_ca_xyz_A") or {}
    coords: dict[str, tuple[float, float, float]] = {}
    for key, xyz in raw.items():
        if not isinstance(xyz, (list, tuple)) or len(xyz) != 3:
            continue
        vals = tuple(float(v) for v in xyz)
        if all(isfinite(v) for v in vals):
            coords[str(key)] = vals
    return coords


def annotate_regions(
    analysis: Mapping[str, Any],
    selected_maps: Sequence[str],
    regions: list[dict[str, Any]],
) -> dict[str, Any]:
    """Add independent-by-definition DESCRIPTIVE fields without changing region order.

    L_c(r): existing canonical nonnegative local support in [0,1].
    D_c(P): sum(L_c(r) for r in P) / |P|.
    F_c(P): number of L_c(r)>0 divided by |P|.
    Contrasts are differences from the medians across the already-defined
    *overlapping candidate patch set*, not from statistically independent draws.
    """
    names = list(dict.fromkeys(selected_maps))
    channels = {name: _local_support(analysis, name) for name in names}
    coords = _distance_coordinates(analysis)
    reference_densities = {name: [] for name in names}
    profiles: list[dict[str, Any]] = []
    for region in regions:
        members = list(dict.fromkeys(str(k) for k in region["member_keys"]))
        if not members:
            continue
        count = len(members)
        density = {
            name: sum(channels[name].get(k, 0.0) for k in members) / count
            for name in names
        }
        fraction = {
            name: sum(channels[name].get(k, 0.0) > 0.0 for k in members) / count
            for name in names
        }
        supporting = {
            name: [k for k in members if channels[name].get(k, 0.0) > 0.0]
            for name in names
        }
        pairwise: dict[str, float | None] = {}
        for a, b in combinations(names, 2):
            distances = [
                dist(coords[x], coords[y])
                for x in supporting[a] for y in supporting[b]
                if x != y and x in coords and y in coords
            ]
            pairwise[a + "|" + b] = round(min(distances), 3) if distances else None
        p = {
            "method_version": METHOD_VERSION,
            "member_count": count,
            "mean_local_support_by_map": {k: round(v, 6) for k, v in density.items()},
            "supporting_residue_fraction_by_map": {k: round(v, 6) for k, v in fraction.items()},
            "supporting_residue_count_by_map": {k: len(v) for k, v in supporting.items()},
            "weakest_mean_support": round(min(density.values()), 6) if density else None,
            "min_distinct_intermap_CA_distance_A": (
                round(min(v for v in pairwise.values() if v is not None), 3)
                if any(v is not None for v in pairwise.values()) else None
            ),
            "intermap_distinct_CA_distances_A": pairwise,
            "relative_to_candidate_median_by_map": {},
            "reference": "same-protein overlapping combined candidate patch set",
            "uses_experimental_labels": False,
            "material_specific_binding_prediction": False,
        }
        profiles.append(p)
        for name in names:
            reference_densities[name].append(density[name])
    medians = {
        name: median(reference_densities[name]) if reference_densities[name] else None
        for name in names
    }
    for row, p in zip(regions, profiles):
        p["relative_to_candidate_median_by_map"] = {
            name: (
                round(p["mean_local_support_by_map"][name] - medians[name], 6)
                if medians[name] is not None else None
            )
            for name in names
        }
        row["regional_profile"] = p
    return {
        "method_version": METHOD_VERSION,
        "definition": {
            "mean_local_support": "mean over all region members of the unchanged canonical map local_score, with absent chemistry support set to zero",
            "supporting_residue_fraction": "fraction of region members with canonical local_score > 0 for each selected map",
            "weakest_mean_support": "minimum of the per-channel means, descriptive only; NOT used to rerank patches",
            "cross_channel_distance": "minimum C-alpha distance between DISTINCT residues supporting different maps; not interatomic contact or bond length",
            "regional_contrast": "per-channel region mean minus the median of the overlapping candidate region means",
        },
        "reference_candidate_region_count": len(regions),
        "reference_median_local_support_by_map": {
            k: round(v, 6) if v is not None else None for k, v in medians.items()
        },
        "regional_proximity_available": bool(coords),
        "ranking_unchanged": True,
        "no_fitted_weights": True,
        "experimental_labels_used": False,
        "inference_warning": (
            "Candidate patches overlap and are not independent observations. "
            "Descriptive support, proximity and contrast are not adsorption affinity, "
            "unique binding orientation, experimental agreement, or calibrated confidence."
        ),
    }
