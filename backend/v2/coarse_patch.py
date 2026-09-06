"""Weight-free coarse interface patch construction for InterfaceScout V2.

V2 predicts *surface regions*, not residue affinities. Experimental adsorption
labels never enter this module.

Construction logic
------------------
1. Frozen V1 chemistry maps provide local patch-persistence values.
2. Candidate centres are local maxima on the exposed protein surface at the
   already-frozen 8 A V1 patch scale.
3. A V2 patch is the non-transitive 8 A surface neighbourhood of one local
   maximum, restricted to the same coarse outward-facing hemisphere.
4. Chemistry, accessibility, dynamics and orientation are reported separately.
5. Patches are compared by Pareto dominance; no empirical weighted sum is used.

The non-transitive definition is deliberate. Connected-component growth can
percolate through a protein surface and turn a local biointerface hypothesis
into an unrealistically large fraction of the protein.
"""

from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np

from .geometry import (
    build_surface_geometry,
    ca_distance,
    same_face,
    patch_orientation_coherence,
    patch_diameter_A,
)

PATCH_SCALE_A = 8.0  # inherited from the frozen V1 multiscale patch definition


def _pair_values(keys: List[str], matrix: np.ndarray, index: Dict[str, int]) -> List[float]:
    vals: List[float] = []
    for i in range(len(keys)):
        if keys[i] not in index:
            continue
        ii = index[keys[i]]
        for j in range(i + 1, len(keys)):
            if keys[j] not in index:
                continue
            jj = index[keys[j]]
            vals.append(float(matrix[ii, jj]))
    return vals


def _local_maxima(channel: dict, geometry: dict) -> List[dict]:
    """Select nonredundant V1 patch maxima without a benchmark-fitted threshold.

    Centres are processed by the frozen V1 persistence ranking. A lower-ranked
    centre is suppressed only when it lies within the same 8 A local patch and
    on the same coarse face as an already retained maximum.
    """
    rows = [
        r for r in channel.get("patch_centers", [])
        if r.get("center_key") in geometry["coords"]
        and float(r.get("multiscale_persistence", 0.0)) > 0.0
    ]
    rows.sort(
        key=lambda r: (
            -float(r.get("multiscale_persistence", 0.0)),
            -float(r.get("multiscale_geomean", 0.0)),
            str(r.get("center_key")),
        )
    )

    kept: List[dict] = []
    for row in rows:
        key = str(row["center_key"])
        redundant = False
        for prev in kept:
            pkey = str(prev["center_key"])
            if ca_distance(key, pkey, geometry) <= PATCH_SCALE_A and same_face(key, pkey, geometry):
                redundant = True
                break
        if not redundant:
            kept.append(row)
    return kept


def _patch_around_center(center: dict, chemistry_rows: Dict[str, dict], geometry: dict) -> dict:
    ckey = str(center["center_key"])
    members = []
    for key in geometry["coords"]:
        if ca_distance(ckey, key, geometry) <= PATCH_SCALE_A and same_face(ckey, key, geometry):
            members.append(key)
    if ckey not in members:
        members.append(ckey)

    seeds = sorted(k for k in members if k in chemistry_rows)
    return {
        "center": ckey,
        "center_row": center,
        "members": sorted(set(members)),
        "seeds": seeds,
    }


def _patch_descriptors(patch: dict, geometry: dict, gnm: dict) -> dict:
    members = list(patch["members"])
    seeds = list(patch["seeds"])
    center = patch["center_row"]

    # Chemistry is represented as composition (fraction of patch residues that
    # are chemically compatible in the frozen V1 channel), while patch
    # coherence is the V1 multiscale local persistence at the patch centre.
    chemistry_fraction = float(len(seeds) / len(members)) if members else 0.0
    patch_coherence = float(center.get("multiscale_persistence", 0.0)) / 100.0

    access_vals = [float(geometry["scrsa"].get(k, 0.0)) for k in members]
    accessibility = float(np.mean(access_vals)) if access_vals else 0.0

    pair_corr = _pair_values(members, gnm["correlation_matrix"], gnm["index"])
    dynamic_coupling = float(np.mean(np.abs(pair_corr))) if pair_corr else 0.0
    dynamic_signed = float(np.mean(pair_corr)) if pair_corr else 0.0
    orientation = patch_orientation_coherence(members, geometry)

    meta = geometry["meta"]
    rows = [meta[k] for k in members if k in meta]
    seed_rows = [meta[k] for k in seeds if k in meta]
    chains = sorted({r["chain"] for r in rows})
    resseqs = [int(r["res_seq"]) for r in rows]

    return {
        "center_key": patch["center"],
        "center_residue": meta.get(patch["center"]),
        "members": members,
        "seed_members": seeds,
        "n_members": len(members),
        "n_chemistry_seeds": len(seeds),
        "chains": chains,
        "residue_span_min": min(resseqs) if resseqs else None,
        "residue_span_max": max(resseqs) if resseqs else None,
        "member_residues": rows,
        "seed_residues": seed_rows,
        "diameter_A": patch_diameter_A(members, geometry),
        "chemistry_support": chemistry_fraction,
        "mean_accessibility": accessibility,
        "patch_coherence": patch_coherence,
        "dynamic_coupling_abs": dynamic_coupling,
        "dynamic_coupling_signed": dynamic_signed,
        "orientation_coherence": orientation,
        "spatial_coherence": "fixed 8 A local neighbourhood around a V1 patch maximum",
    }


def _dominates(a: dict, b: dict) -> bool:
    fields = (
        "chemistry_support",
        "mean_accessibility",
        "patch_coherence",
        "dynamic_coupling_abs",
        "orientation_coherence",
    )
    av = [float(a[f]) for f in fields]
    bv = [float(b[f]) for f in fields]
    return all(x >= y for x, y in zip(av, bv)) and any(x > y for x, y in zip(av, bv))


def assign_pareto_fronts(patches: List[dict]) -> List[dict]:
    remaining = set(range(len(patches)))
    front = 1
    while remaining:
        current = []
        for i in remaining:
            if not any(_dominates(patches[j], patches[i]) for j in remaining if j != i):
                current.append(i)
        for i in current:
            patches[i]["pareto_front"] = front
        remaining -= set(current)
        front += 1

    patches.sort(
        key=lambda p: (
            int(p["pareto_front"]),
            -float(p["patch_coherence"]),
            -float(p["chemistry_support"]),
            str(p.get("center_key", "")),
        )
    )
    for rank, p in enumerate(patches, start=1):
        p["display_rank"] = rank
        p["classification"] = (
            "primary plausible interface" if p["pareto_front"] == 1 else "alternative interface"
        )
    return patches


def build_coarse_patches(*, v1_result: dict, chemistry: str, gnm: dict) -> List[dict]:
    channel = v1_result.get("chemistries", {}).get(chemistry)
    if not channel:
        raise ValueError(f"Unknown or unavailable chemistry channel: {chemistry}")

    geometry = build_surface_geometry(v1_result, gnm)
    chemistry_rows = {
        str(r["key"]): r
        for r in channel.get("residues", [])
        if r.get("key") and float(r.get("local_score", 0.0)) > 0.0
    }
    if not chemistry_rows:
        return []

    maxima = _local_maxima(channel, geometry)
    raw_patches = [_patch_around_center(c, chemistry_rows, geometry) for c in maxima]
    patches = [_patch_descriptors(p, geometry, gnm) for p in raw_patches]
    return assign_pareto_fronts(patches)
