"""Weight-free coarse interface patch construction for InterfaceScout V2.

V2 is intentionally a *region* predictor, not a residue-affinity predictor.
Experimental adsorption labels never enter this module.

Construction logic
------------------
1. Chemistry seeds come from the frozen V1 chemistry channel.
2. Only V1 solvent-exposed residues are eligible.
3. Seeds are grouped into spatially connected, co-facing components using the
   already-frozen 8 A V1 patch scale.
4. Each seed component is expanded by one 8 A shell of exposed residues that
   lie on the same coarse protein face.
5. Dynamics (GNM) and orientation are descriptors of the resulting patch; they
   do not alter membership through fitted cutoffs.
6. Patches are compared by Pareto dominance rather than a weighted sum.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Set

import numpy as np

from .geometry import build_surface_geometry, ca_distance, same_face, patch_orientation_coherence, patch_diameter_A

PATCH_SCALE_A = 8.0  # inherited from the frozen V1 multiscale patch definition


def _connected_components(keys: Iterable[str], geometry: dict) -> List[Set[str]]:
    keys = [k for k in keys if k in geometry["coords"]]
    unseen = set(keys)
    comps: List[Set[str]] = []
    while unseen:
        root = unseen.pop()
        comp = {root}
        stack = [root]
        while stack:
            a = stack.pop()
            neighbors = []
            for b in list(unseen):
                if ca_distance(a, b, geometry) <= PATCH_SCALE_A and same_face(a, b, geometry):
                    neighbors.append(b)
            for b in neighbors:
                unseen.remove(b)
                comp.add(b)
                stack.append(b)
        comps.append(comp)
    return comps


def _expand_component(seed_comp: Set[str], geometry: dict) -> Set[str]:
    patch = set(seed_comp)
    for key in geometry["coords"]:
        if key in patch:
            continue
        for seed in seed_comp:
            if ca_distance(key, seed, geometry) <= PATCH_SCALE_A and same_face(key, seed, geometry):
                patch.add(key)
                break
    return patch


def _merge_overlapping(patches: List[dict]) -> List[dict]:
    """Merge patches when they share at least one physical surface residue.

    This avoids an arbitrary Jaccard threshold.  Overlap means that two seed
    components expand into the same coarse contact region.
    """
    work = [dict(p) for p in patches]
    changed = True
    while changed:
        changed = False
        out: List[dict] = []
        used = [False] * len(work)
        for i, p in enumerate(work):
            if used[i]:
                continue
            members = set(p["members"])
            seeds = set(p["seeds"])
            used[i] = True
            for j in range(i + 1, len(work)):
                if used[j]:
                    continue
                qmem = set(work[j]["members"])
                if members & qmem:
                    members |= qmem
                    seeds |= set(work[j]["seeds"])
                    used[j] = True
                    changed = True
            out.append({"members": members, "seeds": seeds})
        work = out
    return work


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


def _patch_descriptors(patch: dict, chemistry_rows: Dict[str, dict], geometry: dict, gnm: dict) -> dict:
    members = sorted(patch["members"])
    seeds = sorted(patch["seeds"])

    chem_vals = [float(chemistry_rows[k].get("multiscale_persistence", 0.0)) / 100.0 for k in seeds if k in chemistry_rows]
    chemistry_support = float(np.mean(chem_vals)) if chem_vals else 0.0

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
        "chemistry_support": chemistry_support,
        "mean_accessibility": accessibility,
        "dynamic_coupling_abs": dynamic_coupling,
        "dynamic_coupling_signed": dynamic_signed,
        "orientation_coherence": orientation,
        "spatial_coherence": True,
    }


def _dominates(a: dict, b: dict) -> bool:
    fields = (
        "chemistry_support",
        "mean_accessibility",
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
    patches.sort(key=lambda p: (int(p["pareto_front"]), -float(p["chemistry_support"]), -int(p["n_chemistry_seeds"])))
    for rank, p in enumerate(patches, start=1):
        p["display_rank"] = rank
        p["classification"] = "primary plausible interface" if p["pareto_front"] == 1 else "alternative interface"
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
    seeds = [k for k in chemistry_rows if k in geometry["coords"]]
    if not seeds:
        return []

    seed_components = _connected_components(seeds, geometry)
    raw = [
        {"seeds": set(comp), "members": _expand_component(set(comp), geometry)}
        for comp in seed_components
    ]
    merged = _merge_overlapping(raw)
    patches = [_patch_descriptors(p, chemistry_rows, geometry, gnm) for p in merged]
    return assign_pareto_fronts(patches)
