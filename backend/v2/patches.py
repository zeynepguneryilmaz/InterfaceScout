"""Patch assembly for InterfaceScout V2-alpha.

V2-alpha reuses V1's frozen per-chemistry 5/8 Å persistence maps but changes the
prediction unit: a candidate is a spatial patch centre, scored jointly across a
material's mechanism profile rather than by one chemistry channel at a time.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List


def _centers_for_channel(channel: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    return {
        str(row["center_key"]): row
        for row in channel.get("patch_centers", [])
        if row.get("center_key")
    }


def assemble_candidate_patches(
    chemistries: Dict[str, Dict[str, Any]],
    profile_weights: Dict[str, float],
) -> List[Dict[str, Any]]:
    """Create a union of patch centres across all active profile channels."""
    by_channel = {
        chemistry: _centers_for_channel(chemistries.get(chemistry, {}))
        for chemistry in profile_weights
    }
    all_keys = sorted({k for d in by_channel.values() for k in d})
    candidates: List[Dict[str, Any]] = []

    for key in all_keys:
        template = next((d[key] for d in by_channel.values() if key in d), None)
        if template is None:
            continue
        channels: Dict[str, Any] = {}
        member_union = set()
        for chemistry, weight in profile_weights.items():
            row = by_channel[chemistry].get(key)
            persistence = float(row.get("multiscale_persistence", 0.0)) if row else 0.0
            geomean = float(row.get("multiscale_geomean", 0.0)) if row else 0.0
            members = list(row.get("compatible_members_8A", [])) if row else []
            member_union.update(members)
            channels[chemistry] = {
                "weight": float(weight),
                "persistence": persistence,
                "geomean": geomean,
                "compatible_members_8A": members,
            }
        candidates.append({
            "center_key": key,
            "res_name": template.get("res_name"),
            "res_seq": template.get("res_seq"),
            "icode": template.get("icode", ""),
            "chain": template.get("chain"),
            "channels": channels,
            "compatible_member_union_8A": sorted(member_union),
        })
    return candidates


def nonredundant_top_patches(
    ranked: Iterable[Dict[str, Any]],
    top_n: int = 3,
    max_member_jaccard: float = 0.70,
) -> List[Dict[str, Any]]:
    """Keep top patches while suppressing near-duplicate member sets.

    This is deliberately conservative: two high-scoring centres representing
    effectively the same local patch should not occupy two of the Top-3 slots.
    """
    selected: List[Dict[str, Any]] = []
    for candidate in ranked:
        members = set(candidate.get("compatible_member_union_8A", []))
        duplicate = False
        for kept in selected:
            other = set(kept.get("compatible_member_union_8A", []))
            union = members | other
            jaccard = (len(members & other) / len(union)) if union else 1.0
            if jaccard >= max_member_jaccard:
                duplicate = True
                break
        if not duplicate:
            selected.append(candidate)
            if len(selected) >= top_n:
                break
    return selected
