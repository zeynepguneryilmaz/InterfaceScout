"""Transparent V2-alpha patch scoring."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List


def score_patch(candidate: Dict[str, Any]) -> Dict[str, Any]:
    """Score one multi-mechanism patch on a 0-100 scale.

    Alpha score = weighted mean of V1 multiscale persistence channels.
    No fitted coefficients, Ebase values, APBS energies, desolvation terms, or
    orientation energies are used at this stage.
    """
    contributions: Dict[str, float] = {}
    total = 0.0
    for chemistry, channel in candidate.get("channels", {}).items():
        w = float(channel.get("weight", 0.0))
        p = float(channel.get("persistence", 0.0))
        contribution = w * p
        contributions[chemistry] = round(contribution, 4)
        total += contribution

    out = dict(candidate)
    out["score"] = round(total, 4)
    out["mechanism_contributions"] = contributions
    out["dominant_mechanisms"] = [
        k for k, _ in sorted(contributions.items(), key=lambda kv: (-kv[1], kv[0]))
        if contributions[k] > 0
    ]
    return out


def rank_patches(candidates: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    scored = [score_patch(c) for c in candidates]
    scored.sort(
        key=lambda x: (
            -float(x.get("score", 0.0)),
            -len(x.get("compatible_member_union_8A", [])),
            str(x.get("center_key", "")),
        )
    )
    for rank, row in enumerate(scored, start=1):
        row["rank"] = rank
    return scored


def summarize_competition(ranked: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not ranked:
        return {"top_score": 0.0, "second_score": 0.0, "margin": 0.0, "interpretation": "no_patch"}
    s1 = float(ranked[0].get("score", 0.0))
    s2 = float(ranked[1].get("score", 0.0)) if len(ranked) > 1 else 0.0
    margin = s1 - s2
    if s1 <= 0:
        interpretation = "no_positive_signal"
    elif margin < 5:
        interpretation = "multiple_competing_patches"
    elif margin < 15:
        interpretation = "moderate_separation"
    else:
        interpretation = "clear_top_patch"
    return {
        "top_score": round(s1, 4),
        "second_score": round(s2, 4),
        "margin": round(margin, 4),
        "interpretation": interpretation,
    }
