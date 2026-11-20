"""Sequence-aware, one-to-one residue correspondence across separately numbered structures.

No adsorption features or learned parameters enter mapping.  The mapping is used
only for descriptive comparisons/consensus, never for the 11 chemistry scores.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any, Mapping

import numpy as np
from Bio.Align import PairwiseAligner


def _chains(rows: list[Mapping[str, Any]]) -> dict[str, list[Mapping[str, Any]]]:
    by_chain: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    seen: set[str] = set()
    for row in rows:
        key = str(row["key"])
        if key in seen:
            raise ValueError(f"Duplicate residue identifier in a structure: {key}")
        seen.add(key)
        by_chain[str(row["chain"])].append(row)
    for chain_rows in by_chain.values():
        chain_rows.sort(key=lambda r: (int(r["res_seq"]), str(r.get("icode") or "")))
    return dict(by_chain)


def _pairs(a: list[Mapping[str, Any]], b: list[Mapping[str, Any]]) -> list[tuple[int, int]]:
    # Use a high mismatch/gap penalty to avoid sliding short repeated motifs.
    # Residue identities (not one-letter code aliases) determine valid matches.
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    # A stable, one-to-one alignment; ambiguous homomer chain assignment is handled below.
    labels = sorted({str(x["res_name"]) for x in a + b})
    alphabet = {r: chr(0x100 + i) for i, r in enumerate(labels)}
    aa = "".join(alphabet[str(r["res_name"])] for r in a)
    bb = "".join(alphabet[str(r["res_name"])] for r in b)
    al = aligner.align(aa, bb)[0]
    return [(int(i), int(j)) for i, j in zip(*al.indices)
            if i >= 0 and j >= 0 and a[int(i)]["res_name"] == b[int(j)]["res_name"]]


def _coordinate_cost(
    aa: list[Mapping[str, Any]], bb: list[Mapping[str, Any]], pairs: list[tuple[int, int]],
    origin_a: np.ndarray | None, origin_b: np.ndarray | None,
) -> float | None:
    if not pairs or origin_a is None or origin_b is None:
        return None
    try:
        pa = np.array([[float(aa[i][s]) for s in ("x", "y", "z")] for i, _ in pairs])
        pb = np.array([[float(bb[j][s]) for s in ("x", "y", "z")] for _, j in pairs])
        # Only used as an ambiguity resolver when structures share an orientation.
        # It is not a conformational similarity or RMSD metric.
        return float(np.sqrt(np.mean(np.sum(
            ((pa - origin_a) - (pb - origin_b)) ** 2, axis=1))))
    except (KeyError, ValueError, TypeError):
        return None


def _origin(rows: list[Mapping[str, Any]]) -> np.ndarray | None:
    try:
        return np.mean(np.asarray(
            [[float(r[k]) for k in ("x", "y", "z")] for r in rows], dtype=float), axis=0)
    except (KeyError, ValueError, TypeError):
        return None


def match_residues(
    left_rows: list[Mapping[str, Any]], right_rows: list[Mapping[str, Any]],
    min_coverage: float = 0.80,
) -> dict[str, Any]:
    """Return reference (left) residue keys matched to right keys.

    Nonidentical residues and alignment gaps are excluded.  Repeated identical
    chains are never arbitrarily paired when sequence and coordinate evidence
    cannot distinguish their assignment.
    """
    if not left_rows or not right_rows:
        raise ValueError("Both structures must contain standard amino-acid residues")
    left, right = _chains(left_rows), _chains(right_rows)
    oa, ob = _origin(left_rows), _origin(right_rows)
    candidates = []
    for ca, ar in left.items():
        for cb, br in right.items():
            pairs = _pairs(ar, br)
            if not pairs:
                continue
            n = len(pairs)
            cover = n / max(len(ar), len(br))
            candidates.append({
                "left": ca, "right": cb, "pairs": pairs, "count": n,
                "fraction": cover,
                "cost": _coordinate_cost(ar, br, pairs, oa, ob),
                "same_label": ca == cb,
            })

    selected = []
    used_a: set[str] = set()
    used_b: set[str] = set()
    while True:
        options = [c for c in candidates if c["left"] not in used_a
                   and c["right"] not in used_b
                   and c["fraction"] >= min_coverage]
        if not options:
            break
        options.sort(key=lambda c: (-c["count"], -c["fraction"]))
        head = options[0]
        # Only treat comparable chain candidates as alternatives.  Do not choose
        # a homomer subunit by alphabetical chain name after renaming.
        tied = [c for c in options if c["count"] == head["count"]
                and abs(c["fraction"] - head["fraction"]) < 1e-9]
        if len(tied) > 1:
            # Chain-centred coordinates can resolve re-labels in a shared
            # coordinate frame, but not rotations or genuinely symmetric copies.
            with_coordinates = [c for c in tied if c["cost"] is not None]
            if len(with_coordinates) == len(tied):
                with_coordinates.sort(key=lambda c: c["cost"])
                # Only use raw coordinate proximity when the two inputs share
                # nearly the same global frame.  After a rotation, a small
                # difference between two large costs is not valid evidence.
                best = with_coordinates[0]
                competing = [
                    candidate for candidate in with_coordinates[1:]
                    if candidate["left"] == best["left"] or
                       candidate["right"] == best["right"]
                ]
                # Two noncompeting matches can both have zero cost; they should
                # be accepted one by one.  Only alternative pairings of the
                # SAME chain make a match ambiguous.
                if (best["cost"] < 2.0 and
                        (not competing or
                         min(c["cost"] for c in competing) - best["cost"] > 1e-3)):
                    head = best
                    head["assignment_basis"] = "coordinate_consistent"
                else:
                    same = [c for c in tied if c["same_label"]]
                    if len(same) == 1:
                        head = same[0]
                        head["assignment_basis"] = "identical_chain_label_assumed"
                    else:
                        raise ValueError(
                            "Ambiguous repeated protein-chain assignment after "
                            "sequence alignment. Orient structures consistently, "
                            "choose corresponding chain subsets, or supply an "
                            "explicit structural alignment.")
            else:
                same = [c for c in tied if c["same_label"]]
                if len(same) == 1:
                    head = same[0]
                    head["assignment_basis"] = "identical_chain_label_assumed"
                else:
                    raise ValueError(
                        "Ambiguous equivalent protein-chain assignment; provide "
                        "matching chain identities for repeated homologous chains.")
        head.setdefault("assignment_basis", "unique_sequence")
        selected.append(head)
        used_a.add(head["left"])
        used_b.add(head["right"])

    a_to_b: dict[str, str] = {}
    for c in selected:
        aa, bb = left[c["left"]], right[c["right"]]
        for i, j in c["pairs"]:
            a_to_b[str(aa[i]["key"])] = str(bb[j]["key"])
    coverage_a = len(a_to_b) / len(left_rows)
    coverage_b = len(a_to_b) / len(right_rows)
    if min(coverage_a, coverage_b) < min_coverage:
        raise ValueError(
            "Insufficient identity-matched residue coverage after sequence alignment "
            f"({coverage_a:.1%} of A; {coverage_b:.1%} of B; require "
            f"{min_coverage:.0%} for each). Check that structures represent the "
            "same protein and assembly, or select corresponding chains.")

    return {
        "a_to_b": a_to_b,
        "b_to_a": {v: k for k, v in a_to_b.items()},
        "chain_matches": [{"chain_a": c["left"], "chain_b": c["right"],
                           "n_identity_matches": c["count"],
                           "assignment_basis": c["assignment_basis"]} for c in selected],
        "warnings": [
            "Repeated-chain assignment assumed matching chain labels because "
            "coordinates did not establish a unique correspondence. Verify "
            "chain assignments before interpreting conformational changes."
        ] if any(c["assignment_basis"] == "identical_chain_label_assumed"
                 for c in selected) else [],
        "coverage_a": coverage_a, "coverage_b": coverage_b,
        "method": "chain-aware sequence alignment with identity-only residue correspondence",
        "n_matched": len(a_to_b),
    }
