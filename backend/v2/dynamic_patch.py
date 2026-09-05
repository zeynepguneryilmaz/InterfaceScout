"""Dynamic interpretation of V1 anchors using GNM.

The dynamic layer does not create adsorption anchors on its own. It asks whether
V1-defined native-state surface anchors sit in mobile/coherent neighborhoods.
"""

from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np


def _surface_key_set(v1_result: dict) -> set[str]:
    return {str(r.get("key")) for r in v1_result.get("surface_residues", []) if r.get("key")}


def analyze_anchor_dynamics(
    *,
    anchor_rows: Iterable[dict],
    gnm: dict,
    v1_result: dict,
    radius_A: float = 12.0,
    min_abs_corr: float = 0.35,
    max_extensions: int = 12,
) -> List[dict]:
    """Annotate V1 anchor centers with native-state dynamic neighborhoods.

    Candidate extension residues must satisfy all of the following:
    - solvent-exposed according to V1;
    - within ``radius_A`` C-alpha distance of the anchor;
    - absolute GNM cross-correlation >= ``min_abs_corr``.

    This is a descriptive expansion, not a fitted score.
    """
    nodes = gnm["nodes"]
    index: Dict[str, int] = gnm["index"]
    corr: np.ndarray = gnm["correlation_matrix"]
    residue_metrics: Dict[str, dict] = gnm["residue_metrics"]
    surface_keys = _surface_key_set(v1_result)

    out: List[dict] = []
    for rank, anchor in enumerate(anchor_rows, start=1):
        key = str(anchor.get("center_key") or anchor.get("key") or "")
        if key not in index:
            out.append({
                "rank": rank,
                "anchor_key": key,
                "status": "anchor_not_in_gnm_nodes",
                "dynamic_extension": [],
            })
            continue

        ai = index[key]
        acoord = nodes[ai]["coord"]
        am = residue_metrics[key]
        candidates: List[dict] = []

        for node in nodes:
            nkey = node["key"]
            if nkey == key or nkey not in surface_keys:
                continue
            ni = index[nkey]
            dist = float(np.linalg.norm(node["coord"] - acoord))
            c = float(corr[ai, ni])
            if dist > float(radius_A) or abs(c) < float(min_abs_corr):
                continue
            nm = residue_metrics[nkey]
            candidates.append({
                "key": nkey,
                "chain": node["chain"],
                "res_seq": node["res_seq"],
                "icode": node["icode"],
                "res_name": node["res_name"],
                "ca_distance_A": dist,
                "gnm_correlation": c,
                "abs_gnm_correlation": abs(c),
                "normalized_fluctuation": float(nm["normalized_fluctuation"]),
                "contact_degree": int(nm["contact_degree"]),
            })

        candidates.sort(key=lambda x: (x["abs_gnm_correlation"], -x["ca_distance_A"]), reverse=True)
        ext = candidates[: int(max_extensions)]
        positive = [x for x in ext if x["gnm_correlation"] > 0]
        coherent_abs = float(np.mean([x["abs_gnm_correlation"] for x in ext])) if ext else 0.0
        coherent_signed = float(np.mean([x["gnm_correlation"] for x in ext])) if ext else 0.0

        out.append({
            "rank": rank,
            "anchor_key": key,
            "anchor_chain": anchor.get("chain"),
            "anchor_res_seq": anchor.get("res_seq"),
            "anchor_res_name": anchor.get("res_name"),
            "v1_multiscale_persistence": float(anchor.get("multiscale_persistence", 0.0)),
            "anchor_normalized_fluctuation": float(am["normalized_fluctuation"]),
            "anchor_contact_degree": int(am["contact_degree"]),
            "search_radius_A": float(radius_A),
            "min_abs_corr": float(min_abs_corr),
            "n_dynamic_neighbors": len(candidates),
            "n_reported_extensions": len(ext),
            "n_positive_correlated_extensions": len(positive),
            "mean_abs_corr_reported": coherent_abs,
            "mean_signed_corr_reported": coherent_signed,
            "dynamic_extension": ext,
        })
    return out
