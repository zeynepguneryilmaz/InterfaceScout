"""Decision-oriented, additive analyses for InterfaceScout.

This module does not modify the frozen 11-map chemistry kernel, its patch
ranking, or claim to predict adsorption states, contact energies, or kinetics.
"""
from __future__ import annotations

from typing import Any, Mapping, Sequence

RESEARCH_GOALS: dict[str, dict[str, Any]] = {
    "explore": {
        "label": "Explore a protein surface",
        "question": "Which chemically distinct and accessible regions does this structure expose?",
        "recommended_maps": [],
        "essentials": ["residue accessibility", "individual chemistry maps", "candidate regions"],
        "conditional": ["map intersections when a research question involves two chemistries"],
        "not_required": ["GNM", "ANM", "RINpy", "MD trajectory"],
        "caution": "Surface descriptions are structure-dependent; no adsorption prediction is made.",
    },
    "noncovalent": {
        "label": "Plan a noncovalent protein–material experiment",
        "question": "Which accessible protein regions have chemistry relevant to my chosen material surface?",
        "recommended_maps": ["anionic", "cationic", "hydrophobic", "pi_carbon",
                             "hbond_donor", "hbond_acceptor"],
        "essentials": ["material-relevant maps", "spatial regions", "residue and chemistry lists",
                       "overlap between relevant channels"],
        "conditional": ["compare conformations if the protein is flexible or structures disagree",
                        "inspect experimental constraints, pH and functional sites externally"],
        "not_required": ["GNM", "ANM", "RINpy", "MD trajectory"],
        "caution": "Selecting a candidate region is not proof of preferred adsorption pose or affinity.",
    },
    "immobilization": {
        "label": "Design an immobilization or residue-labeling experiment",
        "question": "Which chemically reactive side chains are accessible, and where are they located?",
        "recommended_maps": ["metal_coord", "soft_metal_sulfur", "cationic", "anionic"],
        "essentials": ["exposed Cys/Lys/Asp/Glu/His residues", "side-chain accessibility",
                       "nearby residues and spatial position", "functional-site exclusion by user data"],
        "conditional": ["compare experimental conformations if a candidate residue is buried in one",
                        "inspect residue-level reaction chemistry and experimental conditions"],
        "not_required": ["GNM", "ANM", "RINpy", "MD trajectory"],
        "caution": "Compatibility channels are NOT reaction-rate or covalent-coupling predictions. "
                   "For example, the soft-metal/sulfur channel must not be used for all thiol chemistry.",
    },
    "simulation": {
        "label": "Prepare computational interface simulations",
        "question": "Which nonredundant candidate surfaces should I include in an explicit orientation/MD study?",
        "recommended_maps": ["hydrophobic", "pi_carbon", "anionic", "cationic",
                             "hbond_donor", "hbond_acceptor"],
        "essentials": ["candidate region alternatives", "coverage at a fixed residue budget",
                       "map intersections", "residue lists and reproducible input settings"],
        "conditional": ["multi-PDB comparison before choosing an initial conformation",
                        "MD-trajectory analysis only after simulation results exist"],
        "not_required": ["RINpy", "GNM", "ANM"],
        "caution": "No unique pose, favorable free energy, or adsorption pathway is inferred.",
    },
    "conformations": {
        "label": "Compare experimental or simulated conformations",
        "question": "Which surface properties and regions persist across provided structures?",
        "recommended_maps": [],
        "essentials": ["consistent residue identity", "side-chain accessibility changes",
                       "Top-1/3/5 region overlap by map", "changed/retained residues"],
        "conditional": ["trajectory frame selection if user has MD data",
                        "GNM or ANM only for a separately motivated dynamics question"],
        "not_required": ["RINpy"],
        "caution": "Two PDB structures are not an equilibrium ensemble and do not yield occupancies.",
    },
}


def list_research_goals() -> list[dict[str, Any]]:
    return [{"id": key, **value} for key, value in RESEARCH_GOALS.items()]


def _patch_members(map_result: Mapping[str, Any], k: int,
                   eligible: set[str] | None = None) -> set[str]:
    members: set[str] = set()
    for patch in map_result.get("patches", [])[:k]:
        members.update(str(v) for v in patch.get("members", []))
    return members if eligible is None else members & eligible


def _jaccard(a: set[str], b: set[str]) -> float | None:
    union = a | b
    return len(a & b) / len(union) if union else None


def compare_results(a: Mapping[str, Any], b: Mapping[str, Any]) -> dict[str, Any]:
    """Compare two runs of the canonical engine, without learning new scores.

    Residue keys must correspond. We do not silently align sequences or
    interpret different residue names at the same position as identical.
    """
    left = {r["key"]: r for r in a.get("residue_profiles", [])}
    right = {r["key"]: r for r in b.get("residue_profiles", [])}
    if not left or not right:
        raise ValueError("Both analyses must contain residue_profiles from the canonical engine")
    shared_keys = set(left) & set(right)
    eligible = {key for key in shared_keys
                if left[key].get("res_name") == right[key].get("res_name")}
    coverage = len(eligible) / max(1, min(len(left), len(right)))
    if coverage < 0.8:
        raise ValueError("Residue mapping is insufficient (<80% identity-matched coverage). "
                         "Supply matched chains and residue numbering, or align/renumber structures first.")

    maps_a, maps_b = a.get("maps", {}), b.get("maps", {})
    if set(maps_a) != set(maps_b):
        raise ValueError("Map sets differ between structures")
    if a.get("input", {}).get("pH") != b.get("input", {}).get("pH"):
        raise ValueError("Use the same pH in both runs to isolate structural changes")

    residue_changes = []
    exposed_a, exposed_b = 0, 0
    for key in sorted(eligible):
        ra, rb = left[key], right[key]
        xa, xb = bool(ra.get("surface_exposed")), bool(rb.get("surface_exposed"))
        exposed_a += xa
        exposed_b += xb
        residue_changes.append({
            "key": key,
            "res_name": ra["res_name"],
            "exposed_a": xa, "exposed_b": xb,
            "changed_exposure_class": xa != xb,
            "scRSA_a": ra.get("scrsa_raw"),
            "scRSA_b": rb.get("scrsa_raw"),
            "delta_scRSA": round(float(rb.get("scrsa_raw", 0.0)) -
                                 float(ra.get("scrsa_raw", 0.0)), 5),
            "sidechain_sasa_a": ra.get("sidechain_sasa"),
            "sidechain_sasa_b": rb.get("sidechain_sasa"),
            "delta_sidechain_sasa": round(float(rb.get("sidechain_sasa", 0.0)) -
                                          float(ra.get("sidechain_sasa", 0.0)), 3),
        })

    compared_maps = []
    for name in a.get("map_order", maps_a.keys()):
        ma, mb = maps_a[name], maps_b[name]
        overlap = []
        for k in (1, 3, 5):
            x, y = _patch_members(ma, k, eligible), _patch_members(mb, k, eligible)
            overlap.append({
                "top_k": k, "jaccard": _jaccard(x, y),
                "n_a": len(x), "n_b": len(y), "n_shared": len(x & y),
                "retained": sorted(x & y),
                "lost_in_b": sorted(x - y),
                "new_in_b": sorted(y - x),
                "surface_coverage_a": round(len(x) / max(exposed_a, 1), 4),
                "surface_coverage_b": round(len(y) / max(exposed_b, 1), 4),
            })
        pa, pb = ma.get("patches", []), mb.get("patches", [])
        compared_maps.append({
            "key": name, "label": ma.get("label", name),
            "top1_center_a": pa[0].get("center_key") if pa else None,
            "top1_center_b": pb[0].get("center_key") if pb else None,
            "n_candidate_regions_a": len(pa), "n_candidate_regions_b": len(pb),
            "top_k": overlap,
        })

    return {
        "purpose": "Descriptive comparison of protein-surface candidates across two conformations",
        "structure_a": a.get("input", {}),
        "structure_b": b.get("input", {}),
        "mapping": {
            "n_positions_a": len(left), "n_positions_b": len(right),
            "n_key_matches": len(shared_keys),
            "n_identity_matches": len(eligible),
            "coverage_of_shorter_structure": round(coverage, 4),
            "not_compared": sorted((set(left) | set(right)) - eligible),
        },
        "accessibility": {
            "n_exposed_matched_a": exposed_a,
            "n_exposed_matched_b": exposed_b,
            "n_exposure_status_changes": sum(r["changed_exposure_class"] for r in residue_changes),
        },
        "map_comparisons": compared_maps,
        "residue_changes": residue_changes,
        "interpretation": "Jaccard measures overlap of ranked region MEMBERS at fixed Top-K, "
                          "not adsorption prediction accuracy, energy, or conformation probability.",
    }


def combine_maps(analysis: Mapping[str, Any], selected: Sequence[str],
                 max_regions: int = 25) -> dict[str, Any]:
    """Show exact seed overlap separately from spatial co-localization.

    Candidate neighborhoods are reused from the canonical single-map patches;
    no new learned score, material mixture assumption, or pose prediction.
    """
    names = list(dict.fromkeys(selected))
    if not 2 <= len(names) <= 11:
        raise ValueError("Select 2 to 11 distinct maps")
    if not 1 <= max_regions <= 100:
        raise ValueError("max_regions must be between 1 and 100")
    available = analysis.get("maps", {})
    unknown = [name for name in names if name not in available]
    if unknown:
        raise ValueError("Unknown map(s): " + ", ".join(unknown))

    seeds = {
        name: {str(row["key"]) for row in available[name].get("residues", [])
               if float(row.get("local_score") or 0) > 0}
        for name in names
    }
    exact = set.intersection(*(seeds[name] for name in names))
    candidates = {}
    for name in names:
        for patch in available[name].get("patches", []):
            members = {str(key) for key in patch.get("members", [])}
            if not members:
                continue
            hits = {key: sorted(members & seeds[key]) for key in names}
            represented = [key for key, values in hits.items() if values]
            if len(represented) < 2:
                continue
            shared_in_patch = sorted(members & exact)
            region_key = tuple(sorted(members))
            row = {
                "member_keys": list(region_key),
                "n_members": len(members),
                "represented_maps": represented,
                "n_represented_maps": len(represented),
                "chemistry_seed_members": hits,
                "exact_all_map_seed_overlap": shared_in_patch,
                "has_exact_all_map_seed_overlap": bool(shared_in_patch),
                "all_maps_colocated": len(represented) == len(names),
                "distinct_residue_colocation_only": (
                    len(represented) == len(names) and not shared_in_patch),
                "source_map": name,
                "source_center_key": patch.get("center_key"),
            }
            old = candidates.get(region_key)
            if old is None or (row["n_represented_maps"], len(shared_in_patch)) > (
                    old["n_represented_maps"],
                    len(old["exact_all_map_seed_overlap"])):
                candidates[region_key] = row

    regions = sorted(
        candidates.values(),
        key=lambda r: (-r["n_represented_maps"],
                       -len(r["exact_all_map_seed_overlap"]),
                       r["n_members"], r["source_center_key"] or ""),
    )
    exposed_n = len(analysis.get("diagnostics", {}).get("surface_residue_keys", []))
    return {
        "selected_maps": names,
        "exact_same_residue_overlap": sorted(exact),
        "n_exact_overlap": len(exact),
        "n_exposed_residues": exposed_n,
        "regions_total": len(regions),
        "regions": regions[:max_regions],
        "region_rule": "Reuse single-map 8 A same-face patches; require distinct chemistry "
                       "seeds within the SAME patch to call spatial co-localization.",
        "important": "Distinct-residue co-localization is not exact residue intersection, "
                     "and neither is evidence of adsorption or a calculated binding energy.",
    }
