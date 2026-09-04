"""Protein-derived target interface profile for InterfaceScout 2.0.

The profile is intentionally material-agnostic.  It does not recommend or rank
named materials and it does not use a material library.  Instead, each canonical
InterfaceScout chemistry channel is converted into an independently reported
surface-property requirement together with the protein patch(es) that could
engage that property.

No weighted cross-channel score is produced.  Raw scRSA/state-derived quantities
are exposed so users can inspect the protein-side evidence directly.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple


def _residue_key(chain: str, res_seq: int, icode: str = "") -> str:
    return f"{chain}:{int(res_seq)}:{str(icode or '').strip()}"


def _coords_by_key(all_residues: Iterable[Dict[str, Any]]) -> Dict[str, Tuple[float, float, float]]:
    out: Dict[str, Tuple[float, float, float]] = {}
    for r in all_residues:
        try:
            out[r["key"]] = (float(r["x"]), float(r["y"]), float(r["z"]))
        except Exception:
            continue
    return out


def _distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))


def _site_keys(site_annotations: List[Dict[str, Any]]) -> Set[str]:
    keys: Set[str] = set()
    for site in site_annotations:
        for r in site.get("residues", []):
            try:
                keys.add(_residue_key(r["chain"], int(r["res_seq"]), r.get("icode", "")))
            except Exception:
                continue
    return keys


def _nearest_site_distance(center_key: str, coords: Dict[str, Tuple[float, float, float]], site_keys: Set[str]) -> Optional[float]:
    c = coords.get(center_key)
    if c is None:
        return None
    vals = [_distance(c, coords[k]) for k in site_keys if k in coords]
    return min(vals) if vals else None


def _patch_annotation(
    patch: Dict[str, Any],
    coords: Dict[str, Tuple[float, float, float]],
    pdb_site_keys: Set[str],
    protected_keys: Set[str],
) -> Dict[str, Any]:
    members = list(patch.get("compatible_members_8A", []) or [])
    member_set = set(members)
    pdb_overlap = sorted(member_set & pdb_site_keys)
    protected_overlap = sorted(member_set & protected_keys)
    d_pdb = _nearest_site_distance(str(patch.get("center_key", "")), coords, pdb_site_keys)
    d_protected = _nearest_site_distance(str(patch.get("center_key", "")), coords, protected_keys)

    if protected_overlap:
        relation = "overlaps_user_protected_functional_residue"
    elif pdb_overlap:
        relation = "overlaps_PDB_SITE_annotation"
    elif d_protected is not None and d_protected <= 8.0:
        relation = "within_8A_of_user_protected_functional_residue"
    elif d_pdb is not None and d_pdb <= 8.0:
        relation = "within_8A_of_PDB_SITE_annotation"
    else:
        relation = "no_annotated_functional_site_overlap_detected"

    out = dict(patch)
    out.update({
        "pdb_site_overlap_residues": pdb_overlap,
        "user_protected_overlap_residues": protected_overlap,
        "nearest_PDB_SITE_center_distance_A": round(float(d_pdb), 3) if d_pdb is not None else None,
        "nearest_user_protected_center_distance_A": round(float(d_protected), 3) if d_protected is not None else None,
        "functional_site_relation": relation,
    })
    return out


def build_target_interface_profile(
    chemistries: Dict[str, Dict[str, Any]],
    surface_residues: List[Dict[str, Any]],
    all_residues: List[Dict[str, Any]],
    site_annotations: Optional[List[Dict[str, Any]]] = None,
    protected_residue_keys: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Build an unweighted, protein-only interface requirement profile."""
    site_annotations = site_annotations or []
    protected = {str(k).strip() for k in (protected_residue_keys or []) if str(k).strip()}
    pdb_sites = _site_keys(site_annotations)
    coords = _coords_by_key(all_residues)
    n_surface = max(len(surface_residues), 1)

    channels: List[Dict[str, Any]] = []
    for key, group in chemistries.items():
        residues = list(group.get("residues", []) or [])
        top_patches = list(group.get("top_patches", []) or [])
        local_scores = [float(r.get("local_score", 0.0) or 0.0) for r in residues]
        raw5 = [float(p.get("density_5A_raw", 0.0) or 0.0) for p in top_patches]
        raw8 = [float(p.get("density_8A_raw", 0.0) or 0.0) for p in top_patches]
        annotated_patches = [
            _patch_annotation(p, coords, pdb_sites, protected) for p in top_patches[:5]
        ]

        channels.append({
            "key": key,
            "target_surface_property": group.get("surface_group"),
            "label": group.get("label"),
            "protein_side_interpretation": group.get("description"),
            "n_compatible_surface_residues": len(residues),
            "compatible_surface_fraction": round(len(residues) / n_surface, 6),
            "total_accessible_compatibility_raw": round(sum(local_scores), 6),
            "max_local_compatibility_raw": round(max(local_scores), 6) if local_scores else 0.0,
            "max_patch_density_5A_raw": round(max(raw5), 6) if raw5 else 0.0,
            "max_patch_density_8A_raw": round(max(raw8), 6) if raw8 else 0.0,
            "top_patches": annotated_patches,
            "compatible_residues": [
                {
                    "key": r.get("key"),
                    "res_name": r.get("res_name"),
                    "chain": r.get("chain"),
                    "res_seq": r.get("res_seq"),
                    "local_score": r.get("local_score"),
                    "propensity": r.get("propensity"),
                    "multiscale_persistence": r.get("multiscale_persistence"),
                    "mechanism": r.get("mechanism"),
                }
                for r in residues
            ],
        })

    return {
        "basis": "protein_only",
        "material_library_used": False,
        "named_material_recommendation": False,
        "cross_channel_weighted_score": False,
        "interpretation": (
            "Each channel states a surface property that could engage one or more protein patches. "
            "Channels are independent and are not combined into a named-material ranking."
        ),
        "functional_site_annotation_caution": (
            "PDB SITE records are generic structural-site annotations and are not assumed to be catalytic active sites. "
            "User-provided protected residues take precedence when assessing possible functional-site overlap."
        ),
        "n_surface_residues": len(surface_residues),
        "pdb_site_annotations": site_annotations,
        "user_protected_residue_keys": sorted(protected),
        "interface_channels": channels,
    }
