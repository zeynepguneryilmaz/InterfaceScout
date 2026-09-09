"""Chemistry corrections and selected model settings applied by InterfaceScout V2.

The legacy /analyze_surface endpoint is retained for backward compatibility.
V2 applies the publication chemistry corrections and the independently selected
numerical/spatial settings before generating chemistry maps. Experimental
protein-material benchmark labels are not used in these settings.
"""
from __future__ import annotations

import numpy as np

from .model_settings import (
    SASA_POINTS,
    SASA_PROBE_A,
    SC_RSA_THRESHOLD,
    MULTISCALE_RADII_A,
    COARSE_PATCH_RADIUS_A,
)


def _replace_state(entry, state):
    """Preserve mechanism metadata while replacing the required ionization state."""
    _, mechanism, _ = entry
    return (0.0, mechanism, state)


def _neutralize_ebase(channel):
    """Remove numerical energy magnitude from V2 while preserving mechanism/state."""
    for group in ("favorable", "repulsive"):
        rows = channel.get(group, {})
        for residue, entry in list(rows.items()):
            _, mechanism, state = entry
            rows[residue] = (0.0, mechanism, state)


def _selected_chemistry_map(v1, surface, distances, key, pH):
    """V2 chemistry map using the independently selected 6/9 A aggregation pair."""
    meta = v1.CHEMISTRIES[key]
    fav_defs = meta["favorable"]
    rep_defs = meta["repulsive"]
    n = len(surface)
    local = np.zeros(n, dtype=float)
    rep_local = np.zeros(n, dtype=float)

    for i, r in enumerate(surface):
        rn = r["res_name"]
        if rn in fav_defs:
            _, _, state = fav_defs[rn]
            local[i] = r["scrsa"] * v1.state_availability(rn, state, pH)
        if rn in rep_defs:
            _, _, state = rep_defs[rn]
            rep_local[i] = r["scrsa"] * v1.state_availability(rn, state, pH)

    prop = v1.normalize_to_100(local)
    rep_prop = v1.normalize_to_100(rep_local)

    inner, outer = MULTISCALE_RADII_A
    dens = {}
    dens_norm = {}
    for radius in MULTISCALE_RADII_A:
        if n:
            values = (distances <= radius).astype(float) @ local
        else:
            values = np.zeros(0, dtype=float)
        dens[radius] = values
        maximum = float(np.max(values)) if values.size else 0.0
        dens_norm[radius] = values / maximum if maximum > 0 else np.zeros_like(values)

    if n:
        persistence = 100.0 * np.minimum(dens_norm[inner], dens_norm[outer])
        geomean = 100.0 * np.sqrt(dens_norm[inner] * dens_norm[outer])
    else:
        persistence = np.zeros(0, dtype=float)
        geomean = np.zeros(0, dtype=float)

    inner_tag = f"{int(inner)}A"
    outer_tag = f"{int(outer)}A"
    coarse_tag = f"{int(COARSE_PATCH_RADIUS_A)}A"

    patch_centers = []
    for i, r in enumerate(surface):
        if persistence[i] <= 0:
            continue
        coarse_members = [
            surface[j]["key"] for j in range(n)
            if distances[i, j] <= COARSE_PATCH_RADIUS_A and local[j] > 0
        ]
        row = {
            "center_key": r["key"], "res_name": r["res_name"], "res_seq": r["res_seq"],
            "icode": r["icode"], "chain": r["chain"],
            f"density_{inner_tag}_raw": round(float(dens[inner][i]), 6),
            f"density_{outer_tag}_raw": round(float(dens[outer][i]), 6),
            f"density_{inner_tag}_norm": round(float(dens_norm[inner][i] * 100.0), 3),
            f"density_{outer_tag}_norm": round(float(dens_norm[outer][i] * 100.0), 3),
            "multiscale_persistence": round(float(persistence[i]), 3),
            "multiscale_geomean": round(float(geomean[i]), 3),
            f"compatible_members_{coarse_tag}": coarse_members,
        }
        patch_centers.append(row)
    patch_centers.sort(
        key=lambda x: (-x["multiscale_persistence"], -x["multiscale_geomean"], x["center_key"])
    )

    members = []
    repulsive = []
    for i, r in enumerate(surface):
        rn = r["res_name"]
        if rn in fav_defs:
            ebase, mechanism, state = fav_defs[rn]
            fstate = v1.state_availability(rn, state, pH)
            aux_state = v1.charged_fraction(rn, pH) if state == "auxiliary_only" else None
            members.append({
                "key": r["key"], "res_name": rn, "res_seq": r["res_seq"], "icode": r["icode"], "chain": r["chain"],
                "sidechain_sasa": r["sidechain_sasa"], "scrsa_raw": r["scrsa_raw"], "scrsa": r["scrsa"],
                "state_requirement": state or "none", "state_availability": round(float(fstate), 5),
                "auxiliary_charged_fraction": round(float(aux_state), 5) if aux_state is not None else None,
                "local_score": round(float(local[i]), 6), "propensity": round(float(prop[i]), 3),
                f"patch_density_{inner_tag}": round(float(dens_norm[inner][i] * 100.0), 3),
                f"patch_density_{outer_tag}": round(float(dens_norm[outer][i] * 100.0), 3),
                "multiscale_persistence": round(float(persistence[i]), 3),
                "multiscale_geomean": round(float(geomean[i]), 3),
                "phi": r["phi"],
                "electrostatic_relation": v1.electrostatic_relation(r["phi"], meta.get("expected_phi_sign")),
                "ebase_metadata_kcal_mol": ebase, "mechanism": mechanism,
                "ionization_sensitive": r["ionization_sensitive"],
                "ss": r["ss"], "bfactor": r["bfactor"], "n_neighbors_8A": r["n_neighbors_8A"],
            })
        if rn in rep_defs:
            ebase, mechanism, state = rep_defs[rn]
            fstate = v1.state_availability(rn, state, pH)
            repulsive.append({
                "key": r["key"], "res_name": rn, "res_seq": r["res_seq"], "icode": r["icode"], "chain": r["chain"],
                "scrsa": r["scrsa"], "state_requirement": state or "none",
                "state_availability": round(float(fstate), 5),
                "repulsion_score": round(float(rep_local[i]), 6),
                "repulsion_propensity": round(float(rep_prop[i]), 3),
                "phi": r["phi"], "ebase_metadata_kcal_mol": ebase, "mechanism": mechanism,
            })

    members.sort(key=lambda x: (-x["propensity"], -x["multiscale_persistence"], x["key"]))
    repulsive.sort(key=lambda x: (-x["repulsion_propensity"], x["key"]))

    return {
        "kind": "compatibility", "label": meta["label"], "surface_group": meta["surface_group"],
        "description": meta["description"], "n_favorable_residues": len(members),
        "n_repulsive_residues": len(repulsive), "top_patch": patch_centers[0] if patch_centers else None,
        "top_patches": patch_centers[:10], "patch_centers": patch_centers,
        "residues": members, "repulsive_residues": repulsive,
        "notes": {
            "primary_score": "L = scRSA x binary membership x required-state availability",
            "normalization": "within-map maximum = 100",
            "patch": f"multiscale persistence = 100 x min(normalized {int(inner)} A density, normalized {int(outer)} A density)",
            "radius_selection": "6/9 A selected on an independent label-free development panel before external validation",
            "ebase": "metadata only; excluded from ranking",
            "apbs": "auxiliary descriptor only; excluded from ranking",
        },
    }


def apply_publication_chemistry(v1_module) -> None:
    chem = v1_module.CHEMISTRIES

    donor_surface = chem["hbond_donor"]["favorable"]
    donor_surface.pop("LYS", None)
    donor_surface.pop("ARG", None)
    chem["hbond_donor"]["description"] = (
        "Protein-side compatibility with surfaces capable of donating hydrogen bonds; "
        "favorable residues are restricted to acceptor-capable side chains."
    )

    acceptor_surface = chem["hbond_acceptor"]["favorable"]
    if "LYS" in acceptor_surface:
        acceptor_surface["LYS"] = _replace_state(acceptor_surface["LYS"], "protonated")
    if "ARG" in acceptor_surface:
        acceptor_surface["ARG"] = _replace_state(acceptor_surface["ARG"], "protonated")
    chem["hbond_acceptor"]["description"] = (
        "Protein-side compatibility with surfaces capable of accepting hydrogen bonds; "
        "basic side-chain donors are conditioned on their protonated state."
    )

    oxide = chem["oxide"]["favorable"]
    for residue in ("SER", "THR", "TYR"):
        oxide.pop(residue, None)
    chem["oxide"]["description"] = (
        "Primary protein-side compatibility with exposed oxide metal sites through "
        "deprotonated Asp/Glu carboxylate coordination."
    )

    hap = chem["hydroxyapatite"]["favorable"]
    for residue in ("SER", "THR", "TYR", "PHE", "TRP"):
        hap.pop(residue, None)
    hap["ARG"] = (0.0, "charged-site electrostatic complementarity", "protonated")
    hap["LYS"] = (0.0, "charged-site electrostatic complementarity", "protonated")
    hap["HIS"] = (0.0, "charged-site electrostatic complementarity", "protonated")
    chem["hydroxyapatite"]["description"] = (
        "Primary charged-site compatibility with hydroxyapatite/calcium-phosphate "
        "interfaces: deprotonated Asp/Glu carboxylates and protonated basic residues."
    )

    for channel in chem.values():
        _neutralize_ebase(channel)

    # Apply the independently selected settings to the V2 path. These assignments
    # occur before v1.analyze() is called by interface_engine.
    v1_module.SASA_POINTS = SASA_POINTS
    v1_module.SASA_PROBE_A = SASA_PROBE_A
    v1_module.SC_RSA_THRESHOLD = SC_RSA_THRESHOLD
    v1_module.PATCH_RADII_A = MULTISCALE_RADII_A
    v1_module.chemistry_map = lambda surface, distances, key, pH: _selected_chemistry_map(
        v1_module, surface, distances, key, pH
    )
