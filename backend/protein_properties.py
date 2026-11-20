"""Direct protein-side residue properties for decision-oriented surface analysis.

These are transparent residue classifications, not material-compatibility scores,
nor direct predictions of pKa, reaction kinetics or hydrophobic free energy.
"""
from __future__ import annotations

from typing import Any, Mapping

PROPERTY_CLASSES = {
    "nonpolar": {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"},
    "polar_uncharged": {"ASN", "GLN", "SER", "THR", "TYR", "CYS"},
    "acidic_sidechain": {"ASP", "GLU"},
    "basic_sidechain": {"ARG", "LYS", "HIS"},
    "aromatic": {"PHE", "TYR", "TRP", "HIS"},
    "sulfur_containing": {"CYS", "MET"},
    "cysteine": {"CYS"},
    "lysine": {"LYS"},
    "carboxylate_bearing": {"ASP", "GLU"},
}

def profile_rows(core_result: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Retain experimentally interpretable direct descriptors for ALL residues."""
    rows = []
    for r in core_result.get("all_residues", []):
        rn = str(r["res_name"])
        groups = [name for name, residues in PROPERTY_CLASSES.items() if rn in residues]
        rows.append({
            "key": str(r["key"]),
            "chain": str(r["chain"]),
            "x": float(r["x"]), "y": float(r["y"]), "z": float(r["z"]),
            "res_seq": int(r["res_seq"]),
            "icode": str(r.get("icode", "")),
            "res_name": rn,
            "sidechain_sasa": float(r["sidechain_sasa"]),
            "scrsa_raw": float(r["scrsa_raw"]),
            "scrsa": float(r["scrsa"]),
            "surface_exposed": bool(r["surface_exposed"]),
            "charge_descriptor": float(r["charge_descriptor"]),
            "groups": groups,
        })
    return rows

def summarize_profiles(rows: list[dict[str, Any]]) -> dict[str, Any]:
    exposed = [r for r in rows if r["surface_exposed"]]
    return {
        "n_residues": len(rows),
        "n_surface_residues": len(exposed),
        "exposed_fraction": round(len(exposed) / max(len(rows), 1), 4),
        "total_sidechain_sasa_A2": round(sum(r["sidechain_sasa"] for r in rows), 2),
        "groups": {
            name: {
                "n_exposed_residues": sum(name in r["groups"] for r in exposed),
                "exposed_residues": [r["key"] for r in exposed if name in r["groups"]],
                "exposed_sidechain_sasa_A2": round(sum(
                    r["sidechain_sasa"] for r in exposed if name in r["groups"]), 2),
            }
            for name in PROPERTY_CLASSES
        },
        "notes": [
            "Residue groups are descriptive and can overlap; counts should not be added.",
            "The nonpolar group is a coarse amino-acid class, not a molecular hydropathy measurement.",
            "Acidic/basic and cysteine labels refer to residue identity, not proven ionization "
            "or reaction state under experimental conditions.",
            "Material-compatibility channels are separate from these protein-side properties.",
        ],
    }
