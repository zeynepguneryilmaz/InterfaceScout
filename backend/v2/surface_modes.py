"""Validation-only material-to-chemistry mapping.

The public InterfaceScout application does not ask for or expose material names.
This module exists only so archived and reproducible experimental-validation
manifests can map the experimentally studied interface to the appropriate
canonical surface-chemistry channel. No numerical material weights are used.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SurfaceMode:
    key: str
    label: str
    chemistry: str
    description: str


SURFACE_MODES = {
    "silica": SurfaceMode("silica", "Silica / deprotonated silanol-rich", "anionic", "Validation mapping to the anionic surface-chemistry channel."),
    "citrate_au": SurfaceMode("citrate_au", "Citrate-coated gold", "anionic", "Validation mapping to the anionic coating chemistry."),
    "paa_fe3o4": SurfaceMode("paa_fe3o4", "PAA-coated iron oxide", "anionic", "Validation mapping to exposed anionic PAA chemistry."),
    "mpa_au": SurfaceMode("mpa_au", "MPA-coated gold", "anionic", "Validation mapping to the anionic mercaptopropionate coating chemistry."),
    "polystyrene": SurfaceMode("polystyrene", "Polystyrene", "hydrophobic", "Validation mapping to the hydrophobic surface-chemistry channel."),
    "alumoh": SurfaceMode("alumoh", "Aluminum hydroxide / alum", "cationic", "Validation mapping to the cationic surface-chemistry channel for the reported condition."),
    "fe3o4": SurfaceMode("fe3o4", "Iron oxide", "oxide", "Validation mapping to the metal-oxide surface-chemistry channel."),
    "hydroxyapatite": SurfaceMode("hydroxyapatite", "Hydroxyapatite", "hydroxyapatite", "Validation mapping to calcium/phosphate charged-site chemistry."),
    "calcium_fluoride": SurfaceMode("calcium_fluoride", "Calcium fluoride nanoparticle", "cationic", "Validation mapping to the cationic channel for experiments implicating acidic protein regions contacting exposed calcium-rich sites."),
    "hydrophobic": SurfaceMode("hydrophobic", "Generic hydrophobic surface", "hydrophobic", "Generic validation mapping."),
    "anionic": SurfaceMode("anionic", "Generic anionic surface", "anionic", "Generic validation mapping."),
    "cationic": SurfaceMode("cationic", "Generic cationic surface", "cationic", "Generic validation mapping."),
}


def get_surface_mode(key: str) -> SurfaceMode:
    k = (key or "").strip().lower()
    if k not in SURFACE_MODES:
        raise ValueError(f"Unknown validation surface mapping {key!r}. Available: {', '.join(sorted(SURFACE_MODES))}")
    return SURFACE_MODES[k]
