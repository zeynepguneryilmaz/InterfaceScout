"""Categorical material-to-chemistry mapping for InterfaceScout V2.

Unlike the retired V2-alpha surface profiles, these mappings contain no fitted
or hand-tuned numeric weights. They only select the primary frozen V1 chemistry
channel used to seed coarse interface patches.
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
    "silica": SurfaceMode("silica", "Silica / deprotonated silanol-rich", "anionic", "Primary electrostatic compatibility with a negatively charged silica surface under neutral/basic conditions."),
    "citrate_au": SurfaceMode("citrate_au", "Citrate-coated gold", "anionic", "Primary protein-side compatibility with the negatively charged citrate corona."),
    "paa_fe3o4": SurfaceMode("paa_fe3o4", "PAA-coated iron oxide", "anionic", "Primary compatibility with exposed anionic carboxylate groups of PAA."),
    "mpa_au": SurfaceMode("mpa_au", "MPA-coated gold", "anionic", "Primary compatibility with the negatively charged mercaptopropionate coating."),
    "polystyrene": SurfaceMode("polystyrene", "Polystyrene", "hydrophobic", "Primary compatibility with a nonpolar polymer surface."),
    "alumoh": SurfaceMode("alumoh", "Aluminum hydroxide / alum", "cationic", "Primary compatibility with a positively charged hydroxide-rich surface in the benchmark condition."),
    "fe3o4": SurfaceMode("fe3o4", "Iron oxide", "oxide", "Primary compatibility with a metal-oxide interface."),
    "hydroxyapatite": SurfaceMode("hydroxyapatite", "Hydroxyapatite", "hydroxyapatite", "Primary compatibility with calcium-phosphate-rich hydroxyapatite."),
    "calcium_fluoride": SurfaceMode(
        "calcium_fluoride",
        "Calcium fluoride nanoparticle",
        "cationic",
        "Primary protein-side compatibility with exposed calcium-rich sites; used categorically for CaF2 systems where experiments identify acidic protein regions contacting surface Ca2+ sites. No numeric material weight is assigned.",
    ),
    "hydrophobic": SurfaceMode("hydrophobic", "Generic hydrophobic surface", "hydrophobic", "Generic nonpolar surface mode."),
    "anionic": SurfaceMode("anionic", "Generic anionic surface", "anionic", "Generic negatively charged surface mode."),
    "cationic": SurfaceMode("cationic", "Generic cationic surface", "cationic", "Generic positively charged surface mode."),
}


def get_surface_mode(key: str) -> SurfaceMode:
    k = (key or "").strip().lower()
    if k not in SURFACE_MODES:
        raise ValueError(f"Unknown surface mode {key!r}. Available: {', '.join(sorted(SURFACE_MODES))}")
    return SURFACE_MODES[k]
