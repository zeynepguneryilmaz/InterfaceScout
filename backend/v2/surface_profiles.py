"""Material surface profiles for InterfaceScout V2-alpha.

Profiles are mechanism mixtures, not fitted adsorption energies.  They provide
transparent weights for patch-level compatibility features.  Values are frozen
for V2-alpha and must not be retuned on the locked external-test set.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


FEATURES = (
    "cationic",
    "anionic",
    "hbond_donor",
    "hbond_acceptor",
    "pi_carbon",
    "hydrophobic",
    "oxide",
    "hydroxyapatite",
    "metal_coord",
    "gold",
    "phosphate",
)


@dataclass(frozen=True)
class SurfaceProfile:
    key: str
    label: str
    description: str
    weights: Dict[str, float]

    def normalized_weights(self) -> Dict[str, float]:
        w = {k: max(float(self.weights.get(k, 0.0)), 0.0) for k in FEATURES}
        s = sum(w.values())
        if s <= 0:
            raise ValueError(f"Surface profile {self.key!r} has no positive weights")
        return {k: v / s for k, v in w.items() if v > 0}


SURFACE_PROFILES: Dict[str, SurfaceProfile] = {
    "silica": SurfaceProfile(
        key="silica",
        label="Silica / silanol-rich surface",
        description="Negative silanolate character with strong hydrogen-bonding/polar contributions.",
        weights={"anionic": 0.50, "hbond_donor": 0.25, "hbond_acceptor": 0.25},
    ),
    "citrate_au": SurfaceProfile(
        key="citrate_au",
        label="Citrate-coated gold",
        description="Protein primarily encounters citrate carboxylates; direct Au contribution is weak in the alpha model.",
        weights={"anionic": 0.65, "hbond_acceptor": 0.20, "gold": 0.15},
    ),
    "paa_fe3o4": SurfaceProfile(
        key="paa_fe3o4",
        label="PAA-coated Fe3O4",
        description="Carboxylate-rich polymer shell with hydrogen-bonding support.",
        weights={"anionic": 0.70, "hbond_acceptor": 0.20, "oxide": 0.10},
    ),
    "fe3o4": SurfaceProfile(
        key="fe3o4",
        label="Iron-oxide surface",
        description="Oxide compatibility with additional transition-metal coordination character.",
        weights={"oxide": 0.60, "metal_coord": 0.30, "hbond_donor": 0.10},
    ),
    "polystyrene": SurfaceProfile(
        key="polystyrene",
        label="Polystyrene",
        description="Hydrophobic/aromatic interface dominated by nonpolar and pi-associated contacts.",
        weights={"hydrophobic": 0.70, "pi_carbon": 0.30},
    ),
    "alumoh": SurfaceProfile(
        key="alumoh",
        label="Al(OH)3 / AlOOH-like positive hydroxylated surface",
        description="Positive hydroxylated alumina-like interface at near-neutral pH.",
        weights={"cationic": 0.65, "hbond_donor": 0.20, "hbond_acceptor": 0.15},
    ),
    "mpa_au": SurfaceProfile(
        key="mpa_au",
        label="MPA-coated gold",
        description="Carboxylate-terminated mercaptopropionic-acid shell with minor Au contribution.",
        weights={"anionic": 0.70, "hbond_acceptor": 0.20, "gold": 0.10},
    ),
    "hydroxyapatite": SurfaceProfile(
        key="hydroxyapatite",
        label="Hydroxyapatite",
        description="Calcium-phosphate interface represented by HAp and phosphate channels.",
        weights={"hydroxyapatite": 0.70, "phosphate": 0.30},
    ),
}


def get_surface_profile(key: str) -> SurfaceProfile:
    k = (key or "").strip().lower()
    if k not in SURFACE_PROFILES:
        available = ", ".join(sorted(SURFACE_PROFILES))
        raise KeyError(f"Unknown surface profile {key!r}. Available: {available}")
    return SURFACE_PROFILES[k]
