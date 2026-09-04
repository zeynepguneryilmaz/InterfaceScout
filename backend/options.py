"""Optional analysis controls and material mechanism profiles for InterfaceScout v5.2 development.

The publication core remains fast and deterministic. Expensive or scope-expanding
analyses are opt-in and report separate outputs; they do not modify the canonical
compatibility or 5/8 Å persistence scores.
"""

from __future__ import annotations

from typing import Dict, List

# Lightweight profiles map a material label to pre-declared mechanistic channels.
# Channels remain separate; no weighted sum is produced.
MATERIAL_PROFILES: Dict[str, Dict[str, object]] = {
    "citrate_au": {
        "label": "Citrate-coated gold",
        "channels": ["anionic", "hbond_acceptor"],
        "notes": "Citrate-facing interface; gold/soft-metal map is not assumed exposed through the coating.",
    },
    "bare_au": {
        "label": "Bare gold / nonpolar Au-like interface",
        "channels": ["hydrophobic", "pi_carbon", "gold"],
        "notes": "Dispersion/hydrophobic and sulfur/soft-metal channels are reported separately.",
    },
    "oxidized_silica": {
        "label": "Oxidized / siloxide-rich silica",
        "channels": ["anionic", "hbond_acceptor", "hbond_donor"],
        "notes": "Electrostatic and hydrogen-bond channels are not combined numerically.",
    },
    "positive_silica": {
        "label": "Positively charged silica-like interface",
        "channels": ["cationic", "hbond_donor", "hbond_acceptor"],
        "notes": "Use optional electrostatic steering when whole-protein orientation may be charge-driven.",
    },
    "graphitic_carbon": {
        "label": "Graphitic / aromatic carbon",
        "channels": ["pi_carbon", "hydrophobic"],
        "notes": "Aromatic and nonpolar compatibility are reported as distinct channels.",
    },
    "oxidized_carbon": {
        "label": "Oxidized carbonaceous interface",
        "channels": ["anionic", "hbond_acceptor", "hbond_donor", "pi_carbon", "hydrophobic"],
        "notes": "Mixed chemistry profile; interpret channels independently with material characterization.",
    },
    "metal_oxide": {
        "label": "Metal-oxide interface",
        "channels": ["oxide", "cationic", "anionic", "hbond_donor", "hbond_acceptor"],
        "notes": "Select physically relevant channels using oxide identity, pH and surface charge.",
    },
    "hydroxyapatite": {
        "label": "Hydroxyapatite / calcium-phosphate",
        "channels": ["hydroxyapatite", "phosphate", "cationic"],
        "notes": "Calcium-rich and phosphate-rich channels remain separate.",
    },
}


def material_profile(name: str | None) -> Dict[str, object] | None:
    if not name:
        return None
    return MATERIAL_PROFILES.get(str(name).strip().lower())


def available_material_profiles() -> List[Dict[str, object]]:
    return [{"key": k, **v} for k, v in MATERIAL_PROFILES.items()]
