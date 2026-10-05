"""Frozen numerical settings for the InterfaceScout publication model.

The 200-point Shrake-Rupley sampling level and the 6/9 Å multiscale radius pair
were selected on an adsorption-label-free eight-protein development panel.
The 8 Å coarse patch radius is a predefined geometric contact-region scale and
was not fitted to experimental adsorption labels.
"""
from __future__ import annotations

PUBLIC_VERSION = "1.0-publication"
MODEL_VERSION = "2.4.0-selected-parameters"

SASA_PROBE_A = 1.40
SASA_POINTS = 200
SC_RSA_THRESHOLD = 0.05

MULTISCALE_RADII_A = (6.0, 9.0)
RADIUS_SELECTION_CANDIDATES_A = (
    (5.0, 7.0), (5.0, 8.0), (5.0, 9.0), (5.0, 10.0),
    (6.0, 8.0), (6.0, 9.0), (6.0, 10.0),
    (7.0, 9.0), (7.0, 10.0), (8.0, 10.0),
)

COARSE_PATCH_RADIUS_A = 8.0
SAME_FACE_DOT_THRESHOLD = 0.0
PKA_SENSITIVITY_WINDOW = 1.0

PARAMETER_SELECTION = {
    "purpose": "independent adsorption-label-free core parameter selection",
    "development_pdbs": ["1CRN", "1R69", "1SHG", "2PPN", "4AKE", "1OMP", "1TIM", "5CSC"],
    "sasa_reference_points_per_atom": 1000,
    "selected_sasa_points_per_atom": SASA_POINTS,
    "selected_multiscale_radii_A": list(MULTISCALE_RADII_A),
    "coarse_patch_radius_A": COARSE_PATCH_RADIUS_A,
    "coarse_patch_radius_basis": "predefined geometric contact-region scale; not fitted to adsorption labels",
    "uses_experimental_adsorption_labels": False,
}
