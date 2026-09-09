"""InterfaceScout V2 model settings fixed before external experimental validation."""
from __future__ import annotations

MODEL_VERSION = "2.4.0-selected-parameters"

# Surface exposure
SASA_PROBE_A = 1.40
SASA_POINTS = 200
SC_RSA_THRESHOLD = 0.05

# Multiscale chemistry aggregation. The pair was selected on an independent,
# adsorption-label-free development panel by equal consensus comparison across
# predefined local C-alpha radius pairs.
MULTISCALE_RADII_A = (6.0, 9.0)
RADIUS_SELECTION_CANDIDATES_A = (
    (5.0, 7.0), (5.0, 8.0), (5.0, 9.0), (5.0, 10.0),
    (6.0, 8.0), (6.0, 9.0), (6.0, 10.0),
    (7.0, 9.0), (7.0, 10.0), (8.0, 10.0),
)

# Coarse patch construction remains a separate structural-neighborhood rule.
# It is intentionally not equated with the outer multiscale aggregation radius.
COARSE_PATCH_RADIUS_A = 8.0
SAME_FACE_DOT_THRESHOLD = 0.0

PARAMETER_SELECTION = {
    "development_pdbs": ["1CRN", "1R69", "1SHG", "2PPN", "4AKE", "1OMP", "1TIM", "5CSC"],
    "sasa_reference_points_per_atom": 1000,
    "selected_sasa_points_per_atom": SASA_POINTS,
    "selected_multiscale_radii_A": list(MULTISCALE_RADII_A),
    "selection_workflow_run": 34329264963,
    "selection_commit": "ba0e63b7b82a222856ff75a6dd39409b89b664dd",
    "uses_experimental_adsorption_labels": False,
}
