"""Publication-frozen InterfaceScout settings.

The canonical chemistry definitions now live directly in ``backend/main.py``.
This helper remains for reproducibility scripts written before the repository was
flattened to one public model; calling it simply enforces the selected numerical
settings and is otherwise idempotent.
"""
from __future__ import annotations

from .model_settings import SASA_POINTS, SASA_PROBE_A, SC_RSA_THRESHOLD, MULTISCALE_RADII_A


def apply_publication_chemistry(core_module) -> None:
    core_module.SASA_POINTS = SASA_POINTS
    core_module.SASA_PROBE_A = SASA_PROBE_A
    core_module.SC_RSA_THRESHOLD = SC_RSA_THRESHOLD
    core_module.PATCH_RADII_A = tuple(float(x) for x in MULTISCALE_RADII_A)
