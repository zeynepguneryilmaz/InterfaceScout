"""Optional interpreter bootstrap for canonical-only validation runs.

Python imports sitecustomize automatically when it is on PYTHONPATH. The switch
below is inactive in normal InterfaceScout use. CI validation enables it so
PDB2PQR/APBS/DSSP annotation layers, which are excluded from prediction and
ranking, do not add runtime or availability dependence to canonical tests.
"""
from __future__ import annotations

import os

if os.environ.get("INTERFACESCOUT_VALIDATION_CANONICAL_ONLY") == "1":
    try:
        import main as _main
        _main.PDB2PQR = None
        _main.APBS = None
        _main.MKDSSP = None
    except Exception:
        pass
