"""InterfaceScout V2-alpha.

This package is intentionally isolated from the publication-frozen V1 backend.
V2-alpha predicts and ranks protein surface patches using material surface profiles.
It does not yet perform rigid-body orientation search or adsorption free-energy prediction.
"""

V2_VERSION = "2.0.0-alpha.1"

from .engine import analyze_v2

__all__ = ["V2_VERSION", "analyze_v2"]
