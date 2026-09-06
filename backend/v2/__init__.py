"""InterfaceScout V2 coarse protein-material biointerface predictor.

V2 is isolated from the publication-frozen V1 backend. It predicts plausible
coarse protein surface regions using frozen V1 chemistry/accessibility plus
local geometry, native-state GNM dynamics, and coarse orientation descriptors.
It does not predict adsorption free energy, adsorption amount, or a unique
atomistic orientation.
"""

V2_VERSION = "2.2.1-coarse-prototype"

from .interface_engine import analyze_interface_v2

# Public alias retained for callers that expect an ``analyze_v2`` entry point.
analyze_v2 = analyze_interface_v2

__all__ = ["V2_VERSION", "analyze_interface_v2", "analyze_v2"]
