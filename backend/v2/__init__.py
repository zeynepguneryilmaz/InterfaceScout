"""InterfaceScout V2 coarse protein-material biointerface predictor.

V2 predicts plausible coarse protein surface regions from material-interaction
chemistry, solvent exposure, pH-dependent state availability, multiscale spatial
aggregation, and coarse surface geometry. It does not predict adsorption free
energy, adsorption amount, or a unique atomistic orientation.
"""

from .model_settings import MODEL_VERSION as V2_VERSION
from .interface_engine import analyze_interface_v2

analyze_v2 = analyze_interface_v2

__all__ = ["V2_VERSION", "analyze_interface_v2", "analyze_v2"]
