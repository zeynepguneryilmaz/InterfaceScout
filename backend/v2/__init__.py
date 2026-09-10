"""InterfaceScout protein-material biointerface predictor.

The public package exposes one canonical InterfaceScout model. Internal modules
remain separated only for implementation and reproducibility.
"""

from .model_settings import MODEL_VERSION
from .interface_engine import analyze_interface_v2 as analyze_interface

__all__ = ["MODEL_VERSION", "analyze_interface"]
