"""Canonical InterfaceScout implementation modules.

The ``v2`` directory name is retained only for provenance of the publication
implementation; users interact with one application, InterfaceScout.
"""

from .model_settings import MODEL_VERSION
from .interface_engine import analyze_all_maps

__all__ = ["MODEL_VERSION", "analyze_all_maps"]
