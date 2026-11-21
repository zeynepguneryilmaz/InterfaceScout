from __future__ import annotations
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import core
import coarse_patch
import model_settings as settings
from interface_engine import analyze_all_maps

class PublicationSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pdb_text = (ROOT / "examples" / "1CRN" / "InterfaceScout_1CRN_anionic_Bfactor.pdb").read_text(encoding="utf-8")
        cls.result = analyze_all_maps(
            pdb_text=pdb_text,
            chain="A",
            pH=7.4,
            ionic_mM=150.0,
            temp_K=298.0,
        )

    def test_single_source_settings(self):
        self.assertEqual(tuple(core.PATCH_RADII_A), settings.MULTISCALE_RADII_A)
        self.assertEqual(core.SASA_POINTS, settings.SASA_POINTS)
        self.assertEqual(core.SASA_PROBE_A, settings.SASA_PROBE_A)
        self.assertEqual(core.SC_RSA_THRESHOLD, settings.SC_RSA_THRESHOLD)
        self.assertEqual(core.COARSE_PATCH_RADIUS_A, settings.COARSE_PATCH_RADIUS_A)
        self.assertEqual(coarse_patch.PATCH_SCALE_A, settings.COARSE_PATCH_RADIUS_A)

    def test_frozen_publication_parameters(self):
        self.assertEqual(settings.SASA_POINTS, 200)
        self.assertEqual(settings.MULTISCALE_RADII_A, (6.0, 9.0))
        self.assertEqual(settings.COARSE_PATCH_RADIUS_A, 8.0)
        self.assertEqual(settings.SASA_PROBE_A, 1.40)
        self.assertEqual(settings.SC_RSA_THRESHOLD, 0.05)

    def test_all_maps_smoke(self):
        self.assertEqual(self.result["n_maps"], 11)
        self.assertEqual(self.result["method"]["multiscale_radii_A"], [6.0, 9.0])
        self.assertEqual(self.result["method"]["patch_radius_A"], 8.0)
        self.assertTrue(self.result["maps"]["anionic"]["patches"])

    def test_scope_is_explicit(self):
        scope = self.result["scope"]
        for key in (
            "predicts_absolute_adsorption_free_energy",
            "predicts_adsorption_amount",
            "predicts_unique_orientation",
            "models_adsorption_induced_unfolding",
            "benchmark_fitted_weights",
        ):
            self.assertFalse(scope[key])

if __name__ == "__main__":
    unittest.main()
