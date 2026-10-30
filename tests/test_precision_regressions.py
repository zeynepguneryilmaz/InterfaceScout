"""Regression tests for precision loss before patch construction/ranking."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
import core
from coarse_patch import build_coarse_patches


def surface_row(name="HIS", scrsa=0.1234567890123):
    return dict(key="A:1:", res_name=name, res_seq=1, icode="", chain="A",
                x=1.0, y=0.0, z=0.0, sidechain_sasa=12.34567890123,
                scrsa_raw=scrsa, scrsa=scrsa, ionization_sensitive=False, bfactor=20.0)


class PrecisionRegressionTests(unittest.TestCase):
    def test_extreme_ph_positive_scores_still_produce_patches(self):
        for name, chemistry, ph in [("HIS", "anionic", 14.0), ("ASP", "oxide", 0.0)]:
            with self.subTest(name=name, pH=ph):
                row = surface_row(name, 0.051234567891)
                channel = core.chemistry_map([row], np.zeros((1, 1)), chemistry, ph)
                nodes = [dict(key=row["key"], chain="A", res_seq=1, icode="",
                              res_name=name, coord=np.array([1.0, 0.0, 0.0]))]
                patches = build_coarse_patches(
                    v1_result={"surface_residues": [row], "chemistries": {chemistry: channel}},
                    chemistry=chemistry, ca_nodes=nodes)
                self.assertEqual(len(patches), 1)
                self.assertEqual(patches[0]["seed_members"], ["A:1:"])
                self.assertEqual(patches[0]["chemistry_support"], 1.0)

    def test_local_score_and_state_keep_full_precision(self):
        row = surface_row()
        result = core.chemistry_map([row], np.zeros((1, 1)), "anionic", 14.0)
        member = result["residues"][0]
        availability = 1.0 / (1.0 + 10.0 ** (14.0 - 6.5))
        self.assertEqual(member["local_score"], row["scrsa"] * availability)
        self.assertEqual(member["state_availability"], availability)

    def test_density_and_persistence_keep_full_precision(self):
        rows = [surface_row(scrsa=0.1234567890123), surface_row(scrsa=0.9876543210987)]
        rows[1].update(key="A:2:", res_seq=2)
        result = core.chemistry_map(rows, np.array([[0., 20.], [20., 0.]]), "anionic", 6.5)
        center = next(c for c in result["patch_centers"] if c["center_key"] == "A:1:")
        self.assertEqual(center["density_6A_raw"], rows[0]["scrsa"] * 0.5)
        expected = 100.0 * (rows[0]["scrsa"] / rows[1]["scrsa"])
        self.assertAlmostEqual(center["multiscale_persistence"], expected, places=13)

    def test_accessibility_is_not_rounded_before_scoring(self):
        pdb = ("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
               "ATOM      2  CB  ALA A   1       1.500   0.000   0.000  1.00 20.00           C  \nEND\n")
        asa = 12.345678901234
        def compute(structure, level):
            for atom in structure.get_atoms():
                atom.sasa = asa if atom.name == "CB" else 0.0
        with patch.object(core, "ShrakeRupley") as sr:
            sr.return_value.compute.side_effect = compute
            row = core._surface_residues(pdb, 7.4)[0]
        self.assertEqual(row["sidechain_sasa"], asa)
        self.assertEqual(row["scrsa"], asa / core.SIDECHAIN_REF_ASA["ALA"])
        self.assertEqual(row["scrsa_raw"], row["scrsa"])


if __name__ == "__main__":
    unittest.main()
