"""Regression tests for opt-in PDB2PQR/APBS integration; no binaries required."""
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
import electrostatics as es


class ElectrostaticsTests(unittest.TestCase):
    def test_missing_binary_fails_closed(self):
        with patch.object(es, "availability", return_value={"pdb2pqr": False, "apbs": False}):
            with self.assertRaisesRegex(es.ElectrostaticsError, "missing executable"):
                es.compute("ATOM", [], pH=7.4, ionic_mM=150, temp_K=298)

    def test_pqr_chain_and_charge(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "sample.pqr"
            p.write_text("ATOM 1 CA GLY A 1 1.0 2.0 3.0 0.2 1.7\n"
                         "ATOM 2 CA ALA B 2 1.5 2.5 3.5 -0.2 1.7\n")
            atoms = es._pqr_atoms(p)
            self.assertEqual(atoms[0]["chain"], "A")
            self.assertEqual(atoms[1]["seq"], 2)
            self.assertAlmostEqual(sum(a["charge"] for a in atoms), 0)

    def test_pqr_packed_chain_and_four_digit_residue_index(self):
        # PDB2PQR emits compact A1001 / A1368-style residue tokens when
        # the chain and four-digit number occupy adjacent PQR columns.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "compact.pqr"
            path.write_text(
                "ATOM 1 CA GLY A1001 1.0 2.0 3.0 0.2 1.7\n"
                "ATOM 2 CB ALA A1368 1.5 2.5 3.5 -0.2 1.7\n"
                "ATOM 3 CA SER B2080A 1.4 2.4 3.4 0.1 1.7\n"
            )
            parsed = es._pqr_atoms(path)
            self.assertEqual(
                [(a["chain"], a["seq"], a["icode"]) for a in parsed],
                [("A", 1001, ""), ("A", 1368, ""), ("B", 2080, "A")],
            )

    def test_zero_radius_hydrogen_is_valid_but_heavy_atom_is_not(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "h.pqr"
            path.write_text(
                "ATOM 127 HG1 THR A 7 31.255 40.309 14.948 0.4102 0.0000\n"
            )
            atoms = es._pqr_atoms(path)
            self.assertEqual(atoms[0]["radius"], 0.0)
            path.write_text(
                "ATOM 128 OG1 THR A 7 31.255 40.309 14.948 -0.4102 0.0000\n"
            )
            with self.assertRaisesRegex(es.ElectrostaticsError, "invalid radius"):
                es._pqr_atoms(path)

    def test_dx_interpolation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "sample.dx"
            p.write_text("object 1 class gridpositions counts 2 2 2\n"
                         "origin 0 0 0\n"
                         "delta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\n"
                         "object 2 class gridconnections counts 2 2 2\n"
                         "object 3 class array type double rank 0 items 8 data follows\n"
                         "0 1 2\n3 4 5\n6 7\n")
            grid, origin, spacing = es._read_dx(p)
            self.assertEqual(grid.shape, (2, 2, 2))
            self.assertAlmostEqual(es._interpolate(grid, origin, spacing, (.5, .5, .5)), 3.5)
            self.assertIsNone(es._interpolate(grid, origin, spacing, (-1, 0, 0)))

    def test_generated_apbs_settings(self):
        atoms = [{"xyz": (0., 0., 0.)}, {"xyz": (10., 10., 10.)}]
        content = es._apbs_input(atoms, 150, 310)
        self.assertIn("lpbe", content)
        self.assertIn("ion charge 1 conc 0.15", content)
        self.assertIn("ion charge -1 conc 0.15", content)
        self.assertIn("temp 310.0000", content)
        self.assertIn("write pot dx potential", content)


if __name__ == "__main__":
    unittest.main()
