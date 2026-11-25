"""Regression tests for opt-in PDB2PQR/APBS integration; no binaries required."""
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
import electrostatics as es
import surface_potential as sp
import numpy as np


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

    def test_sas_surface_sampling_includes_glycine_and_preserves_units(self):
        # Constant genuine-format grid: area-weighted shell means must
        # reproduce the constant for GLY as well as another residue.
        atoms = []
        for i in range(12):
            atoms.append({
                "name": "CA", "radius": 1.5, "xyz": (8.0 * i + 15.0, 20.0, 20.0),
                "chain": "A", "seq": 1 if i < 6 else 2, "icode": "",
            })
        grid = np.full((130, 60, 60), -2.75, dtype=float)
        values, meta = sp.sample_residue_sas_potential(
            atoms, grid, np.zeros(3), np.ones(3), {"A:1:", "A:2:"}
        )
        self.assertEqual(set(values), {"A:1:", "A:2:"})
        self.assertAlmostEqual(values["A:1:"]["surface_potential_kT_per_e"], -2.75)
        self.assertAlmostEqual(values["A:2:"]["surface_potential_kT_per_e"], -2.75)
        self.assertEqual(meta["n_exposed_residues_without_shell_samples"], 0)
        self.assertEqual(meta["directions_per_heavy_atom"], 42)
        self.assertEqual(meta["surface_type"], "approximate SAS shell, not molecular SES")

    def test_sas_grid_out_of_bounds_fails_closed(self):
        atoms = [{
            "name": "CA", "radius": 1.5, "xyz": (i * 8., 0., 0.),
            "chain": "A", "seq": i, "icode": "",
        } for i in range(10)]
        with self.assertRaisesRegex(ValueError, "grid excludes"):
            sp.sample_residue_sas_potential(
                atoms, np.zeros((3, 3, 3)), np.zeros(3), np.ones(3),
                {f"A:{i}:" for i in range(10)}
            )

    def test_strict_orphan_n_only_apbs_cleanup_preserves_original(self):
        # Two isolated N-only serines; no CA/C/O atom is synthesized.
        original = (
            "ATOM      1  N   SER A  11      10.000  10.000  10.000  1.00 20.00           N\n"
            "ATOM      2  N   SER E  11      20.000  10.000  10.000  1.00 20.00           N\n"
            "ATOM      3  N   ALA A  12      11.000  10.000  10.000  1.00 20.00           N\n"
            "ATOM      4  CA  ALA A  12      12.000  10.000  10.000  1.00 20.00           C\n"
            "ATOM      5  C   ALA A  12      13.000  10.000  10.000  1.00 20.00           C\n"
            "ATOM      6  O   ALA A  12      14.000  10.000  10.000  1.00 20.00           O\n"
        )
        valid = {"A:11:": "SER", "E:11:": "SER", "A:12:": "ALA"}
        # 2/3 exceeds the predeclared 1% limit: fail closed.
        with self.assertRaisesRegex(es.ElectrostaticsError, "cannot safely recover"):
            es._pdb_orphan_n_only_residues(original, valid)
        # A realistic larger protein stays below the conservative 1% cap.
        valid.update({f"A:{i}:": "ALA" for i in range(20, 420)})
        cleaned, omissions = es._pdb_orphan_n_only_residues(original, valid)
        self.assertEqual({o["key"] for o in omissions}, {"A:11:", "E:11:"})
        self.assertEqual(len(omissions), 2)
        self.assertEqual(sum(x.startswith("ATOM ") for x in cleaned.splitlines()), 4)
        self.assertIn("SER A  11", original)
        self.assertNotIn("SER A  11", cleaned)
        self.assertIn("ALA A  12", cleaned)
        self.assertEqual(cleaned.count("ALA A  12"), 4)

    def test_partial_backbone_is_never_silently_deleted(self):
        src = (
            "ATOM      1  N   SER A  11      10.000  10.000  10.000  1.00 20.00           N\n"
            "ATOM      2  CA  SER A  11      11.000  10.000  10.000  1.00 20.00           C\n"
            "ATOM      3  C   SER A  11      12.000  10.000  10.000  1.00 20.00           C\n"
        )
        cleaned, dropped = es._pdb_orphan_n_only_residues(src, {"A:11:": "SER"})
        self.assertEqual(cleaned, src)
        self.assertEqual(dropped, [])

    def test_generated_apbs_settings(self):
        atoms = [{"xyz": (0., 0., 0.)}, {"xyz": (10., 10., 10.)}]
        content = es._apbs_input(atoms, 150, 310)
        self.assertIn("lpbe", content)
        self.assertIn("    npbe\n", es._apbs_input(atoms, 150, 310, solver="npbe"))
        self.assertNotIn("    lpbe\n", es._apbs_input(atoms, 150, 310, solver="npbe"))
        self.assertIn("ion charge 1 conc 0.15", content)
        self.assertIn("ion charge -1 conc 0.15", content)
        self.assertIn("temp 310.0000", content)
        self.assertIn("write pot dx potential", content)


if __name__ == "__main__":
    unittest.main()
