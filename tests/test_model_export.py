"""The exported coordinates must be the structure that was actually analyzed."""
from __future__ import annotations
from io import StringIO
from pathlib import Path
import sys
import unittest

from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
from interface_engine import analyze_all_maps


class ModelExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = (ROOT / "examples/1CRN/InterfaceScout_1CRN_anionic_Bfactor.pdb").read_text()
        cls.atoms = "\n".join(l for l in source.splitlines() if l.startswith("ATOM  ")) + "\n"

    def test_first_model_identity_and_export_scope(self):
        # Serial 7 is deliberately not ordinal 1. The same residue keys recur in 42.
        raw = "MODEL        7\n" + self.atoms + "ENDMDL\nMODEL       42\n" + self.atoms + "ENDMDL\nEND\n"
        result = analyze_all_maps(pdb_text=raw, chain="A")
        report = result["structure_preparation"]
        self.assertEqual(report.get("model_policy"), "first_model_only")
        self.assertEqual(report.get("selected_model_index"), 0)
        self.assertEqual(report.get("selected_model_serial"), 7)
        exported = result["prepared_pdb"]
        models = list(PDBParser(QUIET=True).get_structure("export", StringIO(exported)))
        self.assertEqual(len(models), 1)
        self.assertEqual(models[0].serial_num, 7)
        reference = analyze_all_maps(pdb_text=self.atoms + "END\n", chain="A")
        self.assertEqual(result["maps"], reference["maps"])

    def test_implicit_single_model_export(self):
        result = analyze_all_maps(pdb_text=self.atoms + "END\n", chain="A")
        self.assertIsNone(result["structure_preparation"].get("selected_model_serial"))
        self.assertIn("prepared_pdb", result.keys())
        self.assertNotIn("MODEL ", result["prepared_pdb"])

    def test_explicit_single_model_retains_serial(self):
        result = analyze_all_maps(pdb_text="MODEL        7\n" + self.atoms + "ENDMDL\nEND\n")
        self.assertEqual(result["structure_preparation"].get("selected_model_serial"), 7)
        self.assertIn("MODEL      7", result["prepared_pdb"])

    def test_export_uses_only_prepared_selected_chains(self):
        other_chain = "\n".join(line[:21] + "B" + line[22:]
                                for line in self.atoms.splitlines()) + "\n"
        result = analyze_all_maps(pdb_text=self.atoms + "TER\n" + other_chain + "END\n", chain="A")
        structure = PDBParser(QUIET=True).get_structure("export", StringIO(result["prepared_pdb"]))
        self.assertEqual([c.id for c in structure.get_chains()], ["A"])


if __name__ == "__main__":
    unittest.main()
