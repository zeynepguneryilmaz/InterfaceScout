"""Regression checks for rejecting C-alpha-only PDBs without rejecting modest missing atoms."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
from prepare import prepare_pdb_text

def build_pdb(n=30, ca_only_positions=frozenset()):
    rows=[]
    serial=1
    for i in range(1,n+1):
        atom_names=("CA",) if i in ca_only_positions else ("N","CA","C","O","CB")
        for j,name in enumerate(atom_names):
            rows.append(
                f"ATOM  {serial:5d} {name:<4s} ALA A{i:4d}    "
                f"{float(i+j/10):8.3f}{1.1:8.3f}{2.2:8.3f}"
                f"{1.:6.2f}{20.:6.2f}          {name[0]:>2s}\n"
            )
            serial+=1
    return "".join(rows)+"END\n"

class TestCAOnlyGuard(unittest.TestCase):
    def test_coarse_ca_only_coordinate_trace_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "C-alpha-only"):
            prepare_pdb_text(build_pdb(ca_only_positions=set(range(1,31))))
    def test_regular_all_atom_protein_is_accepted(self):
        text,report=prepare_pdb_text(build_pdb())
        self.assertIn("ATOM",text)
        self.assertEqual(report["ca_only_residues_detected"],0)
        self.assertEqual(report["modeled_residues_screened"],30)
    def test_modest_atom_missing_is_not_rejected(self):
        text,report=prepare_pdb_text(build_pdb(ca_only_positions={2,4,6,8,10,12}))
        self.assertIn("ATOM",text)
        self.assertEqual(report["ca_only_residues_detected"],6)
        self.assertEqual(report["modeled_residues_screened"],30)

if __name__=="__main__":
    unittest.main()
