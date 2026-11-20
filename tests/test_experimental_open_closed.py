"""Experimental matched-but-different PDB cross-check.

Adenylate kinase experimental open (4AKE) and closed (1AKE) structures have
identical 214-residue chain A numbering and are used as a functional test of
structure-specific results. No claim is made to predict protein motion.
"""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from interface_engine import analyze_all_maps
from ensemble import analyze_structure_set

FIXTURE = ROOT / "tests" / "fixtures"


class MatchedExperimentalStructuresTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pdb_open = (FIXTURE / "4AKE_A.pdb").read_text("utf-8")
        cls.pdb_closed = (FIXTURE / "1AKE_A.pdb").read_text("utf-8")
        cls.open = analyze_all_maps(pdb_text=cls.pdb_open, chain="A", pH=7.4)
        cls.closed = analyze_all_maps(pdb_text=cls.pdb_closed, chain="A", pH=7.4)

    def test_individual_results_are_structure_dependent(self):
        self.assertEqual(self.open["n_maps"], 11)
        self.assertEqual(self.closed["n_maps"], 11)
        self.assertEqual(len(self.open["residue_profiles"]), 214)
        self.assertEqual(len(self.closed["residue_profiles"]), 214)
        a = {r["key"]: r for r in self.open["residue_profiles"]}
        b = {r["key"]: r for r in self.closed["residue_profiles"]}
        self.assertEqual(a.keys(), b.keys())
        self.assertTrue(all(a[key]["res_name"] == b[key]["res_name"] for key in a))
        delta = sum(abs(a[key]["sidechain_sasa"] - b[key]["sidechain_sasa"]) for key in a)
        self.assertGreater(delta, 50.0)
        self.assertNotEqual(self.open["input"]["source_pdb_sha256"],
                            self.closed["input"]["source_pdb_sha256"])

    def test_individual_multi_pdb_parity_and_real_consensus(self):
        result = analyze_structure_set([
            {"name": "open_4AKE", "pdb_text": self.pdb_open, "chain": "A"},
            {"name": "closed_1AKE", "pdb_text": self.pdb_closed, "chain": "A"},
        ], combined_maps=["hydrophobic", "pi_carbon"])
        self.assertEqual(result["n_structures"], 2)
        self.assertEqual(result["structures"][0]["analysis"]["maps"], self.open["maps"])
        self.assertEqual(result["structures"][1]["analysis"]["maps"], self.closed["maps"])
        self.assertEqual(len(result["comparison"]), 1)
        self.assertFalse(result["comparison_skipped"])
        self.assertIsNone(result["consensus_unavailable_reason"])
        consensus = result["consensus"]
        self.assertEqual(consensus["n_matched_residues"], 214)
        self.assertIsNotNone(consensus["combined_consensus"])
        self.assertGreater(result["comparison"][0]["accessibility"]["n_exposure_status_changes"], 0)
        changed = [
            row["delta_sidechain_sasa"]
            for row in result["comparison"][0]["residue_changes"]
        ]
        self.assertGreater(max(abs(x) for x in changed), 5.0)
        # Inputs are distinct structures; it is expected that at least one
        # selected channel changes its Top-1, Top-3 or Top-5 membership.
        values = [entry["jaccard"] for m in result["comparison"][0]["map_comparisons"]
                  for entry in m["top_k"] if entry["jaccard"] is not None]
        self.assertTrue(any(x < 0.99 for x in values))
        for m in consensus["maps"].values():
            for row in m["residues"]:
                self.assertIn(row["compatible_fraction"], (0.0, 0.5, 1.0))
                for k in (1, 3, 5):
                    self.assertIn(row[f"top{k}_membership_fraction"], (0.0, 0.5, 1.0))


if __name__ == "__main__":
    unittest.main()
