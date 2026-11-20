"""Audit frozen experimental evidence without converting it into a new external validation.

Checks the manifest and CSV are mutually consistent, and preserves the
important distinction between near-hit spatial agreement and direct residue
overlap. This is an integrity test of ALREADY USED data; it is NOT an
independent held-out experimental validation.
"""
from __future__ import annotations

import csv
from pathlib import Path
import json
import re
import statistics
import unittest

ROOT = Path(__file__).resolve().parents[1]
ROOT_DATA = ROOT / "validation"


class IndependentEvidenceAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = json.loads((ROOT_DATA / "benchmark_manifest.json").read_text("utf-8"))
        with (ROOT_DATA / "benchmark_case_results.csv").open(encoding="utf-8-sig", newline="") as file:
            cls.rows = list(csv.DictReader(file))

    def test_provenance_present_and_no_duplicate_condition(self):
        case_ids = [c["id"] for c in self.manifest["cases"]]
        result_ids = [r["case_id"] for r in self.rows]
        self.assertEqual(len(case_ids), len(set(case_ids)))
        self.assertEqual(len(result_ids), len(set(result_ids)))
        self.assertEqual(set(case_ids), set(result_ids))
        self.assertEqual(len(case_ids), 9)
        for c in self.manifest["cases"]:
            self.assertRegex(c["doi"], r"^10\.\d{4,9}/\S+$")
            self.assertTrue(c.get("source_url", "").startswith("https://"))
            self.assertTrue(c.get("evidence", ""))
            self.assertTrue(c.get("pdb_id", ""))
            self.assertTrue(c.get("analysis_set") in ("primary", "secondary"))
            self.assertTrue(c.get("gt_exact") or c.get("gt_ranges"))
            row = next(r for r in self.rows if r["case_id"] == c["id"])
            self.assertEqual(row["analysis_set"], c["analysis_set"])
            self.assertEqual(row["doi"], c["doi"])

    def test_primary_summary_not_misrepresented_as_direct_hit(self):
        primary = [r for r in self.rows if r["analysis_set"] == "primary"]
        secondary = [r for r in self.rows if r["analysis_set"] == "secondary"]
        self.assertEqual((len(primary), len(secondary)), (7, 2))
        def truth(r, col):
            value = r[col].strip().lower()
            self.assertIn(value, ("true", "false"))
            return value == "true"
        self.assertEqual(sum(truth(r, "top3_near_8A_hit") for r in primary), 7)
        self.assertEqual(sum(truth(r, "top5_direct_overlap_hit") for r in primary), 4)
        self.assertEqual(sum(truth(r, "top3_near_8A_hit") for r in secondary), 0)
        self.assertEqual(sum(truth(r, "top5_direct_overlap_hit") for r in secondary), 0)
        median = statistics.median(float(r["top5_near_8A_recall"]) for r in primary)
        self.assertAlmostEqual(median, 0.8, places=10)

    def test_validation_is_not_mislabeled_independent(self):
        policy = self.manifest["policy"]
        self.assertIn("literature", policy["purpose"])
        self.assertIn("before experimental", policy["prediction_ordering"])
        self.assertEqual(policy["model_settings_locked_before_scoring"]["sasa_points_per_atom"], 200)
        self.assertEqual(policy["model_settings_locked_before_scoring"]["multiscale_radii_A"], [6.0, 9.0])


if __name__ == "__main__":
    unittest.main()
