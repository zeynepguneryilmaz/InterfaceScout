"""End-to-end cross-check against a real bundled PDB, not synthetic map rows.

The 3Dmol screenshot itself is browser-dependent. These tests validate the Python
contract with a valid PNG supplied by the client, while the frontend syntax
and browser capture API are checked separately.
"""
from __future__ import annotations

import base64
import csv
from io import BytesIO, StringIO
import json
from pathlib import Path
import sys
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from PIL import Image as PILImage

from interface_engine import analyze_all_maps
from ensemble import analyze_structure_set
from export_bundle import MAP_NAMES, build_export_zip
from app import StructureSetRequest, ExportRequest


def png_uri() -> str:
    im = PILImage.new("RGB", (300, 170), "#dbeafe")
    for x in range(30, 260):
        for y in range(20, 150):
            if ((x - 145) / 104)**2 + ((y - 85) / 58)**2 < 1:
                im.putpixel((x, y), (34 + x // 6, 75, 140))
    out = BytesIO()
    im.save(out, format="PNG")
    return "data:image/png;base64," + base64.b64encode(out.getvalue()).decode("ascii")


class RealPdbIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pdb = (ROOT / "examples" / "1CRN" /
                   "InterfaceScout_1CRN_anionic_Bfactor.pdb").read_text("utf-8")
        cls.one = analyze_all_maps(pdb_text=cls.pdb, chain="A", pH=7.4,
                                   ionic_mM=150.0, temp_K=298.0)

    def test_original_single_map_core_not_changed(self):
        self.assertEqual(self.one["n_maps"], 11)
        self.assertEqual(self.one["map_order"], MAP_NAMES)
        self.assertEqual(self.one["method"]["multiscale_radii_A"], [6.0, 9.0])
        self.assertEqual(self.one["method"]["patch_radius_A"], 8.0)
        self.assertEqual(self.one["method"]["sasa_points_per_atom"], 200)
        self.assertTrue(self.one["residue_profiles"])

    def test_two_identical_real_structures_give_identical_consensus(self):
        result = analyze_structure_set([
            {"name": "crambin_1", "pdb_text": self.pdb, "chain": "A"},
            {"name": "crambin_2", "pdb_text": self.pdb, "chain": "A"},
        ], pH=7.4, combined_maps=["hydrophobic", "pi_carbon"])
        self.assertEqual(result["n_structures"], 2)
        self.assertEqual(len(result["comparison"]), 1)
        self.assertIsNotNone(result["consensus"])
        self.assertIsNone(result["consensus_unavailable_reason"])
        self.assertFalse(result["comparison_skipped"])
        self.assertEqual(result["structures"][0]["analysis"]["maps"],
                         self.one["maps"])
        self.assertEqual(result["structures"][1]["analysis"]["maps"],
                         self.one["maps"])
        consensus = result["consensus"]
        self.assertEqual(consensus["n_matched_residues"], len(self.one["residue_profiles"]))
        for r in consensus["residues"]:
            self.assertIn(r["exposure_fraction"], [0.0, 1.0])
        for name in MAP_NAMES:
            for r in consensus["maps"][name]["residues"]:
                self.assertIn(r["compatible_fraction"], (0.0, 1.0))
                for k in (1, 3, 5):
                    self.assertIn(r[f"top{k}_membership_fraction"], (0.0, 1.0))
        for row in result["comparison"][0]["map_comparisons"]:
            for top in row["top_k"]:
                self.assertIn(top["jaccard"], (1.0, None))
        self.assertIsNotNone(consensus["combined_consensus"])

    def test_unmatched_chains_are_still_independently_analyzed(self):
        renamed = "\n".join(
            line[:21] + "B" + line[22:] if line.startswith("ATOM  ")
            and len(line) > 21 and line[21] == "A" else line
            for line in self.pdb.split("\n")
        )
        result = analyze_structure_set([
            {"name": "A", "pdb_text": self.pdb, "chain": "A"},
            {"name": "B", "pdb_text": renamed, "chain": "B"},
        ])
        self.assertEqual(len(result["structures"]), 2)
        self.assertIsNone(result["consensus"])
        self.assertTrue(result["consensus_unavailable_reason"])
        self.assertEqual(len(result["comparison_skipped"]), 1)
        self.assertFalse(result["comparison"])
        for s in result["structures"]:
            self.assertEqual(s["analysis"]["n_maps"], 11)

    def test_complete_real_PDB_zip_has_full_11_by_N_rows(self):
        payload = {
            "structures": [{
                "name": "real_1CRN", "pdb_text": self.pdb, "analysis": self.one,
                "screenshots": {name: png_uri() for name in MAP_NAMES},
                "combined_maps": None, "combined_screenshot": None,
            }],
            "comparison": [], "consensus": None, "consensus_screenshots": {},
            "comparison_skipped": [], "consensus_unavailable_reason": None,
        }
        # Check the exact FastAPI validation contract used by the browser.
        ExportRequest.model_validate(payload)
        zip_bytes = build_export_zip(payload)
        with zipfile.ZipFile(BytesIO(zip_bytes)) as z:
            names = z.namelist()
            base = "structures/real_1CRN"
            self.assertEqual(
                sum(x.startswith(base + "/pdf_maps/") and x.endswith(".pdf")
                    for x in names), 11)
            self.assertEqual(
                sum(x.startswith(base + "/screenshots/") and x.endswith(".png")
                    for x in names), 11)
            self.assertEqual(
                sum(x.startswith(base + "/pdb_bfactor/") and x.endswith(".pdb")
                    for x in names), 44)
            self.assertTrue(z.read(base + "/full_report.pdf").startswith(b"%PDF"))
            self.assertIn(b"/Subtype /Image", z.read(base + "/full_report.pdf"))
            for i, key in enumerate(MAP_NAMES, start=1):
                b = z.read(f"{base}/pdf_maps/{i:02d}_{key}.pdf")
                self.assertTrue(b.startswith(b"%PDF"), key)
                self.assertIn(b"/Subtype /Image", b, key)
                bfactor = z.read(
                    f"{base}/pdb_bfactor/{i:02d}_{key}_Bfactor.pdb"
                ).decode("utf-8")
                self.assertIn("NOT EXPERIMENTAL", bfactor)
            rows = list(csv.DictReader(StringIO(
                z.read(base + "/csv/Map_Residue_Scores.csv").decode("utf-8-sig"))))
            self.assertEqual(len(rows), 11 * len(self.one["residue_profiles"]))
            self.assertEqual({row["Map"] for row in rows}, set(MAP_NAMES))
            self.assertIn(base + "/csv/Map_Favorable_Residues.csv", names)
            self.assertIn("analysis_warnings.json", names)

    def test_request_validation_for_multi_structure(self):
        request = StructureSetRequest.model_validate({
            "structures": [{"name": "first", "pdb_text": self.pdb}],
            "pH": 7.4,
        })
        self.assertEqual(request.structures[0].name, "first")


if __name__ == "__main__":
    unittest.main()
