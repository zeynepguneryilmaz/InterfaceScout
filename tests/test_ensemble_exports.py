"""Deterministic regression checks for multi-PDB and complete output packages."""
from __future__ import annotations

import base64
from io import BytesIO
import sys
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from PIL import Image
from openpyxl import load_workbook
from ensemble import build_consensus
from export_bundle import MAP_NAMES, build_export_zip, _bfactor_pdb
from research_workflows import combine_maps, compare_results
from export_bundle import structure_tables, regional_profile_tables


def minimal_analysis():
    residue_key = "A:1:"
    patch = {
        "display_rank": 1, "pareto_front": 1,
        "center_key": residue_key, "members": [residue_key],
        "seed_members": [residue_key], "n_members": 1,
        "n_chemistry_seeds": 1, "chemistry_support": 1.0,
        "mean_accessibility": .42, "patch_coherence": .88,
        "orientation_coherence": 1.0,
        "member_residues": [{"key": residue_key, "chain": "A",
                             "res_seq": 1, "icode": "", "res_name": "LYS"}],
        "seed_residues": [{"key": residue_key, "chain": "A",
                           "res_seq": 1, "icode": "", "res_name": "LYS"}],
    }
    maps = {}
    for name in MAP_NAMES:
        maps[name] = {
            "label": name, "description": "Canonical test map",
            "n_patches": 1, "n_pareto_primary_patches": 1,
            "patches": [dict(patch)],
            "residues": [{"key": residue_key, "chain": "A", "res_seq": 1,
                          "res_name": "LYS", "scrsa": .42, "local_score": .42,
                          "propensity": 100.0, "multiscale_persistence": 88.0}],
        }
    return {
        "input": {"pH": 7.4, "chain": "A"},
        "method": {"sasa_points_per_atom": 200, "patch_radius_A": 8},
        "structure_preparation": {"selected_chain": "A"},
        "scope": {"models_adsorption_induced_unfolding": False},
        "map_order": list(MAP_NAMES), "maps": maps,
        "diagnostics": {"surface_residue_keys": [residue_key]},
        "residue_profiles": [{
            "key": residue_key, "res_name": "LYS", "chain": "A",
            "res_seq": 1, "icode": "",
            "scrsa": .42, "scrsa_raw": .42, "sidechain_sasa": 68.0,
            "surface_exposed": True, "charge_descriptor": .95,
            "groups": ["lysine", "basic_sidechain"],
        }],
        "protein_properties": {
            "n_residues": 1, "n_surface_residues": 1, "exposed_fraction": 1.0,
            "total_sidechain_sasa_A2": 68.0,
            "groups": {"lysine": {
                "n_exposed_residues": 1, "exposed_residues": [residue_key],
            }},
        },
    }


def png_data_uri():
    pic = Image.new("RGB", (270, 145), (225, 239, 249))
    for x in range(12, 250):
        for y in range(12, 128):
            if (x - 130)**2 / 9800 + (y - 70)**2 / 2800 < 1:
                pic.putpixel((x, y), (41 + x // 10, 80 + y // 4, 135))
    out = BytesIO()
    pic.save(out, format="PNG")
    return "data:image/png;base64," + base64.b64encode(out.getvalue()).decode()


PDB = (
    "ATOM      1  N   LYS A   1      10.000  10.000  10.000  1.00 20.00           N\n"
    "ATOM      2  CA  LYS A   1      11.000  10.000  10.000  1.00 20.00           C\n"
    "ATOM      3  C   LYS A   1      11.000  11.000  10.000  1.00 20.00           C\n"
    "ATOM      4  O   LYS A   1      11.000  12.000  10.000  1.00 20.00           O\n"
    "END\n"
)


class MultiStructureExportTests(unittest.TestCase):
    def setUp(self):
        self.a = minimal_analysis()
        self.b = minimal_analysis()
        self.c = build_consensus([
            {"name": "one", "analysis": self.a},
            {"name": "two", "analysis": self.b},
        ])
        self.pic = png_data_uri()

    def test_all_eleventh_maps_in_consensus(self):
        self.assertEqual(self.c["n_matched_residues"], 1)
        self.assertEqual(set(self.c["map_order"]), set(MAP_NAMES))
        for key in MAP_NAMES:
            r = self.c["maps"][key]["residues"][0]
            self.assertEqual(r["compatible_fraction"], 1.0)
            self.assertEqual(r["top1_membership_fraction"], 1.0)
            self.assertEqual(self.c["maps"][key]["strict_top5_members"], ["A:1:"])

    def test_consensus_fractions_are_not_probabilities(self):
        self.b["maps"]["anionic"]["residues"] = []
        self.b["maps"]["anionic"]["patches"] = []
        c = build_consensus([
            {"name": "one", "analysis": self.a},
            {"name": "two", "analysis": self.b},
        ])
        r = c["maps"]["anionic"]["residues"][0]
        self.assertEqual(r["compatible_fraction"], 0.5)
        self.assertEqual(r["top5_membership_fraction"], 0.5)
        self.assertEqual(c["maps"]["anionic"]["strict_top5_members"], [])

    def test_reject_unmatched_numbering(self):
        self.b["residue_profiles"][0]["key"] = "B:1:"
        with self.assertRaisesRegex(ValueError, "80%"):
            build_consensus([
                {"name": "one", "analysis": self.a},
                {"name": "two", "analysis": self.b},
            ])

    def test_bfactor_field_is_visualization_only(self):
        binary = _bfactor_pdb(PDB, {"A:1:": 72.5}, "Test")
        text = binary.decode()
        self.assertIn("NOT EXPERIMENTAL", text)
        atom = [line for line in text.splitlines() if line.startswith("ATOM")][0]
        self.assertEqual(atom[60:66].strip(), "72.50")
        self.assertEqual(atom[30:54], PDB.splitlines()[0][30:54])

    def test_named_structure_set_preserves_input_provenance(self):
        from ensemble import analyze_structure_set
        pdb_path=ROOT / "examples" / "1CRN" / "InterfaceScout_1CRN_anionic_Bfactor.pdb"
        pdb=pdb_path.read_text(encoding="utf-8")
        result=analyze_structure_set([
            {"name":"one","pdb_text":pdb,"chain":"A"},
            {"name":"two","pdb_text":pdb,"chain":"A"},
        ], pH=7.4)
        item=result["comparison"][0]
        self.assertEqual(item["structure_a"],"one")
        self.assertEqual(item["structure_b"],"two")
        self.assertIsInstance(item["input_a"],dict)
        self.assertIsInstance(item["input_b"],dict)
        self.assertEqual(item["input_a"]["chain"],"A")
        self.assertEqual(item["input_b"]["chain"],"A")

    def test_regional_profile_numeric_tables_and_complete_zip(self):
        # Even a deliberately empty/stale client combination must be rebuilt
        # from the original unchanged 11 canonical maps during export.
        combined = combine_maps(self.a, ["anionic", "hydrophobic"], max_regions=None)
        self.assertEqual(combined["regions_total"], 1)
        tables = structure_tables("one", self.a, combined)
        self.assertIn("Regional_Channels", tables)
        self.assertIn("Regional_Pair_Distances", tables)
        self.assertIn("Regional_Method", tables)
        channels = tables["Regional_Channels"]
        self.assertEqual(len(channels), 2)
        self.assertEqual({r["Map"] for r in channels}, {"anionic", "hydrophobic"})
        for row in channels:
            self.assertAlmostEqual(row["MeanLocalSupport_0_1"], .42)
            self.assertEqual(row["SupportingResidueCount"], 1)
            self.assertEqual(row["SupportingResidueFraction_0_1"], 1.0)
        self.assertEqual(len(tables["Regional_Pair_Distances"]), 1)
        self.assertIsNone(tables["Regional_Pair_Distances"][0]["MinimumDistinctResidueCA_A"])
        self.assertEqual(
            {r["Item"] for r in tables["Regional_Method"]} &
            {"RankingUnchanged", "UsesExperimentalLabels", "UsesFittedWeights"},
            {"RankingUnchanged", "UsesExperimentalLabels", "UsesFittedWeights"})
        unselected = structure_tables("one", self.a, None)
        self.assertNotIn("Regional_Channels", unselected)
        stale = dict(combined)
        stale["regions"] = []  # simulates a stale or UI-truncated payload
        stale["regions_total"] = 0
        packed = build_export_zip({
            "structures": [{
                "name": "one", "analysis": self.a, "pdb_text": PDB,
                "screenshots": {k: self.pic for k in MAP_NAMES},
                "combined_maps": stale, "combined_screenshot": self.pic
            }]
        })
        with zipfile.ZipFile(BytesIO(packed)) as z:
            base = "structures/one"
            for tab in ("Regional_Channels", "Regional_Pair_Distances", "Regional_Method"):
                self.assertIn(f"{base}/csv/{tab}.csv", z.namelist())
            import csv, json
            from io import StringIO
            rows = list(csv.DictReader(StringIO(
                z.read(base + "/csv/Regional_Channels.csv").decode("utf-8-sig"))))
            self.assertEqual(len(rows), 2)
            self.assertEqual({r["Map"] for r in rows}, {"anionic", "hydrophobic"})
            self.assertEqual(float(rows[0]["MeanLocalSupport_0_1"]), .42)
            self.assertEqual(float(rows[0]["SupportingResidueFraction_0_1"]), 1.0)
            full = json.loads(z.read(base + "/combined/combined_maps.json"))
            self.assertEqual(full["regions_total"], 1)
            self.assertEqual(len(full["regions"]), 1)
            self.assertEqual(full["regional_profile_method"]["method_version"],
                             "regional-chemistry-descriptive-v2")
            self.assertFalse(full["regional_profile_method"]["experimental_labels_used"])
            self.assertNotIn("weakest_mean_support", full["regions"][0]["regional_profile"])
            self.assertNotIn("relative_to_candidate_median_by_map",
                             full["regions"][0]["regional_profile"])
            wb = load_workbook(BytesIO(z.read(base + "/complete_results.xlsx")),
                               read_only=True, data_only=True)
            for tab in ("Regional_Channels", "Regional_Pair_Distances", "Regional_Method"):
                self.assertIn(tab, wb.sheetnames)
            self.assertEqual(wb["Regional_Channels"].max_row, 3)
            wb.close()
            self.assertTrue(z.read(base + "/Complete_Surface_Atlas.pdf").startswith(b"%PDF"))

    def test_full_zip_contains_all_maps_and_reports(self):
        comparison = compare_results(self.a, self.b)
        comparison["structure_a"] = "one"
        comparison["structure_b"] = "two"
        result = build_export_zip({
            "structures": [
                {"name": name, "analysis": analysis, "pdb_text": PDB,
                 "screenshots": {k: self.pic for k in MAP_NAMES},
                 "combined_maps": None}
                for name, analysis in [("one", self.a), ("two", self.b)]
            ],
            "consensus": self.c,
            "comparison": [comparison],
            "consensus_screenshots": {key: self.pic for key in MAP_NAMES},
        })
        with zipfile.ZipFile(BytesIO(result)) as z:
            files = z.namelist()
            for structure in ("one", "two"):
                base = f"structures/{structure}"
                self.assertEqual(len([f for f in files if
                    f.startswith(f"{base}/pdf_maps/") and f.endswith(".pdf")]), 11)
                self.assertEqual(len([f for f in files if
                    f.startswith(f"{base}/screenshots/") and f.endswith(".png")]), 11)
                self.assertTrue(z.read(f"{base}/full_report.pdf").startswith(b"%PDF"))
                self.assertIn(f"{base}/Complete_Surface_Atlas.pdf", files)
                self.assertEqual(z.read(f"{base}/Complete_Surface_Atlas.pdf"),
                                 z.read(f"{base}/full_report.pdf"))
                self.assertIn(f"{base}/csv/Protein_Residues.csv", files)
                self.assertIn(f"{base}/csv/Map_Residue_Scores.csv", files)
                self.assertIn(f"{base}/csv/All_Patches.csv", files)
                self.assertIn(f"{base}/csv/Patch_Members.csv", files)
                self.assertIn(f"{base}/complete_results.xlsx", files)
                self.assertEqual(len([f for f in files if
                    f.startswith(f"{base}/pdb_bfactor/") and f.endswith(".pdb")]), 44)
                for key in MAP_NAMES:
                    pdf = z.read(f"{base}/pdf_maps/{MAP_NAMES.index(key)+1:02d}_{key}.pdf")
                    self.assertTrue(pdf.startswith(b"%PDF"), key)
                wb = load_workbook(BytesIO(z.read(f"{base}/complete_results.xlsx")))
                self.assertIn("Map_Residue_Scores", wb.sheetnames)
            self.assertIn("consensus/consensus_report.pdf", files)
            self.assertIn("comparison/comparison_report.pdf", files)
            self.assertEqual(len([f for f in files if
                f.startswith("consensus/pdb_bfactor/") and f.endswith(".pdb")]), 22)
            self.assertIn("consensus/csv/Consensus_Map_Scores.csv", files)
            self.assertIn("comparison/csv/Comparison_Maps.csv", files)


if __name__ == "__main__":
    unittest.main()
