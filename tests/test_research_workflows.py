from __future__ import annotations
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
from research_workflows import RESEARCH_GOALS, combine_maps, compare_results
from protein_properties import profile_rows, summarize_profiles

def res(key, name="LYS", exposure=True, sasa=12.0, scrsa=0.20):
    chain, seq, icode = key.split(":")
    return {
        "key": key, "chain": chain, "res_seq": int(seq), "icode": icode,
        "res_name": name, "sidechain_sasa": sasa, "scrsa_raw": scrsa,
        "scrsa": min(scrsa, 1.0), "surface_exposed": exposure,
        "charge_descriptor": 1.0,
    }

def analysis(patches1, patches2, profiles, selected=("anionic", "hydrophobic")):
    return {
        "input": {"pH": 7.4, "chain": "A"},
        "residue_profiles": profiles,
        "map_order": list(selected),
        "maps": {
            selected[0]: {
                "key": selected[0], "label": "Anionic",
                "residues": [{"key": "A:1:", "local_score": .5}],
                "patches": patches1,
            },
            selected[1]: {
                "key": selected[1], "label": "Hydrophobic",
                "residues": [{"key": "A:2:", "local_score": .8}],
                "patches": patches2,
            },
        },
        "diagnostics": {"surface_residue_keys": [p["key"] for p in profiles if p["surface_exposed"]]},
    }

class DecisionWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.rows = [res("A:1:"), res("A:2:", "LEU"), res("A:3:", "ALA")]
        p = [{"center_key": "A:1:", "members": ["A:1:", "A:2:"]}]
        q = [{"center_key": "A:2:", "members": ["A:2:", "A:3:"]}]
        self.a = analysis(p, q, self.rows)

    def test_goals_omit_unnecessary_dynamics(self):
        self.assertEqual(len(RESEARCH_GOALS), 5)
        self.assertIn("noncovalent", RESEARCH_GOALS)
        for goal in RESEARCH_GOALS.values():
            self.assertFalse("GNM" in goal["essentials"])
            self.assertFalse("ANM" in goal["essentials"])
            self.assertFalse("RINpy" in goal["essentials"])

    def test_same_residue_vs_colocation(self):
        out = combine_maps(self.a, ["anionic", "hydrophobic"])
        self.assertEqual(out["n_exact_overlap"], 0)
        self.assertTrue(out["regions"])
        self.assertTrue(any(r["distinct_residue_colocation_only"] for r in out["regions"]))
        self.assertTrue(any(r["all_maps_colocated"] for r in out["regions"]))

    def test_selection_validation(self):
        with self.assertRaises(ValueError):
            combine_maps(self.a, ["anionic"])
        with self.assertRaises(ValueError):
            combine_maps(self.a, ["anionic", "unknown"])

    def test_properties(self):
        direct = profile_rows({"all_residues": self.rows})
        report = summarize_profiles(direct)
        self.assertEqual(report["n_surface_residues"], 3)
        self.assertEqual(report["groups"]["lysine"]["n_exposed_residues"], 1)
        self.assertEqual(report["groups"]["nonpolar"]["n_exposed_residues"], 2)

    def test_conformation_changes(self):
        changed = [
            res("A:1:", exposure=False, sasa=0.0, scrsa=0.0),
            res("A:2:", "LEU"), res("A:3:", "ALA"),
        ]
        b = analysis(
            [{"center_key": "A:2:", "members": ["A:2:"]}],
            [{"center_key": "A:2:", "members": ["A:2:", "A:3:"]}],
            changed,
        )
        out = compare_results(self.a, b)
        self.assertEqual(out["mapping"]["n_identity_matches"], 3)
        self.assertEqual(out["accessibility"]["n_exposure_status_changes"], 1)
        self.assertLess(out["map_comparisons"][0]["top_k"][0]["jaccard"], 1.0)
        self.assertIn("adsorption", out["interpretation"])

    def test_compare_renamed_chains(self):
        b = analysis(
            [{"center_key": "B:1:", "members": ["B:1:"]}],
            [{"center_key": "B:2:", "members": ["B:2:"]}],
            [res("B:1:"), res("B:2:", "LEU"), res("B:3:", "ALA")],
        )
        result = compare_results(self.a, b)
        self.assertEqual(result["mapping"]["n_identity_matches"], 3)
        self.assertEqual(result["mapping"]["chain_matches"][0]["chain_b"], "B")

if __name__ == "__main__":
    unittest.main()
