"""Regression tests: structural identity is not a PDB chain/numbering label.

Test the public analysis, 11 maps, pairwise comparison, consensus, selection,
altloc resolution, and checks against silent structural data loss.
"""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from interface_engine import analyze_all_maps
from ensemble import build_consensus
from prepare import prepare_pdb_text
from research_workflows import compare_results
from residue_matching import match_residues


def relabel_pdb(text: str, chain: str = "B", offset: int = 100) -> str:
    modified = []
    for line in text.splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            seq = int(line[22:26])
            line = line[:21] + chain + f"{seq + offset:4d}" + line[26:]
        modified.append(line)
    return "".join(modified)


class StructuralIdentityRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = (ROOT / "examples" / "1CRN" /
                   "InterfaceScout_1CRN_anionic_Bfactor.pdb").read_text("utf-8")
        cls.relabelled = relabel_pdb(cls.raw, "Z", 100)
        cls.a = analyze_all_maps(pdb_text=cls.raw, chain="A")
        cls.b = analyze_all_maps(pdb_text=cls.relabelled, chain="Z")

    def test_all_11_maps_keep_scores_and_region_geometry(self):
        self.assertEqual(self.a["n_maps"], 11)
        self.assertEqual(self.a["map_order"], self.b["map_order"])
        by_position_a = {r["res_seq"]: r for r in self.a["residue_profiles"]}
        by_position_b = {r["res_seq"]-100: r for r in self.b["residue_profiles"]}
        self.assertEqual(set(by_position_a), set(by_position_b))
        for seq in by_position_a:
            self.assertEqual(by_position_a[seq]["res_name"], by_position_b[seq]["res_name"])
            self.assertAlmostEqual(
                by_position_a[seq]["scrsa_raw"], by_position_b[seq]["scrsa_raw"], places=5)
        for name in self.a["map_order"]:
            ma, mb = self.a["maps"][name], self.b["maps"][name]
            self.assertEqual(ma["n_patches"], mb["n_patches"], name)
            r_a = {r["res_seq"]: r for r in ma["residues"]}
            r_b = {r["res_seq"]-100: r for r in mb["residues"]}
            self.assertEqual(set(r_a), set(r_b), name)
            for seq in r_a:
                self.assertEqual(r_a[seq]["local_score"], r_b[seq]["local_score"], name)
                self.assertEqual(r_a[seq]["propensity"], r_b[seq]["propensity"], name)
            for pa, pb in zip(ma["patches"], mb["patches"]):
                self.assertEqual(pa["center_residue"]["res_seq"] + 100,
                                 pb["center_residue"]["res_seq"], name)
                for metric in ("pareto_front", "chemistry_support",
                               "mean_accessibility", "patch_coherence",
                               "orientation_coherence", "center_coordinate_A"):
                    self.assertEqual(pa[metric], pb[metric], (name, metric))
                self.assertEqual(
                    sorted(int(k.split(":")[1]) for k in pa["members"]),
                    sorted(int(k.split(":")[1])-100 for k in pb["members"]), name)

    def test_comparison_remaps_different_numbering(self):
        c = compare_results(self.a, self.b)
        self.assertEqual(c["mapping"]["n_identity_matches"],
                         len(self.a["residue_profiles"]))
        self.assertEqual(c["mapping"]["coverage_a"], 1)
        self.assertEqual(c["mapping"]["coverage_b"], 1)
        self.assertEqual(c["accessibility"]["n_exposure_status_changes"], 0)
        self.assertTrue(all(abs(r["delta_scRSA"]) < 1e-8
                            for r in c["residue_changes"]))
        for m in c["map_comparisons"]:
            for k in m["top_k"]:
                if k["jaccard"] is not None:
                    self.assertAlmostEqual(k["jaccard"], 1.0)

    def test_multi_structure_consensus_uses_reference_numbering(self):
        runs = [{"name": "original", "analysis": self.a, "combined_maps": None},
                {"name": "relabelled", "analysis": self.b, "combined_maps": None}]
        c = build_consensus(runs)
        self.assertEqual(c["n_matched_residues"], len(self.a["residue_profiles"]))
        self.assertTrue(all(r["scRSA_range"] == 0 for r in c["residues"]))
        for name in c["map_order"]:
            self.assertTrue(all(row["compatible_fraction"] in (0.0, 1.0)
                                for row in c["maps"][name]["residues"]))
            self.assertTrue(all(row["top5_membership_fraction"] in (0.0, 1.0)
                                for row in c["maps"][name]["residues"]))

    def test_prep_keeps_correct_chain_after_relabelling(self):
        _, report = prepare_pdb_text(self.relabelled, chain="Z")
        self.assertEqual(report["selected_chain"], "Z")
        self.assertGreater(report["protein_residues_retained"], 0)
        with self.assertRaisesRegex(ValueError, "not found"):
            prepare_pdb_text(self.relabelled, chain="A")

    def test_altloc_occupancy_selection_is_one_per_residue(self):
        lines = self.raw.splitlines(keepends=True)
        candidate = next(i for i, line in enumerate(lines)
                         if line.startswith("ATOM  ") and line[12:16].strip() == "CB")
        old = lines[candidate]
        # First side-chain atom has 40% conformer A and 60% conformer B.
        a = old[:16] + "A" + old[17:54] + "  0.40" + old[60:]
        b = old[:16] + "B" + old[17:30] + f"{float(old[30:38])+2.0:8.3f}" + old[38:54] + "  0.60" + old[60:]
        lines[candidate:candidate+1] = [a, b]
        prepared, report = prepare_pdb_text("".join(lines), chain="A")
        self.assertGreaterEqual(report["alternate_residues_resolved"], 1)
        cb = [line for line in prepared.splitlines() if line.startswith("ATOM  ")
              and line[12:16].strip() == "CB" and line[22:26] == old[22:26]]
        self.assertEqual(len(cb), 1)
        self.assertEqual(cb[0][16], "B")

    def test_mmcif_is_not_misparsed_as_pdb(self):
        with self.assertRaisesRegex(ValueError, "mmCIF"):
            prepare_pdb_text("data_test\nloop_\n_atom_site.Cartn_x\n1.0\n")

    def test_repeated_homomer_chain_names_resolved_by_coordinates(self):
        def row(key, xyz):
            chain, seq = key.split(":")
            return dict(key=key+":", chain=chain, res_seq=int(seq),
                        icode="", res_name="ALA", x=xyz, y=0., z=0.)
        a = [row("A:1", 0), row("A:2", 1),
             row("B:1", 50), row("B:2", 51)]
        b = [row("A:101", 50), row("A:102", 51),
             row("B:101", 0), row("B:102", 1)]
        m = match_residues(a, b)
        self.assertEqual(m["a_to_b"]["A:1:"], "B:101:")
        self.assertEqual(m["a_to_b"]["B:1:"], "A:101:")
        self.assertEqual(len(m["chain_matches"]), 2)

    def test_duplicate_residue_keys_rejected(self):
        row = dict(key="A:1:", chain="A", res_seq=1, icode="",
                   res_name="ALA", x=0., y=0., z=0.)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            match_residues([row, row], [row])


if __name__ == "__main__":
    unittest.main()
