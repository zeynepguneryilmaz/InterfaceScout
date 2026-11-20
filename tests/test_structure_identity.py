"""Structure identity invariance and conservative input-reading regression tests.

No experimental adsorption labels are used. These tests deliberately change PDB
metadata while holding atomic coordinates fixed.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from interface_engine import analyze_all_maps
from research_workflows import compare_results
from ensemble import analyze_structure_set
from prepare import prepare_pdb_text

SOURCE = ROOT / "examples" / "1CRN" / "InterfaceScout_1CRN_anionic_Bfactor.pdb"


def relabel(pdb: str, *, chain: str = "B", offset: int = 100,
            insertion: bool = False) -> str:
    result = []
    first_residue = None
    for line in pdb.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            line = line.ljust(80)
            original = int(line[22:26])
            if first_residue is None:
                first_residue = original
            icode = "A" if insertion and original == first_residue else line[26]
            line = line[:21] + chain + f"{original + offset:4d}" + icode + line[27:]
        result.append(line)
    return "\n".join(result) + "\n"


class StructureIdentityRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = SOURCE.read_text("utf-8")
        cls.left = analyze_all_maps(pdb_text=cls.raw, chain="A")
        cls.changed = relabel(cls.raw)
        cls.right = analyze_all_maps(pdb_text=cls.changed, chain="B")

    def test_chain_and_numbering_do_not_change_chemistry_or_top_regions(self):
        a = self.left
        b = self.right
        self.assertEqual(a["n_maps"], 11)
        self.assertEqual(b["n_maps"], 11)
        comparison = compare_results(a, b)
        mapping = comparison["mapping"]
        self.assertEqual(mapping["n_identity_matches"], len(a["residue_profiles"]))
        self.assertEqual(mapping["coverage_a"], 1)
        self.assertEqual(mapping["coverage_b"], 1)
        self.assertEqual(comparison["accessibility"]["n_exposure_status_changes"], 0)
        self.assertTrue(all(abs(r["delta_scRSA"]) < 1e-5
                            for r in comparison["residue_changes"]))
        self.assertTrue(all(abs(r["delta_sidechain_sasa"]) < 1e-3
                            for r in comparison["residue_changes"]))
        a_to_b = mapping["reference_to_comparison"]
        for name in a["map_order"]:
            ar = a["maps"][name]
            br = b["maps"][name]
            self.assertEqual(len(ar["patches"]), len(br["patches"]), name)
            for x, y in zip(ar["patches"], br["patches"]):
                self.assertEqual(a_to_b[x["center_key"]], y["center_key"], name)
                self.assertEqual(
                    {a_to_b[t] for t in x["members"]},
                    set(y["members"]), name)
                for field in ("chemistry_support", "mean_accessibility",
                              "patch_coherence", "orientation_coherence"):
                    self.assertAlmostEqual(x[field], y[field], places=5)
        for m in comparison["map_comparisons"]:
            for k in m["top_k"]:
                self.assertIn(k["jaccard"], (1.0, None), (m["key"], k["top_k"]))

    def test_insertions_and_renumbering_have_explicit_correspondence(self):
        right = analyze_all_maps(pdb_text=relabel(self.raw, insertion=True), chain="B")
        result = compare_results(self.left, right)
        self.assertEqual(result["mapping"]["n_identity_matches"],
                         len(self.left["residue_profiles"]))
        self.assertTrue(any(x["key"] != x["key_b"] for x in result["residue_changes"]))

    def test_multi_structure_consensus_in_reference_numbering(self):
        out = analyze_structure_set([
            {"name": "original", "pdb_text": self.raw, "chain": "A"},
            {"name": "shifted", "pdb_text": self.changed, "chain": "B"},
        ])
        self.assertEqual(len(out["comparison"]), 1)
        self.assertIsNone(out["consensus_unavailable_reason"])
        c = out["consensus"]
        self.assertEqual(c["n_matched_residues"],
                         len(self.left["residue_profiles"]))
        self.assertTrue(all(r["scRSA_range"] < 1e-4 for r in c["residues"]))
        self.assertEqual(set(c["reference_to_structure_keys"]["original"]),
                         set(c["reference_to_structure_keys"]["shifted"]))
        for name in c["map_order"]:
            for row in c["maps"][name]["residues"]:
                for k in (1, 3, 5):
                    self.assertIn(row[f"top{k}_membership_fraction"], (0.0, 1.0))

    def test_unrelated_sequence_cannot_be_equated_by_numbers(self):
        lines = []
        for line in self.raw.splitlines():
            if line.startswith("ATOM  "):
                line = line[:17] + "GLY" + line[20:]
            lines.append(line)
        # Do not call the SASA engine for impossible modified-residue atoms:
        # residue identity matching can be tested directly on profile inputs.
        from residue_matching import match_residues
        fake = [{**r, "res_name": "GLY"} for r in self.left["residue_profiles"]]
        with self.assertRaises(ValueError):
            match_residues(self.left["residue_profiles"], fake)

    def test_detects_missing_ca_and_non_pdb_input(self):
        with self.assertRaisesRegex(ValueError, "mmCIF"):
            prepare_pdb_text("data_test\nloop_\n_atom_site.id\n_atom_site.Cartn_x\n")
        first = next(l for l in self.raw.splitlines() if l.startswith("ATOM  ") and l[12:16].strip() == "CA")
        prefix = first[21:27]
        truncated = "\n".join(
            l for l in self.raw.splitlines()
            if not (l.startswith("ATOM  ") and l[21:27] == prefix
                    and l[12:16].strip() == "CA")
        )
        with self.assertRaisesRegex(ValueError, "C-alpha"):
            prepare_pdb_text(truncated)

    def test_duplicate_atom_record_is_rejected(self):
        first_atom = next(
            l for l in self.raw.splitlines() if l.startswith("ATOM  "))
        lines = self.raw.splitlines()
        index = lines.index(first_atom)
        lines.insert(index + 1, first_atom)
        with self.assertRaisesRegex(ValueError, "Duplicate first-model atom"):
            prepare_pdb_text("\n".join(lines) + "\n", chain="A")

    def test_warns_on_incomplete_sidechain(self):
        first_sidechain = next(
            l for l in self.raw.splitlines()
            if l.startswith("ATOM  ") and l[12:16].strip() == "CB")
        reduced = self.raw.replace(first_sidechain + "\n", "")
        _, report = prepare_pdb_text(reduced, chain="A")
        self.assertGreaterEqual(len(report["incomplete_sidechain_residues"]), 1)
        self.assertTrue(report["warnings"])


if __name__ == "__main__":
    unittest.main()
