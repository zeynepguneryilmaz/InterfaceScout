"""Evidence-grade scope tests: numbering, assemblies, temperature/ionic strength.

These explicitly distinguish a functioning feature from an unsupported feature.
No misleading pass is assigned to automatic residue alignment or automatic
biological-assembly inference: those currently do not exist.
"""
from __future__ import annotations
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
from interface_engine import analyze_all_maps
from ensemble import analyze_structure_set
from research_workflows import compare_results

PDB=(ROOT/"examples"/"1CRN"/"InterfaceScout_1CRN_anionic_Bfactor.pdb").read_text("utf-8")

def renumber(pdb:str, offset:int, chain:str="B")->str:
    lines=[]
    for line in pdb.splitlines():
        if line.startswith(("ATOM  ","HETATM","ANISOU")) and len(line)>=27:
            old=line[22:26].strip()
            if old.lstrip("-").isdigit():
                line=line[:21]+chain+f"{int(old)+offset:4d}"+line[26:]
        lines.append(line)
    return "\n".join(lines)+"\n"

class ScopeRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a=analyze_all_maps(pdb_text=PDB,chain="A",pH=7.4,ionic_mM=150,temp_K=298)
    def test_renumbered_same_sequence_is_not_silently_miscompared(self):
        shifted=renumber(PDB,120,chain="B")
        other=analyze_all_maps(pdb_text=shifted,chain="B",pH=7.4)
        left={r["res_name"] for r in self.a["residue_profiles"]}
        self.assertEqual([r["res_name"] for r in self.a["residue_profiles"]],
                         [r["res_name"] for r in other["residue_profiles"]])
        self.assertNotEqual(
            set(r["key"] for r in self.a["residue_profiles"]),
            set(r["key"] for r in other["residue_profiles"]))
        result=analyze_structure_set([
            {"name":"original","pdb_text":PDB,"chain":"A"},
            {"name":"renumbered","pdb_text":shifted,"chain":"B"},
        ])
        self.assertEqual(len(result["structures"]),2)
        self.assertIsNone(result["consensus"])
        self.assertTrue(result["consensus_unavailable_reason"])
        self.assertEqual(len(result["comparison_skipped"]),1)

    def test_ionic_strength_is_metadata_not_a_physics_model(self):
        hi=analyze_all_maps(pdb_text=PDB,chain="A",pH=7.4,ionic_mM=1500,temp_K=298)
        lo=analyze_all_maps(pdb_text=PDB,chain="A",pH=7.4,ionic_mM=0,temp_K=298)
        self.assertEqual(self.a["maps"],hi["maps"])
        self.assertEqual(self.a["maps"],lo["maps"])
        self.assertNotEqual(hi["input"]["ionic_mM"],lo["input"]["ionic_mM"])

    def test_temperature_is_metadata_not_thermodynamic_scoring(self):
        cold=analyze_all_maps(pdb_text=PDB,chain="A",pH=7.4,ionic_mM=150,temp_K=277)
        hot=analyze_all_maps(pdb_text=PDB,chain="A",pH=7.4,ionic_mM=150,temp_K=330)
        self.assertEqual(cold["maps"],hot["maps"])
        self.assertNotEqual(cold["input"]["temp_K"],hot["input"]["temp_K"])

    def test_pH_is_active_scoring_parameter(self):
        acidic=analyze_all_maps(pdb_text=PDB,chain="A",pH=2.0,ionic_mM=150,temp_K=298)
        basic=analyze_all_maps(pdb_text=PDB,chain="A",pH=11.0,ionic_mM=150,temp_K=298)
        self.assertTrue(any(
            acidic["maps"][key]["residues"]!=basic["maps"][key]["residues"]
            for key in acidic["map_order"]
        ))
        self.assertEqual(acidic["input"]["pH"],2.0)
        self.assertEqual(basic["input"]["pH"],11.0)

    def test_first_model_selected_without_auto_biological_assembly(self):
        # A two-model protein is *not* interpreted as two chains; preparation
        # explicitly selects the first PDB MODEL.
        a="\n".join(line for line in PDB.splitlines() if line.strip() != "END") + "\n"
        multi="MODEL        1\n"+a+"ENDMDL\nMODEL        2\n"+a+"ENDMDL\nEND\n"
        one=analyze_all_maps(pdb_text=multi,chain="A",pH=7.4)
        self.assertEqual(one["structure_preparation"]["input_models"],2)
        self.assertEqual(one["structure_preparation"]["selected_model"],1)
        self.assertEqual(one["n_maps"],11)
        self.assertEqual(one["maps"],self.a["maps"])

if __name__=="__main__":
    unittest.main()
