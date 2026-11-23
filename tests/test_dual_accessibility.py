"""Regression tests for production dual-accessibility model.

Independent atom aggregation (same Shrake-Rupley kernel) verifies total area,
true side-chain atom selection, GLY zero chemical support, geometry retention,
and 11-map public API metadata. Does not claim independent SASA implementation.
"""
from __future__ import annotations
import sys
import math
import unittest
from io import StringIO
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import core
import geometry
from interface_engine import analyze_all_maps
from test_publication_smoke import PublicationSmokeTests

PDB=(ROOT/"examples"/"1CRN"/"InterfaceScout_1CRN_anionic_Bfactor.pdb").read_text()

class DualAccessibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=core._surface_residues(PDB,7.4)
        cls.result=analyze_all_maps(pdb_text=PDB,chain="A",pH=7.4)

    def test_frozen_references(self):
        self.assertEqual(len(core.TOTAL_RESIDUE_REF_ASA),20)
        self.assertEqual(core.TOTAL_RESIDUE_REF_ASA["GLY"],104.)
        self.assertEqual(core.SC_RSA_THRESHOLD,.05)

    def test_independent_atom_bookkeeping_and_glu_gly(self):
        st=PDBParser(QUIET=True).get_structure("ref",StringIO(PDB))
        ShrakeRupley(probe_radius=core.SASA_PROBE_A,n_points=core.SASA_POINTS).compute(st,level="A")
        reference={}
        for chain in next(st.get_models()):
            for residue in chain:
                if residue.resname not in core.TOTAL_RESIDUE_REF_ASA or "CA" not in residue:continue
                key=f"{chain.id}:{int(residue.id[1])}:{str(residue.id[2]).strip()}"
                atoms=[a for a in residue if a.element!="H"]
                side=[a for a in atoms if a.name not in {"N","CA","C","O","OXT"}]
                reference[key]=(residue.resname,sum(float(a.sasa) for a in atoms),
                                sum(float(a.sasa) for a in side))
        self.assertEqual(set(reference),{r["key"] for r in self.rows})
        gly=0
        for r in self.rows:
            name,whole,side=reference[r["key"]]
            self.assertAlmostEqual(r["residue_sasa"],whole,delta=.00051)
            self.assertAlmostEqual(r["trsa_raw"],whole/core.TOTAL_RESIDUE_REF_ASA[name],delta=.0000051)
            self.assertEqual(r["surface_exposed"], whole/core.TOTAL_RESIDUE_REF_ASA[name]>=.05)
            self.assertAlmostEqual(r["sidechain_sasa"],side,delta=.00051)
            if name=="GLY":
                gly+=1
                self.assertEqual(r["sidechain_sasa"],0.)
                self.assertEqual(r["scrsa_raw"],0.)
                self.assertEqual(r["scrsa"],0.)
        self.assertGreater(gly,0)

    def test_public_api_dual_definition(self):
        r=self.result
        self.assertEqual(r["n_maps"],11)
        self.assertEqual(r["method"]["surface_gate"],"total_residue_RSA")
        self.assertEqual(r["method"]["total_rsa_threshold"],.05)
        self.assertEqual(r["method"]["glycine_policy"],"whole-residue RSA for geometry; no side-chain chemistry")
        self.assertEqual(len(r["residue_profiles"]),len(self.rows))
        self.assertTrue(all("trsa" in x and "scrsa" in x for x in r["residue_profiles"]))
        self.assertEqual(len(r["diagnostics"]["surface_residue_keys"]),
                         sum(x["surface_exposed"] for x in r["residue_profiles"]))

    def test_geometry_accessibility_uses_trsa(self):
        analysis=core.analyze(core.AnalyzeRequest(pdb_text=PDB))
        nodes=geometry.extract_ca_nodes(PDB)
        g=geometry.build_surface_geometry(analysis,nodes)
        source={r["key"]:r for r in analysis["surface_residues"]}
        self.assertEqual(set(source),set(g["trsa"]))
        for key,val in g["trsa"].items():
            self.assertEqual(val,source[key]["trsa"])

    def test_glycine_retained_without_chemical_score(self):
        analysis=core.analyze(core.AnalyzeRequest(pdb_text=PDB))
        gly_surface=[r for r in analysis["surface_residues"] if r["res_name"]=="GLY"]
        self.assertTrue(gly_surface,"Fixture must include solvent-exposed glycine")
        for r in gly_surface:
            self.assertGreater(r["trsa"],0)
            self.assertEqual(r["scrsa"],0)
            self.assertEqual(r["sidechain_sasa"],0)
        for channel in analysis["chemistries"].values():
            self.assertFalse(any(r["res_name"]=="GLY" for r in channel["residues"]))
            self.assertFalse(any(r["res_name"]=="GLY" for r in channel["repulsive_residues"]))

if __name__=="__main__":
    unittest.main()
