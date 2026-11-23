"""Unit tests for opt-in electrostatics; no APBS binary required for pure tests."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from optional_electrostatics import (
  ElectrostaticSettings,ElectrostaticsError,atoms_from_pqr,input_for_apbs,status
)

class TestElectrostaticsParameters(unittest.TestCase):
    def setUp(self):
        self.atoms=[([1.,2.,3.],-.5,1.5),([3.,3.,3.],.5,1.7)]

    def test_default_apbs_conditions(self):
        s=ElectrostaticSettings()
        s.validate()
        script,grid=input_for_apbs(s,self.atoms)
        for t in ("lpbe","ion charge +1 conc 0.150000","ion charge -1 conc 0.150000",
                  "pdie 2.0000","sdie 78.5000","temp 298.150","write pot dx potential"):
            self.assertIn(t,script)
        self.assertEqual(grid["electrolyte_concentration_M"],.15)

    def test_zero_salt(self):
        script,_=input_for_apbs(ElectrostaticSettings(ionic_mM=0),self.atoms)
        self.assertNotIn("ion charge",script)

    def test_invalid_settings(self):
        for bad in (dict(pH=15),dict(ionic_mM=-1),dict(temperature_K=200),
                    dict(electrolyte="CaCl2"),dict(protein_dielectric=0)):
            with self.subTest(bad=bad),self.assertRaises(ElectrostaticsError):
                ElectrostaticSettings(**bad).validate()

    def test_pqr_charge_radius_parsing(self):
        pqr="ATOM      1  N   ALA A   1       1.000   2.000   3.000 -0.3000 1.5000\n"
        atoms=atoms_from_pqr(pqr)
        self.assertAlmostEqual(atoms[0][1],-.3)
        self.assertEqual(atoms[0][2],1.5)
        with self.assertRaises(ElectrostaticsError):
            atoms_from_pqr(pqr.replace("-0.3000","nan"))

    def test_charged_zero_radius_hydrogen_preserved(self):
        pqr="ATOM     47  HZ1 LYS A   3      52.270   4.898  11.784  0.3300 0.0000\\n"
        parsed=atoms_from_pqr(pqr)
        self.assertEqual(len(parsed),1)
        self.assertAlmostEqual(parsed[0][1],0.33)
        from physics_surface import read_pqr,charge_by_residue
        atoms=read_pqr(pqr)
        self.assertAlmostEqual(charge_by_residue(atoms)["A:3:"],0.33)

    def test_zero_radius_heavy_atom_is_rejected(self):
        pqr="ATOM     47  NZ  LYS A   3      52.270   4.898  11.784  0.3300 0.0000\\n"
        with self.assertRaises(ElectrostaticsError):
            atoms_from_pqr(pqr)

    def test_optional_does_not_modify_maps(self):
        self.assertFalse(status()["changes_standard_maps"])

if __name__=="__main__":unittest.main()
