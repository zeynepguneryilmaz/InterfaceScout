"""Identity is grounded in the supplied PDB, not inferred from a filename."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
from pdf_identity import pdb_identity


class PdfIdentityTests(unittest.TestCase):
    def test_header_and_compnd_are_verified_evidence(self):
        header="HEADER    PLANT PROTEIN                           15-MAY-82   1CRN"
        pdb="\n".join((header,
            "COMPND    MOL_ID: 1;",
            "COMPND   2 MOLECULE: CRAMBIN;",
            "COMPND   3 CHAIN: A;",
            "ATOM      1  CA  THR A   1      10.000  10.000  10.000  1.00 20.00           C"))
        result=pdb_identity(pdb,{"input":{"chain":"A","pdb_id":"1CRN"},"version":"1.2"},"sample")
        self.assertEqual(result["name"],"CRAMBIN")
        self.assertEqual(result["accession"],"1CRN")
        self.assertEqual(result["name_evidence"],"original PDB COMPND")
        self.assertEqual(result["accession_evidence"],"original PDB HEADER")
        self.assertEqual(result["chain"],"A")

    def test_no_filename_name_inference(self):
        pdb="ATOM      1  CA  THR A   1      10.000  10.000  10.000  1.00 20.00           C\nEND\n"
        result=pdb_identity(pdb,{"input":{"chain":"A"}},"InterfaceScout_1CRN_anionic_Bfactor.pdb")
        self.assertEqual(result["name"],"Protein name not provided in original PDB")
        self.assertEqual(result["accession"],"not verified")
        self.assertEqual(result["accession_evidence"],"not available")
        self.assertIn("1CRN",result["source_label"])

    def test_analysis_requested_id_is_labeled_not_verified(self):
        result=pdb_identity("",{"input":{"pdb_id":"2XYZ","chain":"B"}},"uploaded")
        self.assertEqual(result["accession"],"2XYZ")
        self.assertEqual(result["accession_evidence"],"requested PDB ID")


if __name__=="__main__":
    unittest.main()
