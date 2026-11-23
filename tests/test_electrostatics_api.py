"""HTTP contract tests for optional physics-aware UI and archive routes."""
from __future__ import annotations
import io,sys,unittest,zipfile
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"backend"))
from fastapi.testclient import TestClient
from app import app

class TestElectrostaticsHttp(unittest.TestCase):
    def test_prepare_download_and_dx(self):
        raw=io.BytesIO()
        with zipfile.ZipFile(raw,"w") as archive:
            archive.writestr("metadata/results.json",'{"status":"completed"}')
            archive.writestr("data/potential.dx","# synthetic OpenDX integration fixture\n")
        pseudo={"physics_aware_view":{"map_order":["anionic"],"maps":{"anionic":{"patches":[]}},"version":"test-physics"},
            "canonical_view":{"version":"test-original","map_order":["anionic"],"maps":{"anionic":{"patches":[]}}},
            "parameterized_atomic_charge_e":1.0}
        with patch("app.compute_optional_electrostatics",return_value=(pseudo,raw.getvalue())):
            client=TestClient(app)
            res=client.post("/electrostatics/prepare",json={
                "pdb_text":"ATOM  SYNTHETIC TEST PDB RECORD FOR HTTP CONTRACT\n","pH":7.4,
                "ionic_mM":150,"temperature_K":298.15,"apbs_enabled":True
            })
            self.assertEqual(res.status_code,200,res.text)
            obj=res.json()
            self.assertEqual(obj["analysis"]["version"],"test-physics")
            self.assertEqual(obj["canonical_analysis"]["version"],"test-original")
            self.assertTrue(obj["download_url"].startswith("/electrostatics/download/"))
            self.assertTrue(obj["potential_url"].startswith("/electrostatics/potential/"))
            archive=client.get(obj["download_url"])
            self.assertEqual(archive.status_code,200)
            with zipfile.ZipFile(io.BytesIO(archive.content)) as z:
                self.assertIn("data/potential.dx",z.namelist())
            potential=client.get(obj["potential_url"])
            self.assertEqual(potential.status_code,200)
            self.assertIn("OpenDX",potential.text)
            self.assertNotIn("physics_aware_view",obj["result"])
            self.assertNotIn("canonical_view",obj["result"])

if __name__=="__main__":unittest.main()
