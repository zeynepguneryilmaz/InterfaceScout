"""Actual solver regression on an experimental structure; no mocked success."""
import sys,json,urllib.request
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"backend"))
from interface_engine import analyze_all_maps
from electrostatics import ElectrostaticsError
from app import AnalyzeRequest,StructureSetRequest
assert AnalyzeRequest().pb_solver=="npbe"
with urllib.request.urlopen("https://files.rcsb.org/download/1UBQ.pdb",timeout=60) as u:
    pdb=u.read().decode()
kw=dict(pdb_text=pdb,chain="A",pH=7.4,ionic_mM=150.,temp_K=298.)
base=analyze_all_maps(**kw)
results={}
for mode in ("npbe","lpbe","both"):
    opts={} if mode=="npbe" else {"pb_solver":mode}
    x=analyze_all_maps(**kw,electrostatics=True,**opts)
    assert x["maps"]==base["maps"] and x["residue_profiles"]==base["residue_profiles"]
    e=x["electrostatics"];assert e["status"]=="computed",e
    assert e["solver"]==("npbe" if mode=="both" else mode)
    assert e["quality"]["status"] in ("checks_passed","review_recommended")
    assert not e["fallback_used"]
    if mode=="both":
        assert e["comparison"]["status"]=="computed",e["comparison"]
        assert e["comparison"]["same_pqr_grid_shell"]
        assert e["comparison"]["n_paired"]>0
        assert any(r["npbe_minus_lpbe_kT_per_e"] is not None for r in e["residues"])
    results[mode]=e
with patch("interface_engine.compute_electrostatics",side_effect=ElectrostaticsError("deliberate test failure")):
    failed=analyze_all_maps(**kw,electrostatics=True)
    assert failed["maps"]==base["maps"] and failed["residue_profiles"]==base["residue_profiles"]
    assert failed["electrostatics"]["status"]=="unavailable" and not failed["electrostatics"]["fallback_used"]
results["failure"]=failed["electrostatics"]
out=ROOT/"dist"/"npbe_regression";out.mkdir(parents=True,exist_ok=True)
(out/"1UBQ_all_modes.json").write_text(json.dumps(results,indent=2,allow_nan=False))
print("PASS actual default NPBE, LPBE, paired modes and failure isolation")
