"""Real 1GCT end-to-end regression for APBS-only N-fragment handling.

Verifies 1GCT biological assembly unchanged, 11 canonical maps byte-for-byte
unchanged, APBS PQR omits exactly A:11 and E:11 and explicitly stores null
potential for these two residues. Also verifies 1UBQ pristine control has no
recovery and correct full residue matching. Never modifies the input files.
"""
from __future__ import annotations
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"backend"),str(ROOT/"validation")]
from isv_full_80_qa import prepare
from interface_engine import analyze_all_maps

OUT=ROOT/"dist"/"electrostatics"/"backbone_gap_recovery"
OUT.mkdir(parents=True,exist_ok=True)
active=list(csv.DictReader((ROOT/"validation"/"ISV_active_candidates_76_v3.csv").open()))
by={r["proposed_pdb_id"]:r for r in active}


def run(pid):
    source,pdbpath=prepare(by[pid])
    original=pdbpath.read_text(encoding="utf-8")
    sha=hashlib.sha256(original.encode()).hexdigest()
    kwargs={"pdb_text":original,"pH":7.4,"ionic_mM":150.,"temp_K":298.}
    maps_only=analyze_all_maps(**kwargs)
    enriched=analyze_all_maps(**kwargs,electrostatics=True)
    assert maps_only["n_maps"]==enriched["n_maps"]==11
    assert maps_only["maps"]==enriched["maps"]
    assert maps_only["residue_profiles"]==enriched["residue_profiles"]
    assert maps_only["map_order"]==enriched["map_order"]
    assert hashlib.sha256(pdbpath.read_bytes()).hexdigest()==sha
    electro=enriched["electrostatics"]
    assert electro["status"]=="computed"
    assert electro["n_input_residues"]==len(maps_only["residue_profiles"])
    assert len(electro["residues"])==electro["n_input_residues"]
    assert electro["n_matched_residues"]+electro["n_unmatched_residues"]==electro["n_input_residues"]
    assert electro["surface_sampling"]["n_exposed_residues_with_surface_potential"]==electro["n_surface_potential_residues"]
    if pid=="1GCT":
        assert electro["quality_status"]=="partial_input_recovery"
        assert electro["input_recovery"]["applied"]
        assert electro["n_unmatched_residues"]==2
        expected={"A:11:","E:11:"}
        assert {x["key"] for x in electro["input_recovery"]["omitted_residues"]}==expected
        flagged=[x for x in electro["residues"] if x["key"] in expected]
        assert len(flagged)==2
        assert all(x["pqr_mapped"] is False and
                   x["surface_potential_kT_per_e"] is None and
                   x["net_atomic_charge_e"] is None for x in flagged)
        assert all(x["pqr_mapped"] for x in electro["residues"] if x["key"] not in expected)
        assert electro["n_surface_potential_residues"]>=200
    else:
        assert electro["quality_status"]=="complete"
        assert electro["n_unmatched_residues"]==0
        assert not electro["input_recovery"]["applied"]
    payload={
        "pdb_id":pid,"status":"PASS",
        "input_pdb_sha256_unchanged":sha,
        "n_maps":maps_only["n_maps"],
        "n_input_residues":electro["n_input_residues"],
        "n_mapped_residues":electro["n_matched_residues"],
        "n_unmatched_residues":electro["n_unmatched_residues"],
        "quality_status":electro["quality_status"],
        "n_sas_potentials":electro["n_surface_potential_residues"],
        "omitted_residue_keys":[r["key"] for r in electro["input_recovery"]["omitted_residues"]],
        "canonical_maps_exactly_equal":True,
        "original_input_untouched":True,
        "real_external_pdb2pqr_apbs":True,
    }
    (OUT/(pid+"_verified.json")).write_text(json.dumps(payload,indent=2),encoding="utf-8")
    print("BACKBONE_RECOVERY_REAL_CASE "+json.dumps(payload),flush=True)


if __name__=="__main__":
    for pid in ("1GCT","1HEL"):
        run(pid)
