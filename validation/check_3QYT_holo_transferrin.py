"""Verify 3QYT mutant-free holo transferrin alternative with actual mmCIF
and same 11-map canonical engine; no premature publication claim.
"""
import json
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import selected_10_new_structures_atom_qa as qa
import audit_literature_33_engine_20261122 as engine

OUT=Path("holo_transferrin_3QYT_validation")
OUT.mkdir(exist_ok=True)
qa.OUT=OUT
engine.OUT=OUT
pid="3QYT"
atom=qa.audit(pid)
print("ATOM_QA",pid,atom["status"],atom["missing_backbone_residues"],atom["missing_sidechain_residues"],atom["altloc_residues"],flush=True)
engine.prepare(pid)
maps=engine.work(pid)
print("MAPS",pid,maps["n_maps"],maps["n_residues"],maps["n_surface"],flush=True)
result={"pdb":pid,"atom_quality":atom,"interface_engine":maps,
        "limitations":"Technical 11-map chemistry; Fe, carbonate and glycan may be removed by protein-only prep. Assay-state correspondence not demonstrated."}
(OUT/"3QYT_comparison.json").write_text(json.dumps(result,indent=2))
assert atom["status"]=="COORDINATES_AUDITED"
assert maps["n_maps"]==11
