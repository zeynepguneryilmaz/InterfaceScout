"""Compare apo alpha-lactalbumin monomer with whole deposited asymmetric-unit chains,
and holo alpha-lactalbumin monomer using the same canonical 11-map engine.
"""
import csv
import json
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import audit_literature_33_engine_20261122 as engine
sys.path.insert(0,str(HERE.parent/"backend"))
from prepare import prepare_pdb_text
from interface_engine import analyze_all_maps

OUT=Path("alpha_lactalbumin_chain_comparison")
OUT.mkdir(exist_ok=True)
engine.OUT=OUT
records=[]
for pid in ["1F6R","1F6S"]:
    meta=engine.prepare(pid)
    raw=(OUT/f"{pid}_protein_only.pdb").read_text()
    for selection in ["ALL","A"]:
        prepared,prepreport=prepare_pdb_text(raw,chain=None if selection=="ALL" else "A")
        result=analyze_all_maps(pdb_text=prepared,chain=None,pH=7.4,ionic_mM=150,temp_K=298)
        maps=result["maps"]
        order=result["map_order"]
        row={"pdb":pid,"chain_selection":selection,"n_maps":len(order),
             "n_chains_input":meta["n_chains"],
             "n_residues":len(result["residue_profiles"]),
             "n_surface":result["diagnostics"]["n_surface_residues"],
             "n_patches_per_map":{x:len(maps[x].get("patches",[])) for x in order},
             "status":"PASS_11_MAPS" if len(order)==11 else "FAIL"}
        records.append(row)
        print("RESULT",json.dumps(row),flush=True)
(OUT/"comparison.json").write_text(json.dumps(records,indent=2))
with (OUT/"comparison.csv").open("w",newline="") as file:
    cols=list(records[0])
    writer=csv.DictWriter(file,cols);writer.writeheader()
    for rec in records:
        writer.writerow({k:json.dumps(v) if isinstance(v,dict) else v for k,v in rec.items()})
assert len(records)==4 and all(r["status"]=="PASS_11_MAPS" for r in records)
