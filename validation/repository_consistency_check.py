"""Repository consistency checks that must not depend on optional servers.

Use after edits and in publication CI. This verifies local bookkeeping and
documented scope; it does NOT certify experimental evidence or software output.
"""
from __future__ import annotations
import csv
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def read(path):return (ROOT/path).read_text(encoding="utf-8")
def check():
    manifest=list(csv.DictReader((ROOT/"validation/ISV_active_candidates_76_v3.csv").open(encoding="utf-8",newline="")))
    ids=[r["proposed_pdb_id"].upper() for r in manifest]
    assert len(ids)==len(set(ids))==76, "ISV candidate manifest not exactly 76 unique structures"
    excluded={"3GHG","1QJS","4DN0","6LQH"}
    assert not excluded.intersection(ids), "ISV includes excluded structures"
    assert "1STP" in ids, "ISV excludes conditionally retained 1STP"
    tv=list(csv.DictReader((ROOT/"validation/strict_apo_48_selected_candidates_v2_final.csv").open(encoding="utf-8",newline="")))
    assert len(tv)==48, "TV48 manifest count changed"
    readme=read("README.md")
    active=read("validation/ACTIVE_DATASETS.md")
    validation=read("validation/README.md")
    for path,s in [("README.md",readme),
                   ("validation/ACTIVE_DATASETS.md",active),
                   ("validation/README.md",validation)]:
        assert "76" in s and "48" in s, path+" does not describe both active panels"
    assert "ISV76" in readme and "datasets/README.md" in readme, "README must identify the application panel and data index"
    assert "not 76 accepted independent experimental cases" in validation, "Validation README scientific scope missing"
    assert "76 active" in active, "Authoritative active ISV membership missing"
    for path in [".github/workflows/surfmap-isv76-native.yml",
                 ".github/workflows/pep-patch-1gct-repair-diagnostic.yml",
                 ".github/workflows/pep-patch-failed-five-retry.yml"]:
        assert not (ROOT/path).exists(), "Retired auto-workflow still active: "+path
    for path in [".github/workflows/established-methods-full-panel.yml",
                 ".github/workflows/external-methods-preflight.yml",
                 ".github/workflows/protein-sol-isv76-full.yml"]:
        body=read(path)
        trigger=body.split("permissions:",1)[0]
        assert "workflow_dispatch:" in trigger and "  push:" not in trigger, "Unintended automatic benchmark: "+path
    assert "1.1.0" in read("CITATION.cff")
    assert "1.1.0" in readme
    for path in ("datasets/TV48/protein_manifest.csv", "datasets/ISV76/protein_manifest.csv"):
        source = "validation/strict_apo_48_selected_candidates_v2_final.csv" if "TV48" in path else "validation/ISV_active_candidates_76_v3.csv"
        assert (ROOT/path).read_bytes() == (ROOT/source).read_bytes(), "Data package manifest differs: "+path
    print("REPOSITORY_CONSISTENCY_PASS TV48=48 ISV76=76 retired workflows disabled; publication panel manifests verified")

if __name__=="__main__":check()
