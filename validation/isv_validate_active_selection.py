"""Check authoritative active ISV candidates against archived 80 screen and TV48."""
from __future__ import annotations
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def read(name):
    with (ROOT/"validation"/name).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))
def main():
    active=read("ISV_active_candidates_76_v3.csv")
    archived=read("ISV_protein_candidate_pool_80_v2.csv")
    technical=read("strict_apo_48_selected_candidates_v2_final.csv")
    excluded={"3GHG","1QJS","4DN0","6LQH"}
    a={x["proposed_pdb_id"] for x in active}
    old={x["proposed_pdb_id"] for x in archived}
    tv={x["pdb_id"] for x in technical}
    assert len(active)==len(a)==76,(len(active),len(a))
    assert len(archived)==len(old)==80
    assert old-a==excluded,(old-a)
    assert "1STP" in a
    assert not (a&tv),sorted(a&tv)
    assert sum(x["category"].startswith("biofilm_") for x in active)==10
    print("ISV76_ACTIVE_SELECTION_PASS: 76 distinct active; four definitive exclusions; 1STP retained; 10 biofilm; no exact TV48 PDB overlap.")
if __name__=="__main__": main()
