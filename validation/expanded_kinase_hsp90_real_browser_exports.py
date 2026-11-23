"""Five-structure evidence-independent InterfaceScout expansion (2026-11-23).

Reuses the unchanged, fully validated 35-structure browser export runner.
Biological assembly selection is mandatory for AMPK and HSP90; the script
retains provenance, modified-residue exclusions, exact export QA and 55
pairwise geometric map combinations for each protein.

The outputs are technical and structural investigations, NOT validated
adsorption-site predictions nor chemically faithful phosphorylation models.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import gemmi

try:
    from validation import curated_biological_units_full_browser_exports as base
except ModuleNotFoundError:
    import curated_biological_units_full_browser_exports as base

PDB_IDS = ["4CFE", "5ISO", "2GS2", "2GS7", "8EOB"]
base.GROUPS = [PDB_IDS]
base.OUT = Path("expanded_kinase_hsp90_outputs")
base.ASSEMBLY.update({"4CFE", "5ISO", "8EOB"})
base.SCOPE_NOTES.update({
    "4CFE": "AMPK alpha2-beta1-gamma1 heterotrimer: phosphorylation TPO and SEP and all cofactors are removed by the unchanged protein-only engine; model cannot quantify phosphorylation electrochemistry.",
    "5ISO": "AMPK alpha2-beta1-gamma1 non-T-loop-phosphorylated heterotrimer; other SEP modified residues removed; compare only mapped, homologous residues with 4CFE.",
    "2GS2": "EGFR intracellular catalytic domain only; active-like crystal conformation; not full receptor or experimental adsorption state.",
    "2GS7": "EGFR inactive-like catalytic domain; 1 reported mutation, AMP-PNP analog omitted; cross-state differences are confounded and not causal activation effects.",
    "8EOB": "Human HSP90-beta CLOSED dimer without AIPL1 partner in this deposited structure; 3.1 A cryo-EM; missing/disordered portions and AMP-PNP omitted from chemistry model.",
})
EXPECTED_CHAINS = {"4CFE": 3, "5ISO": 3, "2GS2": 1, "2GS7": 1, "8EOB": 2}
MIN_RESIDUES = {"4CFE": 700, "5ISO": 700, "2GS2": 200, "2GS7": 200, "8EOB": 900}
old_prepare = base.prepare_case

def checked_prepare(pid: str, folder: Path):
    provenance, prepared = old_prepare(pid, folder)
    n_chain = provenance["n_protein_chains"]
    n_res = provenance["n_modelled_residues"]
    if n_chain != EXPECTED_CHAINS[pid]:
        raise ValueError(f"{pid}: expected {EXPECTED_CHAINS[pid]} biological protein chains; found {n_chain}")
    if n_res < MIN_RESIDUES[pid]:
        raise ValueError(f"{pid}: unexpected residue loss; only {n_res} selected")
    source_path = folder / ("rcsb_assembly1.cif" if pid in base.ASSEMBLY else "rcsb_entry.cif")
    original = gemmi.read_structure(str(source_path))
    modified = {}
    for chain in original[0]:
        for residue in chain.get_polymer():
            if residue.name in ("TPO", "SEP", "PTR", "MSE"):
                modified[residue.name] = modified.get(residue.name, 0) + 1
    provenance["modified_polymer_residues_excluded"] = modified
    provenance["mandatory_biological_chain_QC"] = "PASS"
    provenance["adsorption_validation"] = "NOT_CLAIMED"
    provenance["phosphorylation_chemistry_supported"] = False
    (folder / "source_and_selection.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")
    return provenance, prepared

base.prepare_case = checked_prepare

def main():
    base.main()
    folder = base.OUT / "shard_0"
    pair_rows = []
    for a, b in [("4CFE", "5ISO"), ("2GS2", "2GS7")]:
        row = base.compare_pair(a, b, folder)
        row["interpretation"] = (
            "matched-residue descriptive structural sensitivity only; mutation, "
            "missing modified residues, bound cofactors and model differences "
            "preclude a causal state-chemistry or adsorption claim")
        pair_rows.append(row)
    (folder / "new_state_pair_comparisons.json").write_text(
        json.dumps(pair_rows, indent=2, ensure_ascii=False), encoding="utf-8")
    print("EXPANSION_PAIR_STATUS", json.dumps(pair_rows, ensure_ascii=False), flush=True)

if __name__ == "__main__":
    main()
