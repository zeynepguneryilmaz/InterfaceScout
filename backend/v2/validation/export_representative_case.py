"""Export an illustrative InterfaceScout map for a representative protein.

The case is illustrative, not an external-validation case: equine serum albumin
(PDB 4F5U, chain A), pH 7.4, generic anionic interface chemistry.
"""
from __future__ import annotations

import json
import urllib.request

import main as v1
from v2.chemistry_freeze import apply_publication_chemistry
from v2.interface_engine import analyze_interface_v2
from v2.model_settings import MULTISCALE_RADII_A
from v2.prepare import prepare_pdb_text

PDB_ID = "4F5U"
CHAIN = "A"
PH = 7.4
CHEMISTRY = "anionic"


def main() -> None:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{PDB_ID}.pdb", timeout=30) as r:
        raw = r.read().decode("utf-8", errors="replace")
    prepared, _ = prepare_pdb_text(raw, chain=CHAIN)
    apply_publication_chemistry(v1)
    v1_out = v1.analyze(v1.AnalyzeRequest(
        pdb_text=prepared,
        chain=None,
        env=v1.EnvParams(pH=PH, ionic=150.0, temp=298.0),
    ))
    chem = v1_out["chemistries"][CHEMISTRY]
    residue_by_key = {str(r["key"]): r for r in v1_out["surface_residues"]}
    chem_by_key = {str(r["key"]): r for r in chem["residues"]}
    inner, outer = MULTISCALE_RADII_A
    inner_tag = f"{int(inner)}A"
    outer_tag = f"{int(outer)}A"

    surface_rows = []
    for key, r in residue_by_key.items():
        c = chem_by_key.get(key, {})
        surface_rows.append({
            "key": key,
            "res_name": r["res_name"],
            "res_seq": r["res_seq"],
            "x": r["x"], "y": r["y"], "z": r["z"],
            "scrsa": r["scrsa"],
            "local_score": c.get("local_score", 0.0),
            "propensity": c.get("propensity", 0.0),
            f"density_{inner_tag}": c.get(f"patch_density_{inner_tag}", 0.0),
            f"density_{outer_tag}": c.get(f"patch_density_{outer_tag}", 0.0),
            "persistence": c.get("multiscale_persistence", 0.0),
        })

    v2_out = analyze_interface_v2(surface="anionic", pdb_text=raw, chain=CHAIN, pH=PH)
    payload = {
        "case": {"pdb_id": PDB_ID, "chain": CHAIN, "pH": PH, "surface_mode": "anionic"},
        "method": {"multiscale_radii_A": list(MULTISCALE_RADII_A), "coarse_patch_radius_A": v2_out["method"]["patch_radius_A"]},
        "summary": {
            "n_residues": v1_out["stats"]["n_residues"],
            "n_surface_residues": v1_out["stats"]["n_surface_res"],
            "n_favorable_residues": chem["n_favorable_residues"],
            "n_coarse_patches": v2_out["n_patches"],
            "n_pareto_primary_patches": v2_out["n_pareto_primary_patches"],
        },
        "surface_residues": surface_rows,
        "top_hotspots": chem.get("top_patches", [])[:10],
        "coarse_patches": v2_out["patches"],
    }
    with open("representative_4F5U_anionic.json", "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    print("REPRESENTATIVE_SUMMARY " + json.dumps(payload["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
