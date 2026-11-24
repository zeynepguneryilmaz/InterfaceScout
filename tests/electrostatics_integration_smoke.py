"""Real APBS / PDB2PQR integration smoke; requires both executables on PATH.

Use RCSB 1UBQ as a nontrivial experimental structure. Fail on any
electrostatics error and assert the canonical maps remain identical.
"""
import json
from pathlib import Path
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
from interface_engine import analyze_all_maps
from electrostatics import availability

if __name__ == "__main__":
    status = availability()
    if not status["pdb2pqr"] or not status["apbs"]:
        raise SystemExit("Missing required PDB2PQR or APBS binary")
    with urllib.request.urlopen("https://files.rcsb.org/download/1UBQ.pdb", timeout=40) as response:
        structure = response.read().decode("utf-8")
    opts = dict(pdb_text=structure, chain="A", pH=7.4, ionic_mM=150.0, temp_K=298.0)
    baseline = analyze_all_maps(**opts)
    enriched = analyze_all_maps(**opts, electrostatics=True)
    assert baseline["n_maps"] == enriched["n_maps"] == 11
    assert baseline["maps"] == enriched["maps"], "Canonical maps changed!"
    assert baseline["residue_profiles"] == enriched["residue_profiles"]
    result = enriched["electrostatics"]
    assert result["status"] == "computed"
    assert result["n_matched_residues"] == len(enriched["residue_profiles"])
    assert any(r["sidechain_mean_potential_kT_per_e"] is not None for r in result["residues"])
    surface = [r for r in result["residues"] if r["surface_exposed"]]
    shell = [r for r in surface if r["surface_potential_kT_per_e"] is not None]
    assert len(shell) >= 0.90 * len(surface), "Unexpectedly low SAS coverage"
    assert result["n_surface_potential_residues"] == len(shell)
    assert result["surface_sampling"]["n_exposed_residues_with_surface_potential"] == len(shell)
    assert result["surface_sampling"]["surface_type"] == "approximate SAS shell, not molecular SES"
    assert result["surface_sampling"]["probe_A"] == 1.4
    assert result["surface_sampling"]["outward_offset_A"] == 0.25
    gly = [r for r in shell if r["res_name"] == "GLY"]
    assert gly and all(r["sidechain_mean_potential_kT_per_e"] is None for r in gly), (
        "Exposed glycine must have SAS potential but no heavy-atom side-chain comparator"
    )

    # Retain auditable per-residue evidence from the actual external solvers.
    output = ROOT / "dist" / "electrostatics"
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "protein": "1UBQ", "status": "PASS", "input": enriched["input"],
        "canonical_maps_unchanged": True,
        "n_canonical_maps": enriched["n_maps"],
        "electrostatics": result,
    }
    (output / "1UBQ_APBS_PROPKA_electrostatics.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "PASS", "protein": "1UBQ", "maps": enriched["n_maps"],
                      "residues": result["n_matched_residues"],
                      "pqr_atoms": result["n_pqr_atoms"]}, indent=2))
