#!/usr/bin/env python3
"""Regression validation for InterfaceScout 2.0 development.

Checks that the frozen InterfaceScout 1.0 canonical score is unchanged in
legacy structural mode and that new 2.0 structural descriptors remain auxiliary.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "backend"
sys.path.insert(0, str(BACKEND))

import main as core  # noqa: E402
import structural_context as structural  # noqa: E402

core.PDB2PQR = None
core.APBS = None
core.MKDSSP = None


def amap(rows, field):
    return {r.get("key") or r.get("center_key"): r.get(field) for r in rows}


def compare_chemistry(a, b):
    checks = []
    for chem in a["chemistry_list"]:
        ga, gb = a["chemistries"][chem], b["chemistries"][chem]
        for collection, field in [
            ("residues", "local_score"),
            ("residues", "propensity"),
            ("patch_centers", "multiscale_persistence"),
            ("repulsive_residues", "repulsion_propensity"),
        ]:
            checks.append(amap(ga.get(collection, []), field) == amap(gb.get(collection, []), field))
    return all(checks)


def request_core(pdb_id, chain=None, ph=7.4):
    return core.analyze(core.AnalyzeRequest(pdb_id=pdb_id, chain=chain, env=core.EnvParams(pH=ph, ionic=150, temp=298)))


def request_v2(pdb_id, chain=None, context="auto", protrusion=True, profile=None, ph=7.4):
    return structural.analyze_structural(structural.AnalyzeRequest(
        pdb_id=pdb_id,
        chain=chain,
        structure_context=context,
        protrusion=protrusion,
        material_profile=profile,
        env=structural.EnvParams(pH=ph, ionic=150, temp=298),
    ))


def check(name, ok, detail, records):
    records.append({"check": name, "pass": bool(ok), "detail": detail})
    if not ok:
        raise AssertionError(f"{name}: {detail}")


def main():
    records = []

    for pdb_id, chain in [("1MBN", "A"), ("4F5S", "A")]:
        old = request_core(pdb_id, chain)
        new = request_v2(pdb_id, chain, context="selected_chain_legacy", protrusion=False)
        check(f"1.0 score regression {pdb_id}", compare_chemistry(old, new), "canonical maps identical", records)
        check(
            f"1.0 residue-count regression {pdb_id}",
            old["stats"]["n_residues"] == new["stats"]["n_residues"] and old["stats"]["n_surface_res"] == new["stats"]["n_surface_res"],
            {"v1": [old["stats"]["n_residues"], old["stats"]["n_surface_res"]], "v2_legacy": [new["stats"]["n_residues"], new["stats"]["n_surface_res"]]},
            records,
        )

    cx_off = request_v2("1MBN", "A", context="selected_chain_legacy", protrusion=False)
    cx_on = request_v2("1MBN", "A", context="selected_chain_legacy", protrusion=True)
    check("CX does not alter canonical score", compare_chemistry(cx_off, cx_on), "canonical maps identical", records)
    cx_vals = [r.get("cx_sidechain_mean") for r in cx_on["all_residues"] if r.get("cx_sidechain_mean") is not None]
    check("CX finite", bool(cx_vals) and all(math.isfinite(float(v)) and float(v) >= 0 for v in cx_vals), {"n": len(cx_vals)}, records)

    hb_iso = request_v2("2HHB", "A", context="selected_chain_legacy")
    hb_asm = request_v2("2HHB", "A", context="auto")
    check("2HHB uses biological assembly", hb_asm["settings"]["structure_context"] == "biological_assembly_1", hb_asm["settings"]["structure_context"], records)
    iso_sc = {r["key"]: r["scrsa"] for r in hb_iso["all_residues"]}
    asm_sc = {r["key"]: r["scrsa"] for r in hb_asm["all_residues"]}
    changed = [k for k in set(iso_sc) & set(asm_sc) if abs(float(iso_sc[k]) - float(asm_sc[k])) > 1e-5]
    check("assembly context changes shielded exposure", bool(changed), {"n_changed": len(changed)}, records)

    prof = request_v2("1MBN", "A", context="selected_chain_legacy", protrusion=False, profile="graphitic_carbon")
    p = prof["material_profile"]
    check("material profile channels declared", p and p["channels"] == ["pi_carbon", "hydrophobic"], p.get("channels") if p else None, records)
    identical = all(p["channel_results"][ch] == prof["chemistries"][ch] for ch in p["channels"])
    check("material profile does not mix scores", identical and p.get("combination_rule", "").startswith("none"), p.get("combination_rule"), records)

    report = {
        "status": "PASS",
        "version": structural.STRUCTURAL_LAYER_VERSION,
        "frozen_reference": structural.CORE_RELEASE_VERSION,
        "checks_passed": len(records),
        "checks_total": len(records),
        "checks": records,
    }
    target = Path(__file__).with_name("regression_report.json")
    target.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
