#!/usr/bin/env python3
"""Development validation for InterfaceScout v5.2 structural corrections.

This is deliberately a development/regression test, not held-out adsorption
validation.  It checks that the frozen v5.1 score is unchanged in legacy mode,
that assembly-aware context changes only structural inputs as intended, that CX
is auxiliary, and that material profiles do not combine or alter chemistry
channels.
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
import v52_app as v52  # noqa: E402

# Keep regression tests independent of optional external binaries.
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
            ma = amap(ga.get(collection, []), field)
            mb = amap(gb.get(collection, []), field)
            checks.append((chem, collection, field, ma == mb, len(ma), len(mb)))
    return checks


def request_core(pdb_id, chain=None, ph=7.4):
    return core.analyze(core.AnalyzeRequest(pdb_id=pdb_id, chain=chain, env=core.EnvParams(pH=ph, ionic=150, temp=298)))


def request_v52(pdb_id, chain=None, context="auto", protrusion=True, profile=None, ph=7.4):
    return v52.analyze(v52.AnalyzeRequest(
        pdb_id=pdb_id,
        chain=chain,
        structure_context=context,
        protrusion=protrusion,
        material_profile=profile,
        env=v52.EnvParams(pH=ph, ionic=150, temp=298),
    ))


def assert_true(name, value, detail=None, out=None):
    rec = {"check": name, "pass": bool(value)}
    if detail is not None:
        rec["detail"] = detail
    out.append(rec)
    if not value:
        raise AssertionError(f"{name}: {detail}")


def main():
    out = []

    # 1) Frozen-core regression: legacy selected-chain mode must reproduce v5.1.
    for pdb_id, chain, ph in [("1MBN", "A", 7.4), ("4F5S", "A", 7.4)]:
        old = request_core(pdb_id, chain, ph)
        new = request_v52(pdb_id, chain, context="selected_chain_legacy", protrusion=False, ph=ph)
        checks = compare_chemistry(old, new)
        bad = [x for x in checks if not x[3]]
        assert_true(
            f"legacy regression {pdb_id} chain {chain}",
            not bad,
            {"comparisons": len(checks), "mismatches": bad[:5]},
            out,
        )
        assert_true(
            f"legacy residue counts {pdb_id}",
            old["stats"]["n_residues"] == new["stats"]["n_residues"] and old["stats"]["n_surface_res"] == new["stats"]["n_surface_res"],
            {"old": old["stats"], "new": new["stats"]},
            out,
        )

    # 2) CX is auxiliary: toggling it must not alter any canonical chemistry map.
    cx_off = request_v52("1MBN", "A", context="selected_chain_legacy", protrusion=False)
    cx_on = request_v52("1MBN", "A", context="selected_chain_legacy", protrusion=True)
    cx_checks = compare_chemistry(cx_off, cx_on)
    cx_bad = [x for x in cx_checks if not x[3]]
    assert_true("CX score invariance", not cx_bad, {"mismatches": cx_bad[:5]}, out)
    cx_vals = [r.get("cx_sidechain_mean") for r in cx_on["all_residues"] if r.get("cx_sidechain_mean") is not None]
    assert_true(
        "CX finite nonnegative descriptors",
        bool(cx_vals) and all(math.isfinite(float(v)) and float(v) >= 0 for v in cx_vals),
        {"n": len(cx_vals), "min": min(cx_vals) if cx_vals else None, "max": max(cx_vals) if cx_vals else None},
        out,
    )

    # 3) Biological-assembly context: report chain A but compute structural
    # context on the tetramer.  The assembly must contain more residues and at
    # least some chain-A scRSA values should differ from isolated-chain mode.
    hb_iso = request_v52("2HHB", "A", context="selected_chain_legacy", protrusion=True)
    hb_asm = request_v52("2HHB", "A", context="auto", protrusion=True)
    assert_true(
        "2HHB auto resolves biological assembly",
        hb_asm["settings"]["structure_context"] == "biological_assembly_1",
        hb_asm["settings"],
        out,
    )
    assert_true(
        "2HHB assembly context larger than reported chain",
        hb_asm["stats"]["n_residues_context"] > hb_asm["stats"]["n_residues_reported"],
        hb_asm["stats"],
        out,
    )
    assert_true(
        "2HHB report remains chain A",
        all(r["chain"] == "A" for r in hb_asm["all_residues"]),
        {"reported_chains": sorted(set(r["chain"] for r in hb_asm["all_residues"]))},
        out,
    )
    iso_sc = {r["key"]: r["scrsa"] for r in hb_iso["all_residues"]}
    asm_sc = {r["key"]: r["scrsa"] for r in hb_asm["all_residues"]}
    common = sorted(set(iso_sc) & set(asm_sc))
    changed = [k for k in common if abs(float(iso_sc[k]) - float(asm_sc[k])) > 1e-5]
    assert_true(
        "2HHB assembly changes chain-A exposure where shielding exists",
        bool(changed),
        {"n_common": len(common), "n_changed": len(changed), "examples": changed[:10]},
        out,
    )

    # 4) Material profile must be a view over existing channels, not a score mix.
    prof = request_v52("1MBN", "A", context="selected_chain_legacy", protrusion=False, profile="graphitic_carbon")
    p = prof["material_profile"]
    assert_true(
        "material profile channels declared",
        p is not None and p["channels"] == ["pi_carbon", "hydrophobic"],
        p,
        out,
    )
    identical = True
    for ch in p["channels"]:
        identical = identical and p["channel_results"][ch] == prof["chemistries"][ch]
    assert_true("material profile does not alter channel results", identical, p.get("combination_rule"), out)
    assert_true("material profile has no weighted combination", p.get("combination_rule", "").startswith("none"), p, out)

    # 5) Applicability panel should reflect active/inactive modules.
    assert_true(
        "applicability notes emitted",
        bool(prof.get("applicability", {}).get("included_in_this_run")) and bool(prof.get("applicability", {}).get("not_included_or_interpretation_limits")),
        prof.get("applicability"),
        out,
    )

    report = {
        "status": "PASS",
        "version": v52.APP_VERSION,
        "core_version": core.APP_VERSION,
        "checks_passed": sum(1 for r in out if r["pass"]),
        "checks_total": len(out),
        "checks": out,
        "note": "Development/regression validation only; not held-out adsorption validation.",
    }
    target = Path(__file__).with_name("v52_structural_validation_report.json")
    target.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
