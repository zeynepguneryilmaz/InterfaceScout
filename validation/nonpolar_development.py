#!/usr/bin/env python3
"""Development diagnostic for InterfaceScout 2.0 nonpolar orientation energy.

The systems below have already been inspected during model development and are
therefore NOT held-out validation. They are used only to check that the
established-physics implementation is numerically finite, deterministic, and
physically selective on already-known development systems.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "backend"
sys.path.insert(0, str(BACKEND))

import main as core  # noqa: E402
from structural_context import AnalyzeRequest, EnvParams, prepare_context  # noqa: E402
from nonpolar_energy import scan  # noqa: E402
from structure_repair import repair_for_forcefield  # noqa: E402

core.PDB2PQR = None
core.APBS = None
core.MKDSSP = None


def run_case(pdb_id, chain, anchors, pH=7.4, ionic=150.0, context="selected_chain_legacy"):
    import tempfile, shutil
    work = Path(tempfile.mkdtemp(prefix="is_v2_nonpolar_diag_"))
    try:
        req = AnalyzeRequest(
            pdb_id=pdb_id,
            chain=chain,
            structure_context=context,
            protrusion=False,
            env=EnvParams(pH=pH, ionic=ionic, temp=298.0),
        )
        pdb, _, context_label, _ = prepare_context(req, work)
        fixed = work / "forcefield_ready.pdb"
        repair = repair_for_forcefield(pdb, fixed)
        if repair.get("status") != "ok":
            raise RuntimeError(repair)
        struct, _, _, _ = core.build_surface_residues(fixed, pH)
        result = scan(fixed, struct, pH=pH, n_orientations=512)
        if result.get("status") != "ok":
            raise RuntimeError(result)
        top = result["top_orientations"]
        aset = set(anchors)

        def metrics(o):
            ids = {int(r["res_seq"]) for r in o.get("contact_residues", []) if r.get("chain") == chain}
            hits = sorted(aset & ids)
            return {
                "hits": hits,
                "recall": len(hits) / len(aset),
                "precision": len(hits) / len(ids) if ids else 0.0,
                "n_contact_residues": len(ids),
                "total_kj_mol": o["total_energy_change_kj_mol"],
                "vdw_kj_mol": o["vdw_energy_kj_mol"],
                "solvation_kj_mol": o["solvation_energy_change_kj_mol"],
                "minimum_separation_A": o["minimum_separation_A"],
            }

        top1 = metrics(top[0])
        evaluated = [metrics(o) for o in top]
        best_recovery = max(evaluated, key=lambda x: (x["recall"], x["precision"]))
        union10 = set()
        for o in top[:10]:
            union10.update(int(r["res_seq"]) for r in o.get("contact_residues", []) if r.get("chain") == chain)

        repeat = scan(fixed, struct, pH=pH, n_orientations=512)
        deterministic = (
            repeat.get("status") == "ok"
            and repeat.get("best_energy_change_kj_mol") == result.get("best_energy_change_kj_mol")
            and repeat.get("top_orientations", [{}])[0].get("orientation_index") == top[0].get("orientation_index")
        )

        return {
            "pdb": pdb_id,
            "chain": chain,
            "structure_context": context_label,
            "anchors": anchors,
            "structure_repair": repair,
            "status": result["status"],
            "method": result["method"],
            "best_energy_change_kj_mol": result["best_energy_change_kj_mol"],
            "median_best_energy_change_kj_mol": result["median_best_energy_change_kj_mol"],
            "top1": top1,
            "best_anchor_recovery_among_top20_descriptive": best_recovery,
            "union_top10_recall": len(aset & union10) / len(aset),
            "deterministic": deterministic,
            "finite_energy": all(abs(float(o["total_energy_change_kj_mol"])) < 1e12 for o in top),
        }
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    cases = [
        run_case("1AKI", "A", [67, 68, 69, 70, 71, 81], pH=7.0, ionic=150.0),
        run_case("5H7A", "B", [33, 34, 218, 220, 221], pH=7.0, ionic=20.0),
        run_case("5H7A", "C", [33, 34, 218, 220, 221], pH=7.0, ionic=20.0),
    ]
    ok = all(c["status"] == "ok" and c["deterministic"] and c["finite_energy"] for c in cases)
    report = {
        "status": "PASS" if ok else "FAIL",
        "classification": "development diagnostic; NOT held-out validation",
        "version": "2.0.0-dev",
        "physics": "CHARMM36 vdW + integrated neutral carbon plane + Harrison-style SASA solvation",
        "cases": cases,
    }
    out = Path(__file__).with_name("nonpolar_development_report.json")
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
