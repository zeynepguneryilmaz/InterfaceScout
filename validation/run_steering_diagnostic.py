from __future__ import annotations

import json
import sys
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend import main as IS
from backend.steering import compute_electrostatic_steering

OUT = ROOT / "validation_results"
OUT.mkdir(exist_ok=True)


def fetch(pid: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pid}.pdb", timeout=60) as r:
        return r.read().decode("utf-8")


def analyze(pid: str, chain: str | None, pH: float, ionic: float, temp: float):
    # Diagnostic deliberately bypasses APBS because the steering module uses the
    # independent Debye-Huckel charged-plane formulation rather than the auxiliary
    # APBS descriptor used elsewhere in InterfaceScout.
    req = IS.AnalyzeRequest(pdb_text=fetch(pid), chain=chain,
                            env=IS.EnvParams(pH=pH, ionic=ionic, temp=temp))
    with tempfile.TemporaryDirectory(prefix="interfacescout_steering_") as td:
        pdb, _ = IS.prepare_input_pdb(req, Path(td))
        _, all_residues, _, _ = IS.build_surface_residues(pdb, pH)
    surface = [r for r in all_residues if r["surface_exposed"]]
    steering = compute_electrostatic_steering(
        all_residues, surface, ionic, temp,
        n_orientations=4096, footprint_depths_A=(5.0, 8.0, 10.0),
    )
    return {"all_residues": all_residues, "surface_residues": surface}, steering


def footprint_nums(steering, surfkey, depth):
    rows = steering[surfkey]["surface_facing_footprint"][f"within_{int(depth)}A"]
    return {(r["chain"], int(r["res_seq"])) for r in rows}


def anchor_recovery(anchors, footprint):
    hit = [a for a in anchors if a in footprint]
    return {"n": len(anchors), "n_hit": len(hit), "recall": len(hit)/len(anchors) if anchors else None,
            "hits": hit, "misses": [a for a in anchors if a not in footprint]}


def main():
    report = {"model": "v5.2 electrostatic steering diagnostic", "cases": []}

    # Development-only SpA cases; these were already inspected during v5.1 diagnostics.
    for chain in ["B", "C"]:
        result, st = analyze("5H7A", chain, pH=7.0, ionic=20.0, temp=298.0)
        case = {"case": f"SpA_5H7A_chain_{chain}", "net_charge": st.get("net_charge_descriptor_e"),
                "debye_A": st.get("debye_length_A"), "positive_surface": {}, "negative_surface": {}}
        pos_anchors = [(chain, 219), (chain, 220), (chain, 221)]
        neg_anchors = [(chain, 33), (chain, 34), (chain, 35), (chain, 36), (chain, 37)]
        for d in [5,8,10]:
            case["positive_surface"][f"R{d}"] = anchor_recovery(pos_anchors, footprint_nums(st, "positive_surface", d))
            case["negative_surface"][f"R{d}"] = anchor_recovery(neg_anchors, footprint_nums(st, "negative_surface", d))
        case["positive_best_reduced_energy"] = st["positive_surface"]["best_reduced_energy"]
        case["negative_best_reduced_energy"] = st["negative_surface"]["best_reduced_energy"]
        case["positive_top_orientation"] = st["positive_surface"]["best_normal_plane_to_protein"]
        case["negative_top_orientation"] = st["negative_surface"]["best_normal_plane_to_protein"]
        report["cases"].append(case)

    # Whole 2HHB deposited tetramer diagnostic: no chain isolation.
    result, st = analyze("2HHB", None, pH=7.4, ionic=150.0, temp=310.0)
    hb = {"case": "2HHB_whole_tetramer_citrate_like_negative_surface",
          "net_charge": st.get("net_charge_descriptor_e"), "debye_A": st.get("debye_length_A"),
          "negative_surface": {}}
    hb_anchors = (
        [("A",i) for i in list(range(12,26))+list(range(61,82))] +
        [("C",i) for i in list(range(12,26))+list(range(61,82))] +
        [("B",i) for i in range(51,54)] +
        [("D",i) for i in range(45,54)]
    )
    for d in [5,8,10]:
        hb["negative_surface"][f"R{d}"] = anchor_recovery(hb_anchors, footprint_nums(st, "negative_surface", d))
    hb["negative_best_reduced_energy"] = st["negative_surface"]["best_reduced_energy"]
    hb["negative_top_orientation"] = st["negative_surface"]["best_normal_plane_to_protein"]
    report["cases"].append(hb)

    (OUT / "steering_diagnostic.json").write_text(json.dumps(report, indent=2))

    lines = ["# InterfaceScout v5.2 electrostatic steering diagnostic", "",
             "Development cases only; no held-out validation claim is made.", ""]
    for c in report["cases"]:
        lines.append(f"## {c['case']}")
        lines.append(f"- net charge descriptor: {c.get('net_charge')}")
        lines.append(f"- Debye length (Å): {c.get('debye_A')}")
        for side in ["positive_surface", "negative_surface"]:
            if side in c and isinstance(c[side], dict):
                for k,v in c[side].items():
                    if isinstance(v, dict) and "recall" in v:
                        lines.append(f"- {side} {k}: recall={v['recall']:.3f} ({v['n_hit']}/{v['n']})")
        lines.append("")
    (OUT / "STEERING_DIAGNOSTIC.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
