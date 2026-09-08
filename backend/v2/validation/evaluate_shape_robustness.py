"""Adsorption-label-independent robustness audit across structurally diverse proteins.

This audit does not use experimental interface labels and does not retune model
parameters. It asks whether the frozen scRSA surface definition and top-hotspot
ordering remain numerically stable outside the single developmental albumin.
"""
from __future__ import annotations

import json
import urllib.request
from statistics import median

import main as v1
from v2.prepare import prepare_pdb_text

STRUCTURES = [
    {"id": "gb3", "pdb_id": "2OED", "chain": "A", "class": "small compact"},
    {"id": "hcaii", "pdb_id": "2CBA", "chain": "A", "class": "medium globular"},
    {"id": "esa", "pdb_id": "4F5U", "chain": "A", "class": "large globular"},
    {"id": "transferrin", "pdb_id": "1SUV", "chain": "C,E", "class": "multidomain/multichain structural proxy"},
]
THRESHOLDS = [0.02, 0.05, 0.075, 0.10, 0.15]


def fetch_pdb(pdb_id: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=30) as r:
        return r.read().decode("utf-8", errors="replace")


def run_case(raw: str, chain: str, *, points: int, threshold: float) -> dict:
    prepared, _ = prepare_pdb_text(raw, chain=chain)
    old_points, old_thr = v1.SASA_POINTS, v1.SC_RSA_THRESHOLD
    try:
        v1.SASA_POINTS = int(points)
        v1.SC_RSA_THRESHOLD = float(threshold)
        req = v1.AnalyzeRequest(pdb_text=prepared, chain=None, env=v1.EnvParams(pH=7.4, ionic=150.0, temp=298.0))
        out = v1.analyze(req)
    finally:
        v1.SASA_POINTS, v1.SC_RSA_THRESHOLD = old_points, old_thr
    return out


def jaccard(a: set[str], b: set[str]) -> float:
    u = a | b
    return len(a & b) / len(u) if u else 1.0


def top10_sets(result: dict) -> dict[str, set[str]]:
    out = {}
    for chemistry, payload in result.get("chemistries", {}).items():
        out[chemistry] = {str(x.get("center_key")) for x in payload.get("top_patches", [])[:10] if x.get("center_key")}
    return out


def main() -> None:
    rows = []
    for spec in STRUCTURES:
        raw = fetch_pdb(spec["pdb_id"])
        base = run_case(raw, spec["chain"], points=200, threshold=0.05)
        hi = run_case(raw, spec["chain"], points=500, threshold=0.05)
        surf_base = {str(r["key"]) for r in base.get("surface_residues", [])}
        surf_hi = {str(r["key"]) for r in hi.get("surface_residues", [])}
        base_top = top10_sets(base)

        threshold_medians = {}
        for thr in THRESHOLDS:
            cur = base if thr == 0.05 else run_case(raw, spec["chain"], points=200, threshold=thr)
            cur_top = top10_sets(cur)
            vals = [jaccard(base_top.get(c, set()), cur_top.get(c, set())) for c in sorted(base_top)]
            threshold_medians[str(thr)] = float(median(vals)) if vals else None

        rows.append({
            **spec,
            "n_surface_200": len(surf_base),
            "n_surface_500": len(surf_hi),
            "surface_set_jaccard_200_vs_500": jaccard(surf_base, surf_hi),
            "surface_threshold_crossing_count_200_vs_500": len(surf_base ^ surf_hi),
            "median_top10_hotspot_jaccard_vs_threshold_005": threshold_medians,
        })

    payload = {
        "purpose": "label-independent multi-structure robustness audit; no model tuning",
        "structures": rows,
    }
    with open("shape_robustness.json", "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
