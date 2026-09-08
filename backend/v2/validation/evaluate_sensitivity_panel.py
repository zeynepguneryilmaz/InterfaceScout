"""Label-independent sensitivity analysis across a 14-protein structural panel.

The panel spans plasma/corona-relevant proteins, high-resolution adsorption
models, and proteins with diverse sizes/functions. Experimental adsorption
labels are not used here. The analysis asks whether InterfaceScout outputs are
stable to reasonable numerical/model perturbations.
"""
from __future__ import annotations

import json
import math
import urllib.request
from statistics import median

import numpy as np

import main as v1
from v2.chemistry_freeze import apply_publication_chemistry
from v2.prepare import prepare_pdb_text

STRUCTURES = [
    {"id": "albumin", "label": "Bovine serum albumin", "pdb_id": "4F5S", "chain": "A", "group": "plasma_corona"},
    {"id": "fibrinogen", "label": "Human fibrinogen alpha-beta-gamma half-unit", "pdb_id": "3GHG", "chain": "A,B,C", "group": "plasma_corona"},
    {"id": "igg", "label": "Intact IgG2a", "pdb_id": "1IGT", "chain": "A,B,C,D", "group": "plasma_corona"},
    {"id": "transferrin", "label": "Human serum transferrin", "pdb_id": "3QYT", "chain": "A", "group": "plasma_corona"},
    {"id": "lysozyme", "label": "Hen egg-white lysozyme", "pdb_id": "1AKI", "chain": "A", "group": "high_resolution_adsorption_model"},
    {"id": "rnase_a", "label": "Bovine ribonuclease A", "pdb_id": "7RSA", "chain": "A", "group": "high_resolution_adsorption_model"},
    {"id": "myoglobin", "label": "Sperm whale myoglobin", "pdb_id": "1MBN", "chain": "A", "group": "high_resolution_adsorption_model"},
    {"id": "hcaii", "label": "Human carbonic anhydrase II", "pdb_id": "2CBA", "chain": "A", "group": "high_resolution_adsorption_model"},
    {"id": "ubiquitin", "label": "Human ubiquitin", "pdb_id": "1UBQ", "chain": "A", "group": "high_resolution_adsorption_model"},
    {"id": "acp", "label": "Acylphosphatase", "pdb_id": "1APS", "chain": "A", "group": "high_resolution_adsorption_model"},
    {"id": "catalase", "label": "Bovine liver catalase tetramer", "pdb_id": "1TGU", "chain": "A,B,C,D", "group": "diverse_size_function"},
    {"id": "asparaginase", "label": "E. coli L-asparaginase II tetramer", "pdb_id": "3ECA", "chain": "A,B,C,D", "group": "diverse_size_function"},
    {"id": "beta2m", "label": "Human beta2-microglobulin", "pdb_id": "1JNJ", "chain": "A", "group": "diverse_size_function"},
    {"id": "cytochrome_c", "label": "Horse heart cytochrome c", "pdb_id": "1HRC", "chain": "A", "group": "diverse_size_function"},
]

THRESHOLDS = [0.02, 0.05, 0.075, 0.10, 0.15]
PKA_SHIFTS = [-0.5, 0.5]
STATE_SENSITIVE_CHANNELS = {
    "cationic", "anionic", "hbond_acceptor", "pi_carbon", "oxide",
    "hydroxyapatite", "metal_coord", "phosphate",
}


def fetch_pdb(pdb_id: str) -> str:
    with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=45) as r:
        return r.read().decode("utf-8", errors="replace")


def run_case(raw: str, chain: str, *, points: int = 200, threshold: float = 0.05, pka_shift: float = 0.0) -> dict:
    prepared, _ = prepare_pdb_text(raw, chain=chain)
    old_points, old_thr = v1.SASA_POINTS, v1.SC_RSA_THRESHOLD
    old_pka = dict(v1.PKA)
    try:
        v1.SASA_POINTS = int(points)
        v1.SC_RSA_THRESHOLD = float(threshold)
        if pka_shift:
            for k in v1.PKA:
                v1.PKA[k] = float(old_pka[k]) + float(pka_shift)
        req = v1.AnalyzeRequest(
            pdb_text=prepared,
            chain=None,
            env=v1.EnvParams(pH=7.4, ionic=150.0, temp=298.0),
        )
        return v1.analyze(req)
    finally:
        v1.SASA_POINTS, v1.SC_RSA_THRESHOLD = old_points, old_thr
        v1.PKA.clear(); v1.PKA.update(old_pka)


def jaccard(a: set[str], b: set[str]) -> float:
    u = a | b
    return len(a & b) / len(u) if u else 1.0


def top10_sets(result: dict) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for chemistry, payload in result.get("chemistries", {}).items():
        out[chemistry] = {
            str(x.get("center_key"))
            for x in payload.get("top_patches", [])[:10]
            if x.get("center_key")
        }
    return out


def surface_set(result: dict) -> set[str]:
    return {str(r["key"]) for r in result.get("surface_residues", [])}


def scrsa_map(result: dict) -> dict[str, float]:
    out = {}
    for r in result.get("all_residues", []):
        key = str(r.get("key"))
        val = r.get("scrsa_raw", r.get("scrsa"))
        if key and val is not None:
            out[key] = float(val)
    return out


def scrsa_error(a: dict, b: dict) -> tuple[float | None, float | None]:
    ma, mb = scrsa_map(a), scrsa_map(b)
    keys = sorted(set(ma) & set(mb))
    if not keys:
        return None, None
    dif = np.asarray([ma[k] - mb[k] for k in keys], dtype=float)
    return float(np.mean(np.abs(dif))), float(math.sqrt(np.mean(dif * dif)))


def median_hotspot_jaccard(ref: dict[str, set[str]], cur: dict[str, set[str]], channels=None) -> float | None:
    names = sorted(channels if channels is not None else set(ref) | set(cur))
    vals = [jaccard(ref.get(c, set()), cur.get(c, set())) for c in names if c in ref or c in cur]
    return float(median(vals)) if vals else None


def summarize_rows(rows: list[dict]) -> dict:
    good = [r for r in rows if r.get("status") == "ok"]
    summary = {
        "n_requested": len(rows),
        "n_successful": len(good),
        "n_failed": len(rows) - len(good),
        "median_surface_set_jaccard_200_vs_500": None,
        "median_scrsa_mae_200_vs_500": None,
        "median_scrsa_rmse_200_vs_500": None,
        "median_top10_hotspot_jaccard_by_threshold": {},
        "median_state_sensitive_top10_jaccard_by_pka_shift": {},
        "by_group": {},
    }
    if good:
        summary["median_surface_set_jaccard_200_vs_500"] = float(median(r["surface_set_jaccard_200_vs_500"] for r in good))
        mae = [r["scrsa_mae_200_vs_500"] for r in good if r["scrsa_mae_200_vs_500"] is not None]
        rmse = [r["scrsa_rmse_200_vs_500"] for r in good if r["scrsa_rmse_200_vs_500"] is not None]
        summary["median_scrsa_mae_200_vs_500"] = float(median(mae)) if mae else None
        summary["median_scrsa_rmse_200_vs_500"] = float(median(rmse)) if rmse else None
        for thr in THRESHOLDS:
            k = str(thr)
            vals = [r["median_top10_hotspot_jaccard_vs_threshold_005"].get(k) for r in good]
            vals = [v for v in vals if v is not None]
            summary["median_top10_hotspot_jaccard_by_threshold"][k] = float(median(vals)) if vals else None
        for shift in PKA_SHIFTS:
            k = str(shift)
            vals = [r["state_sensitive_top10_jaccard_vs_pka_shift_0"].get(k) for r in good]
            vals = [v for v in vals if v is not None]
            summary["median_state_sensitive_top10_jaccard_by_pka_shift"][k] = float(median(vals)) if vals else None

    for group in sorted({r["group"] for r in rows}):
        gr = [r for r in good if r["group"] == group]
        if not gr:
            summary["by_group"][group] = {"n": 0}
            continue
        summary["by_group"][group] = {
            "n": len(gr),
            "median_surface_set_jaccard_200_vs_500": float(median(r["surface_set_jaccard_200_vs_500"] for r in gr)),
            "median_threshold_0.10_top10_jaccard": float(median(r["median_top10_hotspot_jaccard_vs_threshold_005"]["0.1"] for r in gr)),
            "median_threshold_0.15_top10_jaccard": float(median(r["median_top10_hotspot_jaccard_vs_threshold_005"]["0.15"] for r in gr)),
        }
    return summary


def main() -> None:
    # Use the corrected V2 chemistry definitions for all sensitivity runs.
    apply_publication_chemistry(v1)

    rows = []
    for spec in STRUCTURES:
        row = {**spec}
        try:
            raw = fetch_pdb(spec["pdb_id"])
            base = run_case(raw, spec["chain"], points=200, threshold=0.05)
            hi = run_case(raw, spec["chain"], points=500, threshold=0.05)

            surf_base, surf_hi = surface_set(base), surface_set(hi)
            mae, rmse = scrsa_error(base, hi)
            base_top = top10_sets(base)

            threshold_hotspots = {}
            threshold_surfaces = {}
            for thr in THRESHOLDS:
                cur = base if thr == 0.05 else run_case(raw, spec["chain"], points=200, threshold=thr)
                threshold_hotspots[str(thr)] = median_hotspot_jaccard(base_top, top10_sets(cur))
                threshold_surfaces[str(thr)] = jaccard(surf_base, surface_set(cur))

            pka_hotspots = {}
            for shift in PKA_SHIFTS:
                cur = run_case(raw, spec["chain"], points=200, threshold=0.05, pka_shift=shift)
                pka_hotspots[str(shift)] = median_hotspot_jaccard(
                    base_top, top10_sets(cur), STATE_SENSITIVE_CHANNELS
                )

            row.update({
                "status": "ok",
                "n_residues": int(base.get("stats", {}).get("n_residues", 0)),
                "n_surface_200": len(surf_base),
                "n_surface_500": len(surf_hi),
                "surface_set_jaccard_200_vs_500": jaccard(surf_base, surf_hi),
                "surface_threshold_crossing_count_200_vs_500": len(surf_base ^ surf_hi),
                "scrsa_mae_200_vs_500": mae,
                "scrsa_rmse_200_vs_500": rmse,
                "median_top10_hotspot_jaccard_vs_threshold_005": threshold_hotspots,
                "surface_set_jaccard_vs_threshold_005": threshold_surfaces,
                "state_sensitive_top10_jaccard_vs_pka_shift_0": pka_hotspots,
            })
        except Exception as exc:
            row.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
        rows.append(row)
        print(json.dumps(row, indent=2, sort_keys=True))

    payload = {
        "purpose": "label-independent sensitivity analysis across a 14-protein structural panel",
        "reference_configuration": {
            "pH": 7.4,
            "ionic_mM": 150.0,
            "sasa_points_per_atom": 200,
            "scrsa_threshold": 0.05,
            "patch_radii_A": [5.0, 8.0],
        },
        "perturbations": {
            "sasa_points_per_atom": [200, 500],
            "scrsa_thresholds": THRESHOLDS,
            "uniform_pka_shift": PKA_SHIFTS,
        },
        "structures": rows,
        "summary": summarize_rows(rows),
    }
    with open("sensitivity_panel.json", "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    print(json.dumps(payload["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
