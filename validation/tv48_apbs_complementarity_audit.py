"""TV48: APBS-aware complementary-information audit on frozen production maps.

Real PDB2PQR/PROPKA/APBS is computed for every structure; canonical maps
are independently run WITHOUT and WITH electrostatics and checked identical.
All 48 structures are the sampling units for uncertainty; residue rows are
within-structure observations, not independent biological replicates.
No adsorption affinity/prediction accuracy or trained weights are inferred.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
import traceback
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "backend"), str(ROOT / "validation")]
from interface_engine import analyze_all_maps
from model_settings import MODEL_VERSION, PUBLIC_VERSION
import strict_apo_compat_48 as prep
from export_bundle import MAP_NAMES

MANIFEST = list(csv.DictReader(
    (ROOT / "validation" / "strict_apo_48_selected_candidates_v2_final.csv").open(encoding="utf-8")
))
assert len(MANIFEST) == 48 and len({x["pdb_id"] for x in MANIFEST}) == 48
OUT = Path("tv48_apbs_metrics")
OUT.mkdir(exist_ok=True)


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        w = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def rho(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 5 or np.ptp(x) < 1e-12 or np.ptp(y) < 1e-12:
        return None
    val = float(spearmanr(x, y).statistic)
    return round(val, 7) if np.isfinite(val) else None


def analyze(item):
    pid = item["pdb_id"].upper()
    prep.prepare(pid)
    cif = prep.OUT / f"{pid}-assembly1.cif"
    checksum = hashlib.sha256(cif.read_bytes()).hexdigest()
    if checksum != item["assembly_sha256"]:
        raise ValueError("Biological assembly SHA-256 differs from frozen TV48 manifest")
    pdb = (prep.OUT / f"{pid}_prepared.pdb").read_text(encoding="utf-8")
    kw = dict(pdb_text=pdb, chain=None, pH=7.4, ionic_mM=150.0, temp_K=298.0)
    vanilla = analyze_all_maps(**kw)
    with_apbs = analyze_all_maps(**kw, electrostatics=True)
    if vanilla["maps"] != with_apbs["maps"] or vanilla["residue_profiles"] != with_apbs["residue_profiles"]:
        raise AssertionError("APBS unexpectedly changed the canonical 11 maps or residue profiles")
    if list(with_apbs["map_order"]) != list(MAP_NAMES) or with_apbs["n_maps"] != 11:
        raise AssertionError("Missing canonical map / unexpected map ordering")
    if with_apbs["version"] != MODEL_VERSION:
        raise AssertionError("Production model is not frozen dual-accessibility")
    residue_rows = with_apbs["residue_profiles"]
    electro = with_apbs["electrostatics"]
    if electro.get("status") != "computed" or electro["n_matched_residues"] != len(residue_rows):
        raise AssertionError("Real APBS/PDB2PQR run or residue correspondence failed")
    profile = {r["key"]: r for r in residue_rows}
    potentials = {r["key"]: r for r in electro["residues"]}
    if set(profile) != set(potentials):
        raise AssertionError("Nonunique/missing residue keys in electrostatics results")
    if any(r["res_name"] == "GLY" and (r["scrsa"] != 0 or r["sidechain_sasa"] != 0)
           for r in residue_rows):
        raise AssertionError("Artificial glycine side-chain contribution")

    valid = [
        {"pdb_id": pid, "key": key,
         "res_name": prof["res_name"],
         "surface_exposed": True,
         "trsa_raw": prof["trsa_raw"], "scrsa_raw": prof["scrsa_raw"],
         "net_atomic_charge_e": potentials[key]["net_atomic_charge_e"],
         "APBS_potential_kT_per_e": potentials[key]["sidechain_mean_potential_kT_per_e"]}
        for key, prof in profile.items()
        if prof["surface_exposed"] and potentials[key]["sidechain_mean_potential_kT_per_e"] is not None
    ]
    if len(valid) < 30:
        raise ValueError("Fewer than 30 exposed side-chain potentials")
    values = np.asarray([float(x["APBS_potential_kT_per_e"]) for x in valid])
    charges = np.asarray([float(x["net_atomic_charge_e"]) for x in valid])
    keys = {r["key"] for r in valid}
    n = len(valid)
    k = max(1, math.ceil(0.10 * n))
    hottest = set(x["key"] for x in sorted(valid,
                 key=lambda x: (-abs(float(x["APBS_potential_kT_per_e"])), x["key"]))[:k])
    pos_hot = set(x["key"] for x in sorted(valid,
                 key=lambda x: (-float(x["APBS_potential_kT_per_e"]), x["key"]))[:k])
    neg_hot = set(x["key"] for x in sorted(valid,
                 key=lambda x: (float(x["APBS_potential_kT_per_e"]), x["key"]))[:k])

    means = defaultdict(list)
    for row in valid:
        means[row["res_name"]].append(float(row["APBS_potential_kT_per_e"]))
    grand = float(values.mean())
    total_ss = float(np.sum((values - grand) ** 2))
    within_ss = sum(sum((v - statistics.mean(group)) ** 2 for v in group)
                    for group in means.values())
    within_aa_fraction = within_ss / total_ss if total_ss > 1e-12 else None

    meta = {
        "pdb_id": pid, "status": "PASS", "assembly_sha256": checksum,
        "prepared_sha256": hashlib.sha256(pdb.encode()).hexdigest(),
        "software": PUBLIC_VERSION, "model": MODEL_VERSION,
        "canonical_maps_unchanged": True, "n_maps": 11,
        "n_residues": len(profile), "n_surface_residues": sum(r["surface_exposed"] for r in profile.values()),
        "n_exposed_APBS_sidechains": n, "n_exposed_without_sidechain_potential":
            sum(r["surface_exposed"] for r in profile.values()) - n,
        "n_PQR_atoms": electro["n_pqr_atoms"],
        "APBS_mean_surface_kT_per_e": round(grand, 7),
        "APBS_mean_absolute_surface_kT_per_e": round(float(np.abs(values).mean()), 7),
        "APBS_min_surface_kT_per_e": round(float(values.min()), 7),
        "APBS_max_surface_kT_per_e": round(float(values.max()), 7),
        "APBS_charge_vs_potential_spearman": rho(charges, values),
        "APBS_potential_variance_fraction_within_residue_type": (
            round(within_aa_fraction, 7) if within_aa_fraction is not None else None),
        "surface_gate": "unrounded_total_residue_tRSA_0.05",
        "conditions": "pH 7.4, 150 mM, 298 K",
    }
    maps, joined = [], []
    for ch in MAP_NAMES:
        m = with_apbs["maps"][ch]
        patches = m["patches"]
        favorable = {row["key"]: row for row in m.get("residues", [])}
        props = np.asarray([float(favorable.get(row["key"], {}).get("propensity") or 0.0) for row in valid])
        local = np.asarray([float(favorable.get(row["key"], {}).get("local_score") or 0.0) for row in valid])
        top5 = set(key for patch in patches[:5] for key in patch["members"]) & keys
        top1 = set(patches[0]["members"]) & keys if patches else set()
        h = hottest & top5
        expected = len(top5) * len(hottest) / n
        # Descriptive null uses a uniform redistribution of hotspot identities
        # over the same exposed residues. Spatial autocorrelation is ignored.
        overlap_ratio = len(h) / expected if expected > 0 else None
        fronts = sum(int(p.get("pareto_front", 0)) == 1 for p in patches)
        maps.append({
            "pdb_id": pid, "channel": ch, "n_exposed_APBS_sidechains": n,
            "n_map_favorable_APBS_residues": int(sum(x > 0 for x in local)),
            "n_patches": len(patches), "n_pareto_front1_patches": fronts,
            "n_top5_unique_exposed_members": len(top5),
            "n_top1_unique_exposed_members": len(top1),
            "n_absolute_APBS_hotspot_residues": len(hottest),
            "n_top5_hotspot_overlap": len(h),
            "top5_hotspot_fraction": round(len(h) / len(hottest), 7),
            "top5_hotspot_expected_uniform": round(expected, 7),
            "top5_hotspot_enrichment_vs_uniform": round(overlap_ratio, 7) if overlap_ratio is not None else None,
            "top5_positive_hotspot_overlap": len(pos_hot & top5),
            "top5_negative_hotspot_overlap": len(neg_hot & top5),
            "signed_APBS_vs_map_propensity_spearman": rho(props, values),
            "absolute_APBS_vs_map_propensity_spearman": rho(props, np.abs(values)),
            "signed_APBS_vs_map_local_score_spearman": rho(local, values),
            "absolute_APBS_vs_map_local_score_spearman": rho(local, np.abs(values)),
            "top1_center": patches[0]["center_key"] if patches else "",
            "top5_centers": ";".join(p["center_key"] for p in patches[:5]),
        })
        for row, p, loc in zip(valid, props, local):
            joined.append({
                "pdb_id": pid, "channel": ch, "key": row["key"],
                "res_name": row["res_name"], "trsa_raw": row["trsa_raw"],
                "scrsa_raw": row["scrsa_raw"], "net_atomic_charge_e": row["net_atomic_charge_e"],
                "APBS_potential_kT_per_e": row["APBS_potential_kT_per_e"],
                "APBS_absolute_potential_kT_per_e": abs(row["APBS_potential_kT_per_e"]),
                "map_propensity": round(float(p), 7), "map_local_score": round(float(loc), 7),
                "is_map_top5_patch_member": int(row["key"] in top5),
                "is_top10pct_absolute_potential": int(row["key"] in hottest),
                "is_top10pct_positive_potential": int(row["key"] in pos_hot),
                "is_top10pct_negative_potential": int(row["key"] in neg_hot),
            })
    return meta, maps, valid, joined


MAP_FIELDS = [
    "pdb_id", "channel", "n_exposed_APBS_sidechains",
    "n_map_favorable_APBS_residues", "n_patches", "n_pareto_front1_patches",
    "n_top5_unique_exposed_members", "n_top1_unique_exposed_members",
    "n_absolute_APBS_hotspot_residues", "n_top5_hotspot_overlap",
    "top5_hotspot_fraction", "top5_hotspot_expected_uniform",
    "top5_hotspot_enrichment_vs_uniform", "top5_positive_hotspot_overlap",
    "top5_negative_hotspot_overlap",
    "signed_APBS_vs_map_propensity_spearman",
    "absolute_APBS_vs_map_propensity_spearman",
    "signed_APBS_vs_map_local_score_spearman",
    "absolute_APBS_vs_map_local_score_spearman",
    "top1_center", "top5_centers",
]
RES_FIELDS = [
    "pdb_id", "key", "res_name", "surface_exposed", "trsa_raw", "scrsa_raw",
    "net_atomic_charge_e", "APBS_potential_kT_per_e",
]
JOINED_FIELDS = [
    "pdb_id", "channel", "key", "res_name", "trsa_raw", "scrsa_raw",
    "net_atomic_charge_e", "APBS_potential_kT_per_e",
    "APBS_absolute_potential_kT_per_e", "map_propensity", "map_local_score",
    "is_map_top5_patch_member", "is_top10pct_absolute_potential",
    "is_top10pct_positive_potential", "is_top10pct_negative_potential",
]
META_FIELDS = [
    "pdb_id", "status", "assembly_sha256", "prepared_sha256", "software", "model",
    "canonical_maps_unchanged", "n_maps", "n_residues", "n_surface_residues",
    "n_exposed_APBS_sidechains", "n_exposed_without_sidechain_potential",
    "n_PQR_atoms", "APBS_mean_surface_kT_per_e",
    "APBS_mean_absolute_surface_kT_per_e", "APBS_min_surface_kT_per_e",
    "APBS_max_surface_kT_per_e", "APBS_charge_vs_potential_spearman",
    "APBS_potential_variance_fraction_within_residue_type", "surface_gate",
    "conditions",
]


def do_shard(index, shards):
    metas, maps, res, joined, errors = [], [], [], [], []
    for i, item in enumerate(MANIFEST):
        if i % shards != index:
            continue
        try:
            m, p, r, j = analyze(item)
            metas.append(m)
            maps.extend(p); res.extend(r); joined.extend(j)
            print("APBS_AUDIT_PASS " + json.dumps({
                "pdb_id": item["pdb_id"], "n_maps": len(p),
                "APBS_surface_residues": len(r),
                "map_APBS_rows": len(j),
            }), flush=True)
        except Exception as exc:
            failure = {"pdb_id": item["pdb_id"], "error": str(exc),
                       "error_type": type(exc).__name__, "traceback": traceback.format_exc()[-4000:]}
            errors.append(failure)
            print("APBS_AUDIT_FAIL " + json.dumps(failure), flush=True)
    write_csv(OUT / f"protein_{index}.csv", metas, META_FIELDS)
    write_csv(OUT / f"maps_{index}.csv", maps, MAP_FIELDS)
    write_csv(OUT / f"residues_{index}.csv", res, RES_FIELDS)
    write_csv(OUT / f"joined_{index}.csv", joined, JOINED_FIELDS)
    write_json(OUT / f"errors_{index}.json", errors)
    print("APBS_AUDIT_SHARD " + json.dumps({"shard": index,
        "pass": len(metas), "fail": len(errors)}), flush=True)


def aggregate(shards):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    needed = [OUT / f"{kind}_{i}.{extension}" for i in range(shards)
              for kind, extension in (("protein", "csv"), ("maps", "csv"),
                                      ("residues", "csv"), ("joined", "csv"),
                                      ("errors", "json"))]
    missing = [x.name for x in needed if not x.exists()]
    if missing:
        raise ValueError("Missing shards/records: " + ", ".join(missing))
    m = pd.concat([pd.read_csv(OUT / f"protein_{i}.csv") for i in range(shards)], ignore_index=True)
    maps = pd.concat([pd.read_csv(OUT / f"maps_{i}.csv") for i in range(shards)], ignore_index=True)
    residues = pd.concat([pd.read_csv(OUT / f"residues_{i}.csv") for i in range(shards)], ignore_index=True)
    joined = pd.concat([pd.read_csv(OUT / f"joined_{i}.csv") for i in range(shards)], ignore_index=True)
    errors = sum((json.loads((OUT / f"errors_{i}.json").read_text()) for i in range(shards)), [])
    expected = {x["pdb_id"] for x in MANIFEST}
    if errors or len(m) != 48 or set(m.pdb_id) != expected or not m.pdb_id.is_unique:
        raise ValueError("APBS 48 panel incomplete; failures: " + json.dumps(errors[:4]))
    if len(maps) != 528 or set(maps.channel) != set(MAP_NAMES):
        raise ValueError("Incorrect 528-map coverage")
    if maps.groupby("pdb_id").channel.nunique().ne(11).any():
        raise ValueError("At least one protein lacks one or more maps")
    if len(residues) != int(m.n_exposed_APBS_sidechains.sum()):
        raise ValueError("Residue/sample denominator mismatch")
    if len(joined) != 11 * len(residues):
        raise ValueError("Channel-by-residue matrix incomplete")
    if int(m.n_residues.sum()) != 12068:
        raise ValueError("The frozen TV48 residue count changed")
    if not m.canonical_maps_unchanged.all():
        raise ValueError("Canonical map invariance failure")

    m.to_csv(OUT / "ALL_48_PROTEIN_APBS_METRICS.csv", index=False)
    maps.to_csv(OUT / "ALL_528_MAP_APBS_METRICS.csv", index=False)
    residues.to_csv(OUT / "ALL_EXPOSED_RESIDUE_APBS_VALUES.csv", index=False)
    joined.to_csv(OUT / "ALL_MAP_RESIDUE_APBS_MATCHES.csv", index=False)

    def boot_ci(values, seed, count=2000):
        vals = np.array([float(v) for v in values if pd.notna(v)], dtype=float)
        if len(vals) < 5:
            return [None, None]
        rng = np.random.default_rng(seed)
        samples = rng.choice(vals, (count, len(vals)), replace=True)
        qs = np.percentile(np.median(samples, axis=1), [2.5, 97.5])
        return [round(float(q), 6) for q in qs]

    summary_rows = []
    for ch, part in maps.groupby("channel"):
        signed = part.signed_APBS_vs_map_propensity_spearman.dropna()
        absolute = part.absolute_APBS_vs_map_propensity_spearman.dropna()
        enrich = part.top5_hotspot_enrichment_vs_uniform.dropna()
        summary_rows.append({
            "channel": ch, "n_proteins": part.pdb_id.nunique(),
            "n_signed_correlations_evaluable": len(signed),
            "median_signed_APBS_propensity_spearman": round(float(signed.median()), 6) if len(signed) else None,
            "signed_correlation_protein_bootstrap_CI": boot_ci(signed, 29171),
            "n_absolute_correlations_evaluable": len(absolute),
            "median_absolute_APBS_propensity_spearman": round(float(absolute.median()), 6) if len(absolute) else None,
            "absolute_correlation_protein_bootstrap_CI": boot_ci(absolute, 29173),
            "median_top5_hotspot_enrichment_uniform": round(float(enrich.median()), 6) if len(enrich) else None,
            "hotspot_enrichment_protein_bootstrap_CI": boot_ci(enrich, 29179),
            "mean_patch_count": round(float(part.n_patches.mean()), 5),
        })
    write_json(OUT / "CHANNEL_PROTEIN_LEVEL_SUMMARY.json", summary_rows)
    pd.DataFrame(summary_rows).drop(columns=[
        "signed_correlation_protein_bootstrap_CI",
        "absolute_correlation_protein_bootstrap_CI",
        "hotspot_enrichment_protein_bootstrap_CI",
    ]).to_csv(OUT / "CHANNEL_PROTEIN_LEVEL_SUMMARY.csv", index=False)

    charts = OUT / "figures"
    charts.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 220})
    def save(name):
        plt.tight_layout()
        plt.savefig(charts / (name + ".png"), bbox_inches="tight")
        plt.savefig(charts / (name + ".pdf"), bbox_inches="tight")
        plt.close()
    plt.figure(figsize=(7.2, 4.8))
    plt.hist(m.APBS_mean_absolute_surface_kT_per_e, bins=12)
    plt.xlabel("Protein mean |APBS side-chain potential| (kT/e)")
    plt.ylabel("Proteins (n=48)")
    plt.title("TV48: absolute surface electrostatic potential")
    save("01_TV48_mean_absolute_potential_distribution")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(m.n_exposed_APBS_sidechains, m.APBS_mean_absolute_surface_kT_per_e)
    ax.set_xlabel("Exposed residues with APBS side-chain potential")
    ax.set_ylabel("Mean |potential| (kT/e)")
    ax.set_title("TV48: protein size vs APBS potential magnitude")
    save("02_TV48_surface_residue_count_vs_potential")

    fig, ax = plt.subplots(figsize=(8, 11))
    ordered = m.sort_values("APBS_mean_surface_kT_per_e")
    y = np.arange(len(ordered))
    ax.hlines(y, ordered.APBS_min_surface_kT_per_e,
              ordered.APBS_max_surface_kT_per_e, lw=1)
    ax.scatter(ordered.APBS_mean_surface_kT_per_e, y, s=12)
    ax.set_yticks(y, ordered.pdb_id)
    ax.set_xlabel("Side-chain APBS potential (kT/e): min—max and protein mean")
    ax.set_title("TV48: per-protein surface potential range")
    save("03_TV48_min_mean_max_potential")

    rho_matrix = maps.pivot(index="pdb_id", columns="channel",
                            values="signed_APBS_vs_map_propensity_spearman")
    rho_matrix = rho_matrix.reindex(index=ordered.pdb_id, columns=MAP_NAMES)
    fig, ax = plt.subplots(figsize=(11, 12))
    image = ax.imshow(rho_matrix.to_numpy(dtype=float), aspect="auto",
                      cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_yticks(np.arange(48), labels=rho_matrix.index, fontsize=7)
    ax.set_xticks(np.arange(11), labels=MAP_NAMES, rotation=52, ha="right")
    fig.colorbar(image, ax=ax, label="Within-protein Spearman rho")
    ax.set_title("TV48: APBS potential vs eleven map propensities")
    save("04_TV48_APBS_map_signed_correlation_heatmap")

    vals = [maps.loc[maps.channel == ch,
                     "signed_APBS_vs_map_propensity_spearman"].dropna().values
            for ch in MAP_NAMES]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.boxplot(vals, tick_labels=MAP_NAMES, showfliers=False)
    ax.axhline(0, lw=0.7, alpha=0.7)
    ax.tick_params(axis="x", rotation=52)
    ax.set_ylabel("Within-protein Spearman rho")
    ax.set_title("TV48: distribution across proteins, not pooled residues")
    save("05_TV48_APBS_map_correlation_boxplots")

    fig, ax = plt.subplots(figsize=(9, 5))
    vals = [maps.loc[maps.channel == ch,
                     "top5_hotspot_enrichment_vs_uniform"].dropna().values
            for ch in MAP_NAMES]
    ax.boxplot(vals, tick_labels=MAP_NAMES, showfliers=False)
    ax.axhline(1, lw=0.8, linestyle="--")
    ax.tick_params(axis="x", rotation=52)
    ax.set_ylabel("Observed / uniform expectation (top 10% |APBS| hotspots)")
    ax.set_title("Top-5 region overlap: descriptive residue-label null only")
    save("06_TV48_top5_absolute_APBS_hotspot_enrichment")

    summary = {
        "status": "PASS", "software": PUBLIC_VERSION,
        "model": MODEL_VERSION, "expected_proteins": 48, "proteins_pass": len(m),
        "n_canonical_map_instances": len(maps),
        "n_protein_residues": int(m.n_residues.sum()),
        "n_surface_residues": int(m.n_surface_residues.sum()),
        "n_exposed_sidechain_APBS_residues": len(residues),
        "n_map_by_residue_APBS_observations": len(joined),
        "canonical_map_invariance": bool(m.canonical_maps_unchanged.all()),
        "number_of_evaluable_signed_channel_correlations":
            int(maps.signed_APBS_vs_map_propensity_spearman.notna().sum()),
        "number_of_evaluable_absolute_channel_correlations":
            int(maps.absolute_APBS_vs_map_propensity_spearman.notna().sum()),
        "APBS_within_AA_type_residual_variance_protein_median":
            float(m.APBS_potential_variance_fraction_within_residue_type.median()),
        "graphs": sorted(x.name for x in charts.glob("*.png")),
        "channels": summary_rows,
        "scope": [
            "48 distinct protein structures, not 12,068 independent experiments",
            "Protein-level bootstrap resamples proteins, not residues or 11 correlated channels",
            "APBS values are potential at side-chain heavy atoms, not binding free energy",
            "Top-decile hotspot and uniform-label enrichment are descriptive, not measured contact validation",
            "In-sample AA-type residual variance indicates within-type structural/contextual variation",
            "APBS is an optional annotation; no Pareto reranking or model calibration was done",
            "The current APBS grid and force-field sensitivity must be evaluated separately",
            "TV48 is technical verification, not independent experimental ISV validation",
        ],
    }
    write_json(OUT / "SUMMARY.json", summary)
    print("APBS_AUDIT_SUMMARY " + json.dumps({
        "status": "PASS", "n_proteins": len(m), "n_maps": len(maps),
        "n_apbs_exposed": len(residues),
        "n_correlations": summary["number_of_evaluable_signed_channel_correlations"],
        "within_AA_residual_variance_fraction_median":
            summary["APBS_within_AA_type_residual_variance_protein_median"],
        "figures": len(summary["graphs"])
    }), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=int)
    parser.add_argument("--shards", type=int, default=8)
    parser.add_argument("--aggregate", action="store_true")
    opts = parser.parse_args()
    if opts.aggregate:
        aggregate(opts.shards)
    elif opts.shard is not None and 0 <= opts.shard < opts.shards:
        do_shard(opts.shard, opts.shards)
    else:
        parser.error("Specify --shard N or --aggregate")


if __name__ == "__main__":
    main()
