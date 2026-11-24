"""Frozen TV48: full 3x3 pH x ionic strength SAS-APBS technical sensitivity.

48 verified assemblies x (pH 6.5,7.4,8.5) x (50,150,300 mM).
Each pH gets independently prepared PQR with PROPKA. Each ionic setting gets
a new, real APBS LPBE solution. The primary InterfaceScout SAS estimator
(42 Fibonacci directions, 1.4 A probe + 0.25 A outward offset) is identical
to production; the original 11 maps and Pareto ranking are unchanged.

This is a finite sensitivity study, not experimental prediction validation.
All protein-level metrics are calculated on paired SAME-residue keys and
all intermediate surface potentials are exported. A failed case fails CI.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "backend"), str(ROOT / "validation")]
from electrostatics import _pqr_atoms
from surface_potential import sample_residue_sas_potential
from tv48_apbs_sas_grid_audit import prepare_case, make_pqr, run_apbs, MANIFEST

PH_VALUES = (6.5, 7.4, 8.5)
IONIC_VALUES = (50.0, 150.0, 300.0)
REF = (7.4, 150.0)
N_SHARDS = 12
OUT = Path("tv48_ph_ionic_full")
OUT.mkdir(exist_ok=True)
assert len(MANIFEST) == 48 and len({r["pdb_id"] for r in MANIFEST}) == 48


def label(ph, ionic):
    return "ph" + str(ph).replace(".", "") + "_i" + str(int(ionic))


def json_write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")


def csv_write(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        raise ValueError("No records to save")
    cols = list(dict.fromkeys(k for row in records for k in row))
    with path.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def get_phi(mapping, key):
    return float(mapping[key]["surface_potential_kT_per_e"])


def paired_stats(left, right, *, expected_min=20):
    common = sorted(set(left) & set(right))
    if len(common) < expected_min:
        raise ValueError(f"Insufficient matched SAS residue count ({len(common)})")
    a = np.asarray([get_phi(left, k) for k in common], dtype=float)
    b = np.asarray([get_phi(right, k) for k in common], dtype=float)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Nonfinite matched surface potential")
    delta = a - b
    absdelta = np.abs(delta)
    rho = float(spearmanr(a, b).statistic) if np.ptp(a) > 1e-12 and np.ptp(b) > 1e-12 else None
    nhot = max(1, math.ceil(len(common) * .1))
    hota = {common[int(i)] for i in np.argsort(-np.abs(a), kind="stable")[:nhot]}
    hotb = {common[int(i)] for i in np.argsort(-np.abs(b), kind="stable")[:nhot]}
    overlap = len(hota & hotb) / len(hota | hotb)
    return {
        "common_residues": len(common),
        "left_only": len(left) - len(common),
        "right_only": len(right) - len(common),
        "coverage_fraction_left": round(len(common) / len(left), 7),
        "coverage_fraction_right": round(len(common) / len(right), 7),
        "median_abs_delta_kT_per_e": round(float(np.median(absdelta)), 7),
        "p95_abs_delta_kT_per_e": round(float(np.quantile(absdelta, .95)), 7),
        "max_abs_delta_kT_per_e": round(float(np.max(absdelta)), 7),
        "mean_signed_delta_kT_per_e": round(float(np.mean(delta)), 7),
        "rmse_delta_kT_per_e": round(float(np.sqrt(np.mean(delta ** 2))), 7),
        "spearman": round(rho, 7) if rho is not None and np.isfinite(rho) else None,
        "top10_absolute_potential_jaccard": round(overlap, 7),
    }


def interaction_stats(fields):
    """Difference-in-differences: joint response beyond additive endpoints."""
    keys = sorted(set.intersection(
        *(set(fields[x]) for x in ((6.5, 50.0), (6.5, 300.0),
                                  (8.5, 50.0), (8.5, 300.0)))))
    if len(keys) < 20:
        raise ValueError("Insufficient complete 2x2 endpoint intersection")
    diffs = []
    for key in keys:
        d = ((get_phi(fields[(8.5, 300.0)], key) - get_phi(fields[(8.5, 50.0)], key))
             - (get_phi(fields[(6.5, 300.0)], key) - get_phi(fields[(6.5, 50.0)], key)))
        diffs.append(d)
    a = np.asarray(diffs, float)
    return {
        "common_2x2_residues": len(keys),
        "median_abs_difference_in_differences_kT_per_e": round(float(np.median(np.abs(a))), 7),
        "p95_abs_difference_in_differences_kT_per_e": round(float(np.quantile(np.abs(a), .95)), 7),
        "mean_signed_difference_in_differences_kT_per_e": round(float(a.mean()), 7),
        "interpretation": "Potential response interaction, not experimental synergy or affinity.",
    }


def one_protein(item):
    pid, pdb, analysis, digest = prepare_case(item)
    allkeys = {row["key"] for row in analysis["residue_profiles"]}
    exposed = {row["key"] for row in analysis["residue_profiles"] if row["surface_exposed"]}
    fields, metas, charges, rows, pqr_checksums = {}, {}, {}, [], {}
    with tempfile.TemporaryDirectory(prefix="TV48_sas_3x3_") as td:
        workdir = Path(td)
        for ph in PH_VALUES:
            atoms = make_pqr(workdir, pdb, ph)
            pqr_path = workdir / "structure.pqr"
            pqr_checksums[str(ph)] = hashlib.sha256(pqr_path.read_bytes()).hexdigest()
            residue_keys = {f'{a["chain"]}:{a["seq"]}:{a["icode"]}' for a in atoms}
            if not allkeys.issubset(residue_keys):
                raise AssertionError(f"{pid} PQR/PDB mismatch at pH {ph}")
            charge_by_key = {}
            for atom in atoms:
                key = f'{atom["chain"]}:{atom["seq"]}:{atom["icode"]}'
                charge_by_key[key] = charge_by_key.get(key, 0.) + atom["charge"]
            charges[ph] = {
                "net_pqr_charge_e": round(float(sum(a["charge"] for a in atoms)), 6),
                "n_pqr_atoms": len(atoms),
                "n_residue_charge_differences_vs_pH74": None,
            }
            if ph != 7.4:
                reference_charges = charges[7.4].get("residue_charges", {})
                charges[ph]["n_residue_charge_differences_vs_pH74"] = sum(
                    abs(charge_by_key.get(k, 0) - reference_charges.get(k, 0)) > .01
                    for k in allkeys
                )
            charges[ph]["residue_charges"] = charge_by_key
            for ionic in IONIC_VALUES:
                state = (ph, ionic)
                alias = label(*state)
                grid, metadata = run_apbs(workdir, atoms, ionic, target_grid=.8, label=alias)
                potentials, surface_meta = sample_residue_sas_potential(
                    atoms, *grid, exposed, n_directions=42
                )
                del grid
                if len(potentials) < .90 * len(exposed):
                    raise AssertionError(f"{pid} {alias}: SAS coverage below 90%")
                fields[state] = potentials
                metas[state] = {
                    **metadata,
                    "pH": ph,
                    "ionic_mM": ionic,
                    "n_surface_residues": len(potentials),
                    "n_exposed_surface_residues": len(exposed),
                    "n_shell_samples": surface_meta["n_surface_samples"],
                    "pqr_sha256": pqr_checksums[str(ph)],
                    "net_pqr_charge_e": charges[ph]["net_pqr_charge_e"],
                    "n_pqr_atoms": len(atoms),
                }
                for key, value in potentials.items():
                    rows.append({
                        "pdb_id": pid, "key": key,
                        "pH": ph, "ionic_mM": ionic,
                        "SAS_potential_kT_per_e": value["surface_potential_kT_per_e"],
                        "SAS_area_A2": value["surface_sampling_area_A2"],
                        "n_surface_sampling_points": value["surface_potential_samples"],
                        "pqr_residue_charge_e": round(charge_by_key[key], 6),
                    })
                print("TV48_PH_IONIC_CONDITION " + json.dumps({
                    "pdb_id": pid, "pH": ph, "ionic_mM": ionic,
                    "n_surface_residues": len(potentials),
                    "grid_dims": metadata["grid_dims"],
                    "charge_e": charges[ph]["net_pqr_charge_e"],
                }), flush=True)

    # Remove nested charge mapping from concise metadata; retain charge per residue in CSV.
    for ph in PH_VALUES:
        del charges[ph]["residue_charges"]

    if len(fields) != 9 or len(rows) < 700:
        raise AssertionError(f"Missing factorial conditions or residue results: {pid}")
    ref = fields[REF]
    comparisons = []
    for ph in PH_VALUES:
        for ionic in IONIC_VALUES:
            state = (ph, ionic)
            comparisons.append({
                "pdb_id": pid, "pH": ph, "ionic_mM": ionic,
                "comparison": "vs_reference_pH7.4_150mM",
                **paired_stats(fields[state], ref),
            })
    for ph in PH_VALUES:
        comparisons.append({
            "pdb_id": pid, "pH": ph, "ionic_mM": None,
            "comparison": "ionic_300_vs_50_at_fixed_pH",
            **paired_stats(fields[(ph,300.)], fields[(ph,50.)]),
        })
    for ionic in IONIC_VALUES:
        comparisons.append({
            "pdb_id": pid, "pH": None, "ionic_mM": ionic,
            "comparison": "pH_8.5_vs_6.5_at_fixed_ionic",
            **paired_stats(fields[(8.5,ionic)], fields[(6.5,ionic)]),
        })
    interaction = interaction_stats(fields)
    charge_rows = [
        {"pdb_id": pid, "pH": ph, "ionic_mM": None,
         **{k: v for k, v in charges[ph].items() if k != "residue_charges"}}
        for ph in PH_VALUES
    ]
    out = {
        "pdb_id":pid,"status":"PASS","assembly_sha256":digest,
        "n_residues":len(allkeys),"n_surface_exposed_residues":len(exposed),
        "n_actual_apbs_solves":len(fields), "n_real_propka_pqr_preparations":len(PH_VALUES),
        "n_surface_sas_rows":len(rows),
        "reference":{"pH":7.4,"ionic_mM":150.0,"temperature_K":298.,
                     "sas_directions_per_atom":42,"probe_A":1.4,"outward_offset_A":.25},
        "grid":{label(*state):value for state,value in metas.items()},
        "pqr_charge_by_pH":charge_rows,
        "interaction":interaction,
        "n_canonical_maps_unchanged":11,
    }
    dumpdir = OUT / "proteins" / pid
    csv_write(dumpdir / "all_9_states_residue_potentials.csv", rows)
    csv_write(dumpdir / "paired_comparisons.csv", comparisons)
    csv_write(dumpdir / "pH_charge_audit.csv", charge_rows)
    json_write(dumpdir / "SUMMARY.json", out)
    return out


def run_shard(index, shards):
    results = []
    for i, item in enumerate(MANIFEST):
        if i % shards != index:
            continue
        pid = item["pdb_id"].upper()
        try:
            result = one_protein(item)
        except Exception as exc:
            result = {"pdb_id":pid,"status":"FAIL",
                      "error_type":type(exc).__name__,"error":str(exc),
                      "traceback_tail":traceback.format_exc()[-6000:]}
            json_write(OUT / "failures" / f"{pid}.json", result)
        results.append(result)
        print("TV48_PH_IONIC_CASE " + json.dumps({
            "pdb_id":pid,"status":result["status"],
            "real_APBS":result.get("n_actual_apbs_solves", 0),
            "error":result.get("error","")[:260]
        }), flush=True)
    json_write(OUT / f"shard_{index}.json", results)
    failed = sum(r["status"]!="PASS" for r in results)
    print("TV48_PH_IONIC_SHARD " + json.dumps({
        "shard":index,"expected":len(results),"passed":len(results)-failed,
        "failed":failed
    }), flush=True)
    if failed:
        raise SystemExit(f"{failed} protein(s) failed in shard {index}; see artifacts")


def aggregate(shards):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    records = []
    missing = []
    for i in range(shards):
        path = OUT / f"shard_{i}.json"
        if path.exists():
            records.extend(json.loads(path.read_text(encoding="utf-8")))
        else:
            missing.append(str(path))
    failed = [r for r in records if r["status"] != "PASS"]
    ids = [r["pdb_id"] for r in records]
    if (missing or len(records)!=48 or len(set(ids))!=48 or
        set(ids)!={m["pdb_id"] for m in MANIFEST} or failed or
        sum(r.get("n_actual_apbs_solves",0) for r in records)!=432):
        json_write(OUT / "INCOMPLETE_AUDIT.json", {
            "missing_shards":missing,"failed_proteins":failed,
            "actual_count":len(records),"n_apbs_solved":
                sum(r.get("n_actual_apbs_solves",0) for r in records)
        })
        raise ValueError("TV48 factorial APBS verification incomplete")
    condition_rows = []
    comparisons = []
    charge_rows = []
    interaction_rows = []
    residue_dfs = []
    for pid in [m["pdb_id"] for m in MANIFEST]:
        root = OUT / "proteins" / pid
        state = next(r for r in records if r["pdb_id"]==pid)
        if not root.is_dir():
            raise ValueError("Missing per-protein 9-state evidence for " + pid)
        df = pd.read_csv(root / "all_9_states_residue_potentials.csv")
        if (len(df.groupby(["pH","ionic_mM"]))!=9 or
            df.groupby(["pH","ionic_mM"]).size().min()<20):
            raise ValueError("Missing APBS factorial residue condition " + pid)
        residue_dfs.append(df)
        comparison = pd.read_csv(root / "paired_comparisons.csv")
        if len(comparison)!=15:raise ValueError(pid+" must have 15 comparisons")
        comparisons.append(comparison)
        charge = pd.read_csv(root / "pH_charge_audit.csv")
        if len(charge)!=3:raise ValueError(pid+" must have 3 PQR pH states")
        charge_rows.append(charge)
        interaction_rows.append({"pdb_id":pid, **{
            k:v for k,v in state["interaction"].items() if k!="interpretation"
        }})
        for k, v in state["grid"].items():
            condition_rows.append({"pdb_id":pid,"condition":k,**v})
    big = pd.concat(residue_dfs, ignore_index=True)
    comp = pd.concat(comparisons, ignore_index=True)
    charge = pd.concat(charge_rows, ignore_index=True)
    interactions = pd.DataFrame(interaction_rows)
    conditions = pd.DataFrame(condition_rows)
    if len(conditions)!=432 or len(big)<48000:
        raise ValueError("Missing 9 condition data or too few residue rows")
    if comp.groupby("pdb_id").size().ne(15).any():
        raise ValueError("Incomplete paired comparisons")
    big.to_csv(OUT / "TV48_ALL_9x48_RESIDUE_SAS_POTENTIALS.csv", index=False)
    comp.to_csv(OUT / "TV48_ALL_PROTEIN_PAIRED_SENSITIVITY.csv", index=False)
    charge.to_csv(OUT / "TV48_PQR_CHARGE_CHANGES_BY_pH.csv", index=False)
    interactions.to_csv(OUT / "TV48_pH_IONIC_INTERACTION.csv", index=False)
    conditions.to_csv(OUT / "TV48_APBS_ACTUAL_GRID_AND_PROVENANCE.csv", index=False)

    ref = comp[comp.comparison=="vs_reference_pH7.4_150mM"]
    def median_ci(vals, seed):
        v = np.asarray(vals, dtype=float)
        v = v[np.isfinite(v)]
        if len(v)<5:return None
        rng = np.random.default_rng(seed)
        sample = rng.choice(v, size=(2000,len(v)), replace=True)
        return [round(float(x), 6) for x in np.quantile(np.median(sample,axis=1),[.025,.975])]
    variants = []
    for ph in PH_VALUES:
        for ionic in IONIC_VALUES:
            sub = ref[(ref.pH==ph)&(ref.ionic_mM==ionic)]
            if len(sub)!=48:raise ValueError(f"Missing 48 matched records for {ph}/{ionic}")
            med = float(sub.median_abs_delta_kT_per_e.median())
            variants.append({
                "pH":ph,"ionic_mM":ionic,"n_proteins":48,
                "protein_median_abs_delta_kT_per_e":round(med,7),
                "protein_median_abs_delta_bootstrap_95pct_CI":median_ci(sub.median_abs_delta_kT_per_e,20261124),
                "protein_median_p95_abs_delta_kT_per_e":round(float(sub.p95_abs_delta_kT_per_e.median()),7),
                "protein_median_spearman":round(float(sub.spearman.median()),7),
                "protein_median_top10_jaccard":round(float(sub.top10_absolute_potential_jaccard.median()),7),
                "n_proteins_with_p95_change_over_0p5":int((sub.p95_abs_delta_kT_per_e>.5).sum()),
                "median_coverage_fraction_vs_reference":round(float(sub.coverage_fraction_left.median()),7),
            })
    import math as _math
    if not all(_math.isfinite(v["protein_median_abs_delta_kT_per_e"]) for v in variants):
        raise ValueError("Invalid summary metric")
    summaries_by_contrast=[]
    for name in ("ionic_300_vs_50_at_fixed_pH","pH_8.5_vs_6.5_at_fixed_ionic"):
        sub=comp[comp.comparison==name]
        levelkey="pH" if name.startswith("ionic") else "ionic_mM"
        for level,group in sub.groupby(levelkey):
            if len(group)!=48:raise ValueError("Unbalanced matched contrast")
            summaries_by_contrast.append({
                "contrast":name,"fixed_level":float(level),"n_proteins":len(group),
                "median_protein_median_abs_change":round(float(group.median_abs_delta_kT_per_e.median()),7),
                "median_protein_p95_abs_change":round(float(group.p95_abs_delta_kT_per_e.median()),7),
                "median_protein_spearman":round(float(group.spearman.median()),7),
            })
    figures=OUT/"figures"
    figures.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                         "font.family":"sans-serif","savefig.dpi":240})

    # Excluding redundant baseline from effect figure (zero by construction).
    fig, ax=plt.subplots(figsize=(8.4,4.5))
    conditions_for_plot = [(ph,ionic) for ph in PH_VALUES for ionic in IONIC_VALUES if (ph,ionic)!=REF]
    plotvals=[ref[(ref.pH==ph)&(ref.ionic_mM==ionic)].median_abs_delta_kT_per_e.to_numpy()
              for ph,ionic in conditions_for_plot]
    ax.boxplot(plotvals, tick_labels=[f"pH {ph:g}\n{int(ionic)} mM" for ph,ionic in conditions_for_plot],
               showfliers=False)
    ax.set_ylabel("Protein-level median |ΔSAS potential| (kT/e)")
    fig.tight_layout()
    for suffix in ("png","pdf"):fig.savefig(figures/("01_factorial_change_from_reference."+suffix))
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.7,4.9))
    mat=np.asarray([[next(v["protein_median_abs_delta_kT_per_e"] for v in variants
                          if v["pH"]==ph and v["ionic_mM"]==ionic)
                     for ionic in IONIC_VALUES] for ph in PH_VALUES])
    im=ax.imshow(mat,aspect="auto")
    for i in range(3):
        for j in range(3):
            ax.text(j,i,f"{mat[i,j]:.3f}",ha="center",va="center",
                    color="white" if mat[i,j]>(mat.min()+mat.max())/2 else "black")
    ax.set_xticks(range(3),["50","150","300"])
    ax.set_yticks(range(3),["6.5","7.4","8.5"])
    ax.set_xlabel("Ionic strength (mM)")
    ax.set_ylabel("pH")
    cb=fig.colorbar(im,ax=ax,label="Median of protein medians |Δφ| (kT/e)")
    fig.tight_layout()
    for suffix in ("png","pdf"):fig.savefig(figures/("02_factorial_sensitivity_heatmap."+suffix))
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.5,4.2))
    x=np.arange(len(interactions))
    sort=interactions.sort_values("median_abs_difference_in_differences_kT_per_e")
    ax.scatter(np.arange(len(sort)),sort.median_abs_difference_in_differences_kT_per_e,s=12)
    ax.set_xlabel("Proteins sorted by interaction magnitude (n=48)")
    ax.set_ylabel("Median |pH×ionic difference-in-differences| (kT/e)")
    fig.tight_layout()
    for suffix in ("png","pdf"):fig.savefig(figures/("03_pH_ionic_interaction_by_protein."+suffix))
    plt.close(fig)

    summary = {
        "status":"PASS_COMPUTATION",
        "n_proteins":48,"n_factorial_conditions":9,
        "n_actual_APBS_solves":432,
        "n_actual_PROPKA_PQR_preparations":144,
        "n_sas_residue_condition_rows":len(big),
        "n_protein_condition_comparisons":len(comp),
        "n_mapped_surface_residues_at_baseline":int(
            big[(big.pH==7.4)&(big.ionic_mM==150.0)].shape[0]),
        "all_11_maps_unchanged":True,
        "conditions":{"pH":list(PH_VALUES),"ionic_mM":list(IONIC_VALUES),
                      "baseline_pH":7.4,"baseline_ionic_mM":150.0,
                      "temperature_K":298.0,"grid_target_A":.8,
                      "SAS_directions":42,"SAS_probe_A":1.4,"offset_A":.25},
        "variants":variants,
        "fixed_other_factor_comparisons":summaries_by_contrast,
        "median_pH_ionic_interaction_abs_kT_per_e":round(float(
            interactions.median_abs_difference_in_differences_kT_per_e.median()),7),
        "max_protein_p95_interaction_abs_kT_per_e":round(float(
            interactions.p95_abs_difference_in_differences_kT_per_e.max()),7),
        "limitations":[
            "These are actual PDB2PQR/PROPKA/APBS numerical sensitivity results, not measured adsorption outcomes.",
            "PQR protonation microstates are discrete and may not differ at all pH values for every protein.",
            "The same 42-direction SAS estimator is used at all conditions; sampling convergence is in a separate TV48 audit.",
            "Potential changes are evaluated at matched residue identifiers; coverage changes are reported.",
            "A difference-in-differences is descriptive nonadditivity; it does not imply experimental synergy.",
            "Grid is the production capped mg-auto setting; actual grid dimensions are retained.",
        ],
    }
    json_write(OUT/"SUMMARY.json",summary)
    print("TV48_PH_IONIC_FACTORIAL_SUMMARY "+json.dumps(summary,ensure_ascii=False),flush=True)


def self_test():
    data1={f"A:{i}:":{"surface_potential_kT_per_e":float(i)}
           for i in range(1,41)}
    data2={f"A:{i}:":{"surface_potential_kT_per_e":float(i)+.25}
           for i in range(1,41)}
    p=paired_stats(data1,data2)
    assert p["common_residues"]==40 and abs(p["median_abs_delta_kT_per_e"]-.25)<1e-9
    assert abs(p["spearman"]-1)<1e-7
    sample = {(ph,i):data1 for ph in PH_VALUES for i in IONIC_VALUES}
    interaction=interaction_stats(sample)
    assert abs(interaction["median_abs_difference_in_differences_kT_per_e"])<1e-10
    print("TV48_PH_IONIC_SELFTEST PASS",flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--shard",type=int)
    parser.add_argument("--shards",type=int,default=N_SHARDS)
    parser.add_argument("--aggregate",action="store_true")
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test:self_test()
    elif args.aggregate:aggregate(args.shards)
    elif args.shard is not None and 0<=args.shard<args.shards:
        run_shard(args.shard,args.shards)
    else:parser.error("Provide --self-test, --shard N, or --aggregate")


if __name__=="__main__":
    main()
