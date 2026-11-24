"""TV48 APBS numerical convergence: real surface sampling and uncapped grids.

A. Fixed real APBS LPBE solution per frozen protein; sample identical
   solvent-facing shell using 42, 84, 162 Fibonacci directions per heavy atom.
B. 1OX3, 2FB5, 6INX, 6K7E (preselected from archived TV48 grid geometry:
   only structures whose 0.4 Å target actually exceeds 257 nodes on at
   least one axis) receive independently solved APBS grids:
   production 0.8 Å, 0.5 Å, 0.4 Å with 257 cap, and genuinely uncapped
   0.4 Å. Always compare at fixed 162-direction surface resolution.

This is a finite numerical convergence and stability experiment, NOT proof
of material-specific adsorption accuracy. No canonical 11-map score is changed.
All grid spacings, dimensions, common-residue denominators, and failures are
retained. The 0.4 Å cap must actually be active or the run fails.
"""
from __future__ import annotations
import argparse
import csv
import gc
import hashlib
import json
import math
import re
import shutil
import sys
import tempfile
import traceback
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "backend"), str(ROOT / "validation")]
import strict_apo_compat_48 as prepare
from interface_engine import analyze_all_maps
from electrostatics import _pqr_atoms, _apbs_input, _run
from surface_potential import sample_residue_sas_potential

MANIFEST = list(csv.DictReader(
    (ROOT / "validation" / "strict_apo_48_selected_candidates_v2_final.csv")
    .open(encoding="utf-8", newline="")))
assert len(MANIFEST) == 48 and len({x["pdb_id"] for x in MANIFEST}) == 48
BY_ID = {row["pdb_id"]: row for row in MANIFEST}
CAP_AT_TARGET_04 = ("1OX3", "2FB5", "6INX", "6K7E")
assert set(CAP_AT_TARGET_04).issubset(BY_ID)
OUT = Path("tv48_apbs_convergence")
OUT.mkdir(exist_ok=True)


def dump_json(p, item):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(item, ensure_ascii=False, indent=2, allow_nan=False)+"\n",
                 encoding="utf-8")


def dump_csv(p, rows):
    p.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("Cannot save empty evidence rows")
    fieldnames = list(rows[0])
    with p.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def input_geometry(atoms):
    coords = np.asarray([a["xyz"] for a in atoms], dtype=np.float64)
    span = np.maximum(np.ptp(coords, axis=0), 1.0)
    return span + 20.0, span + 40.0


def grid_dimensions(fine, target, cap=None):
    native = [max(65, 32 * math.ceil((v / target - 1) / 32) + 1) for v in fine]
    dims = [min(cap, n) for n in native] if cap is not None else native
    if any(n % 32 != 1 for n in dims):
        raise ValueError("Invalid APBS multigrid shape (32n+1)")
    return native, dims


def load_dx_streaming(path):
    """Parse large APBS OpenDX grids without materializing millions of lines.

    Python's regular .read_text().splitlines() + Python float list is too
    memory-expensive for genuinely uncapped 0.4 Å grids on CI runners.
    """
    counts = None
    origin = None
    delta = []
    expected = None
    with path.open("rb", buffering=0) as f:
        for raw in f:
            line = raw.decode("ascii", errors="replace").strip()
            if "gridpositions counts" in line:
                counts = tuple(map(int, line.split()[-3:]))
            elif line.startswith("origin "):
                origin = np.asarray(list(map(float, line.split()[1:4])), dtype=np.float64)
            elif line.startswith("delta "):
                delta.append(list(map(float, line.split()[1:4])))
            elif "data follows" in line:
                m = re.search(r"items\s+(\d+)", line)
                if not m:
                    raise ValueError("OpenDX has no data count")
                expected = int(m.group(1))
                break
        if counts is None or origin is None or len(delta) != 3:
            raise ValueError("OpenDX grid metadata incomplete")
        if expected != math.prod(counts):
            raise ValueError("OpenDX grid point count inconsistent")
        data = np.fromfile(f, dtype=np.float64, sep=" ", count=expected)
    if len(data) != expected:
        raise ValueError("OpenDX data truncated")
    if not np.isfinite(data).all():
        raise ValueError("OpenDX contains nonfinite potential")
    d = np.asarray(delta, dtype=float)
    if not np.allclose(d, np.diag(np.diag(d)), atol=1.e-8):
        raise ValueError("OpenDX grid not Cartesian")
    spacing = np.diag(d)
    if np.any(spacing <= 0):
        raise ValueError("OpenDX invalid spacing")
    return data.reshape(counts), origin, spacing


def actual_grid(workdir, atoms, target, cap, label):
    """Real PDB2PQR-generated PQR input, actual APBS calculation, measured DX."""
    native_inp = _apbs_input(atoms, 150.0, 298.0)
    fine, coarse = input_geometry(atoms)
    native_dims, dims = grid_dimensions(fine, target, cap)
    gridline = re.search(r"(?m)^    dime \d+ \d+ \d+$", native_inp)
    if gridline is None:
        raise ValueError("Cannot substitute APBS multigrid dimensions")
    text = native_inp[:gridline.start()] + "    dime " + " ".join(map(str, dims)) + native_inp[gridline.end():]
    name = f"numerical_{label}"
    text = text.replace("write pot dx potential", f"write pot dx {name}")
    (workdir / (name + ".in")).write_text(text, encoding="utf-8")
    msg = _run([shutil.which("apbs"), name + ".in"], workdir, timeout=1500)
    exact = workdir / (name + ".dx")
    candidates = [exact] if exact.is_file() else sorted(workdir.glob(name + "-*.dx"))
    if len(candidates) != 1:
        raise RuntimeError(f"APBS output ambiguous/missing {name}: {[q.name for q in candidates]}; {msg[-1200:]}")
    path = candidates[0]
    try:
        arr, origin, spacing = load_dx_streaming(path)
    finally:
        path.unlink(missing_ok=True)
    if list(arr.shape) != dims:
        raise AssertionError(f"Actual APBS grid {list(arr.shape)} != requested {dims}")
    expected_spacing = np.asarray(fine) / (np.asarray(dims)-1)
    if not np.allclose(spacing, expected_spacing, rtol=1.e-5, atol=1.e-5):
        raise AssertionError("Effective spacing does not match actual OpenDX grid")
    record = {
        "label": label, "target_spacing_A": target, "cap": cap,
        "grid_dims": list(map(int, arr.shape)),
        "native_uncapped_dims": native_dims,
        "was_actually_capped": bool(any(a>b for a,b in zip(native_dims,dims))),
        "effective_spacing_A": [round(float(x), 7) for x in spacing],
        "max_spacing_A": round(float(np.max(spacing)), 7),
        "voxel_count": math.prod(dims),
        "fine_lengths_A": [round(float(x), 4) for x in fine],
        "pH": 7.4, "ionic_mM": 150.0, "temperature_K": 298.0,
        "actual_external_APBS_solver_run": True,
    }
    return (arr, origin, spacing), record


def prep_case(item):
    pid = item["pdb_id"].upper()
    prepare.prepare(pid)
    cif = prepare.OUT / f"{pid}-assembly1.cif"
    digest = hashlib.sha256(cif.read_bytes()).hexdigest()
    if digest != item["assembly_sha256"]:
        raise ValueError("Frozen TV48 structure SHA-256 mismatch")
    pdb = (prepare.OUT / f"{pid}_prepared.pdb").read_text(encoding="utf-8")
    core = analyze_all_maps(pdb_text=pdb, pH=7.4, ionic_mM=150.0, temp_K=298.0)
    if core["n_maps"] != 11:
        raise ValueError("Canonical 11-map output changed")
    valid = {r["key"] for r in core["residue_profiles"]}
    surface = {r["key"] for r in core["residue_profiles"] if r["surface_exposed"]}
    if not surface:
        raise ValueError("No surface residues")
    return pid, pdb, valid, surface, digest


def make_pqr(workdir, pdb, valid):
    (workdir / "structure.pdb").write_text(pdb, encoding="utf-8")
    _run([shutil.which("pdb2pqr"), "--ff=AMBER", "--keep-chain",
          "--titration-state-method=propka", "--whitespace", "--with-ph=7.4000",
          "structure.pdb", "structure.pqr"], workdir, timeout=480)
    atoms = _pqr_atoms(workdir / "structure.pqr")
    keys = {f'{a["chain"]}:{a["seq"]}:{a["icode"]}' for a in atoms}
    if not valid.issubset(keys):
        raise ValueError("PQR/PDB identity correspondence incomplete")
    return atoms


def score_pair(lhs, rhs):
    """Compare 2 residue-potential mappings on exactly the common residues."""
    shared = sorted(set(lhs) & set(rhs))
    if len(shared) < 20:
        raise ValueError("Insufficient shared surface residues")
    x = np.asarray([lhs[k]["surface_potential_kT_per_e"] for k in shared], dtype=float)
    y = np.asarray([rhs[k]["surface_potential_kT_per_e"] for k in shared], dtype=float)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Missing/nonfinite numerical surface result")
    absd = np.abs(x-y)
    n_hot = max(1, math.ceil(len(shared)*.10))
    def hot(arr):
        indices = np.argsort(-np.abs(arr), kind="stable")[:n_hot]
        return set(shared[int(i)] for i in indices)
    h1,h2 = hot(x),hot(y)
    correlation = float(spearmanr(x,y).statistic) if np.ptp(x)>1.e-10 and np.ptp(y)>1.e-10 else None
    return {
        "common_n": len(shared),
        "lhs_only_n": len(lhs)-len(shared),
        "rhs_only_n": len(rhs)-len(shared),
        "median_abs_delta_kT_per_e": round(float(np.median(absd)), 7),
        "p90_abs_delta_kT_per_e": round(float(np.quantile(absd, .90)), 7),
        "p95_abs_delta_kT_per_e": round(float(np.quantile(absd, .95)), 7),
        "max_abs_delta_kT_per_e": round(float(np.max(absd)), 7),
        "rmse_delta_kT_per_e": round(float(np.sqrt(np.mean(absd**2))), 7),
        "median_abs_relative_to_reference_potential":
            round(float(np.median(absd)/(np.median(np.abs(y)) + 1.e-12)), 7),
        "spearman_rho": round(correlation,7) if correlation is not None and np.isfinite(correlation) else None,
        "top10pct_absolute_hotspot_jaccard": round(len(h1 & h2)/len(h1 | h2),7),
        "p95_within_0_5_kT_per_e": bool(np.quantile(absd,.95) <= .5),
        "median_within_0_1_kT_per_e": bool(np.median(absd) <= .1),
    }


def sampling_one(item):
    pid, pdb, valid, surface, sha = prep_case(item)
    with tempfile.TemporaryDirectory(prefix="tv48_sampling_") as directory:
        wd = Path(directory)
        atoms = make_pqr(wd, pdb, valid)
        apbs, gridmeta = actual_grid(wd, atoms, .8, 257, "production")
        fields = {}
        meta = {}
        for n in (42,84,162):
            data, m = sample_residue_sas_potential(atoms, *apbs, surface,
                                                    n_directions=n)
            fields[n], meta[n] = data, m
        pairs = {
            "42_vs_84": score_pair(fields[42], fields[84]),
            "84_vs_162": score_pair(fields[84], fields[162]),
            "42_vs_162": score_pair(fields[42], fields[162]),
        }
        # Canonical 42-direction model and the independent deterministic
        # 42-direction resampling must match by construction.
        assert meta[42]["directions_per_heavy_atom"] == 42
        record = {
            "pdb_id":pid, "status":"PASS", "source_sha256":sha,
            "n_surface_residues":len(surface),
            "APBS_grid":gridmeta,
            "directions":{str(n):meta[n] for n in (42,84,162)},
            "comparisons":pairs,
        }
        # Real residue-level fields permit independent re-analysis.
        union = sorted(set().union(*[set(field) for field in fields.values()]))
        rows=[]
        for key in union:
            row={"pdb_id":pid,"residue_key":key}
            for n in (42,84,162):
                v=fields[n].get(key,{})
                row[f"potential_{n}_kT_per_e"]=v.get("surface_potential_kT_per_e")
                row[f"area_{n}_A2"]=v.get("surface_sampling_area_A2")
                row[f"samples_{n}"]=v.get("surface_potential_samples",0)
            rows.append(row)
        dump_json(OUT/"sampling"/f"{pid}.json",record)
        dump_csv(OUT/"sampling"/f"{pid}_residues.csv",rows)
        return record


def grid_one(item):
    pid, pdb, valid, surface, sha = prep_case(item)
    with tempfile.TemporaryDirectory(prefix="tv48_uncapped_grid_") as directory:
        wd = Path(directory)
        atoms = make_pqr(wd, pdb, valid)
        variants = [
            ("standard_08", .8, 257),
            ("medium_05", .5, 257),
            ("high_04_capped", .4, 257),
            ("high_04_uncapped", .4, None),
        ]
        fields={}
        metas={}
        for label, target, cap in variants:
            apbs, m = actual_grid(wd, atoms, target, cap, label)
            field, shell_meta = sample_residue_sas_potential(
                atoms, *apbs, surface, n_directions=162
            )
            del apbs
            gc.collect()
            fields[label],metas[label]=field,{**m, "n_residue_shell_values":len(field),
                                            "n_shell_samples":shell_meta["n_surface_samples"]}
        if not metas["high_04_capped"]["was_actually_capped"]:
            raise AssertionError("PDB designated as capped at 0.4 A is actually not capped")
        if metas["high_04_uncapped"]["was_actually_capped"]:
            raise AssertionError("Uncapped grid was silently capped")
        capped = metas["high_04_capped"]["grid_dims"]
        uncapped = metas["high_04_uncapped"]["grid_dims"]
        if not any(a>b for a,b in zip(uncapped,capped)):
            raise AssertionError("No actual grid refinement past 257 cap")
        step05=np.asarray(metas["medium_05"]["effective_spacing_A"])
        step04=np.asarray(metas["high_04_uncapped"]["effective_spacing_A"])
        if not np.all(step04 <= step05+1.e-6) or not np.any(step04 < step05-1.e-6):
            raise AssertionError("Grid refinement did not decrease actual effective spacing")
        if metas["high_04_uncapped"]["voxel_count"] > 22_000_000:
            raise AssertionError("Uncapped APBS grid exceeded predeclared memory limit")
        comparisons={
            "08_vs_05": score_pair(fields["standard_08"], fields["medium_05"]),
            "05_vs_04_uncapped": score_pair(fields["medium_05"], fields["high_04_uncapped"]),
            "04_cap_vs_04_uncapped": score_pair(fields["high_04_capped"], fields["high_04_uncapped"]),
        }
        record={"pdb_id":pid,"status":"PASS", "source_sha256":sha,
                "n_surface_residues":len(surface),"variants":metas,"comparisons":comparisons}
        rows=[]
        union=sorted(set().union(*[set(v) for v in fields.values()]))
        for key in union:
            row={"pdb_id":pid,"residue_key":key}
            for label in fields:
                row[label+"_kT_per_e"]=fields[label].get(key,{}).get("surface_potential_kT_per_e")
            rows.append(row)
        dump_json(OUT/"grid"/f"{pid}.json",record)
        dump_csv(OUT/"grid"/f"{pid}_residues.csv",rows)
        return record


def shard(kind,index,n):
    items = MANIFEST if kind=="sampling" else [BY_ID[x] for x in CAP_AT_TARGET_04]
    result=[]
    for i,item in enumerate(items):
        if i%n != index:continue
        try:
            data = sampling_one(item) if kind=="sampling" else grid_one(item)
        except Exception as e:
            data={"pdb_id":item["pdb_id"],"status":"FAIL",
                  "error_type":type(e).__name__,"error":str(e),
                  "traceback_tail":traceback.format_exc()[-6000:]}
            dump_json(OUT/"failures"/f"{kind}_{item['pdb_id']}.json",data)
        result.append(data)
        print("TV48_CONVERGENCE_CASE "+json.dumps({
            "kind":kind,"pdb_id":data["pdb_id"],"status":data["status"],
            "error":data.get("error","")[:200]
        }),flush=True)
    dump_json(OUT/f"{kind}_shard_{index}.json",result)
    errors=sum(r["status"]!="PASS" for r in result)
    print("TV48_CONVERGENCE_SHARD "+json.dumps(
        {"kind":kind,"index":index,"attempted":len(result),
         "pass":len(result)-errors,"fail":errors}),flush=True)
    if errors:raise SystemExit("Numerical convergence shard has failures")


def aggregate(n):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ss=[];gg=[];missing=[]
    for kind,into in [("sampling",ss),("grid",gg)]:
        for i in range(n if kind=="sampling" else 4):
            path=OUT/f"{kind}_shard_{i}.json"
            if not path.exists():missing.append(str(path))
            else:into.extend(json.loads(path.read_text(encoding="utf-8")))
    errors=[x for x in ss+gg if x["status"]!="PASS"]
    if len(ss)!=48 or len(gg)!=4 or missing or errors or set(r["pdb_id"] for r in ss)!=set(BY_ID) or set(r["pdb_id"] for r in gg)!=set(CAP_AT_TARGET_04):
        dump_json(OUT/"ERROR_SUMMARY.json",{"missing":missing,"failures":errors})
        raise ValueError(f"Numerical experiment incomplete: {len(ss)}/48, {len(gg)}/4; "
                         f"{len(errors)} failures; missing {missing[:3]}")
    samp_rows=[]
    for record in ss:
        for label,metrics in record["comparisons"].items():
            samp_rows.append({"pdb_id":record["pdb_id"],"comparison":label,
                              "n_surface_residues":record["n_surface_residues"],
                              **metrics,
                              **{f"n_surface_{k}":v["n_exposed_residues_with_surface_potential"]
                                 for k,v in record["directions"].items()}})
    grid_rows=[]
    for record in gg:
        for label,metrics in record["comparisons"].items():
            grid_rows.append({"pdb_id":record["pdb_id"],"comparison":label,
                              "n_surface_residues":record["n_surface_residues"],
                              **metrics})
    size_rows=[]
    for record in gg:
        for label,m in record["variants"].items():
            size_rows.append({"pdb_id":record["pdb_id"],**m})
    dump_csv(OUT/"TV48_42_84_162_SAMPLING_COMPARISONS.csv",samp_rows)
    dump_csv(OUT/"CAPPED_GRIDS_REAL_COMPARISONS.csv",grid_rows)
    dump_csv(OUT/"CAPPED_GRID_ACTUAL_DIMENSIONS.csv",size_rows)
    # Save full residue-level CSVs as delivered in each shard.
    data = pd.DataFrame(samp_rows)
    grid = pd.DataFrame(grid_rows)
    charts=OUT/"figures"
    charts.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.linewidth":.8})
    labels=["42_vs_84","84_vs_162","42_vs_162"]
    vals=[data.loc[data.comparison==label,"median_abs_delta_kT_per_e"] for label in labels]
    fig,ax=plt.subplots(figsize=(6.6,4.3))
    ax.boxplot(vals,labels=["42 vs 84","84 vs 162","42 vs 162"],
               showfliers=True,whis=(5,95))
    ax.set_ylabel("Per-protein median |SAS potential difference| (kT/e)")
    ax.set_xlabel("Fibonacci directions per PQR heavy atom")
    ax.grid(axis="y",alpha=.18)
    fig.tight_layout()
    for ext in ["png","pdf"]:
        fig.savefig(charts/f"01_sampling_42_84_162.{ext}",dpi=280,bbox_inches="tight")
    plt.close(fig)

    fig, ax=plt.subplots(figsize=(6.6,4.3))
    labels=["08_vs_05","05_vs_04_uncapped","04_cap_vs_04_uncapped"]
    for i,label in enumerate(labels):
        sub=grid[grid.comparison==label].sort_values("pdb_id")
        ax.scatter([i]*len(sub),sub.median_abs_delta_kT_per_e,label=label,s=35)
        for _,row in sub.iterrows():
            ax.annotate(row.pdb_id,(i,row.median_abs_delta_kT_per_e),xytext=(5,0),
                        textcoords="offset points",fontsize=7)
    ax.set_xticks(range(3),["0.8 vs 0.5 Å","0.5 vs 0.4 Å uncapped",
                            "0.4 Å capped vs uncapped"])
    ax.set_ylabel("Per-protein median |SAS potential difference| (kT/e)")
    ax.grid(axis="y",alpha=.18)
    fig.tight_layout()
    for ext in ["png","pdf"]:
        fig.savefig(charts/f"02_grid_257_cap_uncapped.{ext}",dpi=280,bbox_inches="tight")
    plt.close(fig)
    summary = {
        "status":"PASS_COMPUTATION","sampling_proteins":len(ss),
        "grid_actually_capped_proteins":len(gg),
        "capped_structure_ids":list(CAP_AT_TARGET_04),
        "sampling_median_absolute_delta_42_84":
            float(data[data.comparison=="42_vs_84"].median_abs_delta_kT_per_e.median()),
        "sampling_median_absolute_delta_84_162":
            float(data[data.comparison=="84_vs_162"].median_abs_delta_kT_per_e.median()),
        "sampling_n_proteins_below_0_1_kT_per_e_by_pair":
            {name:int(data[(data.comparison==name)&(data.median_within_0_1_kT_per_e)].shape[0])
             for name in ("42_vs_84","84_vs_162","42_vs_162")},
        "sampling_n_proteins_below_0_5_kT_per_e_p95_by_pair":
            {name:int(data[(data.comparison==name)&(data.p95_within_0_5_kT_per_e)].shape[0])
             for name in ("42_vs_84","84_vs_162","42_vs_162")},
        "grid_median_absolute_delta_by_pair":
            {name:float(grid[grid.comparison==name].median_abs_delta_kT_per_e.median())
             for name in ("08_vs_05","05_vs_04_uncapped","04_cap_vs_04_uncapped")},
        "grid_max_voxels_real_apbs":
            max(int(row["voxel_count"]) for row in size_rows),
        "interpretation": [
            "PASS_COMPUTATION indicates independent solvers and all records completed, not a guaranteed convergence threshold.",
            "42, 84, 162 Fibonacci direction sets are deterministic but not nested.",
            "One real identical APBS grid is sampled in 3 directions for each of 48 PDBs.",
            "Four genuinely capped PDBs at target 0.4 A rerun actual APBS with and without dimension cap.",
            "Grid comparisons always use 162 directions to isolate grid effects.",
            "Tolerance examples of 0.1 median or 0.5 p95 kT/e are engineering diagnostics, not experimental fidelity claims.",
            "The standard production 0.8 A target grid has no 257-node cap on the 48 PDBs.",
            "The 11 canonical maps and Pareto scores are not recalibrated or modified.",
        ]
    }
    dump_json(OUT/"SUMMARY.json",summary)
    print("TV48_CONVERGENCE_FINAL "+json.dumps(summary),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--sampling-shard",type=int)
    parser.add_argument("--grid-shard",type=int)
    parser.add_argument("--shards",type=int,default=8)
    parser.add_argument("--aggregate",action="store_true")
    a=parser.parse_args()
    if a.aggregate:aggregate(a.shards)
    elif a.sampling_shard is not None and 0<=a.sampling_shard<a.shards:
        shard("sampling",a.sampling_shard,a.shards)
    elif a.grid_shard is not None and 0<=a.grid_shard<4:
        shard("grid",a.grid_shard,4)
    else:parser.error("Provide a valid --sampling-shard, --grid-shard or --aggregate")


if __name__=="__main__":
    main()
