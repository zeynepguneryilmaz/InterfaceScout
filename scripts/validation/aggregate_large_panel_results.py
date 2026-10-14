#!/usr/bin/env python3
from __future__ import annotations
import csv, json, statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "validation/large_panel_285_2026-10-13"
SRC = BASE / "downloaded_artifacts"
OUT = BASE / "results"
OUT.mkdir(parents=True, exist_ok=True)

RUN_ID = 38614780417
CHANNELS = [
    "anionic","cationic","hydrophobic","pi_carbon","hbond_donor",
    "hbond_acceptor","oxide","calcium_phosphate","metal_coord",
    "soft_metal_sulfur","phosphate",
]

def read_many(name):
    files = sorted(SRC.rglob(name))
    if len(files) != 12:
        raise RuntimeError(f"expected 12 {name} files, found {len(files)}")
    rows, fields = [], None
    for p in files:
        with p.open(newline="", encoding="utf-8-sig") as f:
            r = csv.DictReader(f)
            fields = fields or r.fieldnames
            if r.fieldnames != fields:
                raise RuntimeError(f"schema mismatch: {p}")
            rows.extend(r)
    return fields, rows

def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

def fnum(x):
    try: return float(x)
    except: return None

def members(s):
    return {x for x in (s or "").split(";") if x}

def jac(a,b):
    u=a|b
    return None if not u else len(a&b)/len(u)

def median(v):
    v=[x for x in v if x is not None]
    return statistics.median(v) if v else None

def mean(v):
    v=[x for x in v if x is not None]
    return statistics.fmean(v) if v else None

def fmt(x):
    return "" if x is None else f"{x:.6f}"

sf,status = read_many("protein_status.csv")
cf,channels = read_many("channel_summary.csv")
status.sort(key=lambda r:r["panel_id"])
channels.sort(key=lambda r:(r["panel_id"], CHANNELS.index(r["chemistry_channel"])))

if len(status) != 285 or len({r["panel_id"] for r in status}) != 285:
    raise RuntimeError("frozen denominator mismatch")
counts=Counter(r["status"] for r in status)
success={r["panel_id"] for r in status if r["status"]=="SUCCESS"}
if len(channels) != len(success)*11:
    raise RuntimeError("protein x channel row count mismatch")

by=defaultdict(dict)
for r in channels:
    by[r["panel_id"]][r["chemistry_channel"]]=r
for pid in success:
    if set(by[pid]) != set(CHANNELS):
        raise RuntimeError(f"incomplete channel set: {pid}")

specific=[]
pairvals=defaultdict(list)
allj=[]
for pid in sorted(success):
    rr=by[pid]
    sets={c:members(rr[c]["top_patch_members"]) for c in CHANNELS}
    centers=[rr[c]["top_patch_center"] for c in CHANNELS if rr[c]["top_patch_center"]]
    vals=[]; zeros=0; npairs=0
    for i,a in enumerate(CHANNELS):
        for b in CHANNELS[i+1:]:
            v=jac(sets[a],sets[b])
            if v is None: continue
            vals.append(v); pairvals[(a,b)].append(v); allj.append(v)
            npairs+=1; zeros += (v==0)
    specific.append({
        "panel_id":pid,
        "protein_name":rr[CHANNELS[0]]["protein_name"],
        "uniprot_accession":rr[CHANNELS[0]]["uniprot_accession"],
        "plddt_qc_flag":rr[CHANNELS[0]]["plddt_qc_flag"],
        "n_unique_top_patch_centers":len(set(centers)),
        "unique_center_fraction":fmt(len(set(centers))/11),
        "mean_pairwise_top_patch_jaccard":fmt(mean(vals)),
        "median_pairwise_top_patch_jaccard":fmt(median(vals)),
        "fraction_channel_pairs_zero_overlap":fmt(zeros/npairs if npairs else None),
    })

pairs=[]
for i,a in enumerate(CHANNELS):
    for b in CHANNELS[i+1:]:
        v=pairvals[(a,b)]
        pairs.append({
            "channel_a":a,"channel_b":b,"n_proteins":len(v),
            "mean_top_patch_jaccard":fmt(mean(v)),
            "median_top_patch_jaccard":fmt(median(v)),
            "fraction_zero_overlap":fmt(sum(x==0 for x in v)/len(v) if v else None),
        })

write_csv(OUT/"protein_status_all.csv",sf,status)
write_csv(OUT/"channel_summary_all.csv",cf,channels)
write_csv(OUT/"chemistry_specificity_by_protein.csv",
    ["panel_id","protein_name","uniprot_accession","plddt_qc_flag",
     "n_unique_top_patch_centers","unique_center_fraction",
     "mean_pairwise_top_patch_jaccard","median_pairwise_top_patch_jaccard",
     "fraction_channel_pairs_zero_overlap"], specific)
write_csv(OUT/"chemistry_pair_overlap_summary.csv",
    ["channel_a","channel_b","n_proteins","mean_top_patch_jaccard",
     "median_top_patch_jaccard","fraction_zero_overlap"], pairs)

qc=Counter(r["plddt_qc_flag"] for r in status if r["status"]=="SUCCESS")
full=sum(r["status"]=="SUCCESS" and (fnum(r["af_coverage_fraction"]) or 0)>=.95 for r in status)
partial=sum(r["status"]=="SUCCESS" and (fnum(r["af_coverage_fraction"]) or 0)<.95 for r in status)
no_model=[{"panel_id":r["panel_id"],"protein_name":r["protein_name"],
           "uniprot_accession":r["uniprot_accession"]}
          for r in status if r["status"]=="NO_ALPHAFOLD_MODEL"]

channel_stats={}
for c in CHANNELS:
    rr=[r for r in channels if r["chemistry_channel"]==c]
    channel_stats[c]={
        "n_proteins":len(rr),
        "median_n_patches":median([fnum(r["n_patches"]) for r in rr]),
        "median_n_pareto_primary_patches":median([fnum(r["n_pareto_primary_patches"]) for r in rr]),
        "median_top_patch_chemistry_support":median([fnum(r["top_patch_chemistry_support"]) for r in rr]),
        "median_top_patch_coherence":median([fnum(r["top_patch_coherence"]) for r in rr]),
    }

summary={
 "batch_run_id":RUN_ID,
 "source_panel_n":285,
 "status_counts":dict(sorted(counts.items())),
 "successful_proteins_n":len(success),
 "full_or_near_full_success_n":full,
 "partial_structure_success_n":partial,
 "no_model_proteins":no_model,
 "chemistry_channels_n":11,
 "protein_by_channel_rows_n":len(channels),
 "plddt_qc_counts_among_successes":dict(sorted(qc.items())),
 "reference_condition":{"pH":7.4,"ionic_mM_metadata":150.0,"temp_K_metadata":298.0},
 "chemistry_specificity":{
   "proteins_n":len(specific),
   "overall_mean_pairwise_top_patch_jaccard":mean(allj),
   "overall_median_pairwise_top_patch_jaccard":median(allj),
   "median_unique_top_patch_centers_per_protein":median([float(r["n_unique_top_patch_centers"]) for r in specific]),
   "median_unique_center_fraction":median([fnum(r["unique_center_fraction"]) for r in specific]),
   "median_fraction_channel_pairs_zero_overlap_per_protein":median([fnum(r["fraction_channel_pairs_zero_overlap"]) for r in specific]),
 },
 "channel_descriptive_statistics":channel_stats,
 "interpretation":"Large-panel breadth/chemistry-specific landscape only; not adsorption propensity or experimental accuracy.",
 "full_json_outputs":f"GitHub Actions run {RUN_ID}, shard artifacts 0-11",
}
(OUT/"panel_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")

chem=summary["chemistry_specificity"]
readme=f"""# InterfaceScout 285-protein large-panel results

Batch run: GitHub Actions {RUN_ID}

- Frozen source denominator: **285**
- Successful frozen V1 runs: **{len(success)}**
- No AlphaFold DB model: **{counts.get("NO_ALPHAFOLD_MODEL",0)}**
- Runtime analysis errors: **{counts.get("ERROR",0)}**
- Full/near-full structure runs: **{full}**
- Partial structure runs: **{partial}**
- Chemistry channels: **11**
- Protein x chemistry landscapes: **{len(channels)}**

This is a breadth/applicability and chemistry-specific surface-landscape analysis.
It is not an adsorption amount, affinity, corona-membership, or experimental
interface-accuracy benchmark. Experimental region-level accuracy remains the
separate Level-3 validation layer.

AlphaFold pLDDT is QC metadata only and never enters V1 ranking.

Chemistry-specificity (rank-1 patch members across 55 channel pairs/protein):
- overall median Jaccard: **{chem["overall_median_pairwise_top_patch_jaccard"]:.3f}**
- overall mean Jaccard: **{chem["overall_mean_pairwise_top_patch_jaccard"]:.3f}**
- median distinct top-patch centers: **{chem["median_unique_top_patch_centers_per_protein"]:.1f}/11**
- median zero-overlap pair fraction/protein: **{chem["median_fraction_channel_pairs_zero_overlap_per_protein"]:.3f}**

Files:
- protein_status_all.csv
- channel_summary_all.csv
- chemistry_specificity_by_protein.csv
- chemistry_pair_overlap_summary.csv
- panel_summary.json

Complete per-protein JSON outputs remain in the 12 artifacts from Actions run
{RUN_ID}. No protein/channel was selected based on favorable InterfaceScout
output or PC-DB adsorption labels.
"""
(OUT/"RESULTS_README.md").write_text(readme)
print(json.dumps(summary,indent=2,sort_keys=True))

# Aggregation workflow trigger marker; no analytical effect.
