#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math, statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
R=ROOT/"validation/large_panel_285_2026-10-13/results"

def read_csv(name):
    with (R/name).open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def rankdata(vals):
    order=sorted(range(len(vals)), key=lambda i: vals[i])
    ranks=[0.0]*len(vals); i=0
    while i<len(order):
        j=i+1
        while j<len(order) and vals[order[j]]==vals[order[i]]:
            j+=1
        r=(i+1+j)/2.0
        for k in range(i,j): ranks[order[k]]=r
        i=j
    return ranks

def pearson(x,y):
    mx=sum(x)/len(x); my=sum(y)/len(y)
    dx=[a-mx for a in x]; dy=[b-my for b in y]
    den=(sum(a*a for a in dx)*sum(b*b for b in dy))**0.5
    return None if den==0 else sum(a*b for a,b in zip(dx,dy))/den

def spearman(x,y):
    return pearson(rankdata(x),rankdata(y))

def median(xs):
    return statistics.median(xs) if xs else None

status=read_csv("protein_status_all.csv")
spec=read_csv("chemistry_specificity_by_protein.csv")
pairs=read_csv("chemistry_pair_overlap_summary.csv")

by_status={r["panel_id"]:r for r in status if r["status"]=="SUCCESS"}
rows=[]
for r in spec:
    s=by_status[r["panel_id"]]
    rows.append({
        "panel_id":r["panel_id"],
        "length":float(s["uniprot_length"]),
        "surface":float(s["n_surface_residues"]),
        "mean_j":float(r["mean_pairwise_top_patch_jaccard"]),
        "zero":float(r["fraction_channel_pairs_zero_overlap"]),
        "centers":float(r["n_unique_top_patch_centers"]),
        "qc":s["plddt_qc_flag"],
    })

size_qc={}
for pred in ("length","surface"):
    size_qc[pred]={}
    for metric in ("mean_j","zero","centers"):
        size_qc[pred][metric]=spearman(
            [r[pred] for r in rows],[r[metric] for r in rows]
        )

groups={}
for label,keep in [
    ("STANDARD",lambda q:q=="STANDARD"),
    ("FLAGGED",lambda q:q!="STANDARD"),
]:
    rr=[r for r in rows if keep(r["qc"])]
    groups[label]={
        "n":len(rr),
        "median_uniprot_length":median([r["length"] for r in rr]),
        "median_surface_residues":median([r["surface"] for r in rr]),
        "median_mean_pairwise_jaccard":median([r["mean_j"] for r in rr]),
        "median_zero_overlap_fraction":median([r["zero"] for r in rr]),
        "median_unique_top_patch_centers":median([r["centers"] for r in rr]),
    }

pair_rows=[{
    "a":r["channel_a"],"b":r["channel_b"],
    "mean":float(r["mean_top_patch_jaccard"]),
    "median":float(r["median_top_patch_jaccard"]),
    "zero":float(r["fraction_zero_overlap"]),
} for r in pairs]
highest=sorted(pair_rows,key=lambda r:r["mean"],reverse=True)[:8]
lowest=sorted(pair_rows,key=lambda r:r["mean"])[:8]

qc={
    "analysis_type":"post-hoc descriptive QC; no thresholds or weights fitted",
    "n_successful_proteins":len(rows),
    "spearman_size_associations":{
        "uniprot_length_vs_mean_pairwise_jaccard":size_qc["length"]["mean_j"],
        "uniprot_length_vs_zero_overlap_fraction":size_qc["length"]["zero"],
        "uniprot_length_vs_unique_top_patch_centers":size_qc["length"]["centers"],
        "surface_residues_vs_mean_pairwise_jaccard":size_qc["surface"]["mean_j"],
        "surface_residues_vs_zero_overlap_fraction":size_qc["surface"]["zero"],
        "surface_residues_vs_unique_top_patch_centers":size_qc["surface"]["centers"],
    },
    "alphafold_confidence_strata":groups,
    "highest_mean_overlap_channel_pairs":highest,
    "lowest_mean_overlap_channel_pairs":lowest,
    "key_interpretation":[
        "Cross-channel zero-overlap increases with protein size/surface area, so raw zero-overlap fraction must not be treated as size-independent validation.",
        "Distinct top-patch-center count is only weakly associated with protein size and is a more conservative breadth descriptor.",
        "Chemically related channels show much higher top-patch overlap than chemically dissimilar channels.",
        "At pH 7.4, cationic and H-bond-donor rank-1 patches are identical across this panel; the 11 outputs should therefore not be described as 11 statistically independent dimensions.",
        "AlphaFold-confidence strata differ in protein size as well as chemistry-specificity metrics; confidence-stratified differences are descriptive and not causal.",
    ],
}
(R/"posthoc_qc.json").write_text(json.dumps(qc,indent=2,sort_keys=True)+"\n",encoding="utf-8")

s=qc["spearman_size_associations"]
std=groups["STANDARD"]; flg=groups["FLAGGED"]
def p(x): return f"{x:.3f}"
md=f"""# Post-hoc QC for the 285-protein large panel

This QC was performed after the frozen 285-protein batch and is **descriptive only**.
No threshold, chemistry weight, protein selection, or V1 parameter was changed based on these results.

## Protein-size associations

Spearman correlations:

- UniProt length vs mean pairwise top-patch Jaccard: **{p(s["uniprot_length_vs_mean_pairwise_jaccard"])}**
- UniProt length vs zero-overlap channel-pair fraction: **{p(s["uniprot_length_vs_zero_overlap_fraction"])}**
- UniProt length vs number of unique top-patch centers: **{p(s["uniprot_length_vs_unique_top_patch_centers"])}**
- Surface-residue count vs mean pairwise top-patch Jaccard: **{p(s["surface_residues_vs_mean_pairwise_jaccard"])}**
- Surface-residue count vs zero-overlap channel-pair fraction: **{p(s["surface_residues_vs_zero_overlap_fraction"])}**
- Surface-residue count vs number of unique top-patch centers: **{p(s["surface_residues_vs_unique_top_patch_centers"])}**

Therefore, the raw cross-channel zero-overlap fraction is partly size-dependent and should **not** be presented as a size-independent validation statistic. The number of distinct top-patch centers has a substantially weaker size association.

## Chemistry-pair structure

Highest mean rank-1-patch overlaps include:

{chr(10).join(f'- {r["a"]} vs {r["b"]}: mean Jaccard **{r["mean"]:.3f}**, median **{r["median"]:.3f}**' for r in highest[:5])}

Lowest mean overlaps include:

{chr(10).join(f'- {r["a"]} vs {r["b"]}: mean Jaccard **{r["mean"]:.3f}**, zero-overlap fraction **{r["zero"]:.3f}**' for r in lowest[:5])}

Notably, cationic and H-bond-donor top patches are identical throughout this pH 7.4 panel. The 11 canonical outputs are therefore **not 11 statistically independent dimensions**; related channels can encode strongly overlapping protein-side chemistry. This is retained transparently rather than removing a channel after viewing results.

## AlphaFold-confidence strata

- STANDARD models: n={std["n"]}, median length={std["median_uniprot_length"]:.0f}, median surface residues={std["median_surface_residues"]:.0f}, median mean-Jaccard={std["median_mean_pairwise_jaccard"]:.3f}
- Flagged models: n={flg["n"]}, median length={flg["median_uniprot_length"]:.0f}, median surface residues={flg["median_surface_residues"]:.0f}, median mean-Jaccard={flg["median_mean_pairwise_jaccard"]:.3f}

Because the flagged set also contains larger proteins on average, these differences are treated as QC context rather than evidence of a confidence effect.

## Reporting consequence

For the manuscript, the safest large-panel claim is:

> Across a frozen 285-protein source panel, InterfaceScout was executable for 282 proteins and generated chemistry-dependent candidate surface landscapes across all 11 canonical outputs. Related chemistry channels showed strong overlap, whereas dissimilar channels generally prioritized different regions. Cross-channel differentiation was partly protein-size dependent and is therefore treated as a descriptive applicability result rather than an accuracy metric.

Experimental accuracy remains supported only by the separate Level-3 interface-localization benchmark.
"""
(R/"POSTHOC_QC.md").write_text(md,encoding="utf-8")
print(json.dumps(qc,indent=2,sort_keys=True))
