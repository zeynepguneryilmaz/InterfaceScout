# InterfaceScout held-out external validation protocol

Date frozen: 2026-09-04

This protocol must be finalized before scoring any newly selected held-out system. Systems used previously during development or diagnostics are excluded from the held-out set: 4F5S, 5H7A, 1MBN, 2PTN, 2HHB.

## Model to be evaluated

The model evaluated in the final held-out benchmark must be the post-structural-correction frozen v5.2 core. No chemistry memberships, scRSA threshold, SASA settings, patch radii, score formula, geometry representation or benchmark-specific parameters may be changed after the held-out residue/contact labels are collected.

Canonical quantities remain:

- `L_i,c = I_i,c * scRSA_i * f_state,i,c(pH)`
- within-map propensity `P_i,c`
- Cα-based 5/8 Å multiscale persistence `M_i,c`

Optional/auxiliary descriptors must not be silently folded into the canonical score.

## Predeclared mapping rule

The InterfaceScout chemistry class for each benchmark is chosen from the physical surface description in the source paper before looking at model performance. Mixed materials may be represented by a predeclared mechanism profile, but individual chemistry maps remain separate; no fitted weighted sum is permitted.

## Ground truth hierarchy

Preference order:

1. exact residue-level persistent/contact anchors reported in the original paper or SI;
2. exact residue sets defined from a published contact-frequency threshold;
3. explicitly reported binding regions when exact residues are unavailable.

Region-level labels are secondary validation and must not be pooled as if they were exact residue anchors.

## Primary metrics

For local residue propensity `P`:

- AUROC
- average precision (AP)
- prevalence-normalized AP (`AP / positive prevalence`)
- precision@k and recall@k for predeclared k values

For spatial persistence `M`:

- anchor-to-nearest-top-patch Cα distance
- fraction of anchors recovered within 5 Å and 8 Å of top-ranked patch centers
- exact-residue AUROC/AP as secondary metrics

## Baselines

At minimum:

- scRSA alone
- chemistry-membership × scRSA without patch persistence

The canonical model should be compared with accessibility alone on the same evaluated surface-residue universe.

## Significance / uncertainty

- exposure-matched permutation test for anchor enrichment / ranking significance
- bootstrap confidence intervals for AUROC, AP and spatial-recovery summaries where sample size permits
- report per-system results; do not hide heterogeneous failures inside a pooled mean

## Applicability annotations

Each benchmark must be annotated for physics not represented by the lightweight core, including as applicable:

- strong global electrostatic steering
- explicit hydration/water-layer penetration
- adsorption-induced conformational rearrangement
- oligomeric shielding
- multi-protein crowding
- heterogeneous/mixed surface chemistry
- covalent or highly specific surface chemistry

These annotations are interpretation boundaries, not post hoc exclusions.

## First eligible held-out benchmark identified

### Hen egg-white lysozyme on a negatively charged hydrophilic surface

Source: Kubiak-Ossowska & Mulheran, Langmuir 2010, DOI 10.1021/la102960m.

Published exact adsorption anchors:

- major site: Lys1, Arg5, Arg125, Arg128 (Arg128 described as crucial)
- minor site: Arg68

Predeclared InterfaceScout surface class: `anionic`.

This system is eligible because it was not used in prior InterfaceScout model/radius/geometry/steering development. Before scoring, the exact protein structure/PDB used by the source study and the experimental/simulation pH and ionic strength must be verified from the original article/SI or an unambiguous methods source. If the exact structure cannot be verified, the system is retained as literature evidence but not used as an exact-structure benchmark.

## Freeze rule

Once the v5.2 structural corrections are completed and regression-tested, create a new immutable freeze branch/tag and record its commit SHA here. Only that frozen SHA may be used for the final held-out runs.
