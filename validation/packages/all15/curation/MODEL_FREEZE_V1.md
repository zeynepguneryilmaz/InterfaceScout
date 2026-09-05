# InterfaceScout primary-score freeze v1

Freeze date: 2026-09-05
Branch: `validation-study`
Status: **development model frozen; prospective validation not yet performed**

## Decision

The primary public InterfaceScout score is frozen to the chemistry-aware static multiscale patch score (the existing Core IS). GNM, APBS, radial prominence, and RIN descriptors remain available only as auxiliary/explanatory outputs and are **not multiplicative factors in the primary IS score**.

This decision is based on the source-corrected curated development benchmark, not on the older 5-protein exploratory aggregate.

## Frozen primary score

For residue i and chemistry channel c:

1. Local chemistry/state/exposure contribution

`L_i,c = I_i,c * scRSA_i * f_state,i,c(pH)`

2. Chemistry-normalized local residue score

`P_i,c = 100 * L_i,c / max_j(L_j,c among favorable residues)`

3. Local patch density at radius R

`D_i,c(R) = sum_j L_j,c for d_ij <= R`

4. Frozen InterfaceScout score

`IS_i,c = 100 * min( normalized D_i,c(5 A), normalized D_i,c(8 A) )`

The 5/8 A minimum is interpreted as a conservative multiscale persistence gate: a residue is favored when its chemistry-compatible exposed neighborhood persists at both local scales.

## Frozen parameters

- side-chain relative accessibility threshold: `scRSA >= 0.05`
- SASA probe radius: `1.40 A`
- canonical SASA sampling: `200 points/atom`
- patch radii: `5 A` and `8 A`
- pH-dependent state weighting: existing InterfaceScout chemistry/state rules
- 8 A context output: auxiliary only, not a separate primary multiplier
- no fitted regression coefficients
- no protein-specific learned weights

## Evidence used for model selection

Primary residue-level development set: 12 conditions across 7 unique proteins.

Proteins:
1. Fibronectin III8-10
2. Cytochrome c
3. Lysozyme
4. Alpha-chymotrypsin
5. Myoglobin
6. WNV E protein domain III
7. RNase A

Repeated conditions on the same protein are not treated as independent statistical units. Protein-macro summaries and protein-clustered bootstrap (10,000 resamples) were used.

Canonical structure corrections include:
- 1FNF chain A residues 1236-1509 only for FNIII8-10.
- 3NWV chain A only for cytochrome c.
- 4CHA chains A+B+C as one proteolytically split alpha-chymotrypsin molecule.
- 2HG0 chain A PDB residues 8-109 only for WNV E-domain III.

## Protein-level model-selection evidence

Protein-macro median metrics:

| Model | AUROC | AP | Recall@5 | Recall@10 | Recall@20 | Spatial@10A (top10) | median nearest top10 (A) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Core | 0.7111 | 0.2820 | 0.1429 | 0.3000 | 0.3636 | 0.6429 | 7.176 |
| Core + GNM | 0.7262 | 0.1432 | 0.0000 | 0.1111 | 0.3333 | 0.5000 | 7.453 |
| Core + APBS | 0.7111 | 0.2116 | 0.1429 | 0.3000 | 0.3636 | 0.6364 | 7.176 |
| Core + GNM + APBS | 0.6981 | 0.1432 | 0.0000 | 0.1538 | 0.3571 | 0.6364 | 6.605 |
| Core + GNM + APBS + radial (old IS candidate) | 0.6399 | 0.1274 | 0.0000 | 0.1538 | 0.3077 | 0.5714 | 7.835 |
| + RIN degree | 0.6607 | 0.1457 | 0.0000 | 0.1538 | 0.3077 | 0.6154 | 8.942 |
| + RIN closeness | 0.6435 | 0.1265 | 0.0000 | 0.1538 | 0.3077 | 0.5000 | 9.189 |
| + RIN betweenness | 0.6458 | 0.1373 | 0.0000 | 0.1818 | 0.3077 | 0.6154 | 7.835 |

GNM increases protein-macro median AUROC slightly, but decreases AP and early/spatial recovery. Relative to Core, GNM improves AUROC in 5/7 proteins and worsens 2/7, but Spatial@10A worsens in 6/7 proteins. Its paired protein-level Spatial@10A median delta is -0.1667 with bootstrap 95% CI [-0.2727, -0.1364]. This is incompatible with making GNM a mandatory primary-score multiplier for a residue hotspot-localization tool.

Conditional APBS has useful case-specific behavior for genuinely charged surfaces, but does not improve protein-macro median AUROC over Core and reduces protein-macro AP in this set. Citrate-capped AuNP is particularly problematic when represented as a generic homogeneous anionic surface. APBS therefore remains an auxiliary surface-conditional descriptor rather than a primary multiplier.

Radial prominence was promising in the earlier 5-protein exploratory benchmark but fails protein-level generalization after structure/source curation. The radial-containing candidate falls from Core AUROC 0.7111 to 0.6399 at the protein-macro median and loses early/spatial recovery. It is rejected from the frozen primary score.

RIN modifiers do not recover the lost top-k performance consistently and are rejected from the primary score. RIN degree/closeness/betweenness remain auxiliary descriptors for mechanistic analysis.

## Chemistry-specific negative control

For alpha-chymotrypsin/CNT, the chemistry-matched hydrophobic Core score gives AUROC 0.7225 and AP 0.1932. Applying the deliberately wrong pi-carbon channel to the same ground truth gives Core AUROC 0.3702 and AP 0.0713 with zero Recall@5/10/20.

Adding dynamics/geometry/topology partially rescues the intentionally wrong chemistry (wrong-channel AUROC rises into approximately 0.54-0.64), which demonstrates an undesirable loss of chemistry specificity. This is an additional reason not to let generic geometry/dynamics dominate the chemistry-aware primary score.

## pH robustness before freeze

For source-unresolved primary conditions, Core was tested at pH 6.5, 7.0 and 7.5:
- alpha-chymotrypsin/CNT: identical AUROC, AP, Recall@10 and identical top-10 across the three pH values.
- WNV/graphene: identical AUROC, AP, Recall@10 and identical top-10 across the three pH values.
- myoglobin/citrate-AuNP: AUROC changes by only -0.0123 at pH 6.5 and +0.0130 at pH 7.5 relative to pH 7.0; Recall@10 is unchanged.

The Core selection is therefore not driven by the pH 7 development assumption in these systems.

## Scope and limitations frozen with v1

The primary score predicts chemistry-compatible exposed residue patches. It does not claim to predict:
- adsorption free energy or adsorption capacity;
- a unique rigid-body orientation;
- nanoparticle curvature/radius effects;
- conformational adaptation upon adsorption;
- ligand-shell-specific physics when a material is mapped only to a generic chemistry class.

Particle-radius stress tests confirm that 4-nm and 11-nm silica cannot be distinguished by the current score because radius is not an input. This is a scope limitation, not a parameter to fit after seeing the benchmark.

## Auxiliary outputs retained

The software may calculate and report, but must not multiply into primary IS v1:
- GNM mobility;
- APBS residue potential and chemistry-gated electrostatic compatibility;
- radial prominence;
- RIN degree, closeness and betweenness;
- other explicitly labelled diagnostic/context descriptors.

These should be described as auxiliary mechanistic descriptors, not as parts of the frozen IS v1 score.

## Prospective-validation lock

All proteins and outcomes used above are development evidence because benchmark results were inspected during feature selection. After this freeze, new residue-level benchmark proteins/material systems used for prospective validation must not alter the frozen formula, parameters, chemistry ontology, or cutoffs. If a future model version changes those rules, it must be named/versioned separately and evaluated on a new untouched holdout.
