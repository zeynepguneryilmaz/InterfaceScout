# ISV76 — Experimental anchor and region-coverage audit (2026-11-24)

This is a **frozen-model, outcome-neutral descriptive analysis** of the selected active ISV76 structures. It is NOT an adsorption-pose classifier or quantitative experimental-site validation. All 11 chemistry channels, Top-k patch rules, Pareto objectives and the APBS–SAS settings are unchanged from TV48.

## Audited experimental systems

### 7KWI — Staphylococcus epidermidis R2ab / polystyrene

Original paper: [Yadav et al., J. Biol. Chem. 2026, 302, 113533](https://doi.org/10.1016/j.jbc.2026.113533) (published online September 8, 2026). PRE measurements and subsequent mutation-dependent adsorption-associated biofilm assay support **S776** and **K779** as functionally significant. The paper tested Y722A, Y800A, Y843A, Y844A too: these **did not show the same increased-biofilm effect** as S776A/K779A. A null biofilm phenotype is NOT proof of absence of contact. Source PDB explicitly identifies 7KWI and its author residue numbering corresponds to frozen chain A.

| Site | Experiment | tRSA | SAS–APBS (kT/e) | Membership of Top-5 candidate patch unions |
|---|---|---:|---:|---|
| S776 | altered biofilm assay after alanine mutation | 0.48615 | +1.30805 | phosphate |
| K779 | altered biofilm assay after alanine mutation | 0.71068 | +2.04016 | anionic; hydrophobic; hbond_acceptor |
| Y722 | assayed, no corresponding biofilm increase | 0.64977 | +1.03302 | 7 of 11 channels |
| Y800 | assayed, no corresponding biofilm increase | 0.72655 | +1.38270 | 0 of 11 channels |
| Y843 | assayed, no corresponding biofilm increase | 0.65240 | +1.43065 | 3 of 11 channels |
| Y844 | assayed, no corresponding biofilm increase | 0.33211 | +1.73667 | 3 of 11 channels |

Frozen atlas Cα(S776)–Cα(K779) separation: **10.06 Å**; this exceeds the program's nominal 8 Å local-neighborhood radius, so separate regions are not necessarily inconsistent with the experimental surface.

Top-1 union of **all eleven maps**: 60/126 exposed residues, 1/2 experimentally significant positions. Top-3: 108/126 exposed and 2/2 positions. Top-5: 120/126 exposed and 2/2 positions, **95.2% surface coverage**. Therefore, the final 2/2 count is not strong spatial selection evidence.

### 2OED — GB3 / citrate-coated AuNPs

Source evidence: [2016 NMR mutation study, PMC5096844](https://pmc.ncbi.nlm.nih.gov/articles/PMC5096844) reports decreased AuNP surface coverage for **K4A, K13A, K50A**. [Nature Communications 2022](https://doi.org/10.1038/s41467-022-34749-w) uses experimental K13X competition NMR and explicitly lists **PDB 2OED**. These are position-specific affinity/function changes, not an atomically resolved surface contact map. Exact residue identities for 2OED chain A were verified in InterfaceScout.

| Site | tRSA | SAS–APBS (kT/e) | Number of Top-5 maps containing residue |
|---|---:|---:|---:|
| K4 | 0.33182 | +0.99687 | 7/11 |
| K13 | 0.62196 | +0.70618 | 8/11 |
| K50 | 0.22114 | +0.80674 | 9/11 |

Frozen atlas Cα separations: K4–K13 13.57 Å, K4–K50 6.52 Å, K13–K50 **18.75 Å**. Distinct local patches can be expected for some positions.

Top-1 11-map union: 32/48 exposed residues and 1/3 experimentally significant positions. Top-3: 46/48 exposed and 3/3 positions. Top-5: **48/48 exposed and 3/3 positions**. At Top-5, the all-map union therefore has **zero spatial selectivity** even though all three reference sites are included. No performance claim can be drawn from the 3/3 count.

## Interpretation rules and coverage controls

- Evaluate all 11 maps, including maps that do not contain a study residue. Do not cherry-pick a winning channel from test outcomes.
- Report **surface fraction covered by the patch union** alongside any site coverage. Also distinguish *patch member* from *chemically contributing seed*.
- An annotation from PRE, HDX, mutation, footprinting, activity, or biofilm assays is not interchangeable with a true atomic contact label.
- Do not treat an unlabeled/unassayed residue as a true negative; thus sensitivity, specificity, precision and AUROC are NOT currently established.
- The model is protein-side chemistry mapping, not a polystyrene- or AuNP-specific adsorption model.
- A local 8 Å patch radius does not imply all experimentally implicated sites must belong to one patch.

## Other six source-relevant proteins

- **2ACY (acylphosphatase/silica):** [Langmuir 2013](https://doi.org/10.1021/la401985d) reports an experimental interaction surface distributed over two alpha-helices; exact species/isoform, original residue numbering and the frozen bovine 2ACY correspondence require auditing.
- **7RSA, 1HEL, 1HRC:** [Nano Letters 2012](https://doi.org/10.1021/nl2044524) reports particle-diameter-dependent chemical modification/proteomics of RNase A, lysozyme and cytochrome c. The full individual residue labels and particle-size conditions were not extracted here; 1HRC also has a heme-omission caveat.
- **1A6M:** HDX-MS structural changes are indirect and should not be converted into direct site ground truth.
- **5B1E:** AuNP-associated functional/activity data are not residue-resolved atomic contact labels.

For these six, all 11 maps have produced actual fixed-model patch results, but they have **not** been assigned experimental site-recovery scores. The 76-member technical run remains distinct from an accepted independent experimental site-label panel.

## Reproducible numerical evidence

The conversation artifact `ISV76_Experimental_Anchor_And_Scientific_Region_Audit.zip` contains:
- `EXPERIMENTAL_RESIDUE_TO_MAP_MEMBERSHIP.csv`: each source-annotated residue, tRSA/scRSA, SAS potential, Top-1/3/5 channel membership and chemistry seed membership.
- `EXPERIMENTAL_TOP1_TOP3_TOP5_CHANNEL_COVERAGE.csv`: all 11 channels, exact exposed-area/residue coverage denominators and observed anchor inclusion.
- `EXPERIMENTAL_TOP1_TOP3_TOP5_ALL_CHANNEL_UNION.csv`: correct coverage comparison for pooled maps.
- `EXPERIMENTAL_SITE_PAIR_CA_DISTANCES.csv`: Euclidean Cα distances from the verified atlas's actual exported original PDB.
- Eight experiment-relevant proteins' 88 Top-1 channel regions, all-76 class-level summary, four scientific plots, source audit queue, full report, and exact standalone Python analysis implementation.

**No InterfaceScout parameters were modified.** These are observational, structurally mapped results; a broad site union should not be reinterpreted as validated prediction accuracy.
