# Regional Chemistry Profile — methodological justification and limitations

**Status:** implemented as **descriptive diagnostics** of existing InterfaceScout combined-map patches, not an adsorption-prediction model. Canonical 11 maps, their pH-dependent support, solvent exposure, 6/9 Å multiscale mapping, 8 Å patch geometry, Pareto fronts and their ranks remain unchanged.

## Scientific rationale (method analogies; not validation claims)

1. Shrake A; Rupley JA (1973), *Journal of Molecular Biology* 79, 351–371. DOI [10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9): solvent-accessible surfaces.
2. Jones S; Thornton JM (1997), *Journal of Molecular Biology* 272, 121–132. DOI [10.1006/jmbi.1997.1234](https://doi.org/10.1006/jmbi.1997.1234): accessibility and physicochemical differences among protein surface patches.
3. Jones S; Thornton JM (1997), *Journal of Molecular Biology* 272, 133–143. DOI [10.1006/jmbi.1997.1233](https://doi.org/10.1006/jmbi.1997.1233): patch-based protein-interface assessment.
4. PEP-Patch (2023), *Journal of Chemical Information and Modeling*. DOI [10.1021/acs.jcim.3c01490](https://doi.org/10.1021/acs.jcim.3c01490): separate, quantitatively characterized surface property patches. **PEP-Patch computes electrostatic potential, which the present InterfaceScout does not.** This is methodological context, not numerical equivalence.
5. RNase A adsorption accessibility study (2014), *Langmuir*. DOI [10.1021/la503854a](https://doi.org/10.1021/la503854a): same protein on silica glass, HDPE and PMMA at varied coverages; labelling readout reports solvent-accessibility changes, NOT directly proven protein-material contact residues.

## Exact calculations

For every existing candidate patch P (residue set) and canonical chemistry channel c, the per-residue input **L_c(r)** is the canonical, nonnegative local score obtained from the same analysis (sidechain relative SASA × compatible residue-class indicator × required-ionization-state availability); if a residue lacks positive support for c its value is 0.

- **Mean local support:** D_c(P) = [Σ_(r in P) L_c(r)] / |P|. Unitless; descriptive, not an energy. No learned channel weighting.
- **Supporting residue fraction:** F_c(P) = |{r in P: L_c(r)>0}| / |P|. Describes fraction of region members, not surface area or binding probability.
- **Weakest channel support:** B(P) = min_(selected c) D_c(P). **Reported only as a supplemental diagnostic**, never used to reorder patches or select user-facing default regions.
- **Distinct-residue cross-channel proximity:** for two channels c and d, min ||CA(r)-CA(s)|| over r≠s in P, L_c(r)>0 and L_d(s)>0. Returns NA if no such pair or coordinates unavailable; CA distance is not an atom–atom contact.
- **Regional contrast:** C_c(P) = D_c(P) − median_{Q in existing combined patch candidates} D_c(Q). **Candidate regions overlap**, so this median is a descriptive spatial reference only, not independent sampling, a z-score, or significance.

Canonical 8 Å neighborhoods use CA distances and radial centroid-based hemisphere filtering, a coarse side-of-protein heuristic. Shared patch membership does NOT imply the two contributing residues are ≤8 Å apart from one another, so their direct pair distance is explicitly recorded. There is no proposed 6 Å or 0.5/0.75 threshold in production diagnostics: any such thresholds in RNase A exploratory tests were descriptive sensitivity checks only.

## Defensible scope

- These features answer: what chemically compatible and exposed residue classes exist together in an *existing protein-side surface candidate region*, how dense are they, and how do they differ from neighboring candidate regions?
- They do NOT answer: will a material bind, which conformation will adsorb, adhesion energy, kinetics, surface-induced unfolding, material-specific affinity, or exact contact residues.
- pH is applied via fixed residue-class protonation rules; actual residue-specific pKa shifts, local metal surface chemistry, ionic shielding, solvent, glycosylation and protein-protein adsorption are not modeled.
- Several map channels are chemically correlated or even partly redundant; supporting-channel count must not be translated into number of independent pieces of evidence. Separate channels should remain inspectable.
- The same static protein structure and pH necessarily produce the same canonical maps regardless of material identity. Only a *user-selected chemistry view* changes.
- Spatial contrast is not an experimental validation statistic. The published RNase A observations informed creation of regional diagnostics, hence they cannot serve as an independent holdout.

## Test protocol and future evidence

**Code-level requirements**: analytic toy systems for hand-checked D and F, fixed coordinates for a known cross-channel distance, no-coordinate NA behavior, residue numbering and chain preservation, no effect of arbitrary experimental-label keys, and exact preservation of combined patch order and existing per-channel Pareto front/rank. Cross-platform CI must pass.

**Label-free robustness**: compare fields across alternate experimental structures of the same protein, sampling of SASA points (200/500/1000), 8 Å neighborhood radius sensitivity at 6/7/9/10 Å only as an explicit sensitivity study (not optimized on adsorption labels), and duplicated/chemically correlated maps.

**Independent biological relevance**: curate experiments involving the same folded protein across two or more different independently characterized material surfaces. Pre-register chemistry map selection from surface chemistry only, before outcome access; make residue-number mapping auditable. Score equal surface-residue budgets; compare to SASA-only, individual chemistry maps, randomized 8 Å same-face patches, and other co-localization descriptors; use within-paper matched comparisons and protein-/study-group bootstrap. Distinguish direct experimental contacts from footprinting, epitope exposure, and adsorption-induced folding. Record negatives and failures.

## Public UI

Combined maps now show **mean local support, supporting fraction and minimal distinct-residue CA distance** per existing region, with explicit warning that these do not affect patch rankings or constitute binding scores. The backend /combine response includes complete per-region regional_profile diagnostics and method metadata. Proximity requires the surface CA coordinates recorded in the analysis response. No external dependence and no fitting were added.
