# InterfaceScout V1 frozen benchmark scoring protocol — 2026-10-07

This scoring protocol was committed **before examining any new benchmark performance results**. It consumes the frozen benchmark on `benchmark-audit-freeze-2026-10-06` without changing eligibility, labels, source rows, or material mappings.

## Model

- InterfaceScout public V1 code and frozen numerical settings are used without modification:
  - Shrake–Rupley probe 1.40 Å
  - 200 points/atom
  - scRSA >= 0.05
  - multiscale radii 6/9 Å
  - coarse patch radius 8 Å
- PRS, harmonic ANM, and future CG/flexibility models are not used.
- No benchmark label is used as an input to V1.
- No V1 result can change benchmark inclusion/exclusion.

## Structure-enrichment policy for LEVEL 1

Structure assignment is a downstream scoring-enablement step and does not change benchmark membership.

1. Resolve source protein identity independently of InterfaceScout output.
2. Payne single-accession protein groups use their source UniProt accession directly.
3. Payne multi-accession groups remain unresolved for primary scoring; no arbitrary member is chosen.
4. PC-DB name-only proteins are queried against UniProt for Homo sapiens. An accession is accepted only when the source protein name has a unique exact normalized match, or a unique reviewed exact normalized match. Ambiguous/no-match entries remain unresolved.
5. Experimental structure is preferred only when PDBe reports a mapped PDB chain with UniProt coverage >= 0.70 and resolution <= 4.0 Å (or resolution absent for a non-diffraction structure). Among eligible candidates, choose maximum coverage, then best available resolution, then lexical PDB/chain order.
6. If no qualifying experimental structure is available, use the AlphaFold Database model returned for the resolved UniProt accession and label it `alphafold_predicted`.
7. Experimental and predicted structures are reported separately in all primary LEVEL-1 aggregate analyses and may additionally be shown as a pooled sensitivity analysis.
8. Failed downloads, non-standard/empty prepared structures, or unresolved identities are retained in the audit but not scored.

## Environmental policy

- Use source-reported pH when available.
- If pH is not reported, use the frozen V1 reference/default pH 7.4 and flag `pH_source=reference_default_missing_source_pH`.
- Ionic strength defaults to 150 mM only when unavailable; V1 retains ionic strength as metadata and does not explicitly model electrostatic screening.
- Temperature defaults to 298 K only when unavailable; V1 retains temperature as metadata and does not explicitly model thermal response.
- These defaults are fixed before scoring and are not tuned to benchmark performance.

## Material chemistry for LEVEL 1

- Only experiments with a pre-frozen material mapping are chemistry-conditioned in the primary LEVEL-1 analysis.
- `ambiguous_unmapped` experiments remain in the benchmark but are not assigned a V1 chemistry channel post hoc.
- Multi-channel surfaces use the public V1 composite policy: **no weights** and no fitted linear combination.
- For each descriptor below, the composite descriptor is the maximum value across the pre-mapped canonical chemistry channels.

Canonical mapping to public V1 keys:
- anionic -> anionic
- cationic -> cationic
- hydrophobic -> hydrophobic
- pi/aromatic -> pi_carbon
- H-bond donor -> hbond_donor
- H-bond acceptor -> hbond_acceptor
- metal-oxide -> oxide
- calcium/phosphate charged sites -> calcium_phosphate
- transition-metal coordination -> metal_coord
- soft-metal sulfur affinity -> soft_metal_sulfur
- phosphate-rich -> phosphate

## LEVEL-1 V1 descriptors

No new weighted adsorption score is created. The following V1-native patch descriptors are tested separately:

1. `max_front1_chemistry_support`: maximum chemistry-support fraction among Pareto-front-1 patches across mapped channels.
2. `max_front1_patch_coherence`: maximum patch coherence among front-1 patches across mapped channels.
3. `max_front1_mean_accessibility`: maximum mean accessibility among front-1 patches across mapped channels.
4. `max_front1_orientation_coherence`: maximum orientation coherence among front-1 patches across mapped channels.
5. `n_front1_patches`: total number of front-1 patches across mapped channels.
6. `max_residue_propensity`: maximum residue propensity across mapped channels.

No descriptor is selected or discarded after viewing its performance.

## LEVEL-1 occurrence analysis

To avoid pseudo-replication from the >900k observation matrix:

- Evaluate discrimination **within each experiment**, not by one naive row-level p-value.
- An experiment is evaluable for a descriptor only if it has both present and non-positive observations and at least 10 scored protein entities.
- Compute within-experiment AUROC for each descriptor.
- Summarize across experiments by median, IQR, number evaluable, fraction AUROC > 0.5, and two-sided Wilcoxon signed-rank test versus 0.5 when defined.
- Report PC-DB and Payne separately, experimental-structure and AlphaFold subsets separately, and a clearly labeled pooled sensitivity analysis.

## LEVEL-1 abundance analysis

Because quantitative units differ among studies/experiments:

- Never pool raw abundance values across experiments.
- Within each experiment, compute Spearman correlation between each descriptor and quantitative abundance.
- Primary abundance correlation uses source-positive/detected observations only to avoid treating heterogeneous zero/non-positive encodings as a common quantitative scale.
- Require at least 10 scored positive observations and at least 3 distinct abundance values.
- Summarize experiment-level Spearman rho by median, IQR, number evaluable, fraction rho > 0, and Wilcoxon signed-rank test versus 0 when defined.
- PC-DB and Payne are reported separately.

## LEVEL-2 structural-response analysis

Predefined chemistry mapping:
- polystyrene, no modification -> hydrophobic + pi/aromatic
- iron oxide, no modification -> metal-oxide

Labels remain direction-specific:
- `protected/contact_candidate` is contact-compatible but not an unqualified direct-contact label.
- `newly_exposed/conformational_response` is **never** treated as a direct-contact label.

For numeric regions, compute:
- first direct-overlap patch rank;
- first patch rank within 5 Å;
- first patch rank within 8 Å;
- top-3 and top-5 near-8 Å hit/recall;
- best Pareto-front-1 near-8 Å recall;
- mean and maximum V1 residue propensity in the annotated region.

Terminus annotations without a numeric endpoint are retained but are not converted to invented residue ranges.

Protected and newly exposed records are summarized separately.

## LEVEL-3 contact/localization analysis

Predefined exposed-surface chemistry:
- PAA-coated Fe3O4 -> anionic + H-bond acceptor
- citrate-coated Au -> anionic + H-bond acceptor
- silica -> metal-oxide
- MPA-coated Au -> anionic + H-bond acceptor
- AlumOH at the reported positive-surface condition -> cationic + metal-oxide
- unmodified silica/BSA -> metal-oxide

For each frozen localization record, compute the same spatial concordance family used in the publication validation:
- first direct-overlap rank;
- first near-5 Å rank;
- first near-8 Å rank;
- top-3 near-8 Å hit/recall;
- top-5 near-8 Å hit/recall;
- top-5 direct-overlap hit/recall;
- best Pareto-front-1 near-8 Å recall.

Residue-anchor, residue-set, peptide-region, and broader-region evidence are never silently treated as the same experimental resolution; results are stratified by resolution.

## Multiple testing and interpretation

- LEVEL 1 is association testing, not a claim that InterfaceScout predicts absolute adsorption amount.
- Report effect sizes and distribution summaries; do not interpret tiny p-values from massive N as sufficient validation.
- False-discovery-rate adjusted p-values may be shown across the six predefined LEVEL-1 descriptors, but unadjusted effect sizes remain primary.
- LEVEL 2 and LEVEL 3 are spatial-concordance analyses.
- All failures, unresolved structures, unmapped materials, and non-evaluable experiments are counted and reported.

## Immutability

Any post-scoring change to eligibility, chemistry mapping, structure-resolution thresholds, tested descriptors, or primary metrics requires a new explicitly labeled sensitivity analysis or a new benchmark version. It cannot replace this frozen primary analysis.


## Deterministic multi-channel patch pooling addendum (frozen before scoring)

For LEVEL 2 and LEVEL 3 surfaces mapped to more than one canonical channel, spatial concordance uses a deterministic pooled candidate list:

1. Collect every patch from every pre-specified channel.
2. Deduplicate patches only when both the patch center key and the complete member-key set are identical.
3. Sort the pooled list by: Pareto front ascending, patch coherence descending, chemistry support descending, canonical public channel key lexical ascending, center key lexical ascending.
4. Assign pooled ranks 1..N from that ordering.
5. No channel is chosen by experimental overlap, and no best-channel result replaces the pooled primary result.
6. Composite residue propensity is the maximum V1 propensity for each residue across the pre-specified channels, consistent with the public V1 derived-overlay policy.
