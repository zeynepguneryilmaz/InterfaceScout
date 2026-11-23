# InterfaceScout v1.1.0 — Extended scientific evaluation

**Reproducible source run:** [GitHub Actions 45252414860](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45252414860) — **PASS**, all four protein-panel shards, the experimental-case job, and final aggregation. Source commit `4952f75992a15b7dc32a1dd8f6fbad2160a1f065` · model `2.5.0-dual-accessibility`.

**Permanent result archive:** [Full verified data bundle](evidence/extended_v110_full_evidence.zip) · [Machine-readable summary](evidence/extended_v110_SUMMARY.json) · [Independent-baseline bootstrap](evidence/extended_v110_protein_bootstrap.csv) · [Case-study map combinations](evidence/extended_v110_case_studies.csv). The bundle was sourced from the successful run without recomputation.

**Script:** [extended_v110_science.py](extended_v110_science.py) · [CI workflow](../.github/workflows/extended-v110-science.yml). Download the **extended-v110-full-scientific-evidence** artifact from the successful run for the complete numerical data: all correlations, 21,000 stability trials, source checks, control results, bootstrap intervals, per-region example data, and script-generated summary.

**Scope:** These analyses test *computational descriptor dependence*, *numerical ranking robustness*, *independent candidate construction*, and *descriptive chemical information*. They do not measure experimentally calibrated adsorption probability, adsorption free energy, binding orientation, or a new independent prospective holdout.

## 1. Criterion dependence: 48 proteins, 528 maps

Of 528 map instances, **525 contained candidate patches**, **three produced no patches**. The empty maps were retained in dataset coverage counts, but omitted from correlation and rank-perturbation computations. Of 3,150 pairwise criterion-map combinations, **3,126** were evaluable; four nonempty maps had insufficient variation for a correlation statistic.

Spearman rank correlations were computed **within each map** over candidate patches and summarized as *medians of evaluable map-level coefficients*, not as independent protein replicates.

| Criterion pair | Median Spearman rho |
| --- | ---: |
| Chemistry support – mean tRSA | 0.061 |
| Chemistry support – patch coherence | **0.561** |
| Chemistry support – orientation coherence | 0.028 |
| Mean tRSA – patch coherence | 0.273 |
| Mean tRSA – orientation coherence | **0.622** |
| Patch coherence – orientation coherence | 0.140 |

**Inference:** The four objectives are **partially correlated, not statistically or physically independent**. The two strongest relationships are chemistry support with patch coherence, and mean accessibility with orientation coherence. Objective omission in the previous 2,112-case ablation changes the selected regions, but cannot establish independent underlying physical mechanisms. Appropriate manuscript wording: *four complementary, partially correlated prioritization descriptors*.

## 2. Ranking stability: numerical perturbation

For each of the 525 nonempty map instances, apply **20 repeat draws** of a deterministic-seed bounded independent multiplicative perturbation to the four recorded objective values at each of two nominal relative perturbation standard deviations (1% and 5%); recompute full nondominated Pareto fronts and display-ranked Top-5 regions. This produces **21,000 rank-perturbation evaluations**. A perturbation here is a **designed numerical stress test**, not a measured model-error distribution.

| Numerical perturbation | Mean front-1 Jaccard | Mean Top-5 center Jaccard | Fraction retaining original Top-1 center |
| --- | ---: | ---: | ---: |
| 1% relative SD | 0.715 | 0.743 | **0.913** |
| 5% relative SD | 0.606 | 0.639 | **0.778** |

**Inference:** The first ranked region is more robust than exact front membership, but no perfect robustness can be claimed. These scores characterize sensitivity of *scores and ordering* to synthetic noise, not protein coordinate uncertainty, conformational dynamics, or chemical model uncertainty. Further comparison under alternative PDB conformations would address a different robustness question.

## 3. Independently constructed reference regions

The earlier *chemistry-only* and *tRSA-only* ranking baselines all reused InterfaceScout candidate regions. This new analysis adds two **independently generated candidate-center methods**:

1. **SASA-only local peaks:** select exposed residues by descending unrounded tRSA with 8 Å same-face nonmaximum suppression, then generate their 8 Å coarse neighborhoods.
2. **Geometry-only farthest-point sampling (FPS):** select surface C-alpha points maximizing successive spatial coverage, with no residue-type, chemistry-map, or InterfaceScout patch input, then form 8 Å coarse neighborhoods.

Both reference procedures reuse **the same tRSA surface gate and 8 Å same-face neighborhood definition** for a controlled geometric region comparison, so they are **not completely independent end-to-end software pipelines**. They generate their own centers and do not read the 11 maps, patch maxima, or Pareto rankings. The baseline budget is **five regions per case**, while the InterfaceScout readout is the **mean of 11 separate five-region map outputs**. The latter is descriptive across channels and should not be interpreted as a single material-tailored five-region prediction. No specific material channel was selected after inspecting evidence.

### Retrospective primary evidence: seven conditions on six distinct proteins

| Readout | InterfaceScout, mean across 11 channels | SASA-only | Geometry FPS |
| --- | ---: | ---: | ---: |
| At least one experimental region within 8 Å of Top-5 centers | 51.95% | 42.86% | **57.14%** |
| Fraction of labeled residues within 8 Å of centers | **27.72%** | 19.03% | 14.20% |
| Any direct patch-member overlap | 48.05% | 42.86% | **57.14%** |
| Fraction of labeled residues covered by patch members | **21.88%** | 11.88% | 9.38% |

The four metrics differ: a broad annotated region is easier to *hit at least once*, while its **fractional recovery** penalizes incomplete coverage. These primary values average over historically curated nonheldout conditions and channels, not independent confirmations of prediction accuracy.

### Protein-cluster bootstrap uncertainty

All comparisons are paired by experimental condition; carbonic anhydrase pH variants are clustered within the same protein. Below, **positive** means higher descriptive average recovery for InterfaceScout. The 4,000-draw 95% percentile intervals include zero.

| Readout: InterfaceScout minus baseline | Mean delta | 95% protein-cluster interval |
| --- | ---: | ---: |
| Label recovery within 8 Å vs SASA-only | +0.0484 | [−0.0959, +0.1985] |
| Direct label recovery vs SASA-only | +0.0636 | [−0.0252, +0.1703] |
| Label recovery within 8 Å vs geometry FPS | +0.1048 | [−0.0346, +0.2563] |
| Direct label recovery vs geometry FPS | +0.0928 | [−0.0334, +0.2191] |

**Inference:** There is a descriptive trend toward larger fractions of recovered evidence regions under InterfaceScout's multi-channel map average, but **no statistically convincing advantage** can be claimed with this small retrospective set. The method budget, channel selection, and broad vs residue-localized evidence remain important confounders. Additional untouched proteins and independent matched-budget designs would be needed for an accuracy-superiority claim.

## 4. Two prespecified structural case studies

**Proteins:** human ubiquitin, `1UBQ`, pH 7.7; and human carbonic anhydrase II, `2CBA`, pH 6.3. The pairs were specified in the executable *before* the data were summarized: hydrophobic + aromatic/π, cationic + hydrogen-bond donor, and anionic + hydrogen-bond acceptor.

The resulting numbers describe the **whole exposed structure and its mapped candidate regions**, not exactly the experimentally contacting residues on a particular substrate.

| Protein | Two-map combination | Exposed residues | Exact same-residue overlap | Co-localized candidate regions | Regions supported by different residues only |
| --- | --- | ---: | ---: | ---: | ---: |
| Ubiquitin (1UBQ) | Hydrophobic + π | 61 | **2** | 17 | **11** |
| Ubiquitin (1UBQ) | Cationic + donor | 61 | 29 | 12 | 0 |
| Ubiquitin (1UBQ) | Anionic + acceptor | 61 | 21 | 13 | 0 |
| Carbonic anhydrase II (2CBA) | Hydrophobic + π | 177 | **8** | 55 | **34** |
| Carbonic anhydrase II (2CBA) | Cationic + donor | 177 | 73 | 40 | 0 |
| Carbonic anhydrase II (2CBA) | Anionic + acceptor | 177 | 55 | 38 | 0 |

**Demonstrated information distinction:** For ubiquitin hydrophobic + π, only **two residues** individually support both channels, but **17** candidate regions contain both chemistries, including **11** where different residues provide the support. For carbonic anhydrase II, **eight** residues show exact hydrophobic/π overlap while **34 of 55** co-localized regions are supported by distinct residues only. The high exact overlap of cationic/donor and anionic/acceptor reflects chemically related channel definitions and must **not** be described as multiple independent physicochemical confirmations.

Co-localized regions **overlap** and are not 17 or 55 independent surface sites. All spatial co-localization results reuse existing canonical 8 Å candidate regions; they do not add fitted binding energy or predict interface contacts. Detailed region members, local channel supports, supporting fractions and distinct-residue C-alpha distances are in the CI artifact.

## 5. Scientific decisions

- **Retain the four Pareto objectives**, describe them as *complementary but correlated* and report objective-omission sensitivity transparently.
- **Retain all 11 maps** without combining their existence into an independent-evidence count. Channel correlations are expected in some chemical categories.
- **Retain exact overlap and spatial co-localization as separate outputs**; the case studies establish an actual difference in what those outputs describe.
- **Do not claim statistical superiority** to independently generated SASA or geometric patch controls; all primary recovery confidence intervals presented here cross zero.
- **Do not present this exercise as a new independent experimental validation**. Evidence sources and PDB numbering require separate original-study audits.
- **Report current numerical stress testing distinctly** from future conformational or experimental perturbation testing.

## 6. Reproducibility and files

Workflow: [v1.1 independent baselines, criterion dependence and case studies](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45252414860) (**PASS**).

The run artifact contains: `all_criterion_correlations.csv`, `all_rank_stability_trials.csv`, `all_patch_counts.csv`, `independent_baselines.csv`, `production_per_channel.csv`, `independent_baseline_protein_bootstrap.csv`, `per_condition_independent_baseline_comparison.csv`, `case_study_map_combinations.csv`, `case_study_region_details.csv`, `case_provenance.json`, `SUMMARY.json`, plus per-shard source provenance.

The high-quality structural development panel is frozen in `strict_apo_48_selected_candidates_v2_final.csv` and the selected assembly SHA-256 hash is checked for each protein. User-facing scoring logic remains unchanged by this validation code. This report is **not** a replacement for independent experimental residue-annotation verification.
