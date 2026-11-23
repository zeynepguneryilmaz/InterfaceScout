# InterfaceScout v1.1.0 — final-production validation results
> **Dataset scope (current plan):** The active datasets are [TV48 technical verification and the planned ISV initial scientific-validation panel](ACTIVE_DATASETS.md). The 48-structure tests in this document are TV48 results. Any nine-condition retrospective literature analyses reported below are **historical exploratory analyses**, **not ISV** and not independent holdouts.


**Software:** 1.1.0 · **Model:** `2.5.0-dual-accessibility`  
**Analysis date:** 23 November 2026  
**Status:** Both final-production numerical validation workflows **PASS**. This does not by itself certify every cross-platform release job or independent biological accuracy.

## 1. Frozen 48-protein production verification

**Exact successful run:** [Final v1.1.0 48-structure all-map verification, run 45244648958](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45244648958)

Inputs are derived from 48 SHA-matched biological assembly 1 structures in `strict_apo_48_selected_candidates_v2_final.csv`. All calculations use the **public production `analyze_all_maps` API** and the 11 real maps.

| Test | Verified outcome |
| --- | ---: |
| Protein structures | **48/48 PASS** |
| Public chemistry maps | **528/528 PASS** |
| Residue records | **12,068** |
| Glycine residues | **838** |
| Geometrically exposed glycines (tRSA gate) | **637** |
| Total ranked candidate patches across all 528 map instances | **15,910** |
| Pareto-front-1 patches across all instances | **4,226** |
| Four-objective Pareto omission comparisons | **2,112** |

The workflow checks unique residue keys, exact frozen assembly hashes, correct tRSA surface selection, zero artificial GLY side-chain chemistry, patch geometric membership, and independently recomputed nondominated Pareto fronts.

### Objective-omission audit

Mean Jaccard similarity of full four-objective Pareto-front-1 centers to a three-objective front after removing one objective (48 proteins × 11 maps):

| Removed objective | Mean Jaccard |
| --- | ---: |
| Chemistry support | 0.8054 |
| Mean tRSA accessibility | 0.7650 |
| Patch coherence | 0.7241 |
| Orientation coherence | 0.7147 |

These values demonstrate that excluding each objective **changes the computed prioritized region sets**. They do **not** prove that each criterion is an independent biological mechanism or that the four-objective method predicts adsorption more accurately.

## 2. Retrospective evidence concordance and constrained baselines

**Exact successful run:** [Final v1.1.0 literature and fixed-candidate baselines, run 45244597881](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45244597881)

**9 condition records, 8 distinct structures, 11 channels, 3 ranking strategies = 297 outputs.** Seven primary condition records correspond to **six proteins** because carbonic anhydrase II was considered at two pH values. The two secondary challenge cases were evaluated separately.

All ranking methods use the **identical production candidate-patch set**. Thus chemistry-only and tRSA-only are **prioritization ablations** and not stand-alone alternative site-generation algorithms.

### Primary set (seven conditions)

The following are descriptive **means across the 11 channels and seven conditions**, not independent biological replicates.

| Experimental-region readout, top five candidate patches | Pareto | Chemistry support only | Geometric tRSA only |
| --- | ---: | ---: | ---: |
| Any region within 8 Å of experimental labels | 51.95% | **59.74%** | 48.05% |
| Fraction of experimentally labeled residues within 8 Å | **27.72%** | 23.69% | 25.86% |
| Any direct patch-member overlap | 48.05% | **55.84%** | 48.05% |
| Fraction of experimentally labeled residues among patch members | **21.88%** | 17.28% | 21.02% |

**Interpretation:** Pareto produces higher average *residue recovery* in these retrospective conditions, but chemistry-only ranking more often identifies **at least one** evidence-associated region. Pareto is **not uniformly superior**. Small apparent differences must not be presented as statistical significance.

### Protein-cluster bootstrap uncertainty

For the six primary proteins, the protein-cluster bootstrap 95% interval for the **Pareto minus comparator** mean 8 Å label-recovery fraction is:

- vs **chemistry only**: mean +0.0167; 95% interval **[−0.0403, +0.0879]**;
- vs **tRSA only**: mean +0.0294; 95% interval **[−0.0147, +0.0890]**.

**Both intervals include zero.** No convincing superiority of Pareto can be claimed from this sample. Protein-cluster uncertainty is appropriate because the 11 channels are correlated and two experimental conditions share a protein.

Secondary evidence consists of only **two conditions** (transferrin and catalase on polystyrene). This subgroup is too small to infer generality, and some proteolysis changes could reflect adsorption-induced conformational rearrangement rather than direct native-state contact.

## 3. What the manuscript can and cannot claim

**Supported technical claims**
- The final v1.1.0 tRSA/scRSA model processes all 48 selected high-quality structures and 11 maps per structure.
- tRSA surface geometry and side-chain-scRSA chemistry are implemented consistently; GLY is retained geometrically without a fictitious side chain.
- Four Pareto descriptors can produce different front-1 candidate sets, as directly measured by ablation.
- Selected known literature regions can be compared with computed candidate patches at a coarse, explicitly defined spatial resolution.

**Not supported as claims**
- Pareto ranks are experimentally validated material-affinity probabilities.
- The four-objective Pareto method outperforms every simpler comparator in all metrics.
- A broad proteolytic or fluorescence-supported region is an exact atomically proven adsorbed contact.
- The nine-case retrospective set is a new independent holdout.
- Adsorption free energy, pose, kinetic mechanisms, or a unique future adsorption site are predicted.

## 4. Reproducibility and remaining gates

- [Executable 48-structure production verification](final_v110_48.py) and [workflow](../.github/workflows/final-v110-48.yml).
- [Executable literature/baseline comparison](final_v110_literature.py) and [workflow](../.github/workflows/final-v110-literature.yml).
- [Protocol, evidence-resolution hierarchy, source DOIs and acceptance gates](FINAL_V110_VALIDATION_PROTOCOL.md).
- [Source/target structural-quality screening and frozen hash manifest](STRICT_48_FINAL_ACCEPTANCE_REPORT.md).

**Before release/manuscript submission:** verify a **successful latest-head cross-platform CI** (all Windows/macOS/Linux jobs), retain all evidence artifacts at a permanent release tag or archive, and complete exact source-figure/sequence-to-PDB numbering checks for every benchmark evidence region. Numerical validation success is not a substitute for checking original experimental evidence.

**Publishing rule:** Results and figures in the paper must use this v1.1.0 production-run provenance. Historical pre-v1.1.0 maps, thresholds and rankings are not numerically interchangeable with the final results.
