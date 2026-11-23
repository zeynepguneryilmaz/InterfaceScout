# InterfaceScout v1.1.0 — production validation and manuscript evidence matrix
> **Dataset scope (current plan):** The active datasets are [TV48 technical verification and the planned ISV initial scientific-validation panel](ACTIVE_DATASETS.md). The 48-structure tests in this document are TV48 results. Any nine-condition retrospective literature analyses reported below are **historical exploratory analyses**, **not ISV** and not independent holdouts.


**Scope:** v1.1.0 / model `2.5.0-dual-accessibility`. This file defines the interpretation of final-production evidence and must be read together with the numerical workflow artifacts. An action that computes but has not successfully aggregated is **not** a completed data package.

## A. Production technical correctness — 48 nonredundant structures

Code: [final_v110_48.py](final_v110_48.py); workflow:
[Final v1.1.0 48-structure all-map verification](../.github/workflows/final-v110-48.yml).

Acceptance criteria:

1. All 48 structure inputs correspond to the SHA-256 hashes in the frozen biological-assembly manifest; no silent file substitutions.
2. The **production** `analyze_all_maps` API returns 11 maps per input (**528 map instances**), rather than only a forked prototype.
3. Every unique chain/residue/insertion-code identifier maps to one residue; all exposure statuses are defined by the unrounded `tRSA_raw >= 0.05` condition.
4. GLY is eligible for geometric surface membership while its `sidechain_sasa` and `scrsa` remain zero.
5. Every coarse patch consists of eligible surface residues. Every patch's `mean_accessibility` is a mean of geometric tRSA, not side-chain scRSA.
6. Independent direct implementation of Pareto nondominance over the **same four recorded objectives** reproduces front 1.
7. An omission test quantifies the change in Pareto-front membership after removing each objective (4 x 528 = **2112 ablation comparisons**). A changing score is evidence of algorithmic nonredundancy, not independent biological meaning.

Output files: `all_528_maps.csv`, `all_residue_tRSA_scRSA.csv`,
`all_pareto_ablations.csv`, SHA-256 provenance, and `SUMMARY.json`. The workflow
must show **PASS in all four compute shards and the aggregation job**.

## B. Matched literature evidence and controlled ranking baselines

Code: [final_v110_literature.py](final_v110_literature.py); workflow:
[Final v1.1.0 literature and fixed-candidate baselines](../.github/workflows/final-v110-literature.yml).

- The **nine** historically curated conditions represent **eight unique protein structures**; seven are primary (including two pH conditions for carbonic anhydrase II), two are secondary challenge cases. Two measurements of one protein do **not** count as independent protein replicates.
- Use identical prepared PDB coordinates, pH, reference residue labels, and all 11 chemical channels for each scoring method.
- Compare production Pareto/display ranking against **chemistry-support-only** and **tRSA-accessibility-only** rankings. All three methods rank the **same candidate-patch pool**. These are **ranking ablations**, NOT independently generated candidate-region baselines.
- Readouts: Top-5 presence within 8 Å of evidence-labeled residues, fraction of labels near candidate centers, exact patch-member overlap, and fraction of labels in candidate patches. An 8 Å C-alpha criterion is a coarse regional measure, not atomic contact.
- The reported evidence varies: residue anchors, local NMR region sets, crosslinked peptide regions, and broad regions from proteolysis/NMR. Do **not** pool these as equally precise atomic contacts.
- Report per-case, per-channel results. All-map means are only **descriptive** because channels are chemically correlated, and broad evidence spans make overlap easier.
- Compute paired differences on the same condition, with protein-cluster bootstrap intervals rather than treating all 11 map scores as independent replicates.
- Do **not** pick the best-performing map after inspecting experimental evidence and then claim unbiased model performance.
- This manifest was used in earlier InterfaceScout work: it is **not a new holdout**. Consequently no statistically robust claim of general experimental adsorption-site prediction is justified by this retrospective analysis alone.

## C. Minimum reference-provenance matrix

Evidence source (identified by original DOI):

| Condition | Evidence resolution | Source | Interpretation |
| --- | --- | --- | --- |
| BSA–silica | broad folded-domain region | [10.1016/j.jcis.2005.03.064](https://doi.org/10.1016/j.jcis.2005.03.064) | Domain II is a coarse region; not atomistically determined contacts |
| HSA–PAA/Fe3O4 | cross-linked peptide ranges | [10.1021/ac201889j](https://doi.org/10.1021/ac201889j) | Proximal peptide evidence; verify exact peptide numbering against original paper |
| CA II–silica, pH 6.3 | site-directed fluorescence anchor | [10.1529/biophysj.104.054809](https://doi.org/10.1529/biophysj.104.054809) | Neighborhood around position 10, not a unique atomic contact |
| CA II–silica, pH 8.3 | site-directed fluorescence anchor | [10.1529/biophysj.104.054809](https://doi.org/10.1529/biophysj.104.054809) | Region around position 37; dependent on same protein/study |
| Ubiquitin–citrate AuNP | NMR-localized region | [10.1021/nl101746v](https://doi.org/10.1021/nl101746v) | Q2–I3 and L15–E18 region evidence |
| β2-microglobulin–MPA AuNP | NMR chemical-shift region | [10.1039/C8TB01129J](https://doi.org/10.1039/C8TB01129J) | Coating-specific localized region; exact residue set requires figure/table verification |
| Asparaginase–aluminum hydroxide | solid-state NMR regional perturbation | [10.1038/s41541-019-0115-7](https://doi.org/10.1038/s41541-019-0115-7) | Structural perturbation, not necessarily direct surface contact at every labeled residue |
| Transferrin–polystyrene | proteolytic peptide ranges | [10.1021/acs.analchem.9b00482](https://doi.org/10.1021/acs.analchem.9b00482) | Secondary challenge: accessibility can reflect post-adsorption conformational change |
| Catalase–polystyrene | proteolytic peptide ranges | [10.1021/acs.analchem.9b00482](https://doi.org/10.1021/acs.analchem.9b00482) | Secondary challenge: same limitation; equivalent chains need careful matching |

The bibliographic existence and broad experimental topics have been cross-checked, but **individual residue ranges and cross-study residue/PDB numbering are not all independently source-audited**. Before final manuscript claims, audit every mapped range against the original article figures and supplements and preserve the original supporting quote/figure number in a separate evidence table. This is distinct from software execution PASS.

## D. Manuscript claim policy

| Claim | Required support |
| --- | --- |
| Numerical SASA implementation and residue accounting | Analytical geometry checks, atom aggregation checks, frozen-input provenance |
| Structurally consistent 11-map computations | Production 48 x 11 PASS, unique residue keys, no silent input transformations |
| Candidate regions contain interpretable, nonredundant computational information | Objective ablation, transparent maps, appropriate baseline comparisons, regional-descriptor information audit |
| Recovery of experimental protein-side regions | Source-matched conditions, explicit label resolution, transparent top-k metrics and cautious exploratory scope |
| General prediction of material adsorption site, binding orientation, energy, or kinetics | **Not an InterfaceScout claim and not demonstrated** |

**Release gate:** current HEAD cross-platform CI and both final-production validation
workflows must conclude PASS. Tests from previous commits do not certify current HEAD.
Publication figures and numerical tables must be regenerated from the accepted run SHA,
not from historical pre-dual-accessibility output files.
