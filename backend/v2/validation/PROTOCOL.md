# InterfaceScout publication validation protocol

This protocol keeps parameter selection separate from literature-derived experimental evaluation and makes the publication results reproducible from the repository.

## 1. Development panel

Numerical settings are selected on an adsorption-label-free structural panel. The final panel is:

- 1CRN, chain A — crambin
- 1R69, chain A — phage 434 repressor N-terminal domain
- 1SHG, chain A — alpha-spectrin SH3 domain
- 2PPN, chain A — native FKBP12
- 4AKE, chain A — unligated adenylate kinase
- 1OMP, chain A — unliganded maltodextrin-binding protein
- 1TIM, chains A,B — triosephosphate isomerase dimer
- 5CSC, chain B — open citrate synthase

No material identity, adsorption amount, experimental contact residue, or experimental orientation is used in this stage.

### Shrake–Rupley sampling

Sampling densities of 100, 200, 500, and 1000 points per atom are compared with 1000 points per atom as the numerical reference. The selected setting is the smallest tested value that satisfies all predefined convergence criteria:

- median exposed-residue Jaccard ≥ 0.98;
- worst-case exposed-residue Jaccard ≥ 0.95;
- median scRSA mean absolute error ≤ 0.01;
- median Top-10 hotspot Jaccard ≥ 0.90.

This procedure selects **200 points per atom**.

### Multiscale radii

Candidate radius pairs are 5/7, 5/8, 5/9, 5/10, 6/8, 6/9, 6/10, 7/9, 7/10, and 8/10 Å. Every candidate is compared with every other candidate across the development proteins and chemistry maps.

Pairs are ordered lexicographically by:

1. highest 10th-percentile Top-5 hotspot-set consensus;
2. highest median Top-5 consensus;
3. highest median Top-10 consensus;
4. lowest median outer-neighborhood fraction;
5. smaller outer and then inner radius as final deterministic tie-breakers.

This procedure selects **6/9 Å**. The separate **8 Å** radius used to construct a candidate patch is a fixed structural neighborhood scale, not a fitted outer aggregation radius.

## 2. Literature-derived benchmark

A benchmark case must provide experimental information about the protein-side region associated with a material interface. Evidence can be reported as:

- an individual residue or anchor;
- a localized residue set;
- a peptide or sequence segment;
- a broader protein region.

The comparison is interpreted at the spatial resolution supported by the original experiment. Broad regional evidence is not treated as an atomically exact contact map.

The final benchmark contains seven primary conditions across six proteins and two secondary challenge cases. The exact benchmark manifest is `benchmark_experimental_final.json`.

## 3. Prediction-first ordering

Only these fields enter prediction:

- case identifier;
- protein/PDB identifier;
- selected chain or chains;
- validation-only surface-to-chemistry mapping;
- pH.

Experimental localization residues or regions, DOI metadata, evidence descriptions, and comparison outcomes are not passed to the predictor.

The workflow is:

1. download and archive the exact PDB input;
2. prepare the selected coordinate model;
3. save the prepared PDB snapshot;
4. generate and save InterfaceScout predictions;
5. only then load the independently encoded experimental annotation;
6. calculate direct-overlap, 5 Å, and 8 Å spatial comparisons.

The prepared PDB snapshot saved during prediction is reused by the comparison stage so that prediction and scoring operate on exactly the same coordinates.

## 4. Comparison metrics

Each metric answers a different practical question.

- **Direct-overlap recall**: whether the same experimentally annotated residues occur in the candidate set.
- **Near-5 Å recall**: whether the candidate set reaches the same immediate local neighborhood.
- **Near-8 Å recall**: whether it reaches the same patch-scale protein surface region. The 8 Å criterion matches the structural scale used to construct an InterfaceScout candidate patch.
- **First-hit position**: the first displayed candidate with nonzero recovery under a stated criterion.
- **Cumulative Top-3 / Top-5 recall**: recovery after combining the residues from a small, reproducible shortlist of candidates.

Pareto-front assignment defines the multiobjective priority class. The deterministic within-front order provides a reproducible sequence for inspection.

## 5. Geometry reference

The matched geometry-only reference applies the same 8 Å same-face construction to solvent-exposed centers without chemistry filtering. It provides structural context for the spatial opportunity available on the exposed protein surface.

## 6. Frozen publication settings

- SASA points per atom: 200
- SASA probe: 1.40 Å
- scRSA exposure threshold: 0.05
- multiscale aggregation: 6 and 9 Å
- candidate-patch radius: 8 Å
- material-specific fitted weights: none

Any future change to a prediction-defining setting should be treated as a new model version and evaluated separately.

## 7. Reproduction

From the repository root:

```bash
cd backend
python -m v2.validation.reproduce_publication
```

The command regenerates the development and benchmark outputs under `publication_data/` and writes `publication_data/verification.json`. A successful publication reproduction reports zero verification mismatches.
