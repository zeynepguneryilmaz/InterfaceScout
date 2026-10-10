# Prospective Level-3 Expansion Eligibility Protocol — frozen before scoring

Date frozen: 2026-10-10

## Purpose

This expansion prospectively tests only the frozen InterfaceScout V1 claim:

> InterfaceScout prioritizes protein surface regions that are structurally and chemically compatible with engagement of a specified material surface.

The prediction unit is a coarse protein-side surface region / interface patch.

This expansion does **not** test or claim absolute adsorption amount, adsorption free energy, Kd, kinetics, corona membership, unique atomistic orientation, or adsorption-induced unfolding.

The original `benchmark_v1_2026-10-06` remains immutable. New examples are stored and reported separately.

## Primary Level-3 eligibility

A condition is eligible for primary Level-3 scoring only when all of the following are satisfied before InterfaceScout outputs are inspected:

1. **Experimental protein-side localization**
   The study provides direct experimental evidence that can be mapped to the protein sequence as one or more of:
   - exact residue anchor;
   - localized residue set;
   - localized short residue region;
   - experimentally identified peptide region;
   - experimentally inferred orientation/contact anchor.

2. **Material-facing interpretation**
   The evidence must identify a protein region that is surface-facing, surface-proximal, preferentially engaged, protected by surface contact, or otherwise experimentally interpretable as an adsorption/contact region.
   A region that is only newly exposed, destabilized, or unfolded after adsorption is Level 2, not Level 3.

3. **Structurally interpretable state**
   The adsorbed protein must remain sufficiently folded/native-like for the experimental localization to be meaningfully mapped onto the input structure, or the study must explicitly isolate an early/contact state before major rearrangement.
   Severe unfolding, precipitation, aggregation, or irreversible denaturation excludes the condition from primary Level 3 unless a pre-rearrangement contact region was independently measured.

4. **Structure mapping**
   An experimental PDB matching the protein/construct is preferred.
   Sequence numbering and construct identity must be verified before scoring.
   AlphaFold structures, if ever required, are retained as a separate predicted-structure stratum and are never silently pooled with experimental PDBs.

5. **Material mapping**
   The exposed material chemistry must map prospectively to one or more frozen InterfaceScout canonical maps.
   Ambiguous composites, buried core materials, proprietary coatings, or chemically unresolved outer surfaces are deferred rather than forced into a map.
   For explicitly functionalized surfaces, the exposed outer chemistry takes precedence over the substrate.

6. **Independent experimental evidence**
   Docking, molecular dynamics, computational orientation models, and residue-energy decomposition are not admissible as gold-standard labels.
   They may be used only to understand a paper or choose a structure after the experimental label has independently been established.

## Exclusions / separate strata

- Intrinsically disordered proteins: exclude from primary folded-protein Level 3; may form a separate IDP benchmark.
- Adsorption-induced newly exposed/deformed regions without contact localization: Level 2 only.
- Engineered covalent immobilization that imposes orientation: separate engineered-immobilization stratum, not primary spontaneous adsorption.
- Strong aggregation/precipitation/unfolding: exclude primary Level 3 unless an independent native-like contact state is resolved.
- Model-derived contact sites: exclude as ground truth.
- Chemistry-unmapped materials: verified-deferred, not scored.

## Multiple conditions and pseudoreplication

Different pH values, particle sizes, or surface chemistries may be retained as separate experimental conditions when the source independently resolves different interfaces.

However, headline performance must also be clustered by **study × protein × material-family** so that multiple conditions from one experimental family do not inflate apparent sample size.

Confirmatory studies that repeat an already represented protein/material family are retained and reported, but are not treated as wholly independent new biological systems.

## Frozen chemistry mappings for this expansion

Existing canonical rules remain unchanged:
- silica / silicon oxide / glass / quartz -> metal-oxide
- citrate-coated Au outer surface -> anionic + H-bond acceptor
- bare Au -> soft-metal sulfur affinity
- MPA-Au -> anionic + H-bond acceptor
- AlumOH positive condition -> cationic + metal-oxide
- PAA-coated Fe3O4 -> anionic + H-bond acceptor
- pristine polystyrene -> hydrophobic + pi/aromatic

Prospectively added before any new scoring:
- **pristine C60 fullerene -> hydrophobic + pi/aromatic**
  Rationale: pristine fullerene is an exposed aromatic carbon cage with strong hydrophobic/aromatic interaction character. C60 molecular-nanomaterial examples are reported as a distinct sensitivity stratum because C60 is not an extended solid surface.

No chemistry rule will be added or removed after observing new expansion scoring results.

## pH rule

Source experimental pH is used whenever available.
If pH cannot be recovered, the condition is **deferred from primary scoring** rather than silently assigned pH 7.4 in this prospective expansion.

## Evidence-resolution strata

Results must remain stratified by:
- residue_anchor
- localized_residue_set
- localized_residue_region
- peptide_region
- broader_region

Broader regions are not treated as equivalent to residue-scale evidence.

## Frozen scoring metrics

The same Level-3 spatial metrics used for the original frozen benchmark are retained without modification:
- rank of first direct overlap;
- rank of first predicted patch within 5 Å of ground truth;
- rank of first predicted patch within 8 Å;
- top-3 near-8-Å hit and recall;
- top-5 near-8-Å hit and recall;
- top-5 direct-overlap hit and recall;
- best Pareto-front-1 near-8-Å recall;
- residue-propensity mean/max for the experimental region.

For multi-channel surfaces, the deterministic ordering remains:
1. Pareto front ascending;
2. patch coherence descending;
3. chemistry support descending;
4. map key;
5. center key.

No weights, post-result channel selection, or benchmark-fitted thresholds are allowed.

## Prospective workflow

1. Literature candidate discovered.
2. Candidate entered into `candidate_registry.csv` with status ADMIT / DEFER / EXCLUDE before scoring.
3. Exact experimental residue/region label, pH, PDB/construct, and material mapping are verified.
4. Only ADMIT conditions are copied into `prospective_level3_ground_truth.csv`.
5. This freeze branch is sealed.
6. A separate scoring branch is created and InterfaceScout V1 is run without changing V1.
7. Original 7-condition benchmark and prospective expansion are reported separately and then, where appropriate, as a clearly labeled combined sensitivity analysis.
