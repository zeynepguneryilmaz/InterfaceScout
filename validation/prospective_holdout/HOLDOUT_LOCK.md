# Prospective holdout lock

Lock date: 2026-09-05
Frozen model: `IS-v1-Core-primary`

## Frozen primary rule

The prospective primary InterfaceScout score is the **Core IS** score only:

1. material-chemistry compatibility;
2. side-chain solvent accessibility (scRSA);
3. pH-dependent residue state/compatibility;
4. multiscale 5/8 A spatial patch persistence.

GNM mobility, APBS electrostatic potential/compatibility, radial prominence and RIN descriptors are **auxiliary/diagnostic outputs only** and are not universal multiplicative factors in the IS-v1 primary score.

The choice was made using the curated development benchmark (12 positive conditions / 7 unique proteins) before any prospective holdout outcomes were calculated. The development benchmark itself is not blind or prospective.

## Holdout rules

1. Register PDB/structure, material chemistry, source conditions and ground-truth provenance before calculating IS performance whenever practicable.
2. Do not change IS-v1-Core after viewing holdout performance.
3. Report all registered holdout systems, including failures.
4. Repeated conditions on one protein are clustered at the protein level.
5. Do not use GNM/APBS/radial/RIN holdout behavior to alter the primary Core formula. They may be reported descriptively as pre-existing auxiliary descriptors.
6. If any primary rule, chemistry ontology, cutoff, exposure definition, pH/state rule or aggregation rule is changed after holdout inspection, increment the model version and establish a new untouched holdout set.
7. A holdout is only considered prospective if its outcome was not used in the development benchmark or feature-selection process.

At the time of this clarified lock there are zero scored prospective holdout outcomes. Candidate systems may be registered separately before scoring.
