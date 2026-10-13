# InterfaceScout 285-protein large structural panel — freeze protocol

Freeze date: 2026-10-13

## Purpose

This panel expands the frozen InterfaceScout V1 claim:

> InterfaceScout prioritizes protein surface regions that are structurally and chemically compatible with engagement of a specified material surface.

The 285-protein panel is an **applicability / structural landscape panel**, not an adsorption-amount, affinity, kinetics, or corona-membership benchmark.

## Source cohort and deterministic reproduction

Source publication: Li et al., Scientific Reports (2026), DOI 10.1038/s41598-026-55258-6.

The publication reports:
- 597 human NP-protein samples after human/feature curation;
- 1,513 human protein columns;
- 515 samples in the adsorption-classification cohort after removal of missing required features;
- a 5% non-zero detection threshold;
- 285 retained proteins.

The source panel is reproduced only from the pre-existing frozen PC-DB files under:
`validation/benchmark_v1_2026-10-06/raw_sources/chou_pcdb/`

Required classification metadata fields are treated as:
- Modification Type
- Modification Charge
- NP Type
- NP Sub-Type
- Size (nm)
- Zeta Potential (mV)
- Zeta Potential Charge
- Incubation time Group
- Shaking or Agitation During Incubation
- Extraction Centrifugation Time (min)
- Number of Washing Steps

Rows containing blank, `N/A`, or `Not Reported` in these fields are excluded. This reproduces exactly 515 source conditions. The use of `Incubation time Group` is required to reproduce the publication's 515-condition cohort because three records retain a reported ordinal incubation-time group while the raw minute field is marked `Not Reported`.

A protein is included when its source abundance is >0 in at least 5% of these 515 conditions, i.e. in at least 26 conditions. This reproduces exactly 285 proteins.

## Anti-cherry-picking rule

Membership in the 285-protein panel is frozen before any InterfaceScout score, patch, propensity, or material-channel output is inspected. No protein may be removed because InterfaceScout performs poorly or because its patch landscape appears inconvenient.

## Structure audit

Every source protein remains in `source_285_master.csv` even if it later proves unsuitable for primary structure-based analysis.

For each protein, audit:
1. protein identity / UniProt mapping;
2. experimental structure availability;
3. construct and sequence coverage;
4. chain choice and oligomeric context;
5. intrinsically disordered status;
6. whether only a domain/fragment structure is available;
7. whether a biologically meaningful experimental structure can support full-protein V1 analysis.

Primary structural stratum:
- experimentally determined structure;
- defensible identity and sequence mapping;
- sufficiently representative construct;
- folded protein state compatible with V1.

Separate sensitivity strata:
- domain-only experimental structures;
- partial experimental coverage;
- predicted structures (if used at all).

Do not silently pool experimental and predicted structures.

## InterfaceScout batch rule

The source detection label is **not** used as a success/failure target.

For structurally eligible proteins, run frozen V1 without fitting new weights or thresholds. Generate the full standard InterfaceScout outputs and preserve per-protein/per-channel outputs.

The 11 frozen canonical chemistry channels remain:
- anionic
- cationic
- hydrophobic
- pi/aromatic
- H-bond donor
- H-bond acceptor
- metal-oxide/oxide-like
- calcium/phosphate charged sites
- transition-metal coordination
- soft-metal sulfur affinity
- phosphate-rich

All-channel analysis is an applicability/landscape analysis. Real PC-DB material conditions may be annotated afterward for context, but they are not used to retrospectively select a favorable channel.

## Relationship to Level-3 validation

The 285-protein panel does not replace the experimental Level-3 gold standard.

- 285-protein panel: breadth, applicability domain, structure coverage, and patch-landscape behavior.
- Level-3 panel: accuracy against experimentally localized protein-side material-engagement regions using direct-overlap, 5 Å, 8 Å, top-k, Pareto-front, and residue-propensity metrics.

## Reporting

Always report the full source denominator:
- 285 source proteins
- number identity-resolved
- number with suitable experimental structures
- number domain-only/partial
- number predicted-structure only
- number IDP/ineligible
- number unresolved

This prevents easy-structure selection from being mistaken for the source panel itself.
