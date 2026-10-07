# Protein-level direct adsorption propensity protocol — frozen before BAD scoring

## Purpose

This analysis tests whether InterfaceScout V1 can support a protein-level claim that is narrower than absolute adsorption prediction:

> Given a folded protein structure, a material surface chemistry, and solution pH, InterfaceScout estimates the **relative material-specific structural compatibility for adsorption**.

The output is not an absolute adsorption amount, binding free energy, Kd, kinetic rate, or probability.

No InterfaceScout parameters, residue rules, chemistry definitions, structure thresholds, or benchmark eligibility rules are refit to BAD 2.0 adsorption values.

## Primary protein-level descriptor: Global Surface Compatibility (GSC)

For protein surface residue i and mapped canonical material chemistry channel c, InterfaceScout V1 already defines

L_ic = scRSA_i × state_availability_i,c

when residue i is favorable for chemistry c, and zero otherwise.

For a material mapped to one or more canonical chemistry channels, define

GSC = sum_i max_c(L_ic) / sum_i scRSA_i

over all standard amino-acid residues with scRSA >= 0.05 in the prepared protein surface.

Properties:
- range: 0 to 1;
- normalized for exposed protein surface size;
- pH-sensitive through the existing frozen V1 ionization-state availability;
- no fitted coefficients;
- no double-counting when a residue is compatible with multiple mapped channels;
- independent of the experimental adsorption amount.

Interpretation:

**GSC is the fraction of the accessible protein side-chain surface that is chemically compatible with the specified material class under the selected pH, weighted by side-chain accessibility and ionization-state availability.**

This is the primary descriptor. It will not be replaced by a secondary descriptor after observing performance.

## Predefined secondary descriptors

These are diagnostic only and do not replace GSC as the primary test:

1. Compatible residue fraction:
   number of exposed residues with max_c(L_ic) > 0 / number of exposed residues.
2. Compatible scRSA sum:
   sum_i max_c(L_ic), unnormalized; used to diagnose protein-size effects.
3. Strongest coherent patch support:
   maximum front-1 chemistry_support across mapped channels.
4. Strongest coherent patch persistence:
   maximum front-1 patch_coherence across mapped channels.
5. Compatible patch density:
   total number of front-1 patches across mapped channels / number of exposed residues.
6. Maximum V1 residue propensity across mapped channels.

No weighted combination of these descriptors is fitted.

## BAD 2.0 source

Primary direct-adsorption validation uses the live Biomolecular Adsorption Database table at https://bionanoinfo.com/bad/ linked to the peer-reviewed BAD 2.0 publication (DOI 10.1021/acsami.4c06759).

The live table is snapshotted before scoring. Source adsorption values are not used to create the material map or the GSC formula.

## Frozen BAD eligibility

A BAD row is primary-score eligible only if:
- protein has a source PDB ID;
- PDB can be downloaded and prepared by frozen InterfaceScout V1;
- surface name matches a predeclared canonical mapping rule below;
- surface concentration is numeric and > 0;
- solution concentration is numeric and > 0;
- pH is numeric;
- experimental row is not explicitly marked uncertain/invalid in the source validation/notes fields.

Rows failing these criteria remain in the audit and are not silently discarded.

Repeated measurements are not treated as independent proteins. Statistical analysis is performed at experimentally comparable strata as described below.

## Frozen BAD surface mapping rules

Surface-name matching is case-insensitive and based on explicit material/functional-group identity, not adsorption outcome.

- gold / Au, when explicitly bare/unmodified -> soft-metal sulfur affinity
- silica / silicon oxide / glass / quartz -> metal-oxide
- titanium oxide / titania / TiO2 -> metal-oxide
- iron oxide / magnetite / Fe3O4 -> metal-oxide
- alumina / aluminum oxide / aluminium oxide -> metal-oxide
- hydroxyapatite / calcium phosphate -> calcium/phosphate charged sites
- polystyrene / PS -> hydrophobic + pi/aromatic
- poly(tetrafluoroethylene) / PTFE -> hydrophobic
- graphite / graphene / carbon nanotube -> hydrophobic + pi/aromatic
- diamond-like carbon / DLC -> hydrophobic + pi/aromatic
- methyl / CH3-terminated explicit surface -> hydrophobic
- carboxyl / COOH-terminated explicit surface -> anionic + H-bond acceptor
- amine / NH2-terminated explicit surface -> cationic + H-bond donor
- hydroxyl / OH-terminated explicit surface -> H-bond donor + H-bond acceptor
- phosphate / phosphonate-terminated explicit surface -> phosphate-rich + anionic + H-bond acceptor
- PEG / polyethylene glycol outer surface -> H-bond acceptor
- mica -> anionic + H-bond acceptor

Ambiguous composites, proprietary surfaces, biological tissues, or surfaces whose exposed chemistry cannot be identified from the source name remain unmapped.

For explicit functionalized surfaces, the exposed functional group takes precedence over the underlying substrate.

## Experimental target

The primary experimental response is the source-reported BAD surface concentration.

Because adsorption measurements span different proteins, concentrations, methods, materials, and studies, raw surface concentration is not naively pooled across all rows.

Primary association analysis:
1. Create experimental strata with identical adsorbing-surface label, solution concentration, pH, ionic strength, temperature, measurement method, and experiment type where at least 3 different scoreable proteins are present.
2. Within each stratum, compute Spearman correlation between GSC and surface concentration.
3. Summarize stratum-level rho values by median, IQR, fraction > 0, and Wilcoxon signed-rank test versus 0 when defined.

Secondary analysis:
- within each protein × surface combination with >=3 solution concentrations, test whether GSC is constant as expected and do **not** claim it predicts concentration-dependent isotherm shape;
- compare GSC with adsorption normalized by solution concentration only as a sensitivity analysis, not as the primary physical quantity;
- report surface-class-specific results when enough strata exist.

## Ranking interpretation

A successful result supports only:

> Proteins with higher InterfaceScout global surface compatibility tend to show higher adsorption onto the same material under comparable experimental conditions.

It does not support:
- absolute adsorption amount prediction;
- adsorption free energy;
- adsorption kinetics;
- extrapolation to unmapped chemistry;
- conformational-change prediction.

## Corona datasets

PC-DB and Payne are retained only as secondary complex-mixture generalization tests. They do not define the direct-adsorption claim.

## Immutability

GSC, secondary descriptors, BAD mapping rules, eligibility criteria, and primary statistical analysis above are frozen before the BAD adsorption values are scored. Any later alternative formula is a separately labeled sensitivity analysis and cannot replace this primary test.


## Pre-scoring source-inventory addendum

Before any GSC values or adsorption-association results were calculated, inspection of the BAD 2.0 source inventory showed that poly(tetrafluoroethylene) (PTFE) is a recurrent explicit surface label. PTFE is therefore mapped to the hydrophobic canonical channel based solely on its exposed fluorocarbon chemistry. This addendum was committed before the first BAD/GSC scoring run. No surface rule will be added or removed in the primary analysis after performance results are viewed.
