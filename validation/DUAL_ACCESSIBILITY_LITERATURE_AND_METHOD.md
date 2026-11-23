# Dual-accessibility surface model: scientific basis and validation policy

## Motivation
Glycine has no non-backbone heavy-atom side chain. Treating its C-alpha atom
as a chemical side chain under the original scRSA descriptor is a category
error. Removing glycine entirely, however, discards genuine solvent-exposed
backbone geometry.

## Definitions
- **tASA**: sum of solvent-accessible areas of all heavy atoms in each
  standard residue (calculated in the complete biological assembly).
- **tRSA**: tASA / residue-type-specific maximum theoretical ASA from Tien
  et al. (2013), Table 1. This is used as a **geometric surface gate**.
- **scASA**: SASA of non-backbone heavy atoms (excluding N, CA, C, O, OXT).
- **scRSA**: scASA / the pre-existing InterfaceScout side-chain reference
  areas, used for chemistry support for the other 19 residues.
- **Glycine**: scASA = 0 and scRSA = 0; glycine is retained in the geometry
  whenever its tRSA passes the surface gate. This does not suppress nearby
  residue chemistry or remove glycines from CA neighborhood calculations.

## Literature
Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO.
Maximum Allowed Solvent Accessibilites of Residues in Proteins.
PLOS ONE 2013;8(11):e80635.
DOI: https://doi.org/10.1371/journal.pone.0080635

The Tien 2013 theoretical maximum values are established reference
denominators for total-residue RSA (Biopython refers to this as Wilke).
They are *not* experimentally calibrated cutoffs for InterfaceScout.

Shrake A, Rupley JA. Environment and exposure to solvent of protein atoms.
J Mol Biol. 1973;79(2):351-371.
DOI: https://doi.org/10.1016/0022-2836(73)90011-9

## Important methodological caveats
1. This dual-accessibility definition is now the **v1.1.0 production model** (`2.5.0-dual-accessibility`); earlier comparison scripts remain historical developmental evidence.
2. The production threshold **tRSA >= 0.05** is a fixed, sensitivity-tested modeling convention, not an experimentally optimized adsorption-site cutoff.
3. Tien maximum ASA and the InterfaceScout Shrake-Rupley radii/discretization
   are not identical methodology; absolute normalized values are therefore
   convention-dependent.
4. A glycine residue may be a geometric patch member but must not receive
   fictitious side-chain chemistry.
5. Technical sensitivity studies are complete for exploratory model variants; for **final v1.1.0** validation use the frozen 48-structure production re-run and the separate paired literature evaluation. The nine historical literature conditions are not an independent holdout.
6. Mixed tRSA geometry and scRSA chemistry are distinct descriptors and
   must never be labeled as interchangeable SASA measures.

## Reproduction
For final production validation run `.github/workflows/final-v110-48.yml` (48 SHA-matched assemblies, all 11 maps and Pareto ablation) and `.github/workflows/final-v110-literature.yml` (nine evidence conditions, fixed-candidate ablations). The historical exploratory workflow `.github/workflows/dual-accessibility-48.yml` is not the final-production result set.
