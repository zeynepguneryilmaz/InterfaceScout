# InterfaceScout v5.1 canonical freeze — 2026-09-04

This branch is an immutable reference point for the pre-v5.2 canonical model used during development and diagnostic benchmarking.

## Frozen core

- `L_i,c = I_i,c * scRSA_i * f_state,i,c(pH)`
- `P_i,c = 100 * L_i,c / max(L_c,fav)`
- Cα-based patch topology
- patch radii = 5 Å and 8 Å
- `M_i,c = 100 * min(D5_norm, D8_norm)`
- side-chain Shrake–Rupley SASA
- probe radius = 1.40 Å
- points per atom = 200
- surface threshold `scRSA >= 0.05`
- standard reference pKa values
- APBS potential, Ebase values and historical 8 Å context remain auxiliary only and are excluded from the canonical score

## Interpretation freeze

The canonical model predicts protein-side residue/patch compatibility with generalized material chemistries. It does not claim absolute adsorption free energy, adsorption capacity, a unique adsorption orientation, adsorption-induced conformational change, explicit hydration/desolvation, or multi-protein corona organization.

## Development/diagnostic systems already seen

The following systems must **not** be treated as untouched held-out validation systems because they were used during model development, geometry diagnostics, or steering diagnostics:

- 4F5S BSA
- 5H7A Protein A
- 1MBN myoglobin
- 2PTN trypsin
- 2HHB hemoglobin

Any final external validation must use systems selected after this freeze and must not alter the core formulation based on their outcomes.
