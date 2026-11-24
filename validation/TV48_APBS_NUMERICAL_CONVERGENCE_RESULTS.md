# TV48 numerical convergence — real APBS / SAS sampling (2026-11-24)

## Scope and audit provenance
- Frozen SHA-256 checked TV48 structures; real PDB2PQR/PROPKA + APBS
  LPBE, pH 7.4, 150 mM, 298 K.
- Standard production target is 0.8 Å, not a guaranteed exact spacing.
- 48/48 proteins successfully sampled with 42, 84 and 162 deterministic
  Fibonacci directions on *the same* real APBS potential field.
- Four structures preselected from real geometric grid scans (1OX3,
  2FB5, 6INX, 6K7E), the only TV48 entries that truly exceed 257 grid
  points on some axis at target 0.4 Å.
- Each of these four proteins received four *independent* APBS solves:
  0.8 Å, 0.5 Å, 0.4 Å capped at 257, and genuinely uncapped 0.4 Å.
  All four grid comparisons use a fixed 162-direction SAS shell.
- All 48 sampling shards and all 4 uncapped-grid shards completed PASS.
  The initial aggregate job failed solely because 'gemmi' had not been
  installed in its reporting environment. The exact raw artifacts were
  downloaded and aggregated independently; CI aggregate dependency was
  corrected in a subsequent commit.

## Surface sampling: protein-level medians (kT/e)

| Directions | Median of protein median absolute difference | Median of protein residue-p95 absolute difference | Median Spearman rho | Median top-decile | p95 ≤0.5 kT/e |
|---|---:|---:|---:|---:|---:|
| 42 vs 84 | 0.06108 | 0.39901 | 0.98862 | Jaccard 0.826 | 37/48 |
| 84 vs 162 | 0.03788 | 0.24932 | 0.99523 | Jaccard 0.854 | 44/48 |
| 42 vs 162 | 0.05913 | 0.38053 | 0.99012 | Jaccard 0.833 | 39/48 |

The medians are modest, but some local errors are substantially larger:
worst protein-specific p95 is 1.55 kT/e for 42 vs 162.
**42 directions are suitable for broad surface trends but not provably
converged residue-by-residue.** At 84 directions the local tail improves,
although 4/48 proteins still exceed 0.5 kT/e at their p95.

## Real grid convergence: four cap-sensitive PDB structures

| APBS comparison (fixed 162 directions) | Median of four protein median differences (kT/e) | Largest protein residue p95 difference (kT/e) | Minimum rho |
|---|---:|---:|---:|
| target 0.8 vs 0.5 Å | 0.02210 | 0.60999 | 0.99948 |
| target 0.5 vs uncapped 0.4 Å | 0.00919 | 0.17589 | 0.99983 |
| capped vs uncapped 0.4 Å | 0.00406 | 0.06181 | 0.99993 |

Grid cap examples: 1OX3 expands from 161×161×257 to 161×161×321;
2FB5 from 225×257×257 to 225×321×289. The other two uncapped
structures similarly have >257 nodes along their long axis. This is real
numerical refinement confirmed from the OpenDX grid, not only a change
in input target.

### Methodological decision

- Do not change the canonical 11 maps/Pareto.
- No evidence here that production 0.8 Å requires an immediate change:
  0.8 Å does not hit the 257-axis cap anywhere in TV48.
- The current 42-direction default is computationally functional but its
  finite sampling contributes non-negligible local variability. If
  surface-level residue ranking or hotspot membership matters, prefer
  84 directions and consider 162 as a reference; explicitly distinguish
  performance costs before changing application defaults.
- The cap effect at 0.4 Å was small **in these four examples**, not a
  universal property of arbitrarily large biomolecules.
- Fibonacci grids are *not nested* and pairwise errors need not decrease
  monotonically. No experimental adsorption endpoint was tested.

## Reproduction

- Script: [`tv48_numerical_convergence.py`](tv48_numerical_convergence.py)
- CI: [`tv48-numerical-convergence.yml`](../.github/workflows/tv48-numerical-convergence.yml)
- Underlying real input manifest:
  [`strict_apo_48_selected_candidates_v2_final.csv`](strict_apo_48_selected_candidates_v2_final.csv)
