# Independent wwPDB QA follow-up — important correction (2026-11-22)

1PDB failed structural geometry screening despite its favorable resolution and R-free:
- Clashscore = **72.38**, Ramachandran outliers = **4.89%**, rotamer outliers = **35.12%**, RSRZ outliers = 1.08%.
- Decision: **REJECT 1PDB** from primary protein-surface parameter-development panel.

Other previously included candidates, from the independently generated wwPDB XML reports:

| PDB | Clashscore | Ramachandran outliers | Rotamer outliers | RSRZ outliers | Decision |
|---|---:|---:|---:|---:|---|
| 1C5E | 1.67 | 0.00% | 1.88% | 0.35% | Retain |
| 1A1X | 6.35 | 0.00% | 4.21% | unavailable | Retain subject to assembly/source limitations |
| 1STN | 3.60 | 0.00% | 1.74% | unavailable | Retain subject to 91.3% Cα completeness |
| 2CBA | 1.70 | 0.00% | 4.04% | 0.39% | Retain subject to Zn/altloc handling |
| 1RWZ | 9.32 | 0.00% | 0.95% | **16.39%** | Hold for density-fit comparison against alternate PCNA |

Quality screen: https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44976725496

**4M6J** is a high-resolution (1.201 Å, R-free 0.1845) monomeric DHFR-NADPH structure with 184/187 residues and zero missing standard sidechain heavy atoms (22 alternate-conformer residues); selected as a **provisional replacement candidate only**, pending independent wwPDB geometry and end-to-end 11-map smoke. Additional candidate 7FQ7 (1.1 Å, R-free 0.1589, 158/159 aa, no missing sidechain, 29 altloc residues) also under examination.

**1GE8** is an alternative ring PCNA with 2.1 Å, R-free 0.2529, 238/249 aa per chain, zero missing standard sidechain atoms and no altlocs; independent wwPDB report shows 9.2% RSRZ outliers and 0.4% Ramachandran/3.4% sidechain outliers, so not automatically preferable to 1RWZ. PCNA choices must be compared on the predeclared same metrics.

The v4 CSV therefore holds **48 candidate identifiers, not 48 final approved structures**. It deliberately flags 4M6J (pending) and 1RWZ (density review). The last automated independent wwPDB comparison run was queued at https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/44976808346 at the time this report was written. Do not assert completed validation for pending items.

The other 42 proteins in the panel still require homogeneous final audit; successful software execution confirms output generation, not chemistry-map physical correctness.
