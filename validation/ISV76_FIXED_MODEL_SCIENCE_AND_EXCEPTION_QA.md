# ISV76 — fixed TV48 parameters: first full result interpretation and exceptions

**Date:** 2026-11-24. **Scope:** 76 active candidate structures, `ISV_active_candidates_76_v3.csv` (four retired entries excluded). These are **76 structural candidates, not 76 experimentally verified adsorption sites**. All map thresholds, local radii, Pareto objectives, pH 7.4, 150 mM ionic strength, 298 K, APBS nominal 0.8 Å and 42 SAS directions are inherited from TV48 unchanged.

## Actual production results

| Quantity | Observed |
| --- | ---: |
| Structures with 11 canonical maps | 76/76 |
| Canonical map instances | 836 |
| Distinct candidate patches over all maps | 64,592 |
| Standard residue profiles | 51,931 |
| tRSA-exposed residues | 36,273 |
| Real production SAS/APBS successes | 75/76 |
| Exposed residues with reported SAS potential | 34,853 |
| Entries with at least one backbone-incomplete residue | 3 |
| Entries with >10% residues missing one or more expected sidechain atoms | 25 |

**Failure 1GCT:** Original RCSB PDB ID is *gamma*-chymotrypsin; the active manifest's `Alpha chymotrypsin` label is a historical misannotation. In original assembly-1 protein-only coordinates, SER A:11 and SER E:11 each have only the N atom (missing CA/C/O). PDB2PQR appropriately rejects source input. An *explicit* deletion-only diagnostic removing precisely those two N-only records (no fabricated coordinates) allows 7,024-atom PQR / APBS and yields 287 SAS potentials, but the canonical map outputs are **not bitwise identical** to original. This edited structure is a separate diagnostic, **not** a silent replacement of the archived ISV76 result. High-resolution bovine **1YPH** (alpha-chymotrypsin) independently runs 11 maps + real APBS (478 residues; 310 SAS potentials) but is a different protease form and must not be substituted without an explicit scientific inclusion decision.

**8C50 curli-like fibril:** Original LPBE returned extremely negative values (median -380.689 kT/e; minimum -3331.666 kT/e across 683 exposed residues). All three chains have approx -73 e PQR charge, and B-chain median is -1519 kT/e. Real coarser and finer LPBE grids reproduce the extreme magnitudes; increasing ionic strength to 300 mM makes LPBE *far more extreme* (median -1917.186 kT/e, minimum -61377.95 kT/e). Real **nonlinear PB** with the *same* PQR and 0.8 Å nominal grid instead gives median -3.083 kT/e, minimum -7.866 kT/e, and similar chain medians approx -3 kT/e. This demonstrates why the original LPBE output **cannot be trusted for ordinary physical interpretation**. NPBE is a useful diagnostic, **not** separately validated adsorption evidence; low-resolution modeling-assisted fibril coordinates and chemical environment further limit interpretability. Quarantine original 8C50 LPBE magnitude values in pooled summary tables.

## Size and class structure

- Across 76 structures, total patch count and tRSA-exposed residue count have Spearman rho ~0.99459. **Do not interpret raw patch count as model success or use it to rank proteins.**
- Categories: 10 biofilm-related, 19 model enzyme, 22 model protein, 25 serum/corona/adhesion candidate structures. Median patch density per 100 exposed residues is ~172.76, 176.63, 177.05, 175.18, respectively, with different size and chemical contexts. These are descriptive density statistics, not experimental utility comparisons.
- All 55 channel pairs' Top-5 patch-member unions were compared in all 76 proteins: the **cationic and hydrogen-bond-donor** Top-5 member union is **identical in all 76**, anionic/acceptor median Jaccard ~0.8555, oxide/metal-coordinate median ~0.6460. These channels are not 11 statistically independent information sources. High overlap is not by itself a correctness failure.
- 3 structures have any incomplete backbone (1TEN, 1GCT, 1GFL). Many incomplete sidechains/cofactor-stripped assemblies require localized suitability checks, especially for material chemical interactions. 1STP retains the biotin-derived ligand-removed state caveat; 1HC1 is an unusually large hexameric assembly.

## Output/provenance

- Production compute: [ISV76 11-map/SAS-APBS GitHub run](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45385665397); combined status **INCOMPLETE** due explicitly retained 1GCT PDB2PQR failure despite 76/76 map success. Every original per-protein output and failure is archived.
- **Full browser atlas COMPLETE / PASS:** [real 3Dmol 76-structure workflow](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45399987636) completed all 76 case checks and strict final aggregate with **0 failures**. The 76 individual ZIPs contain **836 real map screenshots**, **836 map PDF files**, **3,344 B-factor/top-k PDB outputs**, numerical Excel/CSV/JSON, with ~3.297 GB total ZIP bytes. **75 atlases contain actual production APBS**, and **1GCT is explicitly a map-only atlas**, not a fabricated APBS success. GitHub final artifact: `ISV76-Complete-3D-PDF-Atlas-All-76` (artifact ID `14924292223`).
- Focused structure/physical diagnostic: [1GCT/8C50 real-solvers workflow](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45399840880) and subsequent corrected 1GCT diagnostic [run](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45400431654).

**Decisions:** No canonical map/score recalibration, no PDB substitution, no fabricated missing residues; 8C50 original LPBE excluded from interpretive pooling; original 1GCT APBS remains unavailable, with diagnostic variants separately labeled. Experimental residue-contact evidence acceptance and exact construct/source matching remain separate, prospective scientific validation tasks.

## References for exceptional entries
- RCSB `1GCT`: *Is gamma-chymotrypsin a tetrapeptide acyl-enzyme adduct of gamma-chymotrypsin?* https://www.rcsb.org/structure/1GCT
- RCSB `1YPH`: *High resolution structure of bovine alpha-chymotrypsin*. https://www.rcsb.org/structure/1YPH
- RCSB `8C50`: *Pontibacter korlensis curli subunit CsgA*. https://www.rcsb.org/structure/8C50
