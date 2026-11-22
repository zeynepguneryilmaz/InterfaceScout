# InterfaceScout: final strict-quality nonredundant development structure panel
**Selection date:** 2026-11-22  
**Status:** 48 structure inputs accepted by the stated *structural screening, CATH nonredundancy, residue-level coordinate and software-execution* gates. This does **not** certify the physical predictive accuracy of InterfaceScout's chemistry maps or validate future parameter optimization.

## Reproducible evidence
- Frozen final PDB manifest (all experimental quality metrics, CATH H codes, biological assembly 1 SHA-256): [strict_apo_48_selected_candidates_v2_final.csv](strict_apo_48_selected_candidates_v2_final.csv).
- Initial broad structural selection: [GitHub Actions 45000271643](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45000271643) and [selection algorithm](rebuild_strict_apo_48.py).
- Independent replacement of incomplete 1H99 by 2PCN: [GitHub Actions 45001115049](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45001115049).
- Final 48/48 residue-level and 11-map computation: [GitHub Actions 45001342582](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45001342582).
- Final QA executable code: [strict_apo_residue_48.py](strict_apo_residue_48.py), [strict_apo_compat_48.py](strict_apo_compat_48.py).

## Dataset selection process
1. Initial large panel: **11,477 PDB candidates screened**; 65 passed structural-quality and chemical-context gates; 48 were chosen by sequence and CATH fold nonredundancy.
2. The independent residue-to-sequence quality audit of these 48 detected internal unresolved residues in **1H99**, so it was rejected rather than modeled or repaired.
3. An additional **401 previously unscreened candidates** were examined. **2PCN** passed all experimental, structural, chemical-context, sequence-alignment, CATH H/T and internal-gap checks.
4. The 48-entry manifest was frozen as **v2_final**; a new complete verification run returned **48/48 sequence/coordinate passes** and **48/48 successful 11-map InterfaceScout runs**.

Total tested across the two independent searches: **11,878 candidate records** (11,477 + 401). The 401 second-round records were selected from a previously untested segment of the same deterministic CATH-based candidate enumeration.

## Prespecified acceptance criteria (study-defined, not wwPDB official cutoffs)
| Category | Scientific screen |
|---|---|
| Experimental structure | X-ray diffraction, resolution <=2.50 Å |
| Refinement | wwPDB-reported DCC or deposited R-free <=0.280 |
| wwPDB geometry | clashscore <=15, Ramachandran outliers <=1%, rotamer outliers <=5% |
| Density | wwPDB RSRZ outlier percentage <=5% **when available** |
| Coordinate coverage | >=98% of declared protein sequence C-alpha atoms in **every selected biological assembly chain** |
| Side-chain completeness | 0 missing standard heavy atoms in modeled standard residues |
| Altloc treatment | 0 alternate-conformer residues, eliminating need for occupancy-based side-chain choices |
| Biochemical context | no biological-assembly DNA/RNA, metal ions, or nontrivial ligand residues (>=6 heavy atoms), as screened from assembly-1 mmCIF |
| Analysis unit | biologically annotated assembly 1; at most four protein chains and 1,500 modeled residues; protein entities 70–650 aa |
| Homology / fold nonredundancy | no shared CATH homologous superfamily **H** or topology **T** among selected structures |
| Sequence redundancy | exclude local alignments with >=30% identity and >=70% of shorter sequence covered, based on >=40 aligned residues |
| Missing internal sequence | 0 internally unresolved sequence positions in the final residue-label alignment |
| File provenance | source assembly CIF SHA-256 from initial selection must match the independently refetched file |
| InterfaceScout software | full biological assembly translated to protein-only PDB with unique chain IDs; all 11 maps run successfully; no duplicate residue keys |

No missing side chains or backbone residues were invented or geometrically reconstructed. Software format conversion and removal of solvents/hydrogens/heteroatoms are standard deterministic input preparation, not experimental coordinate repair.

## Final executed verification
- **48 of 48**: unique PDB IDs and distinct study structures.
- **48 of 48**: SHA-256 of downloaded biological assembly 1 matches the frozen source.
- **48 of 48**: residue-number/sequence correspondence; no internal residue gaps; no side-chain heavy-atom omissions or altloc residues; correct chain mapping.
- **48 of 48**: actual InterfaceScout computation succeeded, each with **11 of 11 physicochemical maps** (**528 map instances total**).
- **55 distinct CATH H superfamily codes and 55 distinct CATH T topology codes**, with no domain code repeated across different selected entries.
- Protein input sizes span from 70 to 678 modeled residues, with 27 monomeric, 18 dimeric and 3 trimeric prepared assemblies.

## Remaining scientific limitations / transparent caveats
1. **Six structures lack wwPDB RSRZ percentages** in the queried archived validation reports; unavailable does not mean zero or bad fit. They passed the other explicitly available metrics. The affected original IDs are **1A1X, 1KAF, 2FIP, 1K4Z, 1XTY, 1BM8**.
2. The no-ligand criterion excludes clearly modeled metals, DNA/RNA and larger bound ligands from the selected assembly, but small non-metal crystallization additives (e.g., sulfate) are possible. Their coordinates are omitted from protein-only InterfaceScout inputs. This is not equivalent to proving native apo-state dynamics.
3. Some high-quality entries represent experimentally determined **domains, mutants or constructs**, not necessarily complete unmodified proteins. Their biological identities must be described accurately in any manuscript.
4. CATH family/fold independence and computational correctness are **not** experimental proof that chemistry-map intensities predict adsorption free energy, unique orientations or eventual bound sites.
5. The current software parameter settings were historically selected on a different pilot panel. The new 48 are a quality-controlled development cohort; results of prospective parameter re-selection or stability analysis must be **computed and reported separately**, with any independent validation/holdout set kept untouched.

## Final 48 PDB IDs
1. 1A1X
2. 3FLD
3. 3IPJ
4. 3MX7
5. 4CU5
6. 4KT6
7. 6NDQ
8. 4BJQ
9. 7BOL
10. 3EOI
11. 3FFV
12. 7TVC
13. 3M7G
14. 3UAF
15. 3VJZ
16. 5AY6
17. 1KAF
18. 1DZF
19. 2FIP
20. 1O83
21. 1IU8
22. 2FB5
23. 1PXV
24. 6INX
25. 3M66
26. 1BGF
27. 6K7E
28. 1PBV
29. 6NZS
30. 4LE4
31. 2PCN
32. 6IVD
33. 1ZHQ
34. 7E0M
35. 3UCI
36. 6Q42
37. 3PIL
38. 4KMC
39. 6DGK
40. 6K39
41. 1K4Z
42. 2DC4
43. 1XTY
44. 1BM8
45. 4ULW
46. 1OX3
47. 7X4A
48. 6K0L

**Controlled source:** RCSB PDB, deposited biological assembly 1 mmCIF; wwPDB validation XML; current release CATH domain-level H and T assignments. Study-specific threshold choices are disclosed explicitly rather than presented as official wwPDB pass/fail boundaries.
