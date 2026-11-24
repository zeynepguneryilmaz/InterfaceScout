# InterfaceScout — established-methods comparison: tasks, capabilities, constraints, resources

Date: 2026-11-24. Preliminary source-backed functional comparison; external FreeSASA run is tracked separately and **not declared PASS before a complete action**.

## Comparison design

No globally comparable adsorption-site prediction benchmark exists for these tools without validated material/condition-specific region labels. Distinguish **(i)** direct SASA numeric agreement, **(ii)** descriptive protein surface features, **(iii)** geometric pocket detection, and **(iv)** the operational effort of using a tool. Absence of documentation ≠ proven absence of functionality.

| Tool | Actual published purpose | Directly comparable domain | Distinct strengths / limitations for this study |
|---|---|---|---|
| InterfaceScout (our work) | 11 descriptive protein-side physicochemical maps, tRSA/scRSA, 8 Å local regions, Pareto, overlays and multi-structure comparison | Protein side surface chemistry; descriptive candidate regions | One integrated offline workflow with real optional PDB2PQR–PROPKA–APBS; no actual material adsorption affinity, kinetics or final contact prediction |
| SURFMAP | Automated **2D surface projection** of electrostatics, Kyte–Doolittle / Wimley–White hydrophobicity, stickiness, topographic relief and any PDB B-factor descriptor | Strong conceptual overlap in surface feature visualization | Custom descriptors and 2D projection are capabilities InterfaceScout should **not** claim as its own novel invention; installs via CLI/Docker and uses dependencies including APBS for electrostatics |
| Protein-sol patches | Web-based charged and hydrophobic surface patches, finite-difference Poisson–Boltzmann electrostatics, NPP SASA | Direct conceptual comparator for charged/hydrophobic surface patching | pH-dependent web visualization and B-factor exports already available; external server processing time/RSS not measurable from current evidence |
| PEP-Patch | Quantitative protein-surface electrostatic patches and open-source surface outputs (NPZ/PLY) | Electrostatic patch analysis | Strong direct comparator of electrostatic patch semantics, but neither full 11-channel mapping nor material-specific adsorption prediction can be assumed |
| ProteinTools | Hydrophobic clusters and salt-bridge networks | Some structure chemistry and local grouping | Focus on hydrophobic contact networks/bridges rather than material-side chemistry patch Pareto; do not assert absent features without examination |
| FreeSASA | Independent C/Python SASA computation (Lee–Richards or Shrake–Rupley) | Solvent accessibility and exposure threshold agreement | Not a general candidate-region or binding-site method. Default atomic radius classification differs from InterfaceScout/Biopython |
| CASTp 3.0 | Surface topographic pockets, cavities, channels, volumes | Geometry only | Pocket region is NOT equivalent to a physicochemical material-interface patch |
| PrankWeb/P2Rank | Ligand binding-site prediction | Different task | Small molecule site localization is not solid-material adsorption; a common experimental success metric would be misleading |

**Implementation attribution:** InterfaceScout uses Biopython Shrake–Rupley for SASA and separately invokes real PDB2PQR/PROPKA and APBS for optional potentials. Do not claim those algorithms are novel components or that shared APBS potentials outperform APBS itself.

## Recorded InterfaceScout operational results

The original ISV76 fixed-model production run recorded **76 full per-protein processing durations** (11 maps + attempted real APBS):
- Median **9.60 seconds**, 75th percentile **21.435 seconds**, 90th percentile **54.795 seconds**.
- Maximum **388.17 seconds** for 2FHA (4,128 residues).
- Sum of 76 per-protein task elapsed times **2,367 seconds** = **39.45 minutes**, **NOT** the wall time of the 16-way parallel GitHub Actions job.
- Residue count vs task elapsed time Spearman **rho=0.95734**; strong size dependence.
- 76/76 canonical 11-map sets; 75/76 original APBS runs. An APBS-only, N-only orphan recovery was subsequently implemented and passed separate real 1GCT and 1HEL regression; the original recorded outcome is not overwritten.
- Initial CI log does **not** contain measured peak RAM; do not make memory claims from this run.
- No comparison to competitor time is legitimate unless inputs, enabled calculations, worker hardware, timing boundaries and I/O are comparable.

## External FreeSASA experiment

Reproducible script: [`isv76_external_freesasa_runtime_audit.py`](isv76_external_freesasa_runtime_audit.py); workflow [`isv76-external-freesasa-audit.yml`](../.github/workflows/isv76-external-freesasa-audit.yml).

Paired 76-structure FreeSASA 2.x and InterfaceScout runs on the same Linux runner, same PDB/structure hashes; FreeSASA Shrake–Rupley 200 points per atom, 1.4 Å probe and 1 thread. Case-level results include residue SASA Spearman, median/p95 SASA differences in Å², tRSA=0.05 exposed-state agreement and sampled peak-RSS/wall-clock for the two separate programs. The engines have different *atomic radii* and their timed workloads differ (FreeSASA SASA only; InterfaceScout 11 maps + patches without APBS); the recorded times are transparent operational diagnostics, **not a comparable throughput ratio**. Failures remain FAIL.

**Packaging challenge discovered in actual GitHub trial:** Ubuntu apt-provided FreeSASA 2.1.2 rejects the `--output-depth=residue` argument advertised in a version of upstream documentation; the distribution exposes `--depth=residue`. The pinned Ubuntu-command issue is corrected in the workflow script. This first failure is preserved, not represented as numerical failure of either SASA method.

## What should be assessed but is not yet measured externally?

| Dimension | Planned recording | Constraint |
|---|---|---|
| Input friction | chain selection, prepared PDB, heteroatom support, residue identity validation, preparatory work | Must use live tool, cannot infer user's clicks from a publication |
| Installation | reproducible steps, dependent binaries/containers, supported OS, installation failures | Server-based software moves work to service; no local memory benchmarking |
| Runtime | actual task wall clock with scope, hardware, startup and export time separated | No invented server CPU time from browser load duration |
| Peak RAM | per-process or full process tree RSS; polling interval specified | Estimated sampled lower-bound not exact peak memory |
| Output inspection | local JSON/CSV, residue-level PDB, viewer, figures, annotations, PDF/ZIP | Each tool's output has different semantic unit |
| Limits | >2000-residue assemblies, fibrils, missing atoms, cofactor removal, multi-chain and model support | A failure in one test case is a scope limit, not blanket tool invalidity |
| Repeatability | fixed versions, exact PDB SHA256, determinism and provenance | Browser servers can update silently; record access date and version if displayed |

## Evidence and sources

- SURFMAP (Schweke et al., *J. Chem. Inf. Model.* 2022), doi:10.1021/acs.jcim.1c01269; https://github.com/i2bc/SURFMAP
- Protein-sol (Hebditch et al., *Sci. Rep.* 2019): https://pmc.ncbi.nlm.nih.gov/articles/PMC6374528/ ; https://www.protein-sol.manchester.ac.uk/patches
- PEP-Patch: https://pmc.ncbi.nlm.nih.gov/articles/PMC10685443/
- ProteinTools: https://pmc.ncbi.nlm.nih.gov/articles/PMC8262690/
- FreeSASA: https://freesasa.github.io/doxygen/CLI.html and https://manpages.debian.org/trixie/freesasa/freesasa.1.en.html
- CASTp 3.0: https://pmc.ncbi.nlm.nih.gov/articles/PMC6031066/

**Publication interpretation:** technical interoperability, efficient candidate-region inspection, transparent map combinations and export may be defensible benefits; **superior predictive accuracy or overall speed has not been demonstrated**.
