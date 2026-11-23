# ISV80 — Surface-relevant missing atoms, biological assemblies and cofactor-state audit

**Version:** InterfaceScout 1.1.0 / model `2.5.0-dual-accessibility`. **Successful run:** [45276540364](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45276540364) — all eight shards and combined evidence **PASS**. Source code: [isv80_surface_assembly_audit.py](isv80_surface_assembly_audit.py).

This audit distinguishes (i) a residue with missing heavy side-chain atoms anywhere in the structure, (ii) a missing-atom residue in the **geometrically exposed surface**, and (iii) one belonging to the union of **Top-5 mapped candidate regions across all 11 channels**. It additionally inventories **all omitted nonwater components** and compares the analyzed biological-assembly-1 selection with the list of available deposited RCSB assemblies.

**Methodological caveat:** This is a *localization/quality audit*, **not** an empirical quantification of how reconstructed side chains or retained cofactors would change map values. No missing side chain was rebuilt, no chemical scoring was altered, and no experimentally measured adsorption geometry was supplied.

## 1. Full 80-structure outcome

| Verification | Measured result |
| --- | ---: |
| Proteins successfully examined | **80/80** |
| Proteins with at least one exposed incomplete side chain | **80** |
| Proteins with an incomplete side chain intersecting one of the eleven channels' Top-5 regions | **80** |
| Total exposed incomplete-sidechain residue occurrences across structures | **2,925** |
| Total incomplete-sidechain residue occurrences in the union of Top-5 regions within each structure | **835** |
| Proteins with any omitted nonwater/nonstandard component | **60** |
| PDBs listing an alternative biological assembly in addition to assembly 1 | **18** |

Each residue is counted only once **per structure**, even when it belongs to multiple Top-5 map regions. Across distinct structures, values represent observations and not independent biological samples. This audit does not show that every absent atom changes a channel score: chemistry and local exposure can affect whether it matters.

## 2. Five previously flagged PDBs

| PDB | Total protein residues in selected assembly | Exposed residues | Exposed residues with missing sidechain atoms | Missing-sidechain residues in Top-5 region union | Unique Top-5 member residues |
| --- | ---: | ---: | ---: | ---: | ---: |
| **3GHG** | 1,929 | 1,480 | **89 (6.0%)** | **18** | 289 |
| **1QJS** | 408 | 292 | **14 (4.8%)** | **8** | 190 |
| **1STP** | 484 | 335 | **12 (3.6%)** | **1** | 154 |
| **4DN0** | 279 | 211 | **55 (26.1%)** | **26** | 121 |
| **6LQH** | 2,160 | 1,399 | **108 (7.7%)** | **19** | 280 |

No backbone heavy-atom omissions were detected in the prepared protein residues of these five assemblies under this audit. Their previously flagged C-alpha coverage differences are **deposited-sequence versus observed-chain scope differences**, not automatically absent backbone atoms within otherwise modeled standard residues.

### 3GHG — fibrinogen

- **Protein-only assembly 1:** ten prepared protein chains and 1,929 standard residues. RCSB enumerates assemblies **1 and 2**; use of assembly 1 must be reconciled against the biologically relevant fibrinogen unit and the selected study question before final acceptance.
- **Surface chemical uncertainty:** 89 exposed residues with missing sidechain heavy atoms; 18 occur in the Top-5 region union. Missingness may be localized or affect a mapped chemical seed; no corrected-map comparison has been run.
- **Removed components include:** NAG, MAN, GAL, SIA and BMA carbohydrate entities, plus Ca. Glycosylation and calcium chemistry are therefore not explicitly represented in the protein-only maps.
- **Preliminary status:** **CONDITIONAL — assembly/chain choice + glycan and missing-surface-chemistry audit required**. Do not interpret 3GHG as a homogeneous complete soluble monomer.

### 1QJS — hemopexin

- **Protein-only assembly 1:** one chain, 408 modeled standard residues; the RCSB entry lists assemblies **1 and 2**.
- **Surface chemical uncertainty:** 14 exposed missing-sidechain residues, eight among Top-5 region members.
- **Removed components:** HEM, Na, Cl and phosphate. **HEM is the bound heme cofactor** integral to the reported hemopexin state. In the protein-only representation, the modeled protein atom positions still originate from the heme-bound structure, but heme electrostatics, occlusion and contacts are absent.
- **Preliminary status:** **CONDITIONAL — heme-bound structural state and assembly selection must be disclosed**. For heme-dependent pocket claims, prefer explicit ligand-aware analysis or decline the claim.

### 1STP — streptavidin

- **Protein-only assembly 1:** biologically familiar **tetramer**, four protein chains and 484 standard residues.
- **Surface chemical uncertainty:** 12 exposed incomplete side chains; **one** occurs in the union of Top-5 candidate regions.
- **Removed component:** BTN (**biotin**); the starting protein coordinates are from the streptavidin–biotin complex. The ligand-removed model **cannot be called experimentally determined apo streptavidin**.
- **Preliminary status:** **CONDITIONALLY USABLE** for general protein-side surface mapping with the tetramer and ligand-removed state clearly labeled. For biotin-pocket conclusions, explicitly evaluate bound ligand and missing side chains.

### 4DN0 — PelD

- **Protein-only assembly 1:** one protein chain, **279 observed** standard residues; experimental construct is **PelD 156–455** rather than full-length membrane protein.
- **Surface chemical uncertainty:** 55/211 exposed residues (26.1%) and 26/121 unique Top-5 region members have missing sidechain heavy atoms. This is the **highest priority of these five** for alternative structure selection or local atom repair **with a strictly justified and recorded method**.
- **Removed components:** C2E (**c-di-GMP**), Na and sulfate. c-di-GMP binding and the truncation of PelD are core functional-context differences; they do not imply a direct extracellular material-binding site.
- **Preliminary status:** **HOLD for strict region-level interpretation**; keep as a named cytoplasmic-domain biofilm example only in the conditional exploratory stratum until the missing-sidechain and ligand/construct problems are addressed.

### 6LQH — CsgG–CsgF

- **Protein-only assembly 1:** **18 protein chains**, including the CsgG nonamer and associated CsgF chains, 2,160 modeled standard residues.
- **Surface chemical uncertainty:** 108 exposed incomplete side chains; 19 appear among Top-5 region members.
- **Context:** this is an **outer-membrane secretion complex**, with the surrounding bilayer absent from the model; exposed-to-solvent geometry alone does not distinguish membrane-buried surfaces from experimentally accessible aqueous-facing surfaces. The screening did not find an omitted small-molecule cofactor; this does not address membrane omission.
- **Preliminary status:** **SPECIAL STRUCTURAL CHALLENGE**, not an ordinary independent soluble protein surface sample. Retain only under an explicitly membrane-context-aware structural interpretation; exclude from strict soluble-globular performance pooling.

## 3. Scientific admission decisions

| PDB | Decision for ISV working panel | Decision for strict residue-level interpretations |
| --- | --- | --- |
| 3GHG | Keep conditionally | **Hold** until assembly, glycan and local surface completeness confirmed |
| 1QJS | Keep conditionally | **Hold for heme-sensitive pocket claims** |
| 1STP | Keep | **Conditionally admissible** away from biotin pocket after target-region completeness check |
| 4DN0 | Keep as a domain-level biofilm challenge | **Hold**; large missing exposed sidechain fraction and ligand truncation |
| 6LQH | Keep only as special membrane-complex challenge | **Exclude from soluble-globular pooled statistics**; membrane-context analysis required |

These decisions are **scope-specific**: a protein can be usable for visualization, software robustness and domain-level descriptive mapping yet not qualify for claims about a specific experimentally observed adsorption interface.

## 4. Remaining checks — explicitly not claimed as completed

1. **Exact score perturbation:** Rebuild missing atoms or choose a better experimentally resolved structure *with method provenance*, rerun chemistry maps and quantify changes in regional scores and Pareto ranks. Missingness location alone does not quantify score bias.
2. **Correct assembly determination:** Assembly 1 was explicitly selected and all available RCSB assembly IDs were inventoried, but the script does not establish by experiment that assembly 1 is biologically correct for **all 80 contexts**. In particular, **3GHG and 1QJS** have multiple deposited assembly choices.
3. **Cofactor-dependent chemistry:** The model omits ligands/metal/cofactors. It identifies their source and omission but does not measure electrostatic/steric impact. The interpretation of heme, biotin, c-di-GMP, glycans and bilayer environment requires a separate ligand-aware or spatial-proximity analysis.
4. **Validation scope:** Protein ID and structure verification is not experimental evidence for a protein–material adsorption-site prediction. The ISV working panel remains distinct from TV48.

**Archived output:** the successful run's `ISV80_SURFACE_QA.csv`, `ISV80_SURFACE_QA_DETAILS.json`, `SUMMARY.json`, five PDB-specific missing-atom details, per-shard issue tables and complete provenance are available under the original [GitHub Actions run](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45276540364) and permanently copied to `validation/evidence/ISV80_SURFACE_ASSEMBLY_COMPLETE.zip`.
