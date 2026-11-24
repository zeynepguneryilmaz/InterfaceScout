# Validation and reproducibility

## Active scientific validation: two datasets

The current InterfaceScout v1.1.0 scientific work is organized around **only two active datasets**. See **[Active datasets and acceptance rules](ACTIVE_DATASETS.md)**.

1. **TV48 — Technical verification:** the frozen panel of **48 accepted protein structures**; evaluates implementation and consistency of all 11 maps, tRSA/scRSA handling, patch construction, numerical stability and Pareto descriptor behavior. [Manifest](strict_apo_48_selected_candidates_v2_final.csv) · [Final numerical verification](FINAL_V110_RESULTS.md) · [Extended technical audit](EXTENDED_V110_SCIENTIFIC_RESULTS.md).
2. **ISV76 — Initial scientific validation candidates:** a fixed active list of **76 distinct PDB structures**, selected from an 80-entry historical structural screen after four exclusions. The active [ISV76 manifest](ISV_active_candidates_76_v3.csv) establishes *structural panel membership*, not 76 accepted independent experimental systems. All 76 have completed 11-map computational analysis, whereas experimental region/construct mapping and source-level acceptance remain pending. See the [selection decision](ISV_SELECTION_DECISION_76.md).

### Historical analyses (archived, not active datasets)

The earlier eight-protein parameter-development study, nine-condition retrospective literature analysis, and original 80-structure screen are preserved as methodological provenance and historical work. They are **not current validation datasets**, and the previous nine conditions are **not** the ISV panel.

---

# Archived methodology and previous exploratory audits

This directory preserves the frozen evidence used to justify the InterfaceScout publication settings and the literature-localization benchmark.

## Parameter selection

`parameter_selection.json` is the final adsorption-label-free development-panel analysis. The eight development structures are 1CRN, 1R69, 1SHG, 2PPN, 4AKE, 1OMP, 1TIM, and 5CSC.

- **200 points/atom** was the lowest tested Shrake-Rupley sampling level that passed the predefined convergence criteria against the 1000-points/atom reference.
- **6/9 Å** ranked first among ten candidate multiscale radius pairs in the label-free cross-candidate consensus analysis. The former 5/8 Å pair ranked second and is not the publication model.
- **8 Å coarse patch radius** was not fitted to adsorption data; it is a predefined geometric scale used to turn a local chemistry maximum into a coarse candidate contact region.

## Literature benchmark

`benchmark_manifest.json` records literature-derived protein-side localization evidence and provenance. `benchmark_case_results.csv` contains the frozen case-level comparison metrics.

The primary set contains seven conditions across six proteins for which a native/folded-state structural comparison was considered interpretable. Two additional cases (transferrin-polystyrene and catalase-polystyrene) are retained as **secondary challenge cases** because adsorption-associated structural changes make a static native-state comparison more demanding.

Run:

```bash
python validation/verify_publication.py
```

to recompute the manuscript-level aggregate checks from the frozen case-level results and verify the selected development parameters.

This is a spatial-concordance validation of **candidate protein-side contact regions**, not of atomistic docking, adsorption free energy, adsorption amount, unique bound orientation, or adsorption-induced unfolding.

## Regional Chemistry Profile (descriptive; no trained combined score)

The retained Regional Chemistry Profile exposes **three** different explanations of an existing combined-map patch: each selected channel's mean canonical local support, its count/fraction of supporting surface residues, and the minimum *distinct-residue* Cα separation between different channels. These data do not rerank original single-map Pareto patches and do not predict contacts or affinity. All PDF, CSV, Excel, JSON, and complete ZIP outputs use the same backend definitions.

The [full equations, scope and 48-protein information audit](REGIONAL_CHEMISTRY_EVIDENCE_AND_METHOD.md) explain the decision to keep these descriptors but omit regional contrast (affine duplicate) and the minimum-over-channels combined scalar (repackaged data that could look like a fitted rank).

Run the normal test suite to verify invariance and end-to-end exports. The frozen label-free audit is available in the [successful GitHub Actions job](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45054941106). It is not an experimental interface-site performance validation.
