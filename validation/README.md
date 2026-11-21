# Publication validation and audit data

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
