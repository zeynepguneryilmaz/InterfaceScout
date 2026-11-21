# Literature-localization benchmark

`results/benchmark_case_results.csv` is a **verbatim copy** of `validation/benchmark_case_results.csv` in the same repository revision.

It contains 7 primary and 2 secondary cases. The accompanying experiment-level provenance and experimental localization annotations remain in [the canonical manifest](../../validation/benchmark_manifest.json).

To verify the frozen parameter-selection and aggregate checks from the repository root:

```bash
python validation/verify_publication.py
```

The two secondary cases reflect the difficulty of applying a static native-state protein model to adsorption-induced structural changes. The listed values must not be interpreted as binding energies, adsorption rates or atomistic docking accuracy.

For other results, use the built-in application **Complete export** and archive the exact files with the associated input PDB and version information.
