# Archived literature-localization benchmark

The results are maintained only once, at **[`validation/benchmark_case_results.csv`](../../validation/benchmark_case_results.csv)**. Experimental source provenance and regions are in **[`validation/benchmark_manifest.json`](../../validation/benchmark_manifest.json)**.

There are 7 primary and 2 secondary cases. These benchmark records are separate from the descriptive 44-condition literature coverage and the existing 16-PDB archived examples.

Run publication-integrity checks from the repository root:

```bash
python validation/verify_publication.py
```

The secondary cases illustrate limitations of applying a static native-state fold to adsorption-associated structural changes. Do not interpret these values as adsorption free energies or predicted docking orientations.
