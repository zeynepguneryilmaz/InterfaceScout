# Examples and results

This directory provides end-to-end examples **inside the InterfaceScout application repository**. Each example keeps its instructions and corresponding results nearby. These are **not separate applications**, and the repository's normal Windows/macOS/Linux launchers still apply.

| Example | Input | Results | Status |
|---|---|---|---|
| [1LYZ](1LYZ/README.md) | RCSB PDB ID `1LYZ`, chain A, pH 7.4 | Analysis JSON, 11 map residue values and patch membership CSVs | Generated on your machine; **not precomputed** |
| [literature benchmark](benchmark/README.md) | Frozen published validation case manifest | [case-level results](benchmark/results/benchmark_case_results.csv) | Archived, copied unchanged from `validation/` |

For a **complete visual export with actual viewer screenshots, B-factor PDBs, CSV, JSON, Excel and PDF**, use the application's **Complete export** control. The API-only example in `1LYZ/` provides analysis data but does not pretend to generate browser-captured screenshots.

The benchmark CSV is a **frozen literature-concordance validation**, not the 44-condition descriptive literature panel and not a claim of general adsorption predictive accuracy. Do not equate those sets.
