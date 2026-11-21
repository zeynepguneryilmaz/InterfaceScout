# Examples and results

These are **real InterfaceScout example exports**, kept with the application source code. Each archived condition has a complete Excel workbook and one map-labeled B-factor PDB. The B-factor column is a visualization value, not an experimental crystallographic displacement parameter.

## Existing archived example results

| Purpose | Examples | Contents |
|---|---|---|
| Parameter robustness / development | [1CRN](1CRN/), [1R69](1R69/), [1SHG](1SHG/), [2PPN](2PPN/), [4AKE](4AKE/), [1OMP](1OMP/), [1TIM](1TIM/), [5CSC](5CSC/) | Each contains `*_all_results.xlsx` and `*_Bfactor.pdb` |
| Literature-derived cases / structural demonstrations | [4F5S](4F5S/), [2VUF](2VUF/), [2CBA](2CBA/), [1UBQ](1UBQ/), [1JNJ](1JNJ/), [3ECA](3ECA/), [1SUV](1SUV/), [1TGU](1TGU/) | As above; 2CBA contains separate pH 6.3 and pH 8.3 subfolders |
| Reproducible new run | [1LYZ](1LYZ/README.md) | Script, run instructions and `results/` destination; results are **not precomputed** |

**Important:** These are **16 unique PDB IDs and 17 sets of conditions** (2CBA is included at two pH conditions). Do not conflate the archive with the full 44 literature conditions or with all 27 structure–pH analytical states used in manuscript figures.

## Scientific validation

The frozen literature benchmark result CSV is stored **once**, in [`validation/benchmark_case_results.csv`](../validation/benchmark_case_results.csv), and its provenance is documented in [`validation/benchmark_manifest.json`](../validation/benchmark_manifest.json). The [benchmark explanation](benchmark/README.md) links to that source rather than duplicating data.

The frozen parameter-selection evidence is retained at `validation/parameter_selection.json`; it is not a runtime dependency of the examples folder. Never silently delete or relabel the validation evidence merely to simplify the tree.

## Complete export

The web application's **Complete export** creates per-structure map screenshots, detailed CSV/Excel/JSON tables and annotated PDB/PDF outputs. The archived two-file example sets above are a *subset* of the full export, not evidence that each contains every screenshot or every channel's B-factor PDB.

To generate an up-to-date 1LYZ output, use `python examples/1LYZ/run_example.py` with the app running; see the example README for limitations.
