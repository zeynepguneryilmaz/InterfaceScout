# Publication reproducibility

This directory contains the frozen workflows used for the InterfaceScout manuscript.

## Final development panel

`select_core_parameters_final.py` defines the eight-protein adsorption-label-free development panel. `select_core_parameters_consensus.py` performs the numerical and spatial comparison and writes the selected settings and detailed diagnostics.

The final development PDBs are:

| Protein | PDB | Chain(s) |
|---|---|---|
| Crambin | 1CRN | A |
| Phage 434 repressor N-terminal domain | 1R69 | A |
| Alpha-spectrin SH3 domain | 1SHG | A |
| Native FKBP12 | 2PPN | A |
| Unligated adenylate kinase | 4AKE | A |
| Unliganded maltodextrin-binding protein | 1OMP | A |
| Triosephosphate isomerase dimer | 1TIM | A,B |
| Open citrate synthase | 5CSC | B |

## Final literature-derived benchmark

`benchmark_experimental_final.json` is the only benchmark manifest used for the manuscript. It contains seven primary conditions and two secondary challenge cases.

| Set | System | PDB / chains | pH |
|---|---|---|---:|
| Primary | BSA–silica | 4F5S / A | 8.0 |
| Primary | HSA–PAA/Fe3O4 | 2VUF / A | 6.2 |
| Primary | HCAII–silica | 2CBA / A | 6.3 |
| Primary | HCAII–silica | 2CBA / A | 8.3 |
| Primary | Ubiquitin–citrate/Au | 1UBQ / A | 7.7 |
| Primary | beta2-microglobulin–MPA/Au | 1JNJ / A | 7.0 |
| Primary | L-asparaginase II–Al(OH)3 | 3ECA / A,B,C,D | 7.5 |
| Secondary | Transferrin–polystyrene | 1SUV / C,E | 7.4 |
| Secondary | Catalase–polystyrene | 1TGU / A,B,C,D | 7.4 |

## Files

- `PROTOCOL.md` — frozen development and benchmark protocol.
- `benchmark_experimental_final.json` — final experimental benchmark and source provenance.
- `select_core_parameters_consensus.py` — numerical convergence and radius-pair consensus calculations.
- `select_core_parameters_final.py` — final development-panel definition.
- `generate_experimental_predictions.py` — prediction-first benchmark generation and PDB archiving.
- `compare_experimental_predictions_v2.py` — comparison of saved predictions with experimental annotations.
- `summarize_experimental_validation.py` — compact case-level manuscript summary.
- `reproduce_publication.py` — one-command publication reproduction and verification.

## Generated publication package

Running:

```bash
cd backend
python -m v2.validation.reproduce_publication
```

writes `../publication_data/` with:

```text
publication_data/
├── development/
│   ├── core_parameter_selection_consensus.json
│   └── prepared_pdbs/
├── benchmark/
│   ├── manifest.json
│   ├── predictions/
│   │   ├── raw/
│   │   ├── residues/
│   │   ├── patches/
│   │   ├── pdb/full_rcsb/
│   │   ├── pdb/prepared_for_analysis/
│   │   └── pdb/propensity_tracks/
│   └── comparison/
└── verification.json
```

The exact prepared PDBs, residue-level values, patch tables, propensity tracks, and case-level comparison files are therefore available for direct audit.

## Interpretation

The benchmark is a spatial concordance test between InterfaceScout candidate regions and independently reported experimental protein-side localization. The experimental studies differ in spatial resolution, so direct overlap, 5 Å proximity, and 8 Å patch-scale proximity are reported separately.

See `PROTOCOL.md` for the full interpretation rules.
