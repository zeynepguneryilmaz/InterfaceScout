# InterfaceScout research datasets

Version 1.1.0; model `2.5.0-dual-accessibility`. These packages preserve actual source-run outputs, not new calculations. Each ZIP contains its source-run metadata and file-level SHA-256 inventory. The repository software release is unchanged.

TV48 is the 48-structure technical verification panel; ISV76 is the 76-structure application panel. Retrospective experimental cases are separate and do not constitute an independent 76-case experimental benchmark.

The Protein-sol archive preserves the actual server audit, including incomplete/failed cases; the external aggregate also retains historical P2Rank metrics, which are not a manuscript comparator. The original ISV76 APBS records and the later paired solver outputs are distinct. Large browser/PDF atlases are not bundled.

Use the exact source commits recorded in each ZIP for reproduction. Each panel's local protein manifest is a byte-identical snapshot of its authoritative validation manifest, retained with the dataset for reproducibility. Run conclusions (including partial/failed aggregate or external-server jobs) are explicitly recorded in [MANIFEST.json](MANIFEST.json). See [publication evidence](../validation/PUBLICATION_EVIDENCE.md), [active datasets](../validation/ACTIVE_DATASETS.md), and the [validation scripts](../validation).

| Package | Source run | Files | SHA-256 |
|---|---|---:|---|
| [TV48/TV48_results.zip](TV48/TV48_results.zip) | [45244648958](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45244648958) | 25 | `efb2fd32a283cbe35058a889a79f11e4c51476f55c2c2a52f3fe724a26cb87a5` |
| [TV48/TV48_extended_verification.zip](TV48/TV48_extended_verification.zip) | [45252414860](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45252414860) | 34 | `6372cd1fca5049cccdebfbd1f1220275c8ee1a96412a4883aa2391819321f225` |
| [TV48/TV48_numerical_sensitivity.zip](TV48/TV48_numerical_sensitivity.zip) | [45636665184](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45636665184) | 32 | `cd3f8f82b921b78efee8f649db7ea357e8c0de3e6739aa05a9b91da7b474a183` |
| [TV48/TV48_region_robustness.zip](TV48/TV48_region_robustness.zip) | [45663104217](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45663104217) | 32 | `18aade9a14ef4ff904f000a149b5be3e019561abec7e0fa68e46088c56d36954` |
| [ISV76/ISV76_results.zip](ISV76/ISV76_results.zip) | [45385665397](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45385665397) | 557 | `afb1ef85c4c1d59ccbcf6b05c5da8bfad69154f5daf5ae8e70999ab81f81e693` |
| [ISV76/ISV76_FreeSASA.zip](ISV76/ISV76_FreeSASA.zip) | [45445684505](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45445684505) | 247 | `f431cb4d3771a518e0a6342f429d13618ea9d4ed603c9e29635daaf61d357184` |
| [ISV76/ISV76_PEP_Patch.zip](ISV76/ISV76_PEP_Patch.zip) | [45451763460](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45451763460) | 637 | `e9bc5273f4a7fd37916e23f33bd2c44ac47e23d15fc5b30f4f8b6950a11a64ca` |
| [ISV76/ISV76_external_aggregate.zip](ISV76/ISV76_external_aggregate.zip) | [45490436050](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45490436050) | 5 | `dd7f0cc044bfc77b1cdbea92f624a9e976d4dc0febed1f90a199d78e03c41d54` |
| [ISV76/ISV76_Protein_sol_server_audit.zip](ISV76/ISV76_Protein_sol_server_audit.zip) | [45484412581](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45484412581) | 318 | `40bb82f9752fa8788684b127f69a6acee05c76233344def6b1f3f108231c2a94` |
| [ISV76/ISV76_Protein_sol_result_page_audit.zip](ISV76/ISV76_Protein_sol_result_page_audit.zip) | [45461880631](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45461880631) | 4 | `41ae90861d666b66a23b70ef538aa2e6ad3893f3f6278f715d1451273e6e8e1b` |
| [ISV76/ISV76_Protein_sol_server_diagnostic.zip](ISV76/ISV76_Protein_sol_server_diagnostic.zip) | [45484437867](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45484437867) | 15 | `202c725d4737bda2ada25102a9dff6ec551175a63fc07964cf5b173c6aab6631` |
| [ISV76/ISV76_paired_PB.zip](ISV76/ISV76_paired_PB.zip) | [45653251591](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45653251591) | 1021 | `8dd8981c57344ce670f892b08076b7635ef745a34950370039a89e738f3fe6e1` |
| [ISV76/ISV76_PB_grid_8C50.zip](ISV76/ISV76_PB_grid_8C50.zip) | [45657476304](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45657476304) | 72 | `8609608fb3ef7d0399ad1341f308b6559b3d4ae75645a1cc7f12dc1c35d30e84` |
| [ISV76/ISV76_PB_performance.zip](ISV76/ISV76_PB_performance.zip) | [45663636698](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45663636698) | 548 | `022259650c7394b91aead8f915b1cf4aafbae104380d349e923138e950919b91` |
| [ISV76/PB_solver_regression.zip](ISV76/PB_solver_regression.zip) | [45657475939](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45657475939) | 1 | `04303ee2cbd3fc0b1a4c82d98df77c9194bf51fc16e69e690d634ccbd8bdae5c` |
| [case_studies/retrospective_literature_baselines.zip](case_studies/retrospective_literature_baselines.zip) | [45244597881](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45244597881) | 18 | `299cc53ee6de64236c1ff8b6e3d521328378448a2b6f0f8f8a890c78354f9031` |
