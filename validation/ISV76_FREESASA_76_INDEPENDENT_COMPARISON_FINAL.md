# ISV76 external FreeSASA comparison — FINAL 76/76 PASS

Date: 2026-11-24. Source: [real GitHub Actions run 45445684505](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45445684505) (all 16 shard jobs plus aggregate SUCCESS). Full raw artifact: `ISV76-External-FreeSASA-Full-Evidence` (artifact ID 14944119379).

## Design and completeness

The frozen ISV76 protein panel was run through:
- InterfaceScout production Biopython Shrake–Rupley core, **11 fixed chemistry maps and Pareto patch processing**, with optional APBS **off** for this SASA-only comparison.
- Standalone **FreeSASA 2.1.2** using `--shrake-rupley -n 200 --probe-radius 1.4 --n-threads 1 --format=rsa`. FreeSASA uses its **default ProtOr atomic radii**; InterfaceScout uses Biopython elemental radii.
- Both engines received **exactly the same normalized protein-only PDB coordinates**. SHA256 fingerprints and exact chain/sequence/insertion-code/residue-name matches were required.

**76/76 structures PASS; 0 failures; 51,931/51,931 residues paired; 0 missing.** Full per-protein RSA records, paired residue SASA CSV and resource timing/peak-RSS sampling are in the artifact. IDs such as `1E05` remain text, not scientific notation.

## Results: 76 protein macro statistics and 51,931 paired residue micro statistics

| Indicator | Observed |
|---|---:|
| Median within-protein Spearman rho (residue total SASA) | **0.99549** |
| Protein bootstrap 95% CI for median rho (5,000 structure bootstraps) | **0.99519–0.99584** |
| Minimum protein rho (8C50) | 0.98531 |
| Median of 76 protein-specific medians of absolute total SASA difference | **1.7635 Å²** |
| Bootstrap 95% CI on protein-median SASA absolute difference | **1.688–1.918 Å²** |
| Pooled 51,931-residue median absolute difference | 1.767 Å² |
| Pooled 95th-percentile absolute difference | 7.0865 Å² |
| Pooled mean signed difference (InterfaceScout – FreeSASA) | +0.9482 Å² |
| Pooled RMSE | 3.3683 Å² |
| Median within-protein surface-exposed classification agreement | **97.9256%** |
| Protein bootstrap 95% CI, exposed-state agreement | **97.7124%–98.2609%** |
| Pooled all-residue exposed-state agreement | **97.2829%** |
| Pooled side-chain absolute SASA difference, median | 1.230 Å² |

Both absolute residue areas are normalized with the same Tien et al. residue maxima and the same `tRSA >= 0.05` rule: **34,990** residues exposed in both; **15,530** buried in both; **1,283** InterfaceScout-only exposed and **128** FreeSASA-only exposed. An uneven discordance direction is consistent with radius-definition differences and/or implementation effects; do not assert a unique causal decomposition without a controlled radii-matched experiment.

## Timing and memory are NOT direct speed comparisons

These are separately measured Linux GitHub runner subprocesses, not equivalent algorithmic scopes:

| Executed work | Median per-structure wall time | Median sampled peak RSS |
|---|---:|---:|
| InterfaceScout **all 11 maps and patches; APBS off** | **1.628 s** | **143.78 MiB** |
| FreeSASA **SASA only** | **0.078 s** | **5.51 MiB** |

A 75-ms process-tree RSS polling interval imposes effective measurement floors for FreeSASA short jobs; RSS maxima are **sampled lower bounds**. **Never present the ratio as speedup, slowdown, or memory superiority**. Resource requirements rise substantially for large assemblies; e.g. the largest InterfaceScout job (6TAV, 5,234 residues) lasted **78.38 s** and recorded sampled RSS **1,161.95 MiB** in this APBS-free test.

## Scientific conclusion

Independent *implementation* agreement of solvent-exposure characterization is strong, despite small per-residue area offsets consistent with using different default atomic radii. This is NOT an independent physical-model validation (both use Shrake–Rupley), adsorption-site benchmark, binding-affinity prediction test, or material-specific scoring validation. Results can support publication claims about numerical consistency of the SASA layer and transparent surface classification; they do **not** establish InterfaceScout predictive superiority.

## Exact reproduction and failures encountered

Script: [`isv76_external_freesasa_runtime_audit.py`](isv76_external_freesasa_runtime_audit.py); workflow: [`isv76-external-freesasa-audit.yml`](../.github/workflows/isv76-external-freesasa-audit.yml). Previous failed attempts used unsupported CLI `--output-depth`, JSON tree assumptions and whitespace-split RSA residue columns. Final exact fixed-width RSA parser and CLI passed all structures. All earlier failures remain visible in historical GitHub Actions; final PASS does not erase development history.

## Complete figures and report

A publication-oriented independent analysis of the final raw artifact includes: pooled 51,931-residue SASA density, all 76 per-protein Spearman values, absolute area deviation distribution, exposed/buried agreement matrix, and scope-labeled process wall time and sampled RSS versus protein size, plus a five-page PDF and full paired-residue CSV. These are derived from, and can be regenerated against, the original full raw Actions artifact.
