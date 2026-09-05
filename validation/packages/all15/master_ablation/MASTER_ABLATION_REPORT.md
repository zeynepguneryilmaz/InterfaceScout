# InterfaceScout master ablation report

## Interpretation rule

All currently scored conditions are developmental/prefreeze evidence. Do not call them blind prospective validation. The 15-protein package contains Tier-1 residue-level, Tier-2 orientation/spatial, Tier-3 stress-test, and pending-audit conditions; only residue-level conditions with sufficiently explicit ground truth enter AUROC/AP tables.

## Component consistency versus Core IS

| stage           |   n_conditions |   AUROC_improved |   AUROC_worsened |   AP_improved |   AP_worsened |   R10_improved |   R10_worsened |   R20_improved |   R20_worsened |   median_AUROC |   median_AP |   median_R10 |   median_R20 |
|:----------------|---------------:|-----------------:|-----------------:|--------------:|--------------:|---------------:|---------------:|---------------:|---------------:|---------------:|------------:|-------------:|-------------:|
| apbs            |             10 |                2 |                2 |             3 |             1 |              2 |              0 |              1 |              0 |          0.617 |       0.086 |        0.179 |        0.286 |
| gnm             |             10 |                8 |                2 |             6 |             4 |              3 |              3 |              4 |              4 |          0.685 |       0.090 |        0.050 |        0.243 |
| gnm_apbs        |             10 |                9 |                1 |             7 |             3 |              3 |              2 |              5 |              3 |          0.689 |       0.123 |        0.106 |        0.345 |
| radial          |             10 |                9 |                1 |             8 |             2 |              4 |              2 |              4 |              3 |          0.753 |       0.125 |        0.188 |        0.297 |
| rin_betweenness |             10 |                9 |                1 |             8 |             2 |              4 |              2 |              5 |              3 |          0.739 |       0.132 |        0.132 |        0.321 |
| rin_closeness   |             10 |                8 |                2 |             7 |             3 |              4 |              3 |              4 |              3 |          0.757 |       0.132 |        0.077 |        0.297 |
| rin_degree      |             10 |                9 |                1 |             8 |             2 |              4 |              2 |              4 |              3 |          0.755 |       0.138 |        0.148 |        0.304 |

## Per-condition AUROC / AP

| condition_id   | protein             | model_label                                       |   AUROC |    AP |   Recall@10 |   Recall@20 |   delta_AUROC_vs_core |   delta_AP_vs_core |
|:---------------|:--------------------|:--------------------------------------------------|--------:|------:|------------:|------------:|----------------------:|-------------------:|
| FN_COO         | Fibronectin III8-10 | Core IS                                           |   0.968 | 0.100 |       1.000 |       1.000 |                 0.000 |              0.000 |
| FN_NH3         | Fibronectin III8-10 | Core IS                                           |   0.564 | 0.024 |       0.000 |       0.000 |                 0.000 |              0.000 |
| FN_CH3_HEAD    | Fibronectin III8-10 | Core IS                                           |   0.420 | 0.030 |       0.000 |       0.000 |                 0.000 |              0.000 |
| FN_CH3_SIDE    | Fibronectin III8-10 | Core IS                                           |   0.296 | 0.022 |       0.000 |       0.000 |                 0.000 |              0.000 |
| FN_CH3_BETA    | Fibronectin III8-10 | Core IS                                           |   0.630 | 0.073 |       0.000 |       0.273 |                 0.000 |              0.000 |
| CYTC_CH3       | Cytochrome c        | Core IS                                           |   0.538 | 0.088 |       0.000 |       0.000 |                 0.000 |              0.000 |
| CYTC_COOH      | Cytochrome c        | Core IS                                           |   0.407 | 0.076 |       0.000 |       0.286 |                 0.000 |              0.000 |
| LYZ_SWCNT      | Lysozyme            | Core IS                                           |   0.710 | 0.462 |       0.385 |       0.538 |                 0.000 |              0.000 |
| CHT_CNT_HYDRO  | Alpha-chymotrypsin  | Core IS                                           |   0.725 | 0.207 |       0.300 |       0.400 |                 0.000 |              0.000 |
| MB_CIT_AUNP    | Myoglobin           | Core IS                                           |   0.650 | 0.278 |       0.214 |       0.286 |                 0.000 |              0.000 |
| FN_COO         | Fibronectin III8-10 | Core + GNM                                        |   0.930 | 0.048 |       0.000 |       0.000 |                -0.039 |             -0.052 |
| FN_NH3         | Fibronectin III8-10 | Core + GNM                                        |   0.689 | 0.026 |       0.000 |       0.000 |                 0.125 |              0.001 |
| FN_CH3_HEAD    | Fibronectin III8-10 | Core + GNM                                        |   0.680 | 0.082 |       0.111 |       0.333 |                 0.260 |              0.052 |
| FN_CH3_SIDE    | Fibronectin III8-10 | Core + GNM                                        |   0.493 | 0.032 |       0.000 |       0.000 |                 0.197 |              0.011 |
| FN_CH3_BETA    | Fibronectin III8-10 | Core + GNM                                        |   0.816 | 0.335 |       0.273 |       0.364 |                 0.186 |              0.262 |
| CYTC_CH3       | Cytochrome c        | Core + GNM                                        |   0.562 | 0.098 |       0.000 |       0.286 |                 0.024 |              0.010 |
| CYTC_COOH      | Cytochrome c        | Core + GNM                                        |   0.437 | 0.072 |       0.000 |       0.143 |                 0.029 |             -0.004 |
| LYZ_SWCNT      | Lysozyme            | Core + GNM                                        |   0.659 | 0.280 |       0.154 |       0.385 |                -0.052 |             -0.182 |
| CHT_CNT_HYDRO  | Alpha-chymotrypsin  | Core + GNM                                        |   0.773 | 0.148 |       0.100 |       0.200 |                 0.048 |             -0.058 |
| MB_CIT_AUNP    | Myoglobin           | Core + GNM                                        |   0.765 | 0.410 |       0.286 |       0.571 |                 0.116 |              0.132 |
| FN_COO         | Fibronectin III8-10 | Core + conditional APBS                           |   0.996 | 0.500 |       1.000 |       1.000 |                 0.028 |              0.400 |
| FN_NH3         | Fibronectin III8-10 | Core + conditional APBS                           |   0.721 | 0.078 |       0.500 |       0.500 |                 0.157 |              0.053 |
| FN_CH3_HEAD    | Fibronectin III8-10 | Core + conditional APBS                           |   0.420 | 0.030 |       0.000 |       0.000 |                 0.000 |              0.000 |
| FN_CH3_SIDE    | Fibronectin III8-10 | Core + conditional APBS                           |   0.296 | 0.022 |       0.000 |       0.000 |                 0.000 |              0.000 |
| FN_CH3_BETA    | Fibronectin III8-10 | Core + conditional APBS                           |   0.630 | 0.073 |       0.000 |       0.273 |                 0.000 |              0.000 |
| CYTC_CH3       | Cytochrome c        | Core + conditional APBS                           |   0.538 | 0.088 |       0.000 |       0.000 |                 0.000 |              0.000 |
| CYTC_COOH      | Cytochrome c        | Core + conditional APBS                           |   0.400 | 0.083 |       0.143 |       0.286 |                -0.007 |              0.007 |
| LYZ_SWCNT      | Lysozyme            | Core + conditional APBS                           |   0.710 | 0.462 |       0.385 |       0.538 |                 0.000 |              0.000 |
| CHT_CNT_HYDRO  | Alpha-chymotrypsin  | Core + conditional APBS                           |   0.725 | 0.207 |       0.300 |       0.400 |                 0.000 |              0.000 |
| MB_CIT_AUNP    | Myoglobin           | Core + conditional APBS                           |   0.603 | 0.207 |       0.214 |       0.286 |                -0.047 |             -0.072 |
| FN_COO         | Fibronectin III8-10 | Core + GNM + conditional APBS                     |   0.989 | 0.250 |       1.000 |       1.000 |                 0.021 |              0.150 |
| FN_NH3         | Fibronectin III8-10 | Core + GNM + conditional APBS                     |   0.802 | 0.048 |       0.000 |       0.500 |                 0.239 |              0.024 |
| FN_CH3_HEAD    | Fibronectin III8-10 | Core + GNM + conditional APBS                     |   0.680 | 0.082 |       0.111 |       0.333 |                 0.260 |              0.052 |
| FN_CH3_SIDE    | Fibronectin III8-10 | Core + GNM + conditional APBS                     |   0.493 | 0.032 |       0.000 |       0.000 |                 0.197 |              0.011 |
| FN_CH3_BETA    | Fibronectin III8-10 | Core + GNM + conditional APBS                     |   0.816 | 0.335 |       0.273 |       0.364 |                 0.186 |              0.262 |
| CYTC_CH3       | Cytochrome c        | Core + GNM + conditional APBS                     |   0.562 | 0.098 |       0.000 |       0.286 |                 0.024 |              0.010 |
| CYTC_COOH      | Cytochrome c        | Core + GNM + conditional APBS                     |   0.444 | 0.075 |       0.000 |       0.143 |                 0.037 |             -0.001 |
| LYZ_SWCNT      | Lysozyme            | Core + GNM + conditional APBS                     |   0.659 | 0.280 |       0.154 |       0.385 |                -0.052 |             -0.182 |
| CHT_CNT_HYDRO  | Alpha-chymotrypsin  | Core + GNM + conditional APBS                     |   0.773 | 0.148 |       0.100 |       0.200 |                 0.048 |             -0.058 |
| MB_CIT_AUNP    | Myoglobin           | Core + GNM + conditional APBS                     |   0.699 | 0.361 |       0.286 |       0.357 |                 0.049 |              0.083 |
| FN_COO         | Fibronectin III8-10 | Core + GNM + conditional APBS + radial prominence |   0.972 | 0.111 |       1.000 |       1.000 |                 0.004 |              0.011 |
| FN_NH3         | Fibronectin III8-10 | Core + GNM + conditional APBS + radial prominence |   0.816 | 0.029 |       0.000 |       0.000 |                 0.253 |              0.005 |
| FN_CH3_HEAD    | Fibronectin III8-10 | Core + GNM + conditional APBS + radial prominence |   0.787 | 0.121 |       0.222 |       0.333 |                 0.367 |              0.091 |
| FN_CH3_SIDE    | Fibronectin III8-10 | Core + GNM + conditional APBS + radial prominence |   0.452 | 0.030 |       0.000 |       0.000 |                 0.157 |              0.008 |
| FN_CH3_BETA    | Fibronectin III8-10 | Core + GNM + conditional APBS + radial prominence |   0.868 | 0.204 |       0.273 |       0.364 |                 0.238 |              0.131 |
| CYTC_CH3       | Cytochrome c        | Core + GNM + conditional APBS + radial prominence |   0.625 | 0.129 |       0.286 |       0.286 |                 0.087 |              0.041 |
| CYTC_COOH      | Cytochrome c        | Core + GNM + conditional APBS + radial prominence |   0.498 | 0.080 |       0.000 |       0.143 |                 0.091 |              0.004 |
| LYZ_SWCNT      | Lysozyme            | Core + GNM + conditional APBS + radial prominence |   0.631 | 0.212 |       0.154 |       0.308 |                -0.079 |             -0.250 |
| CHT_CNT_HYDRO  | Alpha-chymotrypsin  | Core + GNM + conditional APBS + radial prominence |   0.762 | 0.129 |       0.000 |       0.100 |                 0.037 |             -0.078 |
| MB_CIT_AUNP    | Myoglobin           | Core + GNM + conditional APBS + radial prominence |   0.745 | 0.395 |       0.357 |       0.429 |                 0.095 |              0.117 |
| FN_COO         | Fibronectin III8-10 | Current IS + RIN degree peripheralness            |   0.982 | 0.167 |       1.000 |       1.000 |                 0.014 |              0.067 |
| FN_NH3         | Fibronectin III8-10 | Current IS + RIN degree peripheralness            |   0.887 | 0.049 |       0.000 |       0.500 |                 0.323 |              0.025 |
| FN_CH3_HEAD    | Fibronectin III8-10 | Current IS + RIN degree peripheralness            |   0.834 | 0.137 |       0.222 |       0.333 |                 0.414 |              0.107 |
| FN_CH3_SIDE    | Fibronectin III8-10 | Current IS + RIN degree peripheralness            |   0.562 | 0.038 |       0.000 |       0.000 |                 0.266 |              0.016 |
| FN_CH3_BETA    | Fibronectin III8-10 | Current IS + RIN degree peripheralness            |   0.875 | 0.237 |       0.182 |       0.273 |                 0.245 |              0.164 |
| CYTC_CH3       | Cytochrome c        | Current IS + RIN degree peripheralness            |   0.575 | 0.106 |       0.143 |       0.286 |                 0.037 |              0.018 |
| CYTC_COOH      | Cytochrome c        | Current IS + RIN degree peripheralness            |   0.502 | 0.078 |       0.000 |       0.143 |                 0.094 |              0.002 |
| LYZ_SWCNT      | Lysozyme            | Current IS + RIN degree peripheralness            |   0.567 | 0.178 |       0.154 |       0.308 |                -0.144 |             -0.285 |
| CHT_CNT_HYDRO  | Alpha-chymotrypsin  | Current IS + RIN degree peripheralness            |   0.765 | 0.139 |       0.100 |       0.300 |                 0.040 |             -0.067 |
| MB_CIT_AUNP    | Myoglobin           | Current IS + RIN degree peripheralness            |   0.745 | 0.435 |       0.357 |       0.500 |                 0.096 |              0.157 |
| FN_COO         | Fibronectin III8-10 | Current IS + RIN closeness peripheralness         |   0.965 | 0.091 |       0.000 |       1.000 |                -0.004 |             -0.009 |
| FN_NH3         | Fibronectin III8-10 | Current IS + RIN closeness peripheralness         |   0.823 | 0.028 |       0.000 |       0.000 |                 0.260 |              0.004 |
| FN_CH3_HEAD    | Fibronectin III8-10 | Current IS + RIN closeness peripheralness         |   0.842 | 0.149 |       0.333 |       0.333 |                 0.422 |              0.119 |
| FN_CH3_SIDE    | Fibronectin III8-10 | Current IS + RIN closeness peripheralness         |   0.494 | 0.033 |       0.000 |       0.000 |                 0.198 |              0.012 |
| FN_CH3_BETA    | Fibronectin III8-10 | Current IS + RIN closeness peripheralness         |   0.881 | 0.215 |       0.273 |       0.364 |                 0.251 |              0.141 |
| CYTC_CH3       | Cytochrome c        | Current IS + RIN closeness peripheralness         |   0.633 | 0.140 |       0.286 |       0.286 |                 0.095 |              0.051 |
| CYTC_COOH      | Cytochrome c        | Current IS + RIN closeness peripheralness         |   0.489 | 0.078 |       0.000 |       0.143 |                 0.082 |              0.002 |
| LYZ_SWCNT      | Lysozyme            | Current IS + RIN closeness peripheralness         |   0.644 | 0.212 |       0.154 |       0.308 |                -0.066 |             -0.251 |
| CHT_CNT_HYDRO  | Alpha-chymotrypsin  | Current IS + RIN closeness peripheralness         |   0.743 | 0.123 |       0.000 |       0.100 |                 0.019 |             -0.083 |
| MB_CIT_AUNP    | Myoglobin           | Current IS + RIN closeness peripheralness         |   0.770 | 0.439 |       0.357 |       0.429 |                 0.121 |              0.161 |
| FN_COO         | Fibronectin III8-10 | Current IS + RIN betweenness peripheralness       |   0.975 | 0.125 |       1.000 |       1.000 |                 0.007 |              0.025 |
| FN_NH3         | Fibronectin III8-10 | Current IS + RIN betweenness peripheralness       |   0.839 | 0.044 |       0.000 |       0.500 |                 0.276 |              0.020 |
| FN_CH3_HEAD    | Fibronectin III8-10 | Current IS + RIN betweenness peripheralness       |   0.822 | 0.124 |       0.111 |       0.333 |                 0.402 |              0.095 |
| FN_CH3_SIDE    | Fibronectin III8-10 | Current IS + RIN betweenness peripheralness       |   0.577 | 0.039 |       0.000 |       0.000 |                 0.282 |              0.017 |
| FN_CH3_BETA    | Fibronectin III8-10 | Current IS + RIN betweenness peripheralness       |   0.903 | 0.231 |       0.273 |       0.364 |                 0.273 |              0.158 |
| CYTC_CH3       | Cytochrome c        | Current IS + RIN betweenness peripheralness       |   0.665 | 0.187 |       0.286 |       0.286 |                 0.127 |              0.098 |
| CYTC_COOH      | Cytochrome c        | Current IS + RIN betweenness peripheralness       |   0.538 | 0.084 |       0.000 |       0.143 |                 0.131 |              0.009 |
| LYZ_SWCNT      | Lysozyme            | Current IS + RIN betweenness peripheralness       |   0.608 | 0.203 |       0.154 |       0.308 |                -0.102 |             -0.260 |
| CHT_CNT_HYDRO  | Alpha-chymotrypsin  | Current IS + RIN betweenness peripheralness       |   0.758 | 0.139 |       0.000 |       0.300 |                 0.033 |             -0.068 |
| MB_CIT_AUNP    | Myoglobin           | Current IS + RIN betweenness peripheralness       |   0.721 | 0.445 |       0.357 |       0.500 |                 0.071 |              0.166 |

## All 15 proteins

| protein                  | pdb   | status           | primary_validation_role                                                                                                   |
|:-------------------------|:------|:-----------------|:--------------------------------------------------------------------------------------------------------------------------|
| BSA                      | 4F5S  | development-seen | full-package application/external validation; BSA-specific literature conditions pending final residue-ground-truth audit |
| Fibronectin III8-10      | 1FNF  | unseen           | charged/hydrophobic SAM chemistry discrimination; HAp extension                                                           |
| Cytochrome c             | 3NWV  | unseen           | functionalized silica chemistry discrimination; curvature stress test                                                     |
| Lysozyme                 | 1LYZ  | development-seen | SWCNT residue-level external validation                                                                                   |
| Alpha-chymotrypsin       | 4CHA  | unseen           | CNT hydrophobic residue-level validation plus chemistry negative control                                                  |
| Myoglobin                | 1MBN  | unseen           | citrate-AuNP residue-level/persistent-contact validation                                                                  |
| Hemoglobin               | 2HHB  | unseen           | citrate-AuNP residue-level/persistent-contact validation                                                                  |
| Trypsin                  | 2PTN  | unseen           | citrate-AuNP residue-level/persistent-contact validation                                                                  |
| RNase A                  | 7RSA  | unseen           | silica residue contacts plus charged-SAM orientation                                                                      |
| Acetylcholinesterase     | 1EEA  | unseen           | charged-SAM/citrate-AuNP orientation and active-site accessibility                                                        |
| WNV E protein domain III | 2HG0  | unseen           | graphene residue contacts/contact-frequency validation                                                                    |
| GB3 / Protein G B1       | 2OED  | unseen           | citrate-AuNP orientation validation                                                                                       |
| Proteinase K             | 5B1E  | unseen           | citrate-AuNP orientation/activity validation                                                                              |
| Carbonic anhydrase II    | 1CA2  | unseen           | citrate-AuNP orientation/activity validation                                                                              |
| Fibrinogen               | 3GHG  | unseen           | citrate-AuNP orientation; graphene mechanism/spatial context                                                              |

## Condition registry

| condition_id      | protein                  | surface_or_condition       | tier             | ground_truth_status            | package_status   |
|:------------------|:-------------------------|:---------------------------|:-----------------|:-------------------------------|:-----------------|
| BSA_TIO2          | BSA                      | TiO2                       | Tier1-candidate  | PENDING_AUDIT                  | planned          |
| BSA_SILICA        | BSA                      | silica                     | Tier1-candidate  | PENDING_AUDIT                  | planned          |
| BSA_GRAPHITIC     | BSA                      | graphitic carbon           | Tier1-candidate  | PENDING_AUDIT                  | planned          |
| FN_COO            | Fibronectin III8-10      | COO- SAM                   | Tier1            | READY                          | ready            |
| FN_NH3            | Fibronectin III8-10      | NH3+ SAM                   | Tier1            | READY                          | ready            |
| FN_CH3_HEAD       | Fibronectin III8-10      | CH3 SAM head-on            | Tier1            | READY                          | ready            |
| FN_CH3_SIDE       | Fibronectin III8-10      | CH3 SAM side-on            | Tier1            | READY                          | ready            |
| FN_CH3_BETA       | Fibronectin III8-10      | CH3 SAM beta-on            | Tier1            | READY                          | ready            |
| FN_OH             | Fibronectin III8-10      | OH SAM                     | Tier1-candidate  | READY_CHEMISTRY_AUDIT_REQUIRED | blocked          |
| CYTC_CH3          | Cytochrome c             | CH3-functionalized silica  | Tier1            | READY                          | ready            |
| CYTC_COOH         | Cytochrome c             | COOH-functionalized silica | Tier1            | READY                          | ready            |
| CYTC_OH           | Cytochrome c             | OH-functionalized silica   | Tier1-candidate  | READY_CHEMISTRY_AUDIT_REQUIRED | blocked          |
| CYTC_NH2          | Cytochrome c             | NH2-functionalized silica  | Tier1-candidate  | PENDING_EXHAUSTIVE_GT_AUDIT    | pending          |
| CYTC_SILICA_4NM   | Cytochrome c             | 4-nm silica                | Tier3            | READY                          | ready            |
| CYTC_SILICA_11NM  | Cytochrome c             | 11-nm silica               | Tier3            | READY                          | ready            |
| LYZ_SWCNT         | Lysozyme                 | SWCNT                      | Tier1            | READY                          | computed         |
| CHT_CNT_HYDRO     | Alpha-chymotrypsin       | CNT                        | Tier1            | READY                          | computed         |
| CHT_CNT_PI_NEG    | Alpha-chymotrypsin       | CNT                        | Negative-control | READY                          | computed         |
| MB_CIT_AUNP       | Myoglobin                | citrate-capped AuNP        | Tier1            | READY                          | ready            |
| HB_CIT_AUNP       | Hemoglobin               | citrate-capped AuNP        | Tier1-candidate  | PENDING_EXACT_RESIDUE_AUDIT    | pending          |
| TRP_CIT_AUNP      | Trypsin                  | citrate-capped AuNP        | Tier1-candidate  | PENDING_EXACT_RESIDUE_AUDIT    | pending          |
| RNASE_SILICA_4NM  | RNase A                  | 4-nm silica                | Tier1-candidate  | PENDING_SOURCE_RECHECK         | pending          |
| RNASE_SILICA_11NM | RNase A                  | 11-nm silica               | Tier3            | PENDING_SOURCE_RECHECK         | pending          |
| RNASE_COOH_SAM    | RNase A                  | COOH SAM                   | Tier2            | READY                          | ready            |
| RNASE_NH2_SAM     | RNase A                  | NH2 SAM                    | Tier2            | READY                          | ready            |
| ACHE_COOH_SAM     | Acetylcholinesterase     | COOH SAM                   | Tier2            | READY                          | ready            |
| ACHE_NH2_SAM      | Acetylcholinesterase     | NH2 SAM                    | Tier2            | READY                          | ready            |
| ACHE_CIT_AUNP     | Acetylcholinesterase     | citrate-capped AuNP        | Tier2            | READY                          | ready            |
| WNV_GRAPHENE      | WNV E protein domain III | graphene                   | Tier1            | READY                          | ready            |
| GB3_CIT_AUNP      | GB3 / Protein G B1       | citrate-capped AuNP        | Tier2            | READY                          | ready            |
| PK_CIT_AUNP       | Proteinase K             | citrate-capped AuNP        | Tier2            | READY                          | ready            |
| CA2_CIT_AUNP      | Carbonic anhydrase II    | citrate-capped AuNP        | Tier2            | READY                          | ready            |
| FIB_CIT_AUNP      | Fibrinogen               | citrate-capped AuNP        | Tier2            | READY                          | ready            |
| FIB_GRAPHENE      | Fibrinogen               | graphene                   | Tier2-candidate  | PENDING_GT_AUDIT               | pending          |