# Experimental holdout audit v1

Date: 2026-09-05
Model: `IS-v1-Core-primary` (frozen before scoring)

## Validity summary

- `B2M_CIT_AUNP_WT`: **VALID prospective experimental test.** Experimental NMR perturbation list was frozen before scoring. Core IS: AUROC 0.667, AP 0.388, Recall@10 0.30, Spatial@10A_top10 0.80, median nearest top10 6.20 A. Interpretation: moderate residue discrimination with strong patch-level spatial agreement.
- `HSA_PAA_FE3O4_XLMS`: **VALID prospective experimental region-level test.** Ground truth consists of experimentally crosslinked peptide regions 373-389, 403-410, 414-428. Core IS: Spatial@5A_top10 0.00, Spatial@8A_top10 0.00, Spatial@10A_top10 0.033, median nearest top10 16.95 A. Interpretation: clear failure for this surface/system under the current generic anionic mapping; do not hide or retune against this result.
- `UBQ_AUNP_NMR`: **PROTOCOL-INVALID as a prospective quantitative test.** The pre-score GT freeze expanded the experimentally highlighted NMR sites into contiguous domains (2-3 and 15-18). After source re-check, the safely supported strongest experimental CSP residues are 2, 15, and 18 rather than every residue in those intervals. Because the IS outcome has already been viewed, the GT must not be edited and rescored as if prospective. Preserve this run only as an audit trail; any future ubiquitin analysis must be labeled retrospective/exploratory.

## Reproducibility

Run: GitHub Actions `one-shot experimental holdout v1`, run 33986766559.
Result commit: `358c813`.
Runner: `validation/run_experimental_holdout_v1.py`.
Frozen GT file: `validation/prospective_holdout/EXPERIMENTAL_GT_FREEZE_V1.csv`.

## Scientific consequence

The first true experimental holdout batch does **not** support claiming universal success. It provides one moderate-to-strong success (beta-2-microglobulin), one clear failure (HSA/PAA-Fe3O4), and one invalid protocol case (ubiquitin). This is exactly why prospective experimental validation must be retained separately from development benchmarks.
