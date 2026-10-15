# Original-benchmark overlap and method-stage pH audit — 2026-10-15

Audit only; original frozen benchmark and prospective ground truth/scoring remain unchanged.

## Directly verified in frozen `validation/benchmark_manifest.json`
- `hsa_paa_fe3o4` is **already included** in original primary Level-3 benchmark as HSA PDB 2VUF, peptide ranges 373–389, 403–410, 414–428, PAA–Fe3O4, and pH 6.2. Critical issue: primary article methods incubated HSA and coated Fe3O4 in deionized water overnight, then performed crosslink capture in MES pH 6.2. Thus 6.2 is **capture-stage pH**, not demonstrated adsorption-stage pH. Preserve original frozen score, flag the experiment-stage mismatch in manuscript/benchmark evidence audit, and **do not** promote L3E028 to independent prospective ADMIT under strict source pH protocol.
- `hcaii_silica_ph63` (pH 6.3; fluorescent H10 anchor) and `hcaii_silica_ph83` (pH 8.3; fluorescent residue-37 anchor) already occur in original primary Level-3 benchmark, source DOI `10.1529/biophysj.104.054809`. L3E034 is confirmatory/duplicate, not independent. Existing original benchmark labels must remain frozen.
- `ubiquitin_citrate_au` already occurs in original primary benchmark (pH 7.7, ranges 2–3/15–18). L3E011 NP1 and L3E012 silica-shell NP2 belong to a different study but ubiquitin–Au family results must be clustered/flagged rather than counted as new independent families.
- Original frozen benchmark has seven primary conditions and two secondary challenge conditions per manifest.

## New literature lead
- L3E033: Site-resolved sequence-unassigned solid-state NMR class-based interface description for ubiquitin–MCM41 silica at pH 5, DOI 10.1116/1.4983273: unsuitable for strict site-level prospective scoring without assignment.
- Do not create gold labels from chemical-residue classes or figure coloring where sequence-site mapping is not explicitly supported.

## Consequence
No new eligible independent prospective ADMIT is supported by the present source audit; prospective score remains 2 conditions. Strict prospective eligibility and legacy benchmark's published scoring provenance must be presented distinctly. Do not use retrospective discovery to silently alter legacy frozen pH or scoring.
