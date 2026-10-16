# InterfaceScout V1 — final compact manuscript decisions (2026-10-16)

- Frozen source 285; V1 successes 282/285; 3 unavailable AlphaFold structures.
- Main **quality-screened descriptive** panel: 239 original STANDARD proteins, 2629 complete chemistry rows. The full 282 remain the applicability denominator; 43 flagged cases remain in provenance/sensitivity, not erased.
- Do not adopt the new 97-protein stringent filter as the primary panel: it is post-hoc and local/biological maturity validation is incomplete.
- For the 239: median unique rank-1 centers 7/11; mean 55-pair rank-1 Jaccard 0.136414; median zero-overlap fraction 0.818182; median frozen top-patch coherence 0.89631 (internal score, not independent validation).
- High-overlap related channels: cationic–H-bond donor mean 1.000; anionic–H-bond acceptor 0.869976; oxide–metal coordination 0.548031. Hydrophobic–calcium-phosphate mean 0.000800.
- Signal-peptide/precursor caveat: 612/2629 rank-1 patch centers occur in precursor residues 1–30 in the 239 panel. A precursor patch is not necessarily a mature circulating-protein interface.
- Overlap is partly size-related, and 11 channels are not statistically independent. This is landscape evidence, **not adsorption or experimental accuracy**.
- Legacy experimental primary: 7/7 top-3 near-8 Å hits and 4/7 top-5 direct hits; two separate secondary structural-response cases miss at top 3.
- Separate prospective Level-3: GB3 citrate Au K13 direct rank 2; apoTf iron oxide near-8 top-3 recall 0.075 and no top-5 direct overlap. Do not pool with legacy or hide weak recall.
- Original HSA PAA-Fe3O4 pH 6.2 is crosslinking-stage rather than demonstrated initial adsorption pH. HCAII fluorescent orientation experiments use engineered constructs. Frozen labels/results preserved, claims qualified.
- MD-derived Level-C concordance remains unscored; no external agreement percentage is justified.
- Conclusion: frozen V1 generated chemistry-dependent candidate protein-side structural landscapes over a broad panel, with heterogeneous limited experimental localization concordance. Do not claim adsorption quantities, affinities, corona membership, or independently proven localization accuracy for all 282 proteins.

Data tables in this directory: Table_S1_channel_pair_Jaccard_239.csv (55 chemistry pairs), Table_S2_protein_specificity_239.csv (239 proteins). Complete original channel outputs in frozen `validation/large_panel_285_2026-10-13/results/channel_summary_all.csv` on `large-panel-freeze-2026-10-13`. All parameters and frozen benchmarks unchanged.
