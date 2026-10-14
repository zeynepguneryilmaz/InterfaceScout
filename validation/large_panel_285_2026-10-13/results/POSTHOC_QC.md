# Post-hoc QC for the 285-protein large panel

This QC was performed after the frozen 285-protein batch and is **descriptive only**.
No threshold, chemistry weight, protein selection, or V1 parameter was changed based on these results.

## Protein-size associations

Spearman correlations:

- UniProt length vs mean pairwise top-patch Jaccard: **-0.219**
- UniProt length vs zero-overlap channel-pair fraction: **0.370**
- UniProt length vs number of unique top-patch centers: **0.115**
- Surface-residue count vs mean pairwise top-patch Jaccard: **-0.273**
- Surface-residue count vs zero-overlap channel-pair fraction: **0.422**
- Surface-residue count vs number of unique top-patch centers: **0.174**

Therefore, the raw cross-channel zero-overlap fraction is partly size-dependent and should **not** be presented as a size-independent validation statistic. The number of distinct top-patch centers has a substantially weaker size association.

## Chemistry-pair structure

Highest mean rank-1-patch overlaps include:

- cationic vs hbond_donor: mean Jaccard **1.000**, median **1.000**
- anionic vs hbond_acceptor: mean Jaccard **0.877**, median **1.000**
- anionic vs phosphate: mean Jaccard **0.567**, median **1.000**
- hbond_acceptor vs phosphate: mean Jaccard **0.537**, median **0.667**
- oxide vs metal_coord: mean Jaccard **0.507**, median **0.556**

Lowest mean overlaps include:

- hydrophobic vs calcium_phosphate: mean Jaccard **0.001**, zero-overlap fraction **0.989**
- cationic vs hydrophobic: mean Jaccard **0.003**, zero-overlap fraction **0.986**
- hydrophobic vs hbond_donor: mean Jaccard **0.003**, zero-overlap fraction **0.986**
- hydrophobic vs oxide: mean Jaccard **0.004**, zero-overlap fraction **0.982**
- cationic vs soft_metal_sulfur: mean Jaccard **0.005**, zero-overlap fraction **0.979**

Notably, cationic and H-bond-donor top patches are identical throughout this pH 7.4 panel. The 11 canonical outputs are therefore **not 11 statistically independent dimensions**; related channels can encode strongly overlapping protein-side chemistry. This is retained transparently rather than removing a channel after viewing results.

## AlphaFold-confidence strata

- STANDARD models: n=239, median length=424, median surface residues=321, median mean-Jaccard=0.127
- Flagged models: n=43, median length=578, median surface residues=578, median mean-Jaccard=0.081

Because the flagged set also contains larger proteins on average, these differences are treated as QC context rather than evidence of a confidence effect.

## Reporting consequence

For the manuscript, the safest large-panel claim is:

> Across a frozen 285-protein source panel, InterfaceScout was executable for 282 proteins and generated chemistry-dependent candidate surface landscapes across all 11 canonical outputs. Related chemistry channels showed strong overlap, whereas dissimilar channels generally prioritized different regions. Cross-channel differentiation was partly protein-size dependent and is therefore treated as a descriptive applicability result rather than an accuracy metric.

Experimental accuracy remains supported only by the separate Level-3 interface-localization benchmark.
