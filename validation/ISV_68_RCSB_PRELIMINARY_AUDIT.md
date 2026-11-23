# ISV expanded screening — 68 protein structures, initial RCSB record audit

**Status (23 November 2026):** PROVISIONAL CANDIDATE POOL — **not 68 accepted protein structures**.

[Candidate CSV](ISV_protein_candidate_pool_68_v1.csv) · [Reproducible metadata script](isv_rcsb_screen.py) · [Successful original RCSB screening run](https://github.com/zeynepguneryilmaz/InterfaceScout/actions/runs/45259307414)

The experimental-method and PDB-record screening was run on **68 distinct named proteins** against RCSB structural data. The 68 proposed PDB IDs were compared against the exact frozen TV48 PDB identifiers, and no exact ID overlap was found. **PDB ID non-overlap is not proof of nonhomology or holdout independence.**

| RCSB metadata screen | Count |
| --- | ---: |
| Submitted structure IDs with accessible protein-bearing RCSB metadata | **68 / 68** |
| Experimental X-ray structures | **64** |
| Solution NMR structures | **4** |
| Preliminary metadata checks without threshold flags | **55** |
| Structural metadata requiring review (mainly crystallographic resolution >2.5 Å) | **13** |
| Inaccessible/unverifiable RCSB records | **0** |

**Important:** The status 'METADATA_PASS_ONLY' means the RCSB record exists, a protein polymer entity exists, the method is recognized, and a crystallographic-resolution threshold issue was not detected. It **does not** establish acceptable complete side chains, exact chain selection, full-length construct, correct species, complete domain representation, or InterfaceScout 11-map acceptance. No final ISV membership is claimed.

## Main structural groups

The pool was designed to cover common protein–surface research examples:

- **Serum/corona and biointerface components:** serum albumins (human and bovine), transferrin, lactoferrin, fibrinogen, antibody, apolipoproteins, coagulation and complement proteins, vitronectin, tenascin and fibronectin domains.
- **Common model proteins:** lysozyme, RNase A, cytochrome c, myoglobin, ubiquitin, carbonic anhydrase, hemoglobin, beta-lactoglobulin, ovalbumin and other standard soluble folds.
- **Immobilized enzymes and protein materials:** catalase, peroxidase, proteinase K, lipase, glucose oxidase, alkaline phosphatase, dehydrogenase, urease, beta-galactosidase, cellulase, laccase, avidin and streptavidin.
- **Other distinct biochemical surface examples:** thioredoxin, SOD, ferritin, fluorescent proteins, acylphosphatase, GB3 and R2ab.

**Known cases requiring extra care:**
- **1TEN** is a **tenascin-derived fibronectin-type-III repeat**, **not fibronectin proper**. The candidate-name field was corrected after inspecting RCSB polymer metadata and source annotation.
- **1S0Q** is titled native bovine trypsin in the archive, but the deposited polymer entity is called trypsinogen. Both the structural annotation and experimental activation state must be recorded; the candidate has been explicitly annotated to avoid conflating biological processing states.
- **1AV1 (ApoA-I), 6TAV (alpha-2-macroglobulin), 2A73 (C3), 1KCT (alpha-1 antitrypsin), 3GHG (fibrinogen)** are examples above the conservative crystallographic **2.5 Å review boundary**. This is a review flag, **not automatic exclusion**.
- **1LPE, 3R62, 2OK5, 1NL2, 6O5E and 1TEN** may represent specific domains or partial constructs. Such structures can still be usable for an explicitly domain-level analysis, but cannot represent the entire protein without qualification.
- **2OED, 7KWI, 2FNB and 1JNJ** use **solution NMR**; check structural ensemble/model choice and residue completeness rather than inventing a crystallographic resolution.

## Admission criteria (next stage)

At minimum, each prospective accepted entry needs:
1. Exact deposited protein identity, species, construct and domain span; biological assembly selection and valid chain ID.
2. Experimental structure-quality review, deposited versus observed residue correspondence, missing internal backbone and side-chain atoms, nonstandard residues and engineered mutations; no blanket exclusion of NMR.
3. Explicit glycan/cofactor/ligand treatment and applicability of protein-only InterfaceScout representations.
4. Family/domain redundancy versus TV48 and within the ISV pool; a second similar family member is not an independent replicate.
5. Final 11-map smoke check, at least one interpretable exposed-surface profile, input SHA-256 and reproducible outputs.
6. Optional experimentally supported localization evidence recorded in its **own column**, with method and reference; protein inclusion does **not** require experimental adsorption-site localization.

**Next target:** construct a rigorous accepted ISV manifest of **at least 50 proteins** from the screened pool and, if needed, additional verified replacements. Do not claim at least 50 accepted until these checks pass.

Earlier eight-protein parameter-development and nine-case retrospective evidence panels remain archived and do not constitute separate active datasets.
