# Complete InterfaceScout 299-protein residue atlas

This artifact reconstructs full ranked patch results from the **exact archived 300-PDB selection** used in the first panel analysis. The Cα-only structure 3IYK is excluded by the prior Stage-1 structure QC, leaving 299 structures.

- atlas_manifest.csv: structure-level provenance.
- map_qc.csv: all 3,289 map metrics, compared with the frozen panel.
- top_k_residue_unions.csv: all 9,867 Top-1/3/5 cumulative residue unions, protein and chemistry IDs.
- all_ranked_patches.csv: ranked patch centers, member lists, Pareto fronts and individual descriptors.
- atlas/*.json: individual PDB coordinate data and per-map Top-k residue keys for offline 3D interactive inspection.
- index.html: interactive 3Dmol.js viewer, using the vendored version of the software.
- serve.py: local viewer launch (run `python serve.py`, open http://127.0.0.1:8765). Python 3 only; no external web access needed.
- atlas_summary.json: counts and agreement checks.

**Validity:** This is a candidate-patch atlas, not an experimentally validated orientation/adsorption atlas. Top-3 and Top-5 represent unions of the first three/five ranked patches. No claim of equal performance across chemistries follows from complete map generation. The viewer displays cartoon coordinates from archived PDB structures. Validate the chemistry-conditional biological interpretation separately.

**Audit:** PDB SHA-256s from the original manifest are checked; each map's number of patches, Top-5 union size, and surface residue count are compared to the original 2026-11-16 output. Any mismatch triggers failure rather than silent acceptance.
