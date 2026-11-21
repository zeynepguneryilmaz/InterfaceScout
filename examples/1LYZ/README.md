# Example: lysozyme (1LYZ)

**Target:** RCSB PDB `1LYZ`, chain `A`, pH `7.4`. The application computes all 11 canonical chemistry maps, solvent-accessibility descriptors, and ranked candidate regions. No material-specific score is fitted.

## Run from the web interface

1. Install and launch InterfaceScout using the platform-specific guide in `docs/`.
2. Enter PDB ID `1LYZ`, choose chain `A`, set pH `7.4`, and click **Analyze**.
3. Inspect individual maps, then select two or more chemistry classes for a combined view.
4. Use **Complete export** to save the genuine viewer snapshots, map-labeled B-factor PDBs, CSV/Excel/JSON outputs and PDF reports. Place the exported ZIP in `examples/1LYZ/results/` if you intend to archive it.

## Scripted, reproducible data example

With InterfaceScout already running at http://127.0.0.1:8000, from the repository root:

```bash
python examples/1LYZ/run_example.py
```

The script calls the same local `/analyze` endpoint as the interface and writes:
- `results/1LYZ_analysis.json`: full server result;
- `results/1LYZ_residues.csv`: per-map residue scores;
- `results/1LYZ_patch_members.csv`: actual ranked patch membership;
- `results/run_manifest.json`: source and run parameters.

For offline analysis, supply a local 1LYZ PDB file:
```bash
python examples/1LYZ/run_example.py --pdb /path/to/1LYZ.pdb
```

Without `--pdb`, the script supplies `pdb_id=1LYZ`, so the backend may need network access to retrieve the structure. The script **does not install dependencies** and does not fabricate results when the local service is unavailable.

Do not confuse `1LYZ_residues.csv` (per-residue chemistry scores) with the patch-membership CSV (which residues belong to candidate regions). When using B-factor PDB exports, the B-factor field contains InterfaceScout visualization values, not experimental crystallographic B-factors.
