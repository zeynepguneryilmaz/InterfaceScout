# InterfaceScout v1.1.0-rc1 | reproducible software publication bundle

This is the **publication candidate** for a PDB-based physicochemical protein-surface characterization tool, not an adsorption-site or binding-free-energy predictor. It is **separate** from the immutable v1.0.0 tag and `main`.

## Bundle layout

| Location | Purpose |
|---|---|
| `backend/` | Versioned scientific numerical engine and FastAPI endpoints |
| `frontend/` | Local browser interface and vendored visualization/Excel libraries |
| `tests/` | Real PDB, structural identity, boundary, export, installation and browser tests |
| `examples/` | Frozen demonstration structures and legacy output examples |
| `validation/reproducibility/inputs.json` | 21 pinned analysis conditions, input SHA-256 and provenance |
| `validation/reproducibility/model_spec.json` | Declared **scientific model** parameters and chemistry definitions |
| `validation/reproducibility/TOLERANCES.md` | Numerical acceptance contract, fixed before reproductions |
| `validation/reproduce.py` | Fresh PDB-to-all-11-maps comparison with an independent scalar oracle and real B-factor exporter |
| `validation/verify_publication.py` | Historical parameter/benchmark **aggregate consistency**, not an independent recalculation |
| `publication/VALIDATION_METHODS.md` | Manuscript-ready method descriptions and acceptance methodology |
| `publication/ACCEPTANCE.md` | Release gates and explicit blocking failures |
| `publication/TEST_MATRIX.csv` | Traceability table: each claim, test and evidence |
| `publication/KNOWN_LIMITATIONS.md` | Scientific and technical boundaries |
| `publication/build_package.py` | Strict assembly and SHA-256 inventory of the verified publication ZIP |
| `validation-evidence/numerical/` | CI-computed outputs + verification.json (bundled after tests) |
| `validation-evidence/browser/{Linux,Windows,macOS}/` | Browser-generated exports and checksums (bundled after tests) |

The bundled ZIP is created **only after** the full CI matrix passes; it is a release-candidate evidence archive, not a claim that all possible PDBs or proteins have been tested.

## Re-run from the published code

```bash
python -m pip install -r backend/requirements.txt
python -m unittest discover -s tests -v
node --test tests/test_bfactor_export.cjs
python validation/verify_publication.py
python validation/reproduce.py --suite full --output validation-evidence/numerical
```

Node.js 22 is required for B-factor-export reproduction, not ordinary GUI use. The numerical suite operates on committed PDB snapshots: **no new web download is part of numerical verification**. The full browser and first-launch workflows use GitHub Actions' Windows/macOS/Linux machines.

CI evidence is archived for 90 days. For long-term scholarly preservation, upload the accepted ZIP and tagged source snapshot to a durable repository (e.g. Zenodo), record the DOI and exact Git commit in the eventual manuscript.

**Release gate:** status is ACCEPTED only when `publication/build_package.py` confirms a PASS report for all 21 numerical inputs, all 3 operating-system browser exports, required files/checksums, and all upstream CI jobs. The PDB coverage is a representative, risk-based suite, not an exhaustive proof of accuracy for every biomolecule.
