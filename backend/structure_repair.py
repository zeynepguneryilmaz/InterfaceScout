"""Minimal force-field preparation for optional InterfaceScout 2.0 physics.

Experimental PDB coordinates are often missing terminal or unresolved heavy
atoms.  The residue-level InterfaceScout core deliberately works from the
observed coordinates, but atomistic force-field parameterization requires a
chemically complete standard-residue topology.

This helper uses PDBFixer only for the optional atomistic physics layer:
- missing whole residues are NOT built;
- non-standard chemistry is NOT silently replaced;
- heterogens/waters are removed because the current planar nonpolar model is a
  protein-only continuum-interface calculation;
- missing heavy atoms and terminal heavy atoms within already present standard
  residues are added deterministically by PDBFixer;
- hydrogens are left for the downstream CHARMM/OpenMM pH-aware step.

The repaired coordinates never replace the structure used by the canonical
InterfaceScout 1.0 scRSA/chemistry/persistence calculation.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Any

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
except Exception:
    PDBFixer = None
    PDBFile = None


def repair_for_forcefield(source: Path, target: Path) -> Dict[str, Any]:
    if PDBFixer is None or PDBFile is None:
        return {"status": "unavailable", "reason": "PDBFixer is not installed"}

    fixer = PDBFixer(filename=str(source))

    # Detect missing residues for transparency, then explicitly suppress their
    # construction: model-generated loops would change the experimental
    # structure and are outside this lightweight screening layer.
    fixer.findMissingResidues()
    n_missing_residues = int(sum(len(v) for v in fixer.missingResidues.values()))
    fixer.missingResidues = {}

    # Current atomistic surface model is protein-only.  structural_context.py
    # has already retained only standard AAs, but this is an additional guard.
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    n_missing_atoms = int(sum(len(v) for v in fixer.missingAtoms.values()))
    n_missing_terminal_atoms = int(sum(len(v) for v in fixer.missingTerminals.values()))
    fixer.addMissingAtoms()

    with target.open("w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    return {
        "status": "ok",
        "method": "PDBFixer heavy-atom/terminal repair for optional atomistic physics",
        "missing_whole_residues_detected_but_not_built": n_missing_residues,
        "missing_heavy_atoms_added": n_missing_atoms,
        "missing_terminal_heavy_atoms_added": n_missing_terminal_atoms,
        "hydrogens_added_here": False,
        "canonical_core_coordinates_changed": False,
        "output": str(target),
    }
