"""Structure-preparation helpers for InterfaceScout V2-alpha."""

from __future__ import annotations

from io import StringIO
from typing import Tuple

from Bio.PDB import PDBParser, PDBIO, Select


class FirstModelHeavyAtomSelect(Select):
    """Keep model 0, standard amino acids, and non-hydrogen atoms only."""

    def __init__(self, chain: str | None = None):
        super().__init__()
        self.chain = (chain or "").strip() or None

    def accept_model(self, model):
        return 1 if int(model.id) == 0 else 0

    def accept_chain(self, chain):
        if self.chain is None:
            return 1
        return 1 if str(chain.id) == self.chain else 0

    def accept_residue(self, residue):
        hetflag = str(residue.id[0]).strip()
        return 1 if hetflag == "" else 0

    def accept_atom(self, atom):
        element = (getattr(atom, "element", "") or "").strip().upper()
        name = atom.get_name().strip().upper()
        return 0 if element == "H" or name.startswith("H") else 1


def prepare_pdb_text(pdb_text: str, chain: str | None = None) -> Tuple[str, dict]:
    """Normalize an input PDB for V2-alpha.

    Current policy:
    - first structural model only;
    - selected chain if supplied;
    - remove explicit hydrogens;
    - remove hetero residues for the alpha patch engine.

    Returns normalized PDB text and a compact preparation report.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("v2_input", StringIO(pdb_text))
    models = list(structure.get_models())
    if not models:
        raise ValueError("No structural model found in PDB input")

    available_chains = [str(c.id) for c in models[0].get_chains()]
    requested_chain = (chain or "").strip()
    if requested_chain and requested_chain not in available_chains:
        raise ValueError(
            f"Chain {requested_chain!r} not found in first model. "
            f"Available: {', '.join(available_chains)}"
        )

    io = PDBIO()
    io.set_structure(structure)
    out = StringIO()
    io.save(out, FirstModelHeavyAtomSelect(requested_chain or None))
    normalized = out.getvalue()
    if "ATOM" not in normalized:
        raise ValueError("Prepared structure contains no protein ATOM records")

    report = {
        "policy": "first_model_heavy_atom",
        "input_models": len(models),
        "selected_model": 1,
        "available_chains_first_model": available_chains,
        "selected_chain": requested_chain or "ALL",
        "explicit_hydrogens_removed": True,
        "hetero_residues_removed": True,
    }
    return normalized, report
