"""Structure-preparation helpers for InterfaceScout V2."""

from __future__ import annotations

from io import StringIO
from typing import Tuple

from Bio.PDB import PDBParser, PDBIO, Select


def _parse_chain_selection(chain: str | None) -> set[str] | None:
    text = (chain or "").strip()
    if not text:
        return None
    selected = {x.strip() for x in text.split(",") if x.strip()}
    return selected or None


class FirstModelHeavyAtomSelect(Select):
    """Keep model 0, selected chain(s), standard amino acids, and heavy atoms."""

    def __init__(self, chain: str | None = None):
        super().__init__()
        self.chains = _parse_chain_selection(chain)

    def accept_model(self, model):
        return 1 if int(model.id) == 0 else 0

    def accept_chain(self, chain):
        if self.chains is None:
            return 1
        return 1 if str(chain.id) in self.chains else 0

    def accept_residue(self, residue):
        hetflag = str(residue.id[0]).strip()
        return 1 if hetflag == "" else 0

    def accept_atom(self, atom):
        element = (getattr(atom, "element", "") or "").strip().upper()
        name = atom.get_name().strip().upper()
        return 0 if element == "H" or name.startswith("H") else 1


def prepare_pdb_text(pdb_text: str, chain: str | None = None) -> Tuple[str, dict]:
    """Normalize an input PDB for V2.

    Policy:
    - first structural model only;
    - one or more selected chains if supplied (e.g. ``A`` or ``C,E``);
    - remove explicit hydrogens;
    - remove hetero residues from the current coarse protein patch engine.

    Hetero removal is a structural-preparation simplification, not a claim that
    cofactors are irrelevant to adsorption.  Systems whose native interface
    depends directly on a retained cofactor require a dedicated future policy.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("v2_input", StringIO(pdb_text))
    models = list(structure.get_models())
    if not models:
        raise ValueError("No structural model found in PDB input")

    available_chains = [str(c.id) for c in models[0].get_chains()]
    requested = _parse_chain_selection(chain)
    if requested:
        missing = sorted(requested - set(available_chains))
        if missing:
            raise ValueError(
                f"Chain(s) {', '.join(missing)} not found in first model. "
                f"Available: {', '.join(available_chains)}"
            )

    io = PDBIO()
    io.set_structure(structure)
    out = StringIO()
    io.save(out, FirstModelHeavyAtomSelect(chain))
    normalized = out.getvalue()
    if "ATOM" not in normalized:
        raise ValueError("Prepared structure contains no protein ATOM records")

    selected_label = ",".join(sorted(requested)) if requested else "ALL"
    report = {
        "policy": "first_model_heavy_atom_selected_chains",
        "input_models": len(models),
        "selected_model": 1,
        "available_chains_first_model": available_chains,
        "selected_chain": selected_label,
        "explicit_hydrogens_removed": True,
        "hetero_residues_removed": True,
    }
    return normalized, report
