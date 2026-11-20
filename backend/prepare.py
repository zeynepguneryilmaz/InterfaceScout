"""Validated protein-only PDB preparation for the 11-map InterfaceScout engine.

Keep the original chain and author residue identifiers in exported results.
Select coherent residue-level alternate conformations by occupancy.  Never
silently include unknown heavy atoms in the SASA neighbourhood.
"""

from __future__ import annotations

from collections import defaultdict
from io import StringIO
from typing import Tuple

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

_BACKBONE = {"N", "CA", "C", "O"}


def _parse_chain_selection(chain: str | None) -> set[str] | None:
    text = (chain or "").strip()
    if not text:
        return None
    selected = {x.strip() for x in text.split(",") if x.strip()}
    return selected or None


def _residue_label(residue) -> str:
    chain = residue.get_parent()
    cid = str(chain.id) if chain is not None else "?"
    return f"{cid}:{int(residue.id[1])}:{str(residue.id[2]).strip()}"


def _is_standard_protein(residue) -> bool:
    return str(residue.id[0]).strip() == "" and is_aa(residue, standard=True)


def _heavy(atom) -> bool:
    element = (getattr(atom, "element", "") or "").strip().upper()
    name = atom.get_name().strip().upper()
    return element not in {"H", "D"} and not (
        not element and (name.startswith("H") or name.startswith("D")))


def _choose_altloc(residue) -> tuple[str | None, bool]:
    """Choose one occupancy-supported alternate *residue* conformation.

    Shared blank-altloc atoms are always retained. For genuinely ambiguous
    equal occupancies choose lexicographically but disclose the ambiguity.
    """
    groups: dict[str, float] = defaultdict(float)
    for atom in residue.get_unpacked_list():
        if not _heavy(atom):
            continue
        label = str(atom.get_altloc()).strip()
        if label:
            occupancy = atom.get_occupancy()
            groups[label] += max(0.0, float(occupancy)) if occupancy is not None else 0.0
    if not groups:
        return None, False
    order = sorted(groups.items(), key=lambda pair: (-pair[1], pair[0]))
    ambiguous = len(order) > 1 and abs(order[0][1] - order[1][1]) < 1e-9
    return order[0][0], ambiguous


class FirstModelHeavyAtomSelect(Select):
    """Keep model 0, selected chain(s), canonical residues, and one altloc."""

    def __init__(self, chain: str | None = None,
                 alternate_choices: dict[int, str | None] | None = None):
        super().__init__()
        self.chains = _parse_chain_selection(chain)
        self.alternate_choices = alternate_choices or {}

    def accept_model(self, model):
        return 1 if int(model.id) == 0 else 0

    def accept_chain(self, chain):
        if self.chains is None:
            return 1
        return 1 if str(chain.id) in self.chains else 0

    def accept_residue(self, residue):
        return 1 if _is_standard_protein(residue) else 0

    def accept_atom(self, atom):
        if not _heavy(atom):
            return 0
        altloc = str(atom.get_altloc()).strip()
        choice = self.alternate_choices.get(id(atom.get_parent()))
        return 1 if not altloc or altloc == choice else 0


def prepare_pdb_text(pdb_text: str, chain: str | None = None) -> Tuple[str, dict]:
    """Validate first-model protein coordinates and normalize PDB input.

    Hetero residues, waters, hydrogens, and unknown polymer residues are omitted.
    Such exclusions and missing atoms are reported, not interpreted as absence
    of those atoms in the original experimental structure.
    """
    if not isinstance(pdb_text, str) or not pdb_text.strip():
        raise ValueError("Empty protein-structure input")
    if pdb_text.lstrip().startswith("data_") and "_atom_site." in pdb_text:
        raise ValueError("Input is mmCIF, not PDB. Convert to PDB without "
                         "changing chain IDs/coordinates or use a PDB-format export. "
                         "mmCIF is not supported by this InterfaceScout reader yet.")
    if not any(line.startswith(("ATOM  ", "HETATM")) for line in pdb_text.splitlines()):
        raise ValueError("PDB input has no ATOM/HETATM coordinate records")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("interfacescout_input", StringIO(pdb_text))
    models = list(structure.get_models())
    if not models:
        raise ValueError("No structural model found in PDB input")

    first = models[0]
    available_chains = [str(c.id) for c in first]
    requested = _parse_chain_selection(chain)
    if requested:
        missing = sorted(requested - set(available_chains))
        if missing:
            raise ValueError(
                f"Chain(s) {', '.join(missing)} not found in first model. "
                f"Available: {', '.join(available_chains)}"
            )

    active = [c for c in first if requested is None or str(c.id) in requested]
    retained = [r for c in active for r in c if _is_standard_protein(r)]
    if not retained:
        raise ValueError("Selected chain(s) contain no standard protein residues")
    labels = [_residue_label(r) for r in retained]
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate protein residue chain/sequence/insertion identifiers")

    missing_ca = [r for r in retained if "CA" not in r]
    if missing_ca:
        raise ValueError(
            "Protein residues without C-alpha coordinates cannot be mapped: "
            + ", ".join(_residue_label(r) for r in missing_ca[:12])
            + (" ..." if len(missing_ca) > 12 else "")
        )

    choices: dict[int, str | None] = {}
    alternate_count = 0
    ambiguous_alternates: list[str] = []
    incomplete: list[str] = []
    for residue in retained:
        choice, ambiguous = _choose_altloc(residue)
        choices[id(residue)] = choice
        alternate_count += choice is not None
        if ambiguous:
            ambiguous_alternates.append(_residue_label(residue))
        available = {
            str(atom.get_name()).strip().upper()
            for atom in residue.get_unpacked_list()
            if _heavy(atom) and (not str(atom.get_altloc()).strip()
                                 or str(atom.get_altloc()).strip() == choice)
        }
        missing = _BACKBONE - available
        if missing:
            incomplete.append(f"{_residue_label(residue)} missing {','.join(sorted(missing))}")

    nonstandard = [
        _residue_label(r)
        for c in active for r in c if not _is_standard_protein(r)
    ]
    io = PDBIO()
    io.set_structure(structure)
    out = StringIO()
    io.save(out, FirstModelHeavyAtomSelect(chain, choices))
    normalized = out.getvalue()
    if not any(line.startswith("ATOM  ") for line in normalized.splitlines()):
        raise ValueError("Prepared structure contains no standard protein ATOM records")

    warnings: list[str] = []
    if len(models) > 1:
        warnings.append(f"Only the first of {len(models)} structural models was analyzed")
    if nonstandard:
        warnings.append(
            f"{len(nonstandard)} hetero/noncanonical residues were omitted; "
            "cofactors, glycans, modified residues, and their steric effects "
            "are not represented by the standard preparation policy")
    if incomplete:
        warnings.append(
            f"{len(incomplete)} protein residues have incomplete backbone atoms; "
            "side-chain SASA may be biased by other missing atoms")
    if alternate_count:
        warnings.append(
            f"Selected one residue-level alternate conformation for "
            f"{alternate_count} residues by heavy-atom occupancy")
    if ambiguous_alternates:
        warnings.append(
            f"{len(ambiguous_alternates)} alternate conformers had equal occupancy "
            "support; alphabetical altloc choice is deterministic but ambiguous")

    report = {
        "policy": "first_model_standard_protein_heavy_atoms_one_altloc_per_residue",
        "input_models": len(models),
        "selected_model": 1,
        "available_chains_first_model": available_chains,
        "selected_chain": ",".join(sorted(requested)) if requested else "ALL",
        "protein_residues_retained": len(retained),
        "alternate_residues_resolved": alternate_count,
        "ambiguous_alternate_residues": ambiguous_alternates,
        "incomplete_backbone_residues": incomplete,
        "excluded_nonstandard_residues": nonstandard[:100],
        "n_excluded_nonstandard_residues": len(nonstandard),
        "explicit_hydrogens_removed": True,
        "hetero_residues_removed": True,
        "warnings": warnings,
    }
    return normalized, report
