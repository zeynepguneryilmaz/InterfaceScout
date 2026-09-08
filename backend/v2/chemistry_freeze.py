"""Publication-freeze chemistry corrections applied by InterfaceScout V2.

The legacy V1 table is retained for backward compatibility of the original
/analyze_surface endpoint. V2 applies the conservative publication corrections
below before any chemistry map is generated. These corrections were defined
before the frozen external benchmark is interpreted.

The publication engine intentionally uses binary mechanistic eligibility rather
than heterogeneous literature interaction-energy magnitudes. Consequently,
channels are restricted to direct/primary interaction motifs that can be
represented without introducing an empirical weighting scheme.
"""
from __future__ import annotations


def _replace_state(entry, state):
    """Preserve mechanism metadata while replacing the required ionization state."""
    _, mechanism, _ = entry
    return (0.0, mechanism, state)


def _neutralize_ebase(channel):
    """Remove numerical energy magnitude from V2 while preserving mechanism/state."""
    for group in ("favorable", "repulsive"):
        rows = channel.get(group, {})
        for residue, entry in list(rows.items()):
            _, mechanism, state = entry
            rows[residue] = (0.0, mechanism, state)


def apply_publication_chemistry(v1_module) -> None:
    chem = v1_module.CHEMISTRIES

    # H-bond donor surface: the protein-facing group must be acceptor-capable.
    # Protonated Lys ammonium and Arg guanidinium are H-bond donors rather than
    # conventional acceptors, so they are excluded from this channel.
    donor_surface = chem["hbond_donor"]["favorable"]
    donor_surface.pop("LYS", None)
    donor_surface.pop("ARG", None)
    chem["hbond_donor"]["description"] = (
        "Protein-side compatibility with surfaces capable of donating hydrogen bonds; "
        "favorable residues are restricted to acceptor-capable side chains."
    )

    # H-bond acceptor surface: Lys/Arg contribute through their protonated
    # donor states. Explicit state availability therefore enters these rows.
    acceptor_surface = chem["hbond_acceptor"]["favorable"]
    if "LYS" in acceptor_surface:
        acceptor_surface["LYS"] = _replace_state(acceptor_surface["LYS"], "protonated")
    if "ARG" in acceptor_surface:
        acceptor_surface["ARG"] = _replace_state(acceptor_surface["ARG"], "protonated")
    chem["hbond_acceptor"]["description"] = (
        "Protein-side compatibility with surfaces capable of accepting hydrogen bonds; "
        "basic side-chain donors are conditioned on their protonated state."
    )

    # Oxide channel: because V2 does not fit or apply interaction-strength
    # weights, retain the direct carboxylate/metal-site coordination prior and
    # exclude weaker neutral hydroxyl contacts from the binary primary map.
    oxide = chem["oxide"]["favorable"]
    for residue in ("SER", "THR", "TYR"):
        oxide.pop(residue, None)
    chem["oxide"]["description"] = (
        "Primary protein-side compatibility with exposed oxide metal sites through "
        "deprotonated Asp/Glu carboxylate coordination."
    )

    # Hydroxyapatite channel: charged residues show the clearest experimentally
    # supported affinity contrast. Weak neutral contacts are not treated as
    # binary-equivalent to direct charged-site interactions in the unweighted
    # publication model.
    hap = chem["hydroxyapatite"]["favorable"]
    for residue in ("SER", "THR", "TYR", "PHE", "TRP"):
        hap.pop(residue, None)
    hap["ARG"] = (0.0, "charged-site electrostatic complementarity", "protonated")
    hap["LYS"] = (0.0, "charged-site electrostatic complementarity", "protonated")
    hap["HIS"] = (0.0, "charged-site electrostatic complementarity", "protonated")
    chem["hydroxyapatite"]["description"] = (
        "Primary charged-site compatibility with hydroxyapatite/calcium-phosphate "
        "interfaces: deprotonated Asp/Glu carboxylates and protonated basic residues."
    )

    # The publication V2 model retains literature-derived mechanism/source
    # provenance but deliberately removes heterogeneous numerical Ebase
    # magnitudes. Ranking depends only on eligibility, scRSA, state availability,
    # and the frozen spatial architecture.
    for channel in chem.values():
        _neutralize_ebase(channel)
