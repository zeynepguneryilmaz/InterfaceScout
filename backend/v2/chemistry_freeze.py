"""Publication-freeze chemistry corrections applied by InterfaceScout V2.

The legacy V1 table is retained for backward compatibility of the original
/analyze_surface endpoint. V2 applies the corrections below before any V1
chemistry map is generated. This keeps the publication/validation engine
chemically explicit without silently changing legacy output.
"""
from __future__ import annotations


def apply_publication_chemistry(v1_module) -> None:
    chem = v1_module.CHEMISTRIES

    # A surface labelled H-bond donor is interpreted from the protein side:
    # favorable residues must offer an acceptor-capable side-chain atom.
    # Protonated Lys ammonium and Arg guanidinium are donors rather than
    # conventional H-bond acceptors and are therefore excluded.
    donor_surface = chem["hbond_donor"]["favorable"]
    donor_surface.pop("LYS", None)
    donor_surface.pop("ARG", None)
    chem["hbond_donor"]["description"] = (
        "Protein-side compatibility with surfaces capable of donating hydrogen bonds; "
        "favorable residues are restricted to acceptor-capable side chains."
    )

    # Numerical Ebase values in the inherited tuples remain legacy provenance
    # metadata only. V2 never uses their magnitude in residue or patch ranking.
