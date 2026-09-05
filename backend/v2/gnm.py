"""Minimal, transparent Gaussian Network Model (GNM) utilities.

This module is intentionally independent of the V1 scoring model. It estimates
native-state C-alpha fluctuations and residue-residue correlations from the
prepared protein structure. It does not predict adsorption energy.
"""

from __future__ import annotations

from io import StringIO
from typing import Dict, List, Tuple

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


def _residue_key(chain_id: str, residue) -> str:
    seq = int(residue.id[1])
    icode = str(residue.id[2]).strip()
    return f"{chain_id}:{seq}:{icode}"


def extract_ca_nodes(pdb_text: str) -> List[dict]:
    """Extract standard-amino-acid C-alpha nodes from the first model."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("gnm", StringIO(pdb_text))
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("No model found for GNM")

    nodes: List[dict] = []
    for chain in model:
        for residue in chain:
            if not is_aa(residue, standard=True) or "CA" not in residue:
                continue
            nodes.append(
                {
                    "key": _residue_key(str(chain.id), residue),
                    "chain": str(chain.id),
                    "res_seq": int(residue.id[1]),
                    "icode": str(residue.id[2]).strip(),
                    "res_name": str(residue.resname).strip(),
                    "coord": np.asarray(residue["CA"].coord, dtype=float),
                }
            )
    if len(nodes) < 3:
        raise ValueError("Too few C-alpha nodes for GNM")
    return nodes


def build_kirchhoff(nodes: List[dict], cutoff_A: float = 7.3) -> Tuple[np.ndarray, np.ndarray]:
    """Build the standard unweighted GNM Kirchhoff matrix."""
    coords = np.vstack([n["coord"] for n in nodes])
    delta = coords[:, None, :] - coords[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", delta, delta)
    adjacency = (d2 <= float(cutoff_A) ** 2) & (d2 > 0.0)

    gamma = np.zeros((len(nodes), len(nodes)), dtype=float)
    gamma[adjacency] = -1.0
    np.fill_diagonal(gamma, -gamma.sum(axis=1))
    return gamma, adjacency.astype(int)


def solve_gnm(pdb_text: str, cutoff_A: float = 7.3, zero_tol: float = 1e-8) -> dict:
    """Solve a C-alpha GNM and return normalized fluctuations/correlations.

    All non-zero modes are retained. Reported fluctuations are normalized by
    their mean, so values >1 indicate above-average native-state mobility.
    """
    nodes = extract_ca_nodes(pdb_text)
    gamma, adjacency = build_kirchhoff(nodes, cutoff_A=cutoff_A)

    eigvals, eigvecs = np.linalg.eigh(gamma)
    keep = eigvals > float(zero_tol)
    if not np.any(keep):
        raise ValueError("GNM Kirchhoff matrix has no non-zero modes")

    inv = (eigvecs[:, keep] / eigvals[keep]) @ eigvecs[:, keep].T
    diag = np.clip(np.diag(inv), 0.0, None)
    mean_diag = float(diag.mean()) if len(diag) else 1.0
    norm_fluct = diag / mean_diag if mean_diag > 0 else diag

    denom = np.sqrt(np.outer(diag, diag))
    corr = np.divide(inv, denom, out=np.zeros_like(inv), where=denom > 0)
    np.fill_diagonal(corr, 1.0)

    degree = adjacency.sum(axis=1).astype(int)
    index = {n["key"]: i for i, n in enumerate(nodes)}

    residue_metrics: Dict[str, dict] = {}
    for i, n in enumerate(nodes):
        residue_metrics[n["key"]] = {
            "key": n["key"],
            "chain": n["chain"],
            "res_seq": n["res_seq"],
            "icode": n["icode"],
            "res_name": n["res_name"],
            "normalized_fluctuation": float(norm_fluct[i]),
            "contact_degree": int(degree[i]),
        }

    return {
        "cutoff_A": float(cutoff_A),
        "n_nodes": len(nodes),
        "n_zero_modes": int((~keep).sum()),
        "n_nonzero_modes": int(keep.sum()),
        "nodes": nodes,
        "index": index,
        "correlation_matrix": corr,
        "adjacency": adjacency,
        "residue_metrics": residue_metrics,
    }
