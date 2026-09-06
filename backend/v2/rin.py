"""Residue interaction network (RIN) context for InterfaceScout V2.

RIN is deliberately downstream of interface prediction. It does not influence
patch construction or ranking. The network asks a separate question: where does
a predicted or experimental material-contact region sit in the native protein's
structural interaction network?

Network definition
------------------
* node: one standard amino-acid residue in the first prepared structural model
* edge: at least one non-hydrogen atom pair from two residues within 4.5 A

The 4.5 A heavy-atom contact definition follows common contact-based RIN practice
and is not fitted to adsorption-interface labels. Centralities are descriptive;
no binary 'hub' threshold is imposed. Percentiles are reported instead.
"""

from __future__ import annotations

from collections import deque
from io import StringIO
from typing import Dict, Iterable, List, Set

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

RIN_HEAVY_ATOM_CUTOFF_A = 4.5


def _residue_key(chain_id: str, residue) -> str:
    return f"{chain_id}:{int(residue.id[1])}:{str(residue.id[2]).strip()}"


def _extract_residues(pdb_text: str) -> List[dict]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rin", StringIO(pdb_text))
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("No structural model found for RIN")

    rows: List[dict] = []
    for chain in model:
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            atoms = []
            for atom in residue.get_atoms():
                element = (getattr(atom, "element", "") or "").strip().upper()
                name = atom.get_name().strip().upper()
                if element == "H" or name.startswith("H"):
                    continue
                atoms.append(np.asarray(atom.coord, dtype=float))
            if not atoms:
                continue
            arr = np.vstack(atoms)
            center = arr.mean(axis=0)
            radius = float(np.max(np.linalg.norm(arr - center, axis=1)))
            rows.append({
                "key": _residue_key(str(chain.id), residue),
                "chain": str(chain.id),
                "res_seq": int(residue.id[1]),
                "icode": str(residue.id[2]).strip(),
                "res_name": str(residue.resname).strip().upper(),
                "atoms": arr,
                "center": center,
                "radius": radius,
            })
    if len(rows) < 2:
        raise ValueError("Too few residues for RIN")
    return rows


def _count_contacts(a: np.ndarray, b: np.ndarray, cutoff_A: float) -> int:
    delta = a[:, None, :] - b[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", delta, delta)
    return int(np.count_nonzero(d2 <= float(cutoff_A) ** 2))


def build_rin(pdb_text: str, cutoff_A: float = RIN_HEAVY_ATOM_CUTOFF_A) -> dict:
    """Build an unweighted heavy-atom-contact RIN plus contact-count metadata."""
    residues = _extract_residues(pdb_text)
    n = len(residues)
    adjacency: List[Set[int]] = [set() for _ in range(n)]
    contact_counts: Dict[str, int] = {}

    for i in range(n):
        ri = residues[i]
        for j in range(i + 1, n):
            rj = residues[j]
            center_dist = float(np.linalg.norm(ri["center"] - rj["center"]))
            if center_dist > ri["radius"] + rj["radius"] + float(cutoff_A):
                continue
            contacts = _count_contacts(ri["atoms"], rj["atoms"], float(cutoff_A))
            if contacts <= 0:
                continue
            adjacency[i].add(j)
            adjacency[j].add(i)
            contact_counts[f"{ri['key']}|{rj['key']}"] = contacts

    index = {r["key"]: i for i, r in enumerate(residues)}
    degree = np.asarray([len(adjacency[i]) for i in range(n)], dtype=float)
    closeness = _closeness(adjacency)
    betweenness = _brandes_betweenness(adjacency)

    metrics: Dict[str, dict] = {}
    for i, r in enumerate(residues):
        metrics[r["key"]] = {
            "key": r["key"],
            "chain": r["chain"],
            "res_seq": r["res_seq"],
            "icode": r["icode"],
            "res_name": r["res_name"],
            "degree": int(degree[i]),
            "degree_normalized": float(degree[i] / (n - 1)) if n > 1 else 0.0,
            "closeness": float(closeness[i]),
            "betweenness": float(betweenness[i]),
        }

    return {
        "cutoff_A": float(cutoff_A),
        "n_nodes": n,
        "n_edges": int(sum(len(x) for x in adjacency) // 2),
        "residues": residues,
        "index": index,
        "adjacency": adjacency,
        "contact_counts": contact_counts,
        "residue_metrics": metrics,
    }


def _single_source_distances(adjacency: List[Set[int]], source: int) -> np.ndarray:
    n = len(adjacency)
    dist = np.full(n, -1, dtype=int)
    dist[source] = 0
    q = deque([source])
    while q:
        v = q.popleft()
        for w in adjacency[v]:
            if dist[w] < 0:
                dist[w] = dist[v] + 1
                q.append(w)
    return dist


def _closeness(adjacency: List[Set[int]]) -> np.ndarray:
    n = len(adjacency)
    out = np.zeros(n, dtype=float)
    for s in range(n):
        dist = _single_source_distances(adjacency, s)
        reachable = dist >= 0
        n_reach = int(reachable.sum()) - 1
        total = int(dist[(dist > 0)].sum())
        if n_reach <= 0 or total <= 0:
            continue
        # Wasserman-Faust correction for disconnected graphs.
        raw = n_reach / total
        out[s] = raw * (n_reach / (n - 1)) if n > 1 else 0.0
    return out


def _brandes_betweenness(adjacency: List[Set[int]]) -> np.ndarray:
    """Normalized Brandes betweenness for an undirected, unweighted graph."""
    n = len(adjacency)
    cb = np.zeros(n, dtype=float)
    for s in range(n):
        stack: List[int] = []
        pred: List[List[int]] = [[] for _ in range(n)]
        sigma = np.zeros(n, dtype=float)
        sigma[s] = 1.0
        dist = np.full(n, -1, dtype=int)
        dist[s] = 0
        q = deque([s])
        while q:
            v = q.popleft()
            stack.append(v)
            for w in adjacency[v]:
                if dist[w] < 0:
                    q.append(w)
                    dist[w] = dist[v] + 1
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = np.zeros(n, dtype=float)
        while stack:
            w = stack.pop()
            if sigma[w] > 0:
                coeff = (1.0 + delta[w]) / sigma[w]
                for v in pred[w]:
                    delta[v] += sigma[v] * coeff
            if w != s:
                cb[w] += delta[w]

    cb *= 0.5  # undirected paths were counted from both ends
    if n > 2:
        cb *= 2.0 / ((n - 1) * (n - 2))
    return cb


def _percentile(value: float, reference: List[float]) -> float | None:
    if not reference:
        return None
    arr = np.asarray(reference, dtype=float)
    # Mid-rank empirical percentile; no arbitrary hub threshold.
    less = float(np.count_nonzero(arr < value))
    equal = float(np.count_nonzero(arr == value))
    return 100.0 * (less + 0.5 * equal) / len(arr)


def annotate_rin_percentiles(rin: dict, surface_keys: Iterable[str]) -> dict:
    """Attach all-residue and exposed-surface percentiles to RIN metrics."""
    metrics = rin["residue_metrics"]
    surface = [k for k in surface_keys if k in metrics]
    fields = ("degree_normalized", "closeness", "betweenness")
    all_refs = {f: [float(m[f]) for m in metrics.values()] for f in fields}
    surf_refs = {f: [float(metrics[k][f]) for k in surface] for f in fields}

    for key, m in metrics.items():
        for field in fields:
            m[f"{field}_percentile_all"] = _percentile(float(m[field]), all_refs[field])
            m[f"{field}_percentile_surface"] = _percentile(float(m[field]), surf_refs[field])
    return rin


def summarize_patch_rin(keys: Iterable[str], rin: dict) -> dict:
    """Summarize RIN context for a predicted or experimental patch."""
    metrics = rin["residue_metrics"]
    valid = [k for k in keys if k in metrics]
    if not valid:
        return {
            "n_rin_residues": 0,
            "residue_metrics": [],
        }

    def values(field: str) -> List[float]:
        return [float(metrics[k][field]) for k in valid if metrics[k].get(field) is not None]

    idx = rin["index"]
    adjacency = rin["adjacency"]
    vset = {idx[k] for k in valid}
    internal_edges = 0
    boundary_edges = 0
    for i in vset:
        for j in adjacency[i]:
            if j in vset:
                if i < j:
                    internal_edges += 1
            else:
                boundary_edges += 1

    possible_internal = len(vset) * (len(vset) - 1) / 2.0
    internal_density = internal_edges / possible_internal if possible_internal > 0 else 0.0
    total_incident = 2 * internal_edges + boundary_edges
    boundary_fraction = boundary_edges / total_incident if total_incident > 0 else 0.0

    summary = {
        "n_rin_residues": len(valid),
        "mean_degree_percentile_surface": float(np.mean(values("degree_normalized_percentile_surface"))) if values("degree_normalized_percentile_surface") else None,
        "max_degree_percentile_surface": float(np.max(values("degree_normalized_percentile_surface"))) if values("degree_normalized_percentile_surface") else None,
        "mean_betweenness_percentile_surface": float(np.mean(values("betweenness_percentile_surface"))) if values("betweenness_percentile_surface") else None,
        "max_betweenness_percentile_surface": float(np.max(values("betweenness_percentile_surface"))) if values("betweenness_percentile_surface") else None,
        "mean_closeness_percentile_surface": float(np.mean(values("closeness_percentile_surface"))) if values("closeness_percentile_surface") else None,
        "max_closeness_percentile_surface": float(np.max(values("closeness_percentile_surface"))) if values("closeness_percentile_surface") else None,
        "internal_edge_density": float(internal_density),
        "boundary_edge_fraction": float(boundary_fraction),
        "residue_metrics": [metrics[k] for k in valid],
    }
    return summary
