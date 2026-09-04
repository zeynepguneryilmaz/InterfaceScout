"""Atom-charge screened electrostatic steering for charged planar interfaces.

This developmental module refines the residue-charge steering prototype using
PDB2PQR atomistic partial charges and atomic radii. It remains strictly separate
from InterfaceScout's canonical local compatibility and 5/8 Å patch scores.

For each sampled plane normal n, the plane is tangent to the PQR atomic envelope:

    h(n) = min_a [ n . r_a - R_a ]

and each charge center has depth

    z_a(n) = n . r_a - h(n) >= R_a.

For a homogeneous charged plane in the linearized PB/Debye-Huckel limit,
orientation ranking uses

    U_tilde(n;s) = s * sum_a q_a exp[-z_a(n)/lambda_D]

where q_a is the PDB2PQR partial atomic charge (e), s=+1 for a positive plane
and s=-1 for a negative plane. U_tilde = U/(e|psi0|), so no absolute adsorption
energy is claimed when |psi0| is unknown.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple
import numpy as np

from .steering import debye_length_A, fibonacci_sphere


def parse_pqr_atoms(pqr_text: str) -> List[Dict[str, Any]]:
    atoms: List[Dict[str, Any]] = []
    for line in pqr_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        tok = line.split()
        if len(tok) < 10:
            continue
        try:
            # Coordinates, charge, radius are reliably the final five fields.
            x, y, z, q, radius = map(float, tok[-5:])
            # Before xyz: record, serial, atom, residue, [chain], resseq.
            prefix = tok[:-5]
            atom_name = prefix[2]
            res_name = prefix[3]
            if len(prefix) >= 6:
                chain = prefix[4]
                res_seq = int(prefix[5])
            else:
                chain = ""
                res_seq = int(prefix[4])
            atoms.append({"atom_name": atom_name, "res_name": res_name,
                          "chain": chain, "res_seq": res_seq,
                          "coord": np.asarray([x,y,z],float), "charge_e": q, "radius_A": radius})
        except Exception:
            continue
    if not atoms:
        raise ValueError("No atom charges could be parsed from PQR")
    return atoms


def _orientation_energy(atoms: Sequence[Dict[str, Any]], surface_sign: int,
                        lambda_A: float, normals: np.ndarray):
    xyz=np.vstack([a["coord"] for a in atoms]); q=np.asarray([a["charge_e"] for a in atoms],float)
    rad=np.asarray([max(float(a["radius_A"]),0.0) for a in atoms],float)
    proj=xyz@normals.T
    plane=np.min(proj-rad[:,None],axis=0)
    z=np.maximum(proj-plane[None,:],0.0)
    energy=float(surface_sign)*np.sum(q[:,None]*np.exp(-z/lambda_A),axis=0)
    return energy,plane,proj,xyz,rad


def residue_footprint(atoms: Sequence[Dict[str, Any]], normal: np.ndarray, plane: float,
                      shell_A: float = 5.0) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str,int,str], Dict[str,Any]]={}
    for a in atoms:
        # Hydrogens participate in energy but not in geometric contact footprint.
        if str(a["atom_name"]).upper().startswith("H"):
            continue
        exposed_depth=float(np.dot(a["coord"],normal)-float(a["radius_A"])-plane)
        key=(str(a["chain"]),int(a["res_seq"]),str(a["res_name"]))
        rec=grouped.setdefault(key,{"chain":key[0],"res_seq":key[1],"res_name":key[2],"min_surface_depth_A":1e9})
        rec["min_surface_depth_A"]=min(rec["min_surface_depth_A"],exposed_depth)
    rows=[]
    for rec in grouped.values():
        if rec["min_surface_depth_A"]<=float(shell_A)+1e-12:
            rec=dict(rec); rec["min_surface_depth_A"]=round(float(rec["min_surface_depth_A"]),4); rows.append(rec)
    rows.sort(key=lambda r:(r["min_surface_depth_A"],r["chain"],r["res_seq"]))
    return rows


def compute_pqr_steering(pqr_text: str, ionic_mM: float, temperature_K: float,
                         n_orientations: int = 4096,
                         footprint_shells_A: Tuple[float,...]=(2.0,5.0,8.0)) -> Dict[str,Any]:
    atoms=parse_pqr_atoms(pqr_text)
    lam=debye_length_A(ionic_mM,temperature_K)
    if lam is None:
        return {"available":False,"reason":"ionic strength must be >0"}
    normals=fibonacci_sphere(n_orientations)
    result={"available":True,"kind":"atom_charge_electrostatic_steering",
            "model":"PDB2PQR atom-charge + atomic-envelope screened charged-plane ranking",
            "canonical_compatibility_scores_modified":False,
            "debye_length_A":round(float(lam),5),"n_atoms":len(atoms),
            "net_pqr_charge_e":round(float(sum(a["charge_e"] for a in atoms)),6),
            "energy_definition":"U_tilde=s*sum_a q_a exp(-z_a/lambda_D); lower is more favorable",
            "absolute_energy":"not reported because |psi0| is unspecified","surfaces":{}}
    for sign,label,linked in [(+1,"positive","cationic"),(-1,"negative","anionic")]:
        E,plane,proj,xyz,rad=_orientation_energy(atoms,sign,lam,normals)
        order=np.argsort(E); bi=int(order[0]); n=normals[bi]
        side={"surface_sign":sign,"linked_local_compatibility_map":linked,
              "best_normal_plane_to_protein":[round(float(x),6) for x in n],
              "best_protein_facing_direction_to_plane":[round(float(x),6) for x in -n],
              "best_reduced_energy":round(float(E[bi]),8),
              "worst_reduced_energy":round(float(np.max(E)),8),
              "reduced_energy_span":round(float(np.max(E)-np.min(E)),8),
              "sampled_orientations":int(len(normals)),"surface_facing_footprint":{}}
        for shell in footprint_shells_A:
            side["surface_facing_footprint"][f"within_{int(round(shell))}A"] = residue_footprint(atoms,n,float(plane[bi]),shell)
        result["surfaces"][label]=side
    result["scope"]={"predicts":"rigid-body electrostatic approach preference to a homogeneous charged plane",
                     "does_not_predict":["absolute adsorption free energy","charge regulation","nonlinear PB",
                                         "specific ion adsorption","hydration-layer penetration","final anchoring","conformational change"]}
    return result
