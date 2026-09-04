"""Atom-level nonpolar surface orientation layer for InterfaceScout v5.3.

Based on the hydrophobic-surface Monte Carlo framework of Harrison, Weidner and
Castner (Biointerphases 12, 02D401, 2017; DOI 10.1116/1.4971381), which combines
rigid-body orientation sampling with a SASA-based solvation term.  This module
implements the *solvation/orientation component only* as a lightweight planar
screening descriptor; it does not claim to reproduce their CHARMM22 van der
Waals term or a full Metropolis trajectory.

The published model assigns equal-magnitude opposite-sign surface-tension terms
to hydrophobic and hydrophilic atom groups.  For ranking only, the common
magnitude is normalized to +1 (hydrophobic) / -1 (hydrophilic), which leaves the
orientation ordering of the solvation component unchanged.

For a candidate planar approach normal n, the protein is translated until the
lowest heavy-atom center is tangent to the plane.  Atoms within 6 A of the plane
(the contact distance used for residue-contact analysis in the published study)
form the candidate footprint.  If an exposed hydrophobic atom becomes occluded,
the solvation contribution is favorable; occluding an exposed hydrophilic atom
is unfavorable.  The screening score is therefore

    S(n) = sum_contact s_a * SASA_a

with s_a=+1 hydrophobic and -1 hydrophilic. Larger S means a more favorable
nonpolar-contact solvation component.  No fitted weights are used.
"""
from __future__ import annotations
from typing import Any, Dict, List
import numpy as np
from Bio.PDB.Polypeptide import is_aa

CONTACT_A = 6.0
HYDROPHOBIC_RES = {"GLY","ALA","VAL","LEU","ILE","MET","PRO","PHE","TRP"}
HYDROPHILIC_RES = {"SER","THR","ASN","GLN","CYS","ARG","ASP","HIS","LYS","GLU"}
TYR_RING = {"CB","CG","CD1","CD2","CE1","CE2","CZ"}


def fibonacci_sphere(n: int = 2048) -> np.ndarray:
    i=np.arange(n,dtype=float)
    phi=(1.0+5.0**0.5)/2.0
    z=1.0-2.0*(i+0.5)/n
    r=np.sqrt(np.maximum(0.0,1.0-z*z))
    theta=2.0*np.pi*i/phi
    return np.column_stack([r*np.cos(theta),r*np.sin(theta),z])


def atom_group(resname: str, atom_name: str) -> str:
    rn=resname.upper(); an=atom_name.upper()
    # Heavy-atom representation of the published backbone grouping:
    # C-alpha and carbonyl C are nonpolar; amide N and carbonyl O are polar.
    if an in {"CA","C"}: return "hydrophobic_backbone"
    if an in {"N","O","OXT"}: return "hydrophilic_backbone"
    if rn == "TYR":
        if an in TYR_RING: return "hydrophobic_sidechain"
        return "hydrophilic_sidechain"
    if rn in HYDROPHOBIC_RES: return "hydrophobic_sidechain"
    return "hydrophilic_sidechain"


def extract_atoms(struct) -> List[Dict[str, Any]]:
    rows=[]
    model=next(struct.get_models())
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=True): continue
            rn=res.get_resname().strip().upper()
            for atom in res.get_atoms():
                if str(getattr(atom,"element","")).upper()=="H": continue
                sasa=float(getattr(atom,"sasa",0.0) or 0.0)
                grp=atom_group(rn,atom.get_name().strip())
                sign=1.0 if grp.startswith("hydrophobic") else -1.0
                rows.append({
                    "coord":np.asarray(atom.coord,dtype=float),"sasa":sasa,"sign":sign,
                    "group":grp,"chain":str(chain.id),"res_seq":int(res.id[1]),
                    "icode":str(res.id[2]).strip(),"res_name":rn,"atom":atom.get_name().strip(),
                })
    return rows


def scan(struct, n_orientations: int = 2048, contact_A: float = CONTACT_A) -> Dict[str, Any]:
    atoms=extract_atoms(struct)
    if not atoms: return {"status":"unavailable"}
    xyz=np.vstack([a["coord"] for a in atoms])
    sasa=np.asarray([a["sasa"] for a in atoms],float)
    sign=np.asarray([a["sign"] for a in atoms],float)
    normals=fibonacci_sphere(n_orientations)
    scores=np.empty(n_orientations,float)
    contact_counts=np.empty(n_orientations,int)
    for k,n in enumerate(normals):
        p=xyz@n
        depth=p-p.min()
        mask=depth<=contact_A
        scores[k]=float(np.sum(sign[mask]*sasa[mask]))
        contact_counts[k]=int(mask.sum())
    order=np.argsort(-scores,kind="mergesort")
    top=[]
    for rank,idx in enumerate(order[:20],1):
        n=normals[idx]; p=xyz@n; depth=p-p.min(); mask=depth<=contact_A
        residues={}
        for j in np.where(mask)[0]:
            a=atoms[j]; key=f"{a['chain']}:{a['res_seq']}:{a['icode']}"
            rr=residues.setdefault(key,{"key":key,"chain":a['chain'],"res_seq":a['res_seq'],"icode":a['icode'],"res_name":a['res_name'],"hydrophobic_sasa":0.0,"hydrophilic_sasa":0.0})
            if a['sign']>0: rr['hydrophobic_sasa']+=a['sasa']
            else: rr['hydrophilic_sasa']+=a['sasa']
        rlist=[]
        for rr in residues.values():
            rr['hydrophobic_sasa']=round(rr['hydrophobic_sasa'],3); rr['hydrophilic_sasa']=round(rr['hydrophilic_sasa'],3)
            rr['net_contact_sasa']=round(rr['hydrophobic_sasa']-rr['hydrophilic_sasa'],3); rlist.append(rr)
        rlist.sort(key=lambda x:(-x['net_contact_sasa'],x['chain'],x['res_seq']))
        top.append({
            "rank":rank,"orientation_index":int(idx),"normal":[round(float(x),6) for x in n],
            "solvation_contact_score":round(float(scores[idx]),5),"n_contact_atoms":int(contact_counts[idx]),
            "contact_residues":rlist,
        })
    return {
        "status":"ok","method":"Harrison-style SASA solvation component; planar rigid-body orientation scan",
        "n_orientations":n_orientations,"contact_distance_A":contact_A,
        "score_definition":"sum(contact hydrophobic atom SASA) - sum(contact hydrophilic atom SASA); larger is more favorable",
        "full_MC_reproduced":False,"vdw_term_included":False,"fitted_weights":False,
        "best_score":round(float(scores[order[0]]),5),"median_score":round(float(np.median(scores)),5),
        "top_orientations":top,
    }
