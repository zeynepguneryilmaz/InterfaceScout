"""Opt-in electrostatics-aware re-evaluation of the canonical 11 chemistry channels.

No affinity prediction: only discrete PDB2PQR assigned ionization and the
fraction of solvent-side APBS sample points with the expected potential sign
enter charge-mediated mechanisms. Nonionic mechanisms are unchanged.
"""
from __future__ import annotations
from collections import defaultdict
import core
from coarse_patch import build_coarse_patches
from geometry import extract_ca_nodes
from interface_engine import PUBLIC_MAP_SPECS, _residue_rows
from physics_surface import read_pqr,read_dx,charge_by_residue,sample_exposed_potentials

MODEL_VERSION="2.5.0-physics-context-prototype"

def build_physical_11(*,prepared_protein_pdb,pqr_text,dx_text,pH,temperature_K,
                      baseline_analysis=None):
    atoms=read_pqr(pqr_text)
    charges=charge_by_residue(atoms)
    raw=core._surface_residues(prepared_protein_pdb,pH)
    surface=[r for r in raw if r["surface_exposed"]]
    if not surface:raise ValueError("No surface residues in prepared structure")
    keys={r["key"] for r in surface}
    if dx_text is None:raise ValueError("APBS DX required for physical 11-map model")
    potential=sample_exposed_potentials(atoms,read_dx(dx_text),temperature_K,keys)
    # The biological assembly must match the PQR including chain identifiers.
    if len(keys&charges)<len(keys)*.95:
        raise ValueError(f"PQR residue mapping incomplete: {len(keys&charges)}/{len(keys)}")
    physical={}
    for key in keys:
        if key not in charges:raise ValueError("Missing PQR charge at "+key)
        phi=potential.get(key)
        physical[key]={"assigned_charge_e":charges[key],
                      "fraction_positive":phi.get("fraction_positive") if phi else None,
                      "potential_mean_mV":phi.get("potential_mean_mV") if phi else None}
    distances=core.build_distances(surface)
    parent={"chemistries":{},"surface_residues":surface}
    for public,internal,_,_ in PUBLIC_MAP_SPECS:
        parent["chemistries"][internal]=core.chemistry_map(surface,distances,internal,pH,
                                                          physical_evidence=physical)
    nodes=extract_ca_nodes(prepared_protein_pdb)
    records={}
    for public,internal,label,_ in PUBLIC_MAP_SPECS:
        channel=parent["chemistries"][internal]
        patches=build_coarse_patches(v1_result=parent,chemistry=internal,ca_nodes=nodes)
        records[public]={"key":public,"label":label,"description":next(spec[3] for spec in PUBLIC_MAP_SPECS if spec[0]==public),
                         "internal_chemistry_key":internal,"n_patches":len(patches),
                         "n_pareto_primary_patches":sum(int(p["pareto_front"])==1 for p in patches),
                         "primary_patches":[p for p in patches if int(p["pareto_front"])==1],
                         "patches":patches,"residues":_residue_rows(channel),
                         "n_chemistry_residues":len(channel["residues"])}
    return {
        "model":MODEL_VERSION,
        "note":"Physically informed candidate-region ranks, not adsorption energies or experimentally validated interfaces.",
        "physics":{"surface_residues":len(keys),"potential_sampled_residues":len(potential),
                   "charged_residues_mapped":len(keys & charges),
                   "external_protonation":"PDB2PQR with PROPKA",
                   "external_potential":"APBS linearized PB",
                   "electrostatic_gate":"Fraction of probe-accessible APBS points of compatible sign for explicit ionic mechanisms",
                   "nonionic_channels_unchanged":True,
                   "unmodeled":"ligands/metal ions/PTMs removed by protein-only preparation"},
        "maps":records,
    }

def top_members(mapdata,k=5):
    return set().union(*(set(p["members"]) for p in mapdata["patches"][:k])) if mapdata["patches"] else set()

def comparison_rows(old,new,pdb):
    out=[]
    for channel,changed in new["maps"].items():
        reference=old["maps"][channel]
        old_5=top_members(reference);new_5=top_members(changed)
        union=old_5|new_5
        j=len(old_5&new_5)/len(union) if union else 1.
        out.append({"pdb":pdb,"map":channel,"model_from":old["version"],
                    "model_to":MODEL_VERSION,"baseline_patches":len(reference["patches"]),
                    "new_patches":changed["n_patches"],
                    "baseline_top5_members":len(old_5),
                    "new_top5_members":len(new_5),
                    "top5_jaccard":round(j,6),
                    "top1_same_center":int(bool(reference["patches"] and changed["patches"] and
                           reference["patches"][0]["center_key"]==changed["patches"][0]["center_key"])),
                    "n_surface":new["physics"]["surface_residues"],
                    "n_potential_sampled":new["physics"]["potential_sampled_residues"]})
    return out
