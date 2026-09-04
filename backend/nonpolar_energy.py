"""Physics-based nonpolar-interface orientation screening for InterfaceScout 2.0.

This module replaces the earlier SASA-only orientation prototype with a
unit-consistent energy model built from established components:

1. Protein atom Lennard-Jones parameters are assigned with the CHARMM36 force
   field through OpenMM.
2. A neutral graphitic carbon surface is represented as a continuum plane
   obtained by analytically integrating the standard 12-6 Lennard-Jones
   potential over an infinite graphene-like carbon sheet (the standard 10-4
   planar LJ potential).
3. The hydrophobic solvation term follows the atom-group classification used by
   Harrison et al. (Biointerphases 2017, DOI 10.1116/1.4971381): hydrophobic
   groups have +sigma and hydrophilic groups -sigma, with |sigma|=100
   kJ mol^-1 A^-2 and a 1.4 A solvent probe.  The orientation-dependent change
   in protein SASA is estimated from solvent-accessible spherical-cap burial by
   the plane.  Burial of the hydrophobic carbon surface is approximated by an
   equal interfacial-area term, a continuum-interface approximation.

For each rigid-body approach normal, the minimum protein-plane separation is
optimized on a fixed, predeclared grid.  The reported screening energy is

    Delta E(n,z) = E_vdW(n,z) + Delta G_solv(n,z)

No fitted cross-term weights or adsorption labels are used.

Important: this is a lightweight deterministic adaptation of established
physics, not an exact reproduction of Harrison's explicit-graphene GROMACS
Monte Carlo trajectory.  It uses CHARMM36 rather than historical CHARMM22,
an analytically integrated neutral carbon plane rather than atomically
corrugated graphene, and deterministic orientation/separation scanning rather
than Metropolis sampling.  These distinctions are reported in the output.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from Bio.PDB.Polypeptide import is_aa

try:
    import openmm
    from openmm import unit
    from openmm.app import ForceField, Modeller, NoCutoff, PDBFile
except Exception:  # optional at import time; required for this physics layer
    openmm = None
    unit = None
    ForceField = Modeller = NoCutoff = PDBFile = None

# Harrison atom-group classification.
HYDROPHOBIC_RES = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP"}
TYR_RING = {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"}
BACKBONE_HYDROPHOBIC = {"CA", "C"}
BACKBONE_HYDROPHILIC = {"N", "O", "OXT"}

# Published Harrison surface-tension magnitude and solvent probe.
HARRISON_SIGMA_KJ_MOL_A2 = 100.0
SASA_PROBE_A = 1.4
CONTACT_A = 6.0

# Neutral graphitic/aromatic carbon parameters from CHARMM protein parameter
# convention: epsilon = 0.070 kcal/mol; Rmin/2 = 1.9924 A for aromatic CA.
SURFACE_EPSILON_KJ_MOL = 0.070 * 4.184
SURFACE_RMIN2_A = 1.9924
SURFACE_SIGMA_A = (2.0 * SURFACE_RMIN2_A) / (2.0 ** (1.0 / 6.0))
GRAPHENE_CC_A = 1.42
GRAPHENE_AREAL_DENSITY_A2 = 4.0 / (3.0 * math.sqrt(3.0) * GRAPHENE_CC_A**2)

# Deterministic scan settings; these are model settings, not benchmark-tuned.
DEFAULT_N_ORIENTATIONS = 1024
MIN_SEPARATIONS_A = tuple(float(x) for x in np.arange(2.8, 5.61, 0.15))


def fibonacci_sphere(n: int = DEFAULT_N_ORIENTATIONS) -> np.ndarray:
    i = np.arange(n, dtype=float)
    phi = (1.0 + 5.0**0.5) / 2.0
    z = 1.0 - 2.0 * (i + 0.5) / n
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    theta = 2.0 * np.pi * i / phi
    return np.column_stack([r * np.cos(theta), r * np.sin(theta), z])


def _group_sign(resname: str, atom_name: str) -> float:
    """Harrison group sign: +1 hydrophobic, -1 hydrophilic."""
    rn = resname.upper()
    an = atom_name.upper()
    if an in BACKBONE_HYDROPHOBIC:
        return 1.0
    if an in BACKBONE_HYDROPHILIC:
        return -1.0
    if rn == "TYR":
        return 1.0 if an in TYR_RING else -1.0
    return 1.0 if rn in HYDROPHOBIC_RES else -1.0


def _biopython_sasa_map(struct) -> Dict[Tuple[str, int, str, str], float]:
    out: Dict[Tuple[str, int, str, str], float] = {}
    model = next(struct.get_models())
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=True):
                continue
            key0 = (str(chain.id), int(res.id[1]), str(res.id[2]).strip())
            for atom in res.get_atoms():
                if str(getattr(atom, "element", "")).upper() == "H":
                    continue
                out[(key0[0], key0[1], key0[2], atom.get_name().strip())] = float(getattr(atom, "sasa", 0.0) or 0.0)
    return out


def _openmm_atoms(pdb_path: Path, pH: float, sasa_map: Dict[Tuple[str, int, str, str], float]) -> Tuple[List[Dict[str, Any]], str]:
    if openmm is None:
        raise RuntimeError("OpenMM is not installed; install openmm>=8.2 for the full nonpolar energy model")

    pdb = PDBFile(str(pdb_path))
    ff = ForceField("charmm36.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    # Add hydrogens with force-field aware protonation assignment.  The canonical
    # InterfaceScout score is unaffected; this is only for the optional vdW layer.
    modeller.addHydrogens(ff, pH=float(pH))
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)

    nb = None
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            nb = force
            break
    if nb is None:
        raise RuntimeError("CHARMM36 system did not contain an OpenMM NonbondedForce")

    pos = np.asarray([[p.x, p.y, p.z] for p in modeller.positions.value_in_unit(unit.angstrom)], dtype=float)
    rows: List[Dict[str, Any]] = []
    for idx, atom in enumerate(modeller.topology.atoms()):
        _, sigma, epsilon = nb.getParticleParameters(idx)
        sigma_A = float(sigma.value_in_unit(unit.angstrom))
        eps_kj = float(epsilon.value_in_unit(unit.kilojoule_per_mole))
        residue = atom.residue
        chain = residue.chain
        try:
            seq = int(residue.id)
        except Exception:
            # OpenMM normally preserves integer PDB residue IDs; fail explicitly
            # rather than silently mis-map atom-level descriptors.
            raise RuntimeError(f"Non-integer residue id '{residue.id}' in OpenMM topology")
        icode = str(getattr(residue, "insertionCode", "") or "").strip()
        element = atom.element.symbol.upper() if atom.element is not None else ""
        heavy = element != "H"
        sasa = sasa_map.get((str(chain.id), seq, icode, atom.name), 0.0) if heavy else 0.0
        rmin2_A = sigma_A * (2.0 ** (1.0 / 6.0)) / 2.0
        rows.append({
            "coord": pos[idx],
            "sigma_A": sigma_A,
            "epsilon_kj_mol": eps_kj,
            "rmin2_A": rmin2_A,
            "chain": str(chain.id),
            "res_seq": seq,
            "icode": icode,
            "res_name": residue.name.upper(),
            "atom": atom.name,
            "element": element,
            "heavy": heavy,
            "sasa_A2": float(sasa),
            "solvation_sign": _group_sign(residue.name, atom.name) if heavy else 0.0,
        })
    return rows, "OpenMM CHARMM36"


def _plane_lj_per_atom(z_A: np.ndarray, sigma_A: np.ndarray, epsilon_kj: np.ndarray) -> np.ndarray:
    """Integrated 12-6 LJ over an infinite carbon plane (10-4 potential)."""
    cross_sigma = 0.5 * (sigma_A + SURFACE_SIGMA_A)
    cross_eps = np.sqrt(np.maximum(epsilon_kj, 0.0) * SURFACE_EPSILON_KJ_MOL)
    z = np.maximum(z_A, 1e-6)
    s6 = cross_sigma**6
    return 4.0 * math.pi * GRAPHENE_AREAL_DENSITY_A2 * cross_eps * (
        (cross_sigma**12) / (5.0 * z**10) - s6 / (2.0 * z**4)
    )


def _solvation_delta(
    z_A: np.ndarray,
    atom_sasa_A2: np.ndarray,
    signs: np.ndarray,
    rmin2_A: np.ndarray,
    heavy_mask: np.ndarray,
) -> Tuple[float, float, float, float]:
    """Continuum estimate of orientation-dependent Harrison SASA burial.

    The solvent-accessible sphere of each heavy atom is intersected by the
    solvent-accessible plane of the graphitic surface.  The isolated atom SASA
    is multiplied by the spherical-cap fraction to estimate burial.  The
    surface-side buried area is approximated as the same continuum interfacial
    area as the total protein-side burial.
    """
    if not np.any(heavy_mask):
        return 0.0, 0.0, 0.0, 0.0
    idx = np.where(heavy_mask)[0]
    R = rmin2_A[idx] + SASA_PROBE_A
    surface_accessible_height = SURFACE_RMIN2_A + SASA_PROBE_A
    d = z_A[idx] - surface_accessible_height
    frac = np.clip((R - d) / (2.0 * R), 0.0, 1.0)
    buried = atom_sasa_A2[idx] * frac
    signed_burial = float(np.sum(signs[idx] * buried))
    protein_buried = float(np.sum(buried))
    # Equal-area continuum approximation for the opposing planar surface.
    surface_buried = protein_buried
    delta_g_protein = -HARRISON_SIGMA_KJ_MOL_A2 * signed_burial
    delta_g_surface = -HARRISON_SIGMA_KJ_MOL_A2 * surface_buried
    return delta_g_protein + delta_g_surface, protein_buried, surface_buried, signed_burial


def scan(
    pdb_path: Path,
    struct,
    pH: float = 7.4,
    n_orientations: int = DEFAULT_N_ORIENTATIONS,
    separations_A: Optional[Tuple[float, ...]] = None,
) -> Dict[str, Any]:
    """Scan rigid-body orientation and separation for a neutral nonpolar plane."""
    sep_grid = separations_A or MIN_SEPARATIONS_A
    try:
        sasa_map = _biopython_sasa_map(struct)
        atoms, parameter_source = _openmm_atoms(Path(pdb_path), pH, sasa_map)
    except Exception as exc:
        return {
            "status": "unavailable",
            "method": "vdW + Harrison-style SASA solvation planar screen",
            "reason": str(exc),
        }

    xyz = np.vstack([a["coord"] for a in atoms]).astype(float)
    sigma = np.asarray([a["sigma_A"] for a in atoms], dtype=float)
    epsilon = np.asarray([a["epsilon_kj_mol"] for a in atoms], dtype=float)
    rmin2 = np.asarray([a["rmin2_A"] for a in atoms], dtype=float)
    sasa = np.asarray([a["sasa_A2"] for a in atoms], dtype=float)
    signs = np.asarray([a["solvation_sign"] for a in atoms], dtype=float)
    heavy = np.asarray([bool(a["heavy"]) for a in atoms], dtype=bool)
    normals = fibonacci_sphere(int(n_orientations))

    best_energy = np.full(len(normals), np.inf, dtype=float)
    best_sep = np.zeros(len(normals), dtype=float)
    best_vdw = np.zeros(len(normals), dtype=float)
    best_solv = np.zeros(len(normals), dtype=float)
    best_pb = np.zeros(len(normals), dtype=float)
    best_sb = np.zeros(len(normals), dtype=float)
    best_signed = np.zeros(len(normals), dtype=float)

    for k, n in enumerate(normals):
        p = xyz @ n
        depth = p - float(np.min(p))
        for sep in sep_grid:
            z = depth + float(sep)
            vdw = float(np.sum(_plane_lj_per_atom(z, sigma, epsilon)))
            solv, pb, sb, signed = _solvation_delta(z, sasa, signs, rmin2, heavy)
            total = vdw + solv
            if total < best_energy[k]:
                best_energy[k] = total
                best_sep[k] = float(sep)
                best_vdw[k] = vdw
                best_solv[k] = solv
                best_pb[k] = pb
                best_sb[k] = sb
                best_signed[k] = signed

    order = np.argsort(best_energy, kind="mergesort")
    top: List[Dict[str, Any]] = []
    for rank, idx in enumerate(order[:20], 1):
        n = normals[idx]
        p = xyz @ n
        depth = p - float(np.min(p))
        z = depth + best_sep[idx]
        residues: Dict[str, Dict[str, Any]] = {}
        for j, a in enumerate(atoms):
            if not a["heavy"] or z[j] > CONTACT_A:
                continue
            key = f"{a['chain']}:{a['res_seq']}:{a['icode']}"
            rr = residues.setdefault(key, {
                "key": key,
                "chain": a["chain"],
                "res_seq": a["res_seq"],
                "icode": a["icode"],
                "res_name": a["res_name"],
                "min_atom_height_A": float(z[j]),
                "contact_atom_count": 0,
            })
            rr["min_atom_height_A"] = min(rr["min_atom_height_A"], float(z[j]))
            rr["contact_atom_count"] += 1
        rlist = sorted(residues.values(), key=lambda x: (x["min_atom_height_A"], x["chain"], x["res_seq"]))
        for rr in rlist:
            rr["min_atom_height_A"] = round(float(rr["min_atom_height_A"]), 4)
        top.append({
            "rank": rank,
            "orientation_index": int(idx),
            "normal": [round(float(x), 6) for x in n],
            "minimum_separation_A": round(float(best_sep[idx]), 4),
            "total_energy_change_kj_mol": round(float(best_energy[idx]), 4),
            "vdw_energy_kj_mol": round(float(best_vdw[idx]), 4),
            "solvation_energy_change_kj_mol": round(float(best_solv[idx]), 4),
            "protein_buried_sasa_A2": round(float(best_pb[idx]), 4),
            "surface_buried_sasa_A2_approx": round(float(best_sb[idx]), 4),
            "signed_protein_burial_A2": round(float(best_signed[idx]), 4),
            "contact_residues": rlist,
        })

    return {
        "status": "ok",
        "method": "CHARMM36 vdW + Harrison-style SASA-solvation deterministic planar orientation screen",
        "energy_definition": "DeltaE = integrated carbon-plane Lennard-Jones + orientation-dependent SASA solvation change",
        "parameter_source": parameter_source,
        "n_orientations": int(n_orientations),
        "separation_grid_A": [round(float(x), 3) for x in sep_grid],
        "contact_distance_A": CONTACT_A,
        "graphitic_surface": {
            "neutral": True,
            "carbon_epsilon_kj_mol": SURFACE_EPSILON_KJ_MOL,
            "carbon_Rmin2_A": SURFACE_RMIN2_A,
            "carbon_sigma_A": SURFACE_SIGMA_A,
            "areal_density_atoms_A2": GRAPHENE_AREAL_DENSITY_A2,
            "plane_model": "integrated 12-6 LJ (10-4 planar potential)",
        },
        "solvation": {
            "Harrison_sigma_kj_mol_A2": HARRISON_SIGMA_KJ_MOL_A2,
            "probe_A": SASA_PROBE_A,
            "surface_burial_approximation": "equal protein/surface interfacial buried area in the continuum limit",
        },
        "no_fitted_weights": True,
        "exact_Harrison_MC_reproduction": False,
        "differences_from_Harrison": [
            "CHARMM36 rather than historical CHARMM22",
            "integrated neutral carbon plane rather than atomically corrugated explicit graphene",
            "deterministic orientation/separation scan rather than Metropolis Monte Carlo sampling",
            "continuum spherical-cap/equal-area estimate of SASA burial rather than recomputing explicit-complex SASA in GROMACS",
        ],
        "best_energy_change_kj_mol": round(float(best_energy[order[0]]), 4),
        "median_best_energy_change_kj_mol": round(float(np.median(best_energy)), 4),
        "top_orientations": top,
    }
