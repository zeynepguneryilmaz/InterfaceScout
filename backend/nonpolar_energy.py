"""InterfaceScout 2.0 nonpolar-interface orientation energy.

Lightweight deterministic adaptation of established hydrophobic-surface physics:

    DeltaE(n,z) = E_vdW(n,z) + DeltaG_solv(n,z)

* Protein atom LJ types/parameters: OpenMM CHARMM36.
* Surface vdW: neutral graphitic continuum plane (integrated 12-6 -> 10-4 LJ).
* Solvation: Harrison-style signed SASA burial with 1.4 A probe.
* Sampling: deterministic Fibonacci orientation grid and fixed separation grid.

The implementation is deliberately not fitted to adsorption labels and is not an
absolute adsorption-free-energy calculator. It is a rigid-body orientation ranking
model. The original InterfaceScout 1.0 residue/patch scores remain untouched.
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
except Exception:
    openmm = None
    unit = None
    ForceField = Modeller = NoCutoff = PDBFile = None

HYDROPHOBIC_RES = {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP"}
TYR_RING = {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"}
BACKBONE_HYDROPHOBIC = {"CA", "C"}
BACKBONE_HYDROPHILIC = {"N", "O", "OXT"}

HARRISON_SIGMA_KJ_MOL_A2 = 100.0
SASA_PROBE_A = 1.4
CONTACT_A = 6.0

# CHARMM aromatic carbon convention used as the neutral graphitic surface atom.
SURFACE_EPSILON_KJ_MOL = 0.070 * 4.184
SURFACE_RMIN2_A = 1.9924
SURFACE_SIGMA_A = (2.0 * SURFACE_RMIN2_A) / (2.0 ** (1.0 / 6.0))
GRAPHENE_CC_A = 1.42
GRAPHENE_AREAL_DENSITY_A2 = 4.0 / (3.0 * math.sqrt(3.0) * GRAPHENE_CC_A**2)

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
    rn, an = resname.upper(), atom_name.upper()
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
            k0 = (str(chain.id), int(res.id[1]), str(res.id[2]).strip())
            for atom in res.get_atoms():
                if str(getattr(atom, "element", "")).upper() == "H":
                    continue
                out[(k0[0], k0[1], k0[2], atom.get_name().strip())] = float(getattr(atom, "sasa", 0.0) or 0.0)
    return out


def _table_values_2d(func) -> Tuple[int, int, List[float]]:
    """Read an OpenMM Discrete2DFunction-like coefficient table."""
    params = func.getFunctionParameters()
    if len(params) != 3:
        raise RuntimeError(f"Unsupported CHARMM LJ table representation: {func.__class__.__name__}")
    width, height, values = int(params[0]), int(params[1]), list(params[2])
    return width, height, [float(v) for v in values]


def _diag_value(width: int, height: int, values: List[float], t: int) -> float:
    if t < 0 or t >= width or t >= height:
        raise RuntimeError(f"CHARMM LJ atom type index {t} outside {width}x{height} table")
    # OpenMM Discrete2DFunction stores values with x varying fastest.
    return float(values[t + width * t])


def _extract_lj_parameters(system) -> Tuple[np.ndarray, np.ndarray, str]:
    """Decode CHARMM36 per-particle LJ sigma [A] and epsilon [kJ/mol].

    OpenMM's CHARMM36 port represents LJ as
        acoef(type1,type2)/r^12 - bcoef(type1,type2)/r^6
    in a CustomNonbondedForce. Particle data contain only an integer type index;
    A/B pair coefficients reside in tabulated functions. For each atom type,
    the self-pair coefficients give the equivalent standard 12-6 parameters:
        sigma = (A/B)^(1/6)
        epsilon = B^2/(4A)
    where r and sigma are in nm in OpenMM internal units.
    """
    for force in system.getForces():
        if not isinstance(force, openmm.CustomNonbondedForce):
            continue
        pnames = [str(force.getPerParticleParameterName(i)).lower() for i in range(force.getNumPerParticleParameters())]
        fnames = [str(force.getTabulatedFunctionName(i)).lower() for i in range(force.getNumTabulatedFunctions())]
        if pnames != ["type"] or "acoef" not in fnames or "bcoef" not in fnames:
            continue

        functions = {str(force.getTabulatedFunctionName(i)).lower(): force.getTabulatedFunction(i) for i in range(force.getNumTabulatedFunctions())}
        aw, ah, avals = _table_values_2d(functions["acoef"])
        bw, bh, bvals = _table_values_2d(functions["bcoef"])
        if (aw, ah) != (bw, bh):
            raise RuntimeError("CHARMM36 acoef/bcoef tables have inconsistent dimensions")

        sig, eps = [], []
        cache: Dict[int, Tuple[float, float]] = {}
        for i in range(force.getNumParticles()):
            t = int(round(float(force.getParticleParameters(i)[0])))
            if t not in cache:
                A = _diag_value(aw, ah, avals, t)
                B = _diag_value(bw, bh, bvals, t)
                if A <= 0.0 or B <= 0.0:
                    cache[t] = (0.0, 0.0)
                else:
                    sigma_nm = (A / B) ** (1.0 / 6.0)
                    epsilon_kj = (B * B) / (4.0 * A)
                    cache[t] = (sigma_nm * 10.0, epsilon_kj)
            s, e = cache[t]
            sig.append(s); eps.append(e)
        return np.asarray(sig), np.asarray(eps), "OpenMM CHARMM36 type-indexed acoef/bcoef tables"

    # Generic force-field fallback.
    for force in system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            sig, eps = [], []
            for i in range(force.getNumParticles()):
                _, s, e = force.getParticleParameters(i)
                sig.append(float(s.value_in_unit(unit.angstrom)))
                eps.append(abs(float(e.value_in_unit(unit.kilojoule_per_mole))))
            arr_e = np.asarray(eps)
            if np.max(arr_e, initial=0.0) > 1e-12:
                return np.asarray(sig), arr_e, "OpenMM NonbondedForce"
    raise RuntimeError("No nonzero CHARMM/OpenMM Lennard-Jones parameter source was found")


def _openmm_atoms(pdb_path: Path, pH: float, sasa_map: Dict[Tuple[str, int, str, str], float]):
    if openmm is None:
        raise RuntimeError("OpenMM is not installed")
    pdb = PDBFile(str(pdb_path))
    ff = ForceField("charmm36.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff, pH=float(pH))
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
    lj_sigma_A, lj_eps_kj, source = _extract_lj_parameters(system)
    positions = modeller.positions.value_in_unit(unit.angstrom)
    xyz = np.asarray([[p.x, p.y, p.z] for p in positions], dtype=float)

    rows: List[Dict[str, Any]] = []
    for i, atom in enumerate(modeller.topology.atoms()):
        element = atom.element.symbol.upper() if atom.element is not None else ""
        if element == "H":
            continue
        s, e = float(lj_sigma_A[i]), float(lj_eps_kj[i])
        if s <= 0 or e <= 0:
            raise RuntimeError(f"Invalid LJ parameters for heavy atom {atom}: sigma={s}, epsilon={e}")
        res = atom.residue
        try:
            seq = int(res.id)
        except Exception:
            raise RuntimeError(f"Non-integer residue id '{res.id}' in force-field topology")
        icode = str(getattr(res, "insertionCode", "") or "").strip()
        chain = str(res.chain.id)
        sasa = sasa_map.get((chain, seq, icode, atom.name), 0.0)
        rows.append({
            "coord": xyz[i], "sigma_A": s, "epsilon_kj_mol": e,
            "rmin2_A": s * (2.0 ** (1.0 / 6.0)) / 2.0,
            "chain": chain, "res_seq": seq, "icode": icode,
            "res_name": res.name.upper(), "atom": atom.name,
            "sasa_A2": float(sasa), "solvation_sign": _group_sign(res.name, atom.name),
        })
    if not rows:
        raise RuntimeError("No heavy atoms available after CHARMM36 typing")
    eps = np.asarray([a["epsilon_kj_mol"] for a in rows])
    sig = np.asarray([a["sigma_A"] for a in rows])
    diagnostics = {
        "n_heavy_atoms": len(rows),
        "epsilon_min_kj_mol": float(eps.min()), "epsilon_median_kj_mol": float(np.median(eps)), "epsilon_max_kj_mol": float(eps.max()),
        "sigma_min_A": float(sig.min()), "sigma_median_A": float(np.median(sig)), "sigma_max_A": float(sig.max()),
    }
    return rows, source, diagnostics


def _plane_lj_per_atom(z_A: np.ndarray, sigma_A: np.ndarray, epsilon_kj: np.ndarray) -> np.ndarray:
    cross_sigma = 0.5 * (sigma_A + SURFACE_SIGMA_A)
    cross_eps = np.sqrt(epsilon_kj * SURFACE_EPSILON_KJ_MOL)
    z = np.maximum(z_A, 1e-6)
    return 4.0 * math.pi * GRAPHENE_AREAL_DENSITY_A2 * cross_eps * (
        cross_sigma**12 / (5.0 * z**10) - cross_sigma**6 / (2.0 * z**4)
    )


def _solvation_delta(z_A, sasa_A2, signs, rmin2_A) -> Tuple[float, float, float, float]:
    R = rmin2_A + SASA_PROBE_A
    surface_accessible_height = SURFACE_RMIN2_A + SASA_PROBE_A
    d = z_A - surface_accessible_height
    frac = np.clip((R - d) / (2.0 * R), 0.0, 1.0)
    buried = sasa_A2 * frac
    signed = float(np.sum(signs * buried))
    protein_buried = float(np.sum(buried))
    surface_buried = protein_buried
    dg = -HARRISON_SIGMA_KJ_MOL_A2 * signed - HARRISON_SIGMA_KJ_MOL_A2 * surface_buried
    return dg, protein_buried, surface_buried, signed


def scan(pdb_path: Path, struct, pH: float = 7.4, n_orientations: int = DEFAULT_N_ORIENTATIONS,
         separations_A: Optional[Tuple[float, ...]] = None) -> Dict[str, Any]:
    sep_grid = separations_A or MIN_SEPARATIONS_A
    try:
        atoms, source, ljdiag = _openmm_atoms(Path(pdb_path), pH, _biopython_sasa_map(struct))
    except Exception as exc:
        return {"status": "unavailable", "method": "vdW + Harrison-style SASA solvation planar screen", "reason": str(exc)}

    xyz = np.vstack([a["coord"] for a in atoms])
    sigma = np.asarray([a["sigma_A"] for a in atoms])
    epsilon = np.asarray([a["epsilon_kj_mol"] for a in atoms])
    rmin2 = np.asarray([a["rmin2_A"] for a in atoms])
    sasa = np.asarray([a["sasa_A2"] for a in atoms])
    signs = np.asarray([a["solvation_sign"] for a in atoms])
    normals = fibonacci_sphere(int(n_orientations))

    N = len(normals)
    bestE = np.full(N, np.inf); bestZ = np.zeros(N); bestV = np.zeros(N); bestS = np.zeros(N)
    bestPB = np.zeros(N); bestSB = np.zeros(N); bestSigned = np.zeros(N)
    for k, n in enumerate(normals):
        projection = xyz @ n
        depth = projection - float(projection.min())
        for sep in sep_grid:
            z = depth + float(sep)
            v = float(np.sum(_plane_lj_per_atom(z, sigma, epsilon)))
            s, pb, sb, signed = _solvation_delta(z, sasa, signs, rmin2)
            e = v + s
            if e < bestE[k]:
                bestE[k], bestZ[k], bestV[k], bestS[k] = e, sep, v, s
                bestPB[k], bestSB[k], bestSigned[k] = pb, sb, signed

    order = np.argsort(bestE, kind="mergesort")
    top = []
    for rank, idx in enumerate(order[:20], 1):
        n = normals[idx]
        depth = xyz @ n; depth = depth - float(depth.min())
        z = depth + bestZ[idx]
        residues: Dict[str, Dict[str, Any]] = {}
        for j, a in enumerate(atoms):
            if z[j] > CONTACT_A:
                continue
            key = f"{a['chain']}:{a['res_seq']}:{a['icode']}"
            row = residues.setdefault(key, {"key":key,"chain":a["chain"],"res_seq":a["res_seq"],"icode":a["icode"],"res_name":a["res_name"],"min_atom_height_A":float(z[j]),"contact_atom_count":0})
            row["min_atom_height_A"] = min(row["min_atom_height_A"], float(z[j])); row["contact_atom_count"] += 1
        rlist = sorted(residues.values(), key=lambda r:(r["min_atom_height_A"],r["chain"],r["res_seq"]))
        for r in rlist: r["min_atom_height_A"] = round(r["min_atom_height_A"],4)
        top.append({
            "rank":rank,"orientation_index":int(idx),"normal":[round(float(x),6) for x in n],
            "minimum_separation_A":round(float(bestZ[idx]),4),
            "total_energy_change_kj_mol":round(float(bestE[idx]),4),
            "vdw_energy_kj_mol":round(float(bestV[idx]),4),
            "solvation_energy_change_kj_mol":round(float(bestS[idx]),4),
            "protein_buried_sasa_A2":round(float(bestPB[idx]),4),
            "surface_buried_sasa_A2_approx":round(float(bestSB[idx]),4),
            "signed_protein_burial_A2":round(float(bestSigned[idx]),4),
            "contact_residues":rlist,
        })

    return {
        "status":"ok",
        "method":"CHARMM36 vdW + Harrison-style SASA-solvation deterministic planar orientation screen",
        "energy_definition":"DeltaE = integrated graphitic-plane Lennard-Jones + orientation-dependent SASA solvation change",
        "parameter_source":source,"lj_parameter_diagnostics":ljdiag,
        "n_orientations":int(n_orientations),"separation_grid_A":[round(float(x),3) for x in sep_grid],"contact_distance_A":CONTACT_A,
        "graphitic_surface":{"neutral":True,"carbon_epsilon_kj_mol":SURFACE_EPSILON_KJ_MOL,"carbon_Rmin2_A":SURFACE_RMIN2_A,"carbon_sigma_A":SURFACE_SIGMA_A,"areal_density_atoms_A2":GRAPHENE_AREAL_DENSITY_A2,"plane_model":"integrated 12-6 LJ (10-4 planar potential)"},
        "solvation":{"Harrison_sigma_kj_mol_A2":HARRISON_SIGMA_KJ_MOL_A2,"probe_A":SASA_PROBE_A,"surface_burial_approximation":"equal protein/surface interfacial buried area in continuum limit"},
        "no_fitted_weights":True,"absolute_adsorption_free_energy":False,"exact_Harrison_MC_reproduction":False,
        "differences_from_Harrison":["CHARMM36 instead of historical CHARMM22/CMAP","integrated neutral graphitic plane instead of atomically corrugated explicit graphene","deterministic orientation/separation scan instead of Metropolis Monte Carlo","continuum spherical-cap/equal-area SASA-burial approximation"],
        "best_energy_change_kj_mol":round(float(bestE[order[0]]),4),"median_best_energy_change_kj_mol":round(float(np.median(bestE)),4),"top_orientations":top,
    }
