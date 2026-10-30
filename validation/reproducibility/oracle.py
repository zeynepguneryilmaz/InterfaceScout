"""Independent scalar checks of the declared model, not an adsorption oracle.

No imports from the application. Shared dependencies: PDBParser, the declared
chemistry/ASA tables, and (for large structures only) Biopython ShrakeRupley.
Discrete decisions use verified, unrounded application values; scalar numerical
re-evaluation checks those values without silently changing floating-point ties.
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from io import StringIO
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley, ATOMIC_RADII

ATOL = 1e-10
RTOL = 1e-10
BFACTOR_ATOL = 0.0050000001
SPEC_PATH = Path(__file__).with_name("model_spec.json")


def spec():
    return json.loads(SPEC_PATH.read_text(encoding="utf-8"))


class Checks:
    def __init__(self):
        self.count = 0
        self.max_abs_error = 0.0
        self.near_ties = {}

    def equal(self, actual, expected, label):
        self.count += 1
        if actual != expected:
            raise AssertionError(f"{label}: {actual!r} != {expected!r}")

    def close(self, actual, expected, label):
        self.count += 1
        a, e = float(actual), float(expected)
        error = abs(a-e)
        if not math.isfinite(a) or not math.isfinite(e) or error > ATOL + RTOL*abs(e):
            raise AssertionError(f"{label}: {a:.17g} != {e:.17g}; abs_error={error:.3g}")
        self.max_abs_error = max(self.max_abs_error, error)


def sampled_sasa(atoms, probe=1.4, points=200):
    """All-pairs atom occlusion, no Biopython SASA or KD-tree implementation.

    Uses the declared golden-spiral quadrature and float32 unit-sphere storage,
    matching the numerical sampling convention rather than a different lattice.
    Intended for small fixtures; O(atoms^2 + nearby_atoms * points).
    """
    xyz = np.asarray([a.coord for a in atoms], dtype=np.float64)
    radii = np.asarray([ATOMIC_RADII[a.element]+probe for a in atoms])
    dz = 2.0/points
    angle, z = 0.0, 1.0-dz/2
    unit = []
    for _ in range(points):
        radial = math.sqrt(1-z*z)
        unit.append((math.cos(angle)*radial, math.sin(angle)*radial, z))
        z -= dz
        angle += math.pi*(3-math.sqrt(5))
    unit = np.asarray(unit, dtype=np.float32)
    areas = []
    for i, radius in enumerate(radii):
        samples = unit*radius + xyz[i]
        d = xyz-xyz[i]
        neighbors = np.flatnonzero(np.sum(d*d, axis=1) < (radius+radii)**2)
        visible = np.ones(points, dtype=bool)
        for j in neighbors:
            if j != i:
                delta = samples-xyz[j]
                visible &= np.sum(delta*delta, axis=1) > radii[j]**2
        areas.append(int(visible.sum()) * (radius*radius*(4*math.pi/points)))
    return areas


def accessibility(prepared, independent_sasa=False):
    cfg = spec()
    structure = PDBParser(QUIET=True).get_structure("reference", StringIO(prepared))
    atoms = list(structure.get_atoms())
    if independent_sasa:
        areas = sampled_sasa(atoms, cfg["parameters"]["sasa_probe_A"], cfg["parameters"]["sasa_points_per_atom"])
        for atom, area in zip(atoms, areas):
            atom.sasa = area
    else:
        ShrakeRupley(probe_radius=cfg["parameters"]["sasa_probe_A"],
                    n_points=cfg["parameters"]["sasa_points_per_atom"]).compute(structure, level="A")
    rows = []
    for chain in next(structure.get_models()):
        for residue in chain:
            name = residue.resname.strip()
            if not is_aa(residue, standard=True) or "CA" not in residue:
                continue
            selected = [residue["CA"]] if name == "GLY" else [
                a for a in residue if a.name not in {"N", "CA", "C", "O", "OXT"} and a.element != "H"]
            asa = sum(a.sasa for a in selected)
            raw = asa/cfg["reference_sidechain_asa"][name]
            rows.append(dict(key=f"{chain.id}:{residue.id[1]}:{residue.id[2].strip()}",
                             res_name=name, coord=tuple(float(v) for v in residue["CA"].coord),
                             sidechain_sasa=asa, scrsa_raw=raw, scrsa=min(1.0,max(0.0,raw)),
                             surface_exposed=raw >= cfg["parameters"]["scrsa_threshold"]))
    return rows


def availability(name, state, ph, cfg):
    if state not in {"protonated", "deprotonated"} or name not in cfg["pka"]:
        return 1.0
    exponent = ph-cfg["pka"][name]
    if state == "deprotonated":
        exponent = -exponent
    return 1.0/(1.0+10.0**exponent)


def scalar_channel(rows, chemistry, ph):
    cfg = spec()
    favorable = cfg["chemistries"][chemistry]["favorable"]
    local = {r["key"]: r["scrsa"]*availability(r["res_name"], favorable[r["res_name"]][2], ph, cfg)
             if r["res_name"] in favorable else 0.0 for r in rows}
    distances = {(a["key"], b["key"]): math.dist(a["coord"], b["coord"]) for a in rows for b in rows}
    densities = {}
    for radius in cfg["parameters"]["multiscale_radii_A"]:
        densities[radius] = {a["key"]: sum(local[b["key"]] for b in rows if distances[a["key"],b["key"]] <= radius) for a in rows}
    maxlocal = max(local.values(), default=0.0)
    maxima = {r:max(v.values(), default=0.0) for r,v in densities.items()}
    values = {}
    for row in rows:
        key = row["key"]
        norm = [densities[r][key]/maxima[r] if maxima[r] else 0.0 for r in densities]
        values[key] = dict(local_score=local[key], propensity=100*local[key]/maxlocal if maxlocal else 0.,
                           multiscale_persistence=100*min(norm), multiscale_geomean=100*math.sqrt(math.prod(norm)))
        for r in densities:
            values[key][f"density_{int(r)}A_raw"] = densities[r][key]
            values[key][f"density_{int(r)}A_norm"] = 100*densities[r][key]/maxima[r] if maxima[r] else 0.
    return values, distances


def dominance_depths(patches):
    """Longest path in the strict-dominance DAG, rather than front peeling."""
    fields = ("chemistry_support", "mean_accessibility", "patch_coherence", "orientation_coherence")
    predecessors = {}
    for key, a in patches.items():
        predecessors[key] = [k for k,b in patches.items() if k != key and
                             all(b[f] >= a[f] for f in fields) and any(b[f] > a[f] for f in fields)]
    @lru_cache(None)
    def depth(key):
        return 1+max((depth(k) for k in predecessors[key]), default=0)
    return {key:depth(key) for key in patches}


def verify_map(rows, all_rows, channel, public_map, chemistry, ph, checks):
    values, distances = scalar_channel(rows, chemistry, ph)
    cfg = spec()
    fav = cfg["chemistries"][chemistry]["favorable"]
    bykey = {r["key"]:r for r in rows}
    members = {r["key"]:r for r in channel["residues"]}
    checks.equal(set(members), {r["key"] for r in rows if r["res_name"] in fav}, chemistry+" residue set")
    for key,r in members.items():
        for field in ("local_score", "propensity", "multiscale_persistence", "multiscale_geomean"):
            checks.close(r[field],values[key][field],f"{chemistry}/{key}/{field}")
        if values[key]["local_score"] > 0:
            checks.equal(r["local_score"] > 0,True,chemistry+" positive local score")
        checks.close(r["state_availability"], availability(bykey[key]["res_name"],fav[bykey[key]["res_name"]][2],ph,cfg),"availability")
    exported = {r["key"]:r for r in public_map["residues"]}
    checks.equal(set(exported),set(members),"API residue identities")
    for key,r in exported.items():
        for field in ("scrsa", "local_score", "propensity", "multiscale_persistence"):
            checks.close(r[field],members[key][field],"API "+key+" "+field)
    centers = {c["center_key"]:c for c in channel["patch_centers"]}
    checks.equal(set(centers),{k for k,v in values.items() if v["multiscale_persistence"]>0},"centers")
    for key,c in centers.items():
        for field in ("density_6A_raw","density_9A_raw","density_6A_norm","density_9A_norm","multiscale_persistence","multiscale_geomean"):
            checks.close(c[field],values[key][field],chemistry+"/"+key+"/"+field)
    # Parser coordinates are promoted to float64 before the geometric centroid.
    centroid = np.asarray([r["coord"] for r in all_rows],dtype=np.float64).mean(axis=0)
    normals = {}
    for r in rows:
        vector = [v-float(c) for v,c in zip(r["coord"],centroid)]
        norm = math.sqrt(sum(v*v for v in vector))
        normals[r["key"]] = [v/norm for v in vector] if norm>1e-12 else [0.,0.,0.]
    def face(a,b):
        return sum(x*y for x,y in zip(normals[a],normals[b]))>0
    radius=cfg["parameters"]["patch_radius_A"]
    priority=sorted(centers,key=lambda k:(-centers[k]["multiscale_persistence"],-centers[k]["multiscale_geomean"],k))
    retained=[]
    for key in priority:
        if not any(distances[key,k]<=radius and face(key,k) for k in retained):
            retained.append(key)
    if not any(r["local_score"]>0 for r in members.values()):
        retained=[]
    patches={p["center_key"]:p for p in public_map["patches"]}
    checks.equal(set(patches),set(retained),"selected patch centers")
    for key,p in patches.items():
        ids=sorted({k for k in bykey if distances[key,k]<=radius and face(key,k)}|{key})
        seeds=sorted(k for k in ids if k in members and members[k]["local_score"]>0)
        checks.equal(p["members"],ids,chemistry+"/"+key+" members")
        checks.equal(p["seed_members"],seeds,chemistry+"/"+key+" seeds")
        checks.equal(p["n_members"],len(ids),"member count")
        checks.equal(p["n_chemistry_seeds"],len(seeds),"seed count")
        checks.close(p["chemistry_support"],len(seeds)/len(ids),"chemistry support")
        checks.close(p["mean_accessibility"],sum(bykey[k]["scrsa"] for k in ids)/len(ids),"patch accessibility")
        checks.close(p["patch_coherence"],values[key]["multiscale_persistence"]/100,"patch coherence")
        resultant=[sum(normals[k][j] for k in ids)/len(ids) for j in range(3)]
        checks.close(p["orientation_coherence"],math.sqrt(sum(x*x for x in resultant)),"orientation")
        checks.close(p["diameter_A"],max(distances[a,b] for a in ids for b in ids),"diameter")
    fronts=dominance_depths(patches)
    # Report sensitivity without rounding descriptors or changing dominance.
    near = {}
    for field in ("chemistry_support", "mean_accessibility", "patch_coherence", "orientation_coherence"):
        distinct=sorted({p[field] for p in patches.values()})
        gaps=[b-a for a,b in zip(distinct,distinct[1:])]
        near[field]=dict(min_positive_gap=min(gaps,default=None),
            adjacent_gaps_within_tolerance=sum(b-a<=ATOL+RTOL*abs(b) for a,b in zip(distinct,distinct[1:])))
    checks.near_ties[chemistry]=near
    order=sorted(patches,key=lambda k:(fronts[k],-patches[k]["patch_coherence"],-patches[k]["chemistry_support"],k))
    checks.equal(list(patches),order,"patch ranking")
    for rank,key in enumerate(order,1):
        checks.equal(patches[key]["pareto_front"],fronts[key],"Pareto front")
        checks.equal(patches[key]["display_rank"],rank,"display rank")
    checks.equal([p["center_key"] for p in public_map["primary_patches"]],[k for k in order if fronts[k]==1],"primary patches")
    checks.equal(public_map["n_patches"],len(patches),"patch count")
    checks.equal(public_map["n_pareto_primary_patches"],sum(f==1 for f in fronts.values()),"primary count")


def verify_result(result, core_result, independent_sasa=False):
    checks=Checks()
    all_rows=accessibility(result["prepared_pdb"],independent_sasa)
    checks.equal([r["key"] for r in core_result["all_residues"]],[r["key"] for r in all_rows],"all residue keys")
    for actual,expected in zip(core_result["all_residues"],all_rows):
        for field in ("sidechain_sasa","scrsa_raw","scrsa"):
            checks.close(actual[field],expected[field],expected["key"]+"/"+field)
        checks.equal(actual["surface_exposed"],expected["surface_exposed"],"surface flag")
    rows=[r for r in all_rows if r["surface_exposed"]]
    checks.equal(result["diagnostics"]["surface_residue_keys"],[r["key"] for r in rows],"surface keys")
    checks.equal(result["n_maps"],11,"11 maps")
    checks.equal(result["map_order"],list(spec()["public_to_internal"]),"map order")
    checks.equal(set(result["maps"]),set(spec()["public_to_internal"]),"map set")
    for public,internal in spec()["public_to_internal"].items():
        verify_map(rows,all_rows,core_result["chemistries"][internal],result["maps"][public],internal,result["input"]["pH"],checks)
    return dict(checks=checks.count,max_abs_error=checks.max_abs_error,
                near_ties=checks.near_ties,
                sasa_oracle="independent_all_pairs" if independent_sasa else "shared_biopython_1.85")
