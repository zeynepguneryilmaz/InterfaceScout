"""Postprocess PDB2PQR partial charges and APBS OpenDX on probe-accessible sites.

The probe-sphere point sampling is descriptive: not an exact triangulated SES.
The output is physical electrostatic potential in kBT/e and mV, not affinity.
"""
from __future__ import annotations
from dataclasses import dataclass
import math,re
import numpy as np

class PhysicsInputError(ValueError):pass

@dataclass(frozen=True)
class PQRAtom:
    key: str
    atom: str
    xyz: tuple[float,float,float]
    charge: float
    radius: float

def read_pqr(text):
    out=[]
    for line in text.splitlines():
        if not line.startswith(("ATOM  ","HETATM")):continue
        w=line.split()
        if len(w)<10:raise PhysicsInputError("Invalid PQR atom line")
        xyz=tuple(map(float,w[-5:-2]));q,r=map(float,w[-2:])
        if not all(math.isfinite(x) for x in (*xyz,q,r)) or r<0 or (r==0 and not re.match(r"^(?:[123])?H",w[2].upper())):
            raise PhysicsInputError("Invalid atomic charge/radius")
        # PDB2PQR variants can use chain/resseq/icode as one, two or
        # three whitespace tokens. Never infer chain by token count alone.
        identity=w[4:-5]
        if len(identity)==1:
            chain="";residue_token=identity[0];icode_extra=""
        elif len(identity)==2:
            if re.fullmatch(r"[+-]?\\d+[A-Za-z]?",identity[0]):
                chain="";residue_token=identity[0];icode_extra=identity[1]
            else:
                chain=identity[0];residue_token=identity[1];icode_extra=""
        elif len(identity)==3:
            chain,residue_token,icode_extra=identity
        else:
            raise PhysicsInputError("Unsupported PQR residue fields: "+repr(line[:160]))
        residue=re.fullmatch(r"([+-]?\\d+)([A-Za-z]?)",residue_token)
        if not residue or (icode_extra and not re.fullmatch(r"[A-Za-z]",icode_extra)):
            raise PhysicsInputError("Invalid PQR residue number: "+repr(line[:160]))
        icode=(residue[2] or icode_extra).strip()
        if residue[2] and icode_extra and residue[2]!=icode_extra:
            raise PhysicsInputError("Conflicting PQR insertion codes: "+repr(line[:160]))
        out.append(PQRAtom(f"{chain}:{int(residue[1])}:{icode}",w[2],xyz,q,r))
    if not out:raise PhysicsInputError("No PQR atoms")
    return out

@dataclass(frozen=True)
class OpenDX:
    origin: np.ndarray
    step: np.ndarray
    values: np.ndarray
    def sample(self,points):
        xyz=np.asarray(points,float)
        if xyz.ndim!=2 or xyz.shape[1]!=3:raise PhysicsInputError("Expected Nx3 coordinates")
        fractional=(xyz-self.origin)/self.step
        i=np.floor(fractional).astype(int);f=fractional-i
        ok=np.all((i>=0)&(i<np.array(self.values.shape)-1),axis=1)
        out=np.full(len(xyz),np.nan)
        if ok.any():
            j=i[ok];t=f[ok];v=np.zeros(len(j),float)
            for x in range(2):
                for y in range(2):
                    for z in range(2):
                        weight=(t[:,0] if x else 1-t[:,0])*(t[:,1] if y else 1-t[:,1])*(t[:,2] if z else 1-t[:,2])
                        v+=weight*self.values[j[:,0]+x,j[:,1]+y,j[:,2]+z]
            out[ok]=v
        return out

def read_dx(text):
    counts=None;origin=None;steps=[];items=None;values=[];reading=False
    for line in text.splitlines():
        s=line.strip()
        if "class gridpositions" in s and "counts" in s:
            counts=tuple(int(x) for x in s.rsplit("counts",1)[1].split()[:3])
        elif s.startswith("origin "):origin=np.array(list(map(float,s.split()[1:4])))
        elif s.startswith("delta "):steps.append(np.array(list(map(float,s.split()[1:4]))))
        elif "class array" in s and "data follows" in s:
            match=re.search(r"\bitems\s+(\d+)",s)
            if not match:raise PhysicsInputError("DX array size missing")
            items=int(match[1]);reading=True
        elif reading and len(values)<items:
            if s.startswith(("attribute","object","component")):break
            values.extend(float(x) for x in s.split())
    if counts is None or origin is None or len(steps)!=3 or items is None:
        raise PhysicsInputError("Incomplete OpenDX")
    if math.prod(counts)!=items or len(values)<items:
        raise PhysicsInputError("OpenDX dimensions inconsistent")
    delta=np.asarray(steps)
    if np.any(abs(delta-np.diag(np.diag(delta)))>1e-8) or np.any(np.diag(delta)<=0):
        raise PhysicsInputError("Only axis-aligned positive grids supported")
    array=np.asarray(values[:items],float).reshape(counts)
    if not np.isfinite(array).all():raise PhysicsInputError("Non-finite potential")
    return OpenDX(origin,np.diag(delta),array)

def charge_by_residue(atoms):
    charges={}
    for a in atoms:charges[a.key]=charges.get(a.key,0)+a.charge
    return charges

def sample_exposed_potentials(atoms,grid,temperature_K,eligible_keys,n_directions=96,probe_radius_A=1.4):
    """Sample probe-inflated heavy sidechain atom spheres, excluding sterically occluded points.

    Equal angular sampling is not a true area-weighted SES mesh. APBS native units
    are kBT/e, positive where the electrostatic potential is positive.
    """
    if n_directions<20:raise PhysicsInputError("Insufficient sampling directions")
    from scipy.spatial import cKDTree  # optional heavy dependency, only for APBS sampling
    atoms=list(atoms);eligible=set(eligible_keys)
    xyz=np.asarray([a.xyz for a in atoms],float)
    radii=np.asarray([a.radius for a in atoms],float)
    heavy=np.array([i for i,a in enumerate(atoms) if not a.atom.startswith("H")],dtype=int)
    if not len(heavy):raise PhysicsInputError("No heavy PQR atoms")
    tree=cKDTree(xyz[heavy]);extended=radii[heavy]+probe_radius_A
    t=np.arange(n_directions)+0.5;z=1-2*t/n_directions;angle=np.pi*(3-np.sqrt(5))*t
    direction=np.stack([np.cos(angle)*np.sqrt(1-z*z),np.sin(angle)*np.sqrt(1-z*z),z],axis=1)
    seen={}
    for atomindex in heavy:
        a=atoms[int(atomindex)]
        if a.key not in eligible or a.atom in {"N","CA","C","O","OXT"}:continue
        points=np.asarray(a.xyz)[None,:]+(a.radius+probe_radius_A)*direction
        valid=[]
        for p in points:
            neighbors=tree.query_ball_point(p,float(np.max(extended))+0.1)
            valid.append(all(int(heavy[j])==atomindex or np.sum((xyz[heavy[j]]-p)**2)>=extended[j]**2-1e-5 for j in neighbors))
        if not any(valid):continue
        values=grid.sample(points[np.array(valid,bool)])
        finite=values[np.isfinite(values)]
        if len(finite):seen.setdefault(a.key,[]).extend(finite.tolist())
    mV=temperature_K*0.08617333262145
    return {k:{"potential_mean_kBT_e":float(np.mean(v)),
               "potential_mean_mV":float(np.mean(v)*mV),
               "potential_min_mV":float(np.min(v)*mV),
               "potential_max_mV":float(np.max(v)*mV),
               "fraction_positive":float(np.mean(np.asarray(v)>0)),
               "n_probe_surface_points":len(v)} for k,v in seen.items()}
