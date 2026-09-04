#!/usr/bin/env python3
from __future__ import annotations
import json, sys, tempfile, shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BACKEND=ROOT/'backend'
sys.path.insert(0,str(BACKEND))
import main as core
from structural_context import AnalyzeRequest, EnvParams, prepare_context
from structure_repair import repair_for_forcefield
from openmm.app import PDBFile, ForceField, Modeller, NoCutoff
import openmm

core.PDB2PQR=None; core.APBS=None; core.MKDSSP=None
work=Path(tempfile.mkdtemp(prefix='is_lj_inspect_'))
try:
    req=AnalyzeRequest(pdb_id='1AKI',chain='A',structure_context='selected_chain_legacy',protrusion=False,env=EnvParams(pH=7.0,ionic=150,temp=298))
    pdb,_,_,_=prepare_context(req,work)
    repaired=work/'forcefield_ready.pdb'
    repair=repair_for_forcefield(pdb,repaired)
    if repair.get('status')!='ok':
        raise RuntimeError(repair)
    p=PDBFile(str(repaired))
    ff=ForceField('charmm36.xml')
    m=Modeller(p.topology,p.positions)
    m.addHydrogens(ff,pH=7.0)
    system=ff.createSystem(m.topology,nonbondedMethod=NoCutoff,constraints=None,rigidWater=False)
    out=[]
    for fi,f in enumerate(system.getForces()):
        row={'index':fi,'class':f.__class__.__name__}
        if isinstance(f,openmm.CustomNonbondedForce):
            row['energy_function']=f.getEnergyFunction()
            row['per_particle_parameters']=[f.getPerParticleParameterName(i) for i in range(f.getNumPerParticleParameters())]
            row['global_parameters']=[{'name':f.getGlobalParameterName(i),'default':float(f.getGlobalParameterDefaultValue(i))} for i in range(f.getNumGlobalParameters())]
            row['tabulated_functions']=[f.getTabulatedFunctionName(i) for i in range(f.getNumTabulatedFunctions())]
            row['n_particles']=f.getNumParticles()
            row['first_particle_params']=[[float(x) for x in f.getParticleParameters(i)] for i in range(min(8,f.getNumParticles()))]
        elif isinstance(f,openmm.NonbondedForce):
            row['n_particles']=f.getNumParticles()
            vals=[]
            from openmm import unit
            for i in range(min(8,f.getNumParticles())):
                q,s,e=f.getParticleParameters(i)
                vals.append({'q_e':float(q.value_in_unit(unit.elementary_charge)),'sigma_A':float(s.value_in_unit(unit.angstrom)),'epsilon_kj_mol':float(e.value_in_unit(unit.kilojoule_per_mole))})
            row['first_particle_params']=vals
        out.append(row)
    report={'repair':repair,'forces':out}
    Path('validation/charmm_lj_force_layout.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
finally:
    shutil.rmtree(work,ignore_errors=True)
