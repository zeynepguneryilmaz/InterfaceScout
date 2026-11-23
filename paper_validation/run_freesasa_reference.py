"""Independent FreeSASA side-chain SASA reference for the manuscript.

Runs ONLY as an external evaluation. Does not change InterfaceScout scoring.
FreeSASA default Lee-Richards surface discretization and atomic radii are
different from the Biopython Shrake-Rupley implementation used by InterfaceScout.
Therefore report rank correlation, area differences and exposure agreement
without expecting numerical identity.

Eight pre-specified biological-unit examples provide a manageable independent
sanity check. These are NOT eight independent experimental adsorption outcomes.
"""
from __future__ import annotations
import csv, json, math, sys, warnings
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'backend'))
sys.path.insert(0,str(ROOT/'validation'))
import freesasa
from interface_engine import analyze_all_maps
from curated_biological_units_full_browser_exports import prepare_case
import curated_biological_units_full_browser_exports as base

# Freeze cases before seeing agreement statistics. Include a multimer and
# large dimer; NMR, glycosylation, cofactors and disorder limitations remain.
PDBS=['1AKI','1UBQ','1AO6','1TGU','2HAV','6O5E','2GS2','8EOB']
base.ASSEMBLY.update({'8EOB'})  # correctly selected dimer
OUTPUT=Path('paper_validation_outputs/independent_freesasa')

def key_for(chain, resi):
    token=str(resi).strip()
    digits=''
    for i,x in enumerate(token):
        if x in '+-' and i==0 or x.isdigit():
            digits+=x
        else: break
    if not digits or digits in ('-','+'):return None
    icode=token[len(digits):].strip()
    try: return f'{str(chain).strip()}:{int(digits)}:{icode}'
    except ValueError:return None

def reference_areas(pdb_file):
    struct=freesasa.Structure(str(pdb_file))
    result=freesasa.calc(struct) # FreeSASA default Lee-Richards
    rows={}
    for chain, dictionary in result.residueAreas().items():
        for resno, area in dictionary.items():
            key=key_for(chain,resno)
            if key is None: continue
            value=float(area.relativeSideChain)
            rows[key]={
                'freesasa_sidechain_A2':float(area.sideChain),
                'freesasa_rel_sidechain':value if math.isfinite(value) else None,
                'freesasa_res_type':area.residueType,
            }
    return rows, float(result.totalArea())

def main():
    OUTPUT.mkdir(parents=True,exist_ok=True)
    records=[];cases=[]
    for pdb in PDBS:
        folder=OUTPUT/pdb
        try:
            prov, protein=prepare_case(pdb,folder)
            actual=analyze_all_maps(pdb_text=protein.read_text(),pH=7.4,ionic_mM=150,temp_K=298)
            external, total=reference_areas(protein)
            observed={r['key']:r for r in actual['residue_profiles']}
            common=sorted(set(observed)&set(external))
            numeric=[]; calls=[]
            for key in common:
                a=observed[key]; b=external[key]
                val=b['freesasa_rel_sidechain']
                # GLY sidechain has no equivalent well-defined atom group; skip
                # from sidechain correlation. Keep it in the CSV with flag.
                valid=(a['res_name']!='GLY' and val is not None)
                if valid:
                    numeric.append((float(a['sidechain_sasa']),b['freesasa_sidechain_A2']))
                    calls.append((bool(a['surface_exposed']),val>=0.05))
                records.append({'pdb':pdb,'residue_key':key,'residue_type':a['res_name'],
                   'interfacescout_sidechain_A2':a['sidechain_sasa'],
                   'freesasa_sidechain_A2':b['freesasa_sidechain_A2'],
                   'interfacescout_scRSA':a['scrsa'],
                   'freesasa_relative_sidechain':val,
                   'interfacescout_exposed':int(bool(a['surface_exposed'])),
                   'freesasa_exposed':int(val>=0.05) if val is not None else '',
                   'comparison_eligible':int(valid),
                   'method_note':'Same protein-only coordinates; different SASA algorithms/radii; raw scRSA uses different maximum reference values.'})
            if len(numeric)<25:raise ValueError('Less than 25 common non-glycine residue records')
            A=np.array(numeric)
            rho=float(spearmanr(A[:,0],A[:,1]).statistic)
            C=np.array(calls,dtype=int)
            match=float(np.mean(C[:,0]==C[:,1]))
            cases.append({'pdb':pdb,'status':'PASS','n_residues_matched':len(common),
                  'n_non_glycine_compared':len(numeric),'spearman_absolute_sidechain_A2':rho,
                  'median_absolute_difference_sidechain_A2':float(np.median(np.abs(A[:,0]-A[:,1]))),
                  'exposure_agreement_at_0_05':match,'freesasa_total_A2':total,
                  'prepared_pdb_sha256':prov['protein_only_sha256'],
                  'source_url':prov['source_url'],'n_chains':prov['n_protein_chains'],
                  'n_residues_prepared':prov['n_modelled_residues'],
                  'note':'Agreement does not imply independently validated adsorption regions.'})
            print('FREE_SASA_CASE',pdb,'PASS',len(numeric),'rho',round(rho,4),flush=True)
        except Exception as exc:
            cases.append({'pdb':pdb,'status':'FAILED_OR_HELD','error':type(exc).__name__+': '+str(exc)})
            print('FREE_SASA_CASE',pdb,'FAILED',str(exc),flush=True)
    if records:
        with (OUTPUT/'per_residue_comparison.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    headers=sorted({k for r in cases for k in r})
    with (OUTPUT/'per_structure_summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=headers);w.writeheader();w.writerows(cases)
    (OUTPUT/'method.json').write_text(json.dumps({
       'purpose':'Independent method comparison; not adsorption site prediction validation',
       'freesasa_version':getattr(freesasa,'__version__','module installed'),
       'freesasa_algorithm':'Lee-Richards default',
       'interfacescout_algorithm':'Biopython Shrake-Rupley 200 points per atom, 1.4 A probe',
       'same_protein_coordinates':True,
       'do_not_equate_radii_or_normalizations':True,
       'excluded':'glycine sidechain in quantitative comparisons',
       'preset_cases':PDBS,'n_pass':sum(r['status']=='PASS' for r in cases)
    },indent=2),encoding='utf-8')
    if any(r['status']!='PASS' for r in cases):
        raise SystemExit('One or more cases failed; see per_structure_summary.csv')

if __name__=='__main__':main()
