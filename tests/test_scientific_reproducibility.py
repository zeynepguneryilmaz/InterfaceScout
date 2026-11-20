"""Analytic expectations, boundaries and adversarial checks, not golden outputs."""
from __future__ import annotations
import copy
import json
import math
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'backend'),str(ROOT/'validation'),str(ROOT/'validation/reproducibility')]
import core
from coarse_patch import build_coarse_patches, assign_pareto_fronts
from interface_engine import analyze_all_maps
from oracle import accessibility, Checks, dominance_depths, verify_result, scalar_channel
from reproduce import compute, checked_input


def row(seq,name='ALA',scrsa=0.5,x=0.,y=0.):
    return dict(key=f'A:{seq}:',chain='A',res_seq=seq,icode='',res_name=name,
                x=x,y=y,z=0.,coord=(x,y,0.),scrsa=scrsa,scrsa_raw=scrsa,
                sidechain_sasa=scrsa*core.SIDECHAIN_REF_ASA[name],ionization_sensitive=False,bfactor=20.)


class AnalyticNumericalTests(unittest.TestCase):
    def test_hand_calculated_residue_and_density_scores(self):
        # pH=pKa gives 1/2 availability; isolated residues have only self density.
        rows=[row(1,'HIS',0.2,x=0),row(2,'HIS',0.6,x=20)]
        actual=core.chemistry_map(rows,np.array([[0.,20.],[20.,0.]]),'anionic',6.5)
        values={r['key']:r for r in actual['residues']}
        self.assertAlmostEqual(values['A:1:']['local_score'],0.1,places=14)
        self.assertAlmostEqual(values['A:2:']['local_score'],0.3,places=14)
        self.assertAlmostEqual(values['A:1:']['propensity'],100/3,places=12)
        self.assertAlmostEqual(values['A:1:']['multiscale_persistence'],100/3,places=12)
        self.assertEqual(actual['top_patch']['center_key'],'A:2:')
        oracle,_=scalar_channel(rows,'anionic',6.5)
        self.assertAlmostEqual(oracle['A:1:']['local_score'],0.1,places=14)

    def test_density_radius_boundaries_are_inclusive(self):
        rows=[row(1,x=0),row(2,x=6),row(3,x=9),row(4,x=9.000001)]
        distances=core.build_distances(rows)
        actual=core.chemistry_map(rows,distances,'hydrophobic',7.4)
        center=next(c for c in actual['patch_centers'] if c['center_key']=='A:1:')
        self.assertEqual(center['density_6A_raw'],1.)
        self.assertEqual(center['density_9A_raw'],1.5)

    def test_patch_radius_boundary_and_hand_calculated_descriptors(self):
        rows=[row(1,x=10),row(2,x=18),row(3,x=18.000001)]
        nodes=[dict(key=r['key'],chain='A',res_seq=r['res_seq'],icode='',res_name='ALA',coord=np.array(r['coord'])) for r in rows]
        nodes.append(dict(key='A:99:',chain='A',res_seq=99,icode='',res_name='ALA',coord=np.array([-100.,0.,0.])))
        channel=dict(residues=[dict(r,local_score=0.5) for r in rows],
                     patch_centers=[dict(center_key='A:1:',multiscale_persistence=75.,multiscale_geomean=80.)])
        result=build_coarse_patches(v1_result=dict(surface_residues=rows,chemistries={'hydrophobic':channel}),chemistry='hydrophobic',ca_nodes=nodes)
        p=result[0]
        self.assertEqual(p['members'],['A:1:','A:2:'])
        self.assertEqual(p['seed_members'],p['members'])
        for field,value in [('diameter_A',8.),('chemistry_support',1.),('mean_accessibility',0.5),('patch_coherence',0.75),('orientation_coherence',1.)]:
            self.assertAlmostEqual(p[field],value,places=12)

    def test_perpendicular_faces_are_excluded(self):
        rows=[row(1,x=1),row(2,x=0,y=1)]
        nodes=[dict(key=r['key'],chain='A',res_seq=r['res_seq'],icode='',res_name='ALA',coord=np.array(r['coord'])) for r in rows]
        nodes.append(dict(key='A:99:',chain='A',res_seq=99,icode='',res_name='ALA',coord=np.array([-1.,-1.,0.])))
        channel=dict(residues=[dict(r,local_score=0.5) for r in rows],patch_centers=[dict(center_key='A:1:',multiscale_persistence=100.,multiscale_geomean=100.)])
        result=build_coarse_patches(v1_result=dict(surface_residues=rows,chemistries={'hydrophobic':channel}),chemistry='hydrophobic',ca_nodes=nodes)
        self.assertEqual(result[0]['members'],['A:1:'])

    def test_surface_threshold_is_not_a_verification_tolerance(self):
        pdb='ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nATOM      2  CB  ALA A   1       1.500   0.000   0.000  1.00 20.00           C  \nEND\n'
        for value,expected in [(0.05,True),(0.05-1e-12,False)]:
            with self.subTest(scrsa=value):
                def calculate(structure,level):
                    for atom in structure.get_atoms():atom.sasa=value*core.SIDECHAIN_REF_ASA['ALA'] if atom.name=='CB' else 0.
                with patch.object(core,'ShrakeRupley') as sr:
                    sr.return_value.compute.side_effect=calculate
                    self.assertEqual(core._surface_residues(pdb,7.4)[0]['surface_exposed'],expected)

    def test_analytic_isolated_atom_and_buried_atom_sasa(self):
        isolated='ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n'
        area=4*math.pi*(1.7+1.4)**2
        for independent in (False,True):
            r=accessibility(isolated,independent)[0]
            self.assertAlmostEqual(r['sidechain_sasa'],area,places=11)
            self.assertEqual(r['scrsa'],1.)
        buried='ATOM      1  CA  CYS A   1      20.000   0.000   0.000  1.00 20.00           C  \nATOM      2  CB  CYS A   1       0.000   0.000   0.000  1.00 20.00           C  \nATOM      3  SG  CYS A   1       0.000   0.000   0.000  1.00 20.00           S  \nEND\n'
        self.assertAlmostEqual(accessibility(buried,True)[0]['sidechain_sasa'],4*math.pi*3.2**2,places=11)

    def test_empty_channels_remain_empty_and_finite(self):
        raw='ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n'
        result=analyze_all_maps(pdb_text=raw)
        self.assertEqual(result['n_maps'],11)
        self.assertTrue(all(m['patches']==[] and m['residues']==[] for m in result['maps'].values()))
        json.dumps(result,allow_nan=False)

    def test_invalid_environment_and_structure_fail_explicitly(self):
        raw='ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n'
        for kwargs in [dict(pH=-0.01),dict(pH=14.01),dict(ionic_mM=-1),dict(temp_K=0),dict(chain='Z')]:
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):analyze_all_maps(pdb_text=raw,**kwargs)
        with self.assertRaises(ValueError):analyze_all_maps(pdb_text='END\n')

    def test_pareto_fronts_and_ties_have_analytic_order(self):
        ps=[dict(center_key=k,chemistry_support=a,mean_accessibility=b,patch_coherence=c,orientation_coherence=d)
            for k,a,b,c,d in [('A:1:',1.,1.,1.,1.),('A:2:',.5,.5,.5,.5),('A:3:',1.,.2,.8,.9),('A:4:',.1,.1,.1,.1),('A:0:',1.,1.,1.,1.)]]
        result=assign_pareto_fronts(copy.deepcopy(ps))
        self.assertEqual([(p['center_key'],p['pareto_front']) for p in result],[('A:0:',1),('A:1:',1),('A:3:',2),('A:2:',2),('A:4:',3)])
        self.assertEqual(dominance_depths({p['center_key']:p for p in ps}),{'A:1:':1,'A:2:':2,'A:3:':2,'A:4:':3,'A:0:':1})


class VerifierMutationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cases=json.loads((ROOT/'validation/reproducibility/inputs.json').read_text())['cases']
        cls.case=next(c for c in cases if c['id']=='raw_1CRN')
        cls.result,cls.core_result=compute(cls.case)

    def test_unmodified_raw_input_passes_independent_checks(self):
        self.assertGreater(verify_result(self.result,self.core_result,True)['checks'],1000)

    def test_input_hash_is_enforced(self):
        with self.assertRaisesRegex(ValueError,'SHA-256'):checked_input(dict(self.case,sha256='0'*64))

    def test_every_versioned_input_matches_manifest_bytes(self):
        cases=json.loads((ROOT/'validation/reproducibility/inputs.json').read_text())['cases']
        self.assertEqual(len(cases),21)
        for case in cases:
            with self.subTest(case=case['id']):self.assertTrue(checked_input(case))

    def test_corrupted_scores_membership_and_order_are_rejected(self):
        for kind in ['residue_score','patch_score','membership','ranking','front','nan','map_count']:
            with self.subTest(kind=kind):
                result=copy.deepcopy(self.result); cr=copy.deepcopy(self.core_result)
                m=result['maps']['anionic']
                if kind=='residue_score':m['residues'][0]['local_score']+=1e-3
                elif kind=='patch_score':m['patches'][0]['mean_accessibility']+=1e-3
                elif kind=='membership':m['patches'][0]['members'].pop()
                elif kind=='ranking':m['patches'][0],m['patches'][1]=m['patches'][1],m['patches'][0]
                elif kind=='front':m['patches'][0]['pareto_front']+=1
                elif kind=='nan':m['residues'][0]['propensity']=float('nan')
                else:result['n_maps']=10
                with self.assertRaises(AssertionError):verify_result(result,cr)

    def test_positive_scores_cannot_be_erased_by_absolute_tolerance(self):
        with self.assertRaises(AssertionError):Checks().close(float('nan'),0.,'nonfinite')
        # The separate positive-score invariant is needed below the absolute tolerance.
        from oracle import verify_map
        rows=[row(1,'HIS',1e-5,x=1)]
        channel=core.chemistry_map(rows,np.zeros((1,1)),'anionic',14.)
        channel['residues'][0]['local_score']=0.
        with self.assertRaisesRegex(AssertionError,'positive local score'):
            verify_map(rows,rows,channel,{},'anionic',14.,Checks())


if __name__=='__main__':unittest.main()
