"""Canonical adsorption-label-free parameter-selection run.

The publication development panel contains eight structurally diverse proteins
and is independent of the external experimental adsorption benchmark.
"""
from v2.validation import select_core_parameters_consensus as base

base.STRUCTURES = [
    {"id":"crambin","label":"Crambin","pdb_id":"1CRN","chain":"A","class":"very small compact"},
    {"id":"repressor434","label":"Phage 434 repressor N-terminal domain","pdb_id":"1R69","chain":"A","class":"small all-alpha"},
    {"id":"sh3","label":"Alpha-spectrin SH3 domain","pdb_id":"1SHG","chain":"A","class":"small all-beta"},
    {"id":"fkbp12","label":"Native FKBP12","pdb_id":"2PPN","chain":"A","class":"small alpha-beta"},
    {"id":"adenylate_kinase","label":"Unligated adenylate kinase","pdb_id":"4AKE","chain":"A","class":"medium alpha-beta enzyme"},
    {"id":"maltose_binding_protein","label":"Unliganded maltodextrin-binding protein","pdb_id":"1OMP","chain":"A","class":"large multidomain"},
    {"id":"tim_dimer","label":"Triosephosphate isomerase dimer","pdb_id":"1TIM","chain":"A,B","class":"dimeric enzyme"},
    {"id":"citrate_synthase","label":"Open citrate synthase","pdb_id":"5CSC","chain":"B","class":"large enzyme"},
]

if __name__ == '__main__':
    base.main()
