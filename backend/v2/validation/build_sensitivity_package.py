from __future__ import annotations

import csv
import json
import shutil
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt

from v2.prepare import prepare_pdb_text
from v2.validation.evaluate_sensitivity_panel import STRUCTURES

ROOT = Path('InterfaceScout_sensitivity_package')
PDB_FULL = ROOT / 'PDB' / 'full_rcsb'
PDB_SELECTED = ROOT / 'PDB' / 'selected_for_analysis'
CSV_DIR = ROOT / 'CSV'
FIG_DIR = ROOT / 'Figures'
RAW_DIR = ROOT / 'Raw_JSON'


def ensure_dirs():
    for p in (PDB_FULL, PDB_SELECTED, CSV_DIR, FIG_DIR, RAW_DIR):
        p.mkdir(parents=True, exist_ok=True)


def download_pdbs():
    manifest = []
    for s in STRUCTURES:
        pid = s['pdb_id']
        url = f'https://files.rcsb.org/download/{pid}.pdb'
        with urllib.request.urlopen(url, timeout=60) as r:
            raw = r.read().decode('utf-8', errors='replace')
        full_path = PDB_FULL / f'{pid}.pdb'
        full_path.write_text(raw, encoding='utf-8')
        selected, report = prepare_pdb_text(raw, chain=s['chain'])
        sel_path = PDB_SELECTED / f"{s['id']}_{pid}_chains_{s['chain'].replace(',', '-')}.pdb"
        sel_path.write_text(selected, encoding='utf-8')
        manifest.append({
            'protein_id': s['id'], 'protein_label': s.get('label', s['id']),
            'group': s['group'], 'pdb_id': pid, 'selected_chains': s['chain'],
            'rcsb_url': url, 'full_pdb_file': str(full_path.relative_to(ROOT)),
            'analysis_pdb_file': str(sel_path.relative_to(ROOT)),
            'preparation_policy': report['policy'],
            'hetero_residues_removed': report['hetero_residues_removed'],
            'explicit_hydrogens_removed': report['explicit_hydrogens_removed'],
        })
    return manifest


def write_csv(path, rows, fieldnames):
    with path.open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)


def build_csvs(sens, radius, manifest):
    basic = []
    thresholds = []
    pka = []
    for r in sens['structures']:
        basic.append({
            'protein_id': r['id'], 'protein_label': r.get('label', r['id']), 'group': r['group'],
            'pdb_id': r['pdb_id'], 'chain': r['chain'], 'n_residues': r.get('n_residues'),
            'n_surface_200': r.get('n_surface_200'), 'n_surface_500': r.get('n_surface_500'),
            'surface_set_jaccard_200_vs_500': r.get('surface_set_jaccard_200_vs_500'),
            'surface_threshold_crossing_count_200_vs_500': r.get('surface_threshold_crossing_count_200_vs_500'),
            'scrsa_mae_200_vs_500': r.get('scrsa_mae_200_vs_500'),
            'scrsa_rmse_200_vs_500': r.get('scrsa_rmse_200_vs_500'), 'status': r.get('status'),
        })
        for thr, hv in r.get('median_top10_hotspot_jaccard_vs_threshold_005', {}).items():
            thresholds.append({
                'protein_id': r['id'], 'protein_label': r.get('label', r['id']), 'group': r['group'],
                'pdb_id': r['pdb_id'], 'scrsa_threshold': thr,
                'top10_hotspot_jaccard_vs_0.05': hv,
                'surface_set_jaccard_vs_0.05': r.get('surface_set_jaccard_vs_threshold_005', {}).get(thr),
            })
        for shift, v in r.get('state_sensitive_top10_jaccard_vs_pka_shift_0', {}).items():
            pka.append({
                'protein_id': r['id'], 'protein_label': r.get('label', r['id']), 'group': r['group'],
                'pdb_id': r['pdb_id'], 'uniform_pka_shift': shift,
                'state_sensitive_top10_hotspot_jaccard_vs_default': v,
            })

    radius_rows = []
    labels = {r['id']: r.get('label', r['id']) for r in sens['structures']}
    for r in radius['structures']:
        for pair, v in r.get('median_top10_hotspot_jaccard_vs_5_8', {}).items():
            radius_rows.append({
                'protein_id': r['id'], 'protein_label': labels.get(r['id'], r['id']), 'group': r['group'],
                'pdb_id': r['pdb_id'], 'radius_pair_A': pair,
                'top10_hotspot_jaccard_vs_5_8': v, 'n_surface': r.get('n_surface'), 'status': r.get('status'),
            })

    write_csv(CSV_DIR/'sasa_sampling_by_protein.csv', basic, list(basic[0].keys()))
    write_csv(CSV_DIR/'scrsa_threshold_sensitivity.csv', thresholds, list(thresholds[0].keys()))
    write_csv(CSV_DIR/'pka_sensitivity.csv', pka, list(pka[0].keys()))
    write_csv(CSV_DIR/'radius_pair_sensitivity.csv', radius_rows, list(radius_rows[0].keys()))
    write_csv(CSV_DIR/'pdb_manifest.csv', manifest, list(manifest[0].keys()))

    summary_rows = [
        {'analysis':'SASA sampling','metric':'median surface-set Jaccard, 200 vs 500','value':sens['summary']['median_surface_set_jaccard_200_vs_500']},
        {'analysis':'SASA sampling','metric':'median scRSA MAE, 200 vs 500','value':sens['summary']['median_scrsa_mae_200_vs_500']},
        {'analysis':'SASA sampling','metric':'median scRSA RMSE, 200 vs 500','value':sens['summary']['median_scrsa_rmse_200_vs_500']},
    ]
    for k,v in sens['summary']['median_top10_hotspot_jaccard_by_threshold'].items():
        summary_rows.append({'analysis':'scRSA threshold','metric':f'median top10 hotspot Jaccard at threshold {k} vs 0.05','value':v})
    for k,v in sens['summary']['median_state_sensitive_top10_jaccard_by_pka_shift'].items():
        summary_rows.append({'analysis':'pKa perturbation','metric':f'median state-sensitive top10 Jaccard, shift {k}','value':v})
    for k,v in radius['summary']['panel_median_top10_hotspot_jaccard_vs_5_8'].items():
        summary_rows.append({'analysis':'radius pair','metric':f'median top10 hotspot Jaccard, {k} vs 5/8 Å','value':v})
    write_csv(CSV_DIR/'panel_summary.csv', summary_rows, ['analysis','metric','value'])


def make_figures(sens, radius):
    ok = [r for r in sens['structures'] if r.get('status') == 'ok']
    labels = [r['id'] for r in ok]
    vals = [r['surface_set_jaccard_200_vs_500'] for r in ok]
    plt.figure(figsize=(10,5.5))
    plt.bar(labels, vals)
    plt.ylim(0.94,1.005); plt.ylabel('Surface-set Jaccard (200 vs 500)')
    plt.xlabel('Protein'); plt.title('SASA sampling sensitivity across the 14-protein panel')
    plt.xticks(rotation=55, ha='right'); plt.tight_layout()
    plt.savefig(FIG_DIR/'Figure_S1_SASA_sampling_sensitivity.png', dpi=300); plt.close()

    thr = sens['summary']['median_top10_hotspot_jaccard_by_threshold']
    x = [float(k) for k in thr.keys()]; y = [thr[k] for k in thr.keys()]
    plt.figure(figsize=(7,5))
    plt.plot(x, y, marker='o')
    plt.ylim(0,1.05); plt.xlabel('scRSA threshold'); plt.ylabel('Median top-10 hotspot Jaccard vs 0.05')
    plt.title('Sensitivity to the surface-exposure threshold'); plt.tight_layout()
    plt.savefig(FIG_DIR/'Figure_S2_scRSA_threshold_sensitivity.png', dpi=300); plt.close()

    pr = radius['summary']['panel_median_top10_hotspot_jaccard_vs_5_8']
    pairs = list(pr.keys()); y = [pr[k] for k in pairs]
    plt.figure(figsize=(7.5,5))
    plt.bar(pairs, y)
    plt.ylim(0,1.05); plt.xlabel('Radius pair (Å)'); plt.ylabel('Median top-10 hotspot Jaccard vs 5/8 Å')
    plt.title('Dependence on spatial aggregation scale'); plt.tight_layout()
    plt.savefig(FIG_DIR/'Figure_S3_radius_pair_sensitivity.png', dpi=300); plt.close()

    pka = sens['summary']['median_state_sensitive_top10_jaccard_by_pka_shift']
    shifts = [float(k) for k in pka.keys()]; y = [pka[k] for k in pka.keys()]
    plt.figure(figsize=(6.5,4.8))
    plt.bar([str(x) for x in shifts], y)
    plt.ylim(0,1.05); plt.xlabel('Uniform pKa shift'); plt.ylabel('Median top-10 hotspot Jaccard')
    plt.title('pKa perturbation sensitivity at pH 7.4'); plt.tight_layout()
    plt.savefig(FIG_DIR/'Figure_S4_pKa_sensitivity.png', dpi=300); plt.close()


def write_readme(sens, radius):
    txt = f"""# InterfaceScout 14-protein sensitivity package

This package contains the label-independent sensitivity analysis used to assess numerical and model-parameter robustness of InterfaceScout. No experimental adsorption patch labels were used in these analyses.

## Reference configuration
- pH: 7.4
- ionic strength: 150 mM
- SASA sampling: 200 points/atom
- scRSA threshold: 0.05
- spatial radii: 5 and 8 Å

## Panel
14 proteins: albumin, fibrinogen, IgG, transferrin, lysozyme, RNase A, myoglobin, carbonic anhydrase II, ubiquitin, acylphosphatase, catalase, L-asparaginase II, beta2-microglobulin, and cytochrome c.

Fibrinogen was represented by an alpha-beta-gamma structural half-unit from PDB 3GHG (chains A/B/C). Intact IgG2a used PDB 1IGT chains A/B/C/D. Catalase and asparaginase were analyzed as tetramers.

## Main panel-level results
- Successful structures: {sens['summary']['n_successful']}/{sens['summary']['n_requested']}
- Median surface-set Jaccard, 200 vs 500 SASA points/atom: {sens['summary']['median_surface_set_jaccard_200_vs_500']:.4f}
- Median scRSA MAE, 200 vs 500: {sens['summary']['median_scrsa_mae_200_vs_500']:.5f}
- Median scRSA RMSE, 200 vs 500: {sens['summary']['median_scrsa_rmse_200_vs_500']:.5f}
- Median top-10 hotspot Jaccard at scRSA thresholds 0.02, 0.05, 0.075, 0.10, 0.15: {sens['summary']['median_top10_hotspot_jaccard_by_threshold']}
- Median state-sensitive top-10 hotspot Jaccard after uniform pKa shifts of -0.5/+0.5: {sens['summary']['median_state_sensitive_top10_jaccard_by_pka_shift']}
- Radius-pair median top-10 hotspot Jaccard vs 5/8 Å: {radius['summary']['panel_median_top10_hotspot_jaccard_vs_5_8']}

## Directory contents
- PDB/full_rcsb: original structures downloaded from RCSB PDB.
- PDB/selected_for_analysis: exact chain-selected, first-model, heavy-atom protein structures used by the analysis; hetero residues and explicit hydrogens removed according to InterfaceScout preprocessing.
- CSV: flattened analysis tables and PDB manifest.
- Figures: publication-ready PNG summaries (300 dpi).
- Raw_JSON: complete machine-readable outputs.

## Interpretation
SASA sampling, scRSA threshold over 0.02-0.10, and +/-0.5 pKa perturbations were robust across the panel. Spatial radius choice showed genuine scale dependence and should therefore be treated as part of the model definition rather than as a numerical convergence parameter.
"""
    (ROOT/'README.md').write_text(txt, encoding='utf-8')


def main():
    ensure_dirs()
    sens = json.loads(Path('sensitivity_panel.json').read_text(encoding='utf-8'))
    radius = json.loads(Path('radius_sensitivity_panel.json').read_text(encoding='utf-8'))
    shutil.copy2('sensitivity_panel.json', RAW_DIR/'sensitivity_panel.json')
    shutil.copy2('radius_sensitivity_panel.json', RAW_DIR/'radius_sensitivity_panel.json')
    manifest = download_pdbs()
    build_csvs(sens, radius, manifest)
    make_figures(sens, radius)
    write_readme(sens, radius)
    shutil.make_archive('InterfaceScout_14protein_sensitivity_package', 'zip', ROOT)
    print('Built InterfaceScout_14protein_sensitivity_package.zip')

if __name__ == '__main__':
    main()
