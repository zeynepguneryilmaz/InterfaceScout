"""Generate publication-frozen InterfaceScout outputs for the BSA manuscript cases."""
from pathlib import Path
import json

from v2.validation.generate_full_candidate_outputs import process_case, _write_csv
from v2.chemistry_freeze import apply_publication_chemistry

OUT = Path('bsa_manuscript_outputs')


def main():
    import main as v1
    apply_publication_chemistry(v1)
    OUT.mkdir(parents=True, exist_ok=True)
    cases = [
        {'id':'BSA_pH4p7','protein':'Bovine serum albumin','pdb_id':'4F5U','chain':'A','surface':'hydrophobic','pH':4.7},
        {'id':'BSA_pH7p4','protein':'Bovine serum albumin','pdb_id':'4F5U','chain':'A','surface':'hydrophobic','pH':7.4},
    ]
    # Reuse the full exporter by temporarily pointing its output root here.
    import v2.validation.generate_full_candidate_outputs as exp
    exp.OUTROOT = OUT
    summary=[]
    for case in cases:
        print('BSA_OUTPUT', case['id'], flush=True)
        summary.append(process_case(case, v1))
    _write_csv(OUT/'panel_summary.csv', summary)
    (OUT/'manifest.json').write_text(json.dumps({'cases':cases}, indent=2), encoding='utf-8')
    print('BSA_SUMMARY', json.dumps(summary), flush=True)

if __name__ == '__main__':
    main()
