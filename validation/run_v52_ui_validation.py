#!/usr/bin/env python3
"""Static validation of the lightweight v5.2 UI augmentation."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "backend"
sys.path.insert(0, str(BACKEND))

from v52_ui import inject_v52_ui  # noqa: E402


def main():
    base = (ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
    html = inject_v52_ui(base)
    tests = {
        "structure context control": 'id="structureContext"' in html and 'biological_assembly_1' in html,
        "protrusion control": 'id="protrusion"' in html and 'CX protrusion' in html,
        "material profile control": 'id="materialProfile"' in html and '/material_profiles' in html,
        "analysis sends structure context": 'structure_context:document.getElementById(\'structureContext\').value' in html,
        "analysis sends protrusion": 'protrusion:document.getElementById(\'protrusion\').checked' in html,
        "analysis sends material profile": 'material_profile:document.getElementById(\'materialProfile\').value||null' in html,
        "analysis preserves PDB ID": 'pdb_id:loadedPdbId' in html,
        "applicability panel": 'id="applicabilityPanel"' in html and 'renderApplicability()' in html,
        "CSV context fields": "'structure_context','material_profile'" in html,
        "CSV CX fields": "'cx_residue_mean','cx_sidechain_mean','cx_max','cx_ca'" in html,
        "theory candidate wording": 'v5.2 structural candidate' in html and 'Pintar CX protrusion' in html,
        "frozen equations preserved": 'L<sub>i,c</sub> = I<sub>i,c</sub> × scRSA<sub>i</sub> × f<sub>state,i,c</sub>(pH)' in html and 'M<sub>i,c</sub> = 100 × min' in html,
    }
    report = {
        "status": "PASS" if all(tests.values()) else "FAIL",
        "checks_passed": sum(tests.values()),
        "checks_total": len(tests),
        "checks": tests,
    }
    target = Path(__file__).with_name("v52_ui_validation_report.json")
    target.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
