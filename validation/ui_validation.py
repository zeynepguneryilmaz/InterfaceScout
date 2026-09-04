#!/usr/bin/env python3
"""Static validation for the InterfaceScout 2.0 protein-centered UI."""
from __future__ import annotations
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "backend"
sys.path.insert(0, str(BACKEND))
from ui import inject_ui  # noqa: E402


def main():
    base = (ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
    html = inject_ui(base)
    tests = {
        "structure context control": 'id="structureContext"' in html and 'biological_assembly_1' in html,
        "protrusion control": 'id="protrusion"' in html and 'CX protrusion' in html,
        "protected residues control": 'id="protectedResidues"' in html and 'PROTECTED FUNCTIONAL RESIDUES' in html,
        "no material profile control": 'id="materialProfile"' not in html and '/material_profiles' not in html,
        "analysis sends structure context": "structure_context:document.getElementById('structureContext').value" in html,
        "analysis sends protrusion": "protrusion:document.getElementById('protrusion').checked" in html,
        "analysis sends protected residues": "protected_residue_keys:parseProtectedResidues()" in html,
        "target profile panel": 'id="targetProfilePanel"' in html and 'renderTargetProfile()' in html,
        "applicability panel": 'id="applicabilityPanel"' in html and 'renderApplicability()' in html,
        "protein centered wording": 'protein-centered target interface profiling' in html and 'does not use a named-material library' in html,
        "frozen equations preserved": 'L<sub>i,c</sub> = I<sub>i,c</sub> × scRSA<sub>i</sub> × f<sub>state,i,c</sub>(pH)' in html and 'M<sub>i,c</sub> = 100 × min' in html,
    }
    report = {"status": "PASS" if all(tests.values()) else "FAIL", "checks_passed": sum(tests.values()), "checks_total": len(tests), "checks": tests}
    target = Path(__file__).with_name("ui_validation_report.json")
    target.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
