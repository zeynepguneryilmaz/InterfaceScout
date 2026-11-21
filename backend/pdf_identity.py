"""Conservative PDB metadata extraction for report identity.

Only original PDB HEADER, COMPND and TITLE records are evidence.
File names or B-factor exports are never used to infer a protein name.
"""
from __future__ import annotations
import re
from typing import Mapping, Any
from datetime import datetime, timezone

def pdb_identity(pdb_text: str | None, analysis: Mapping[str, Any], source_name: str) -> dict[str,str]:
    raw=str(pdb_text or "")
    lines=raw.splitlines()
    header_id=""
    for line in lines:
        if line.startswith("HEADER"):
            candidate=line[62:66].strip().upper() if len(line)>=66 else ""
            if re.fullmatch(r"[0-9][A-Z0-9]{3}",candidate):
                header_id=candidate
            break
    comp_parts=[line[10:].strip() for line in lines if line.startswith("COMPND")]
    comp=" ".join(re.sub(r"^\d+\s+","",x) for x in comp_parts)
    match=re.search(r"\bMOLECULE\s*:\s*([^;]+)",comp,re.I)
    molecule=re.sub(r"\s+"," ",match.group(1)).strip(" .") if match else ""
    if len(molecule)>95:
        molecule=molecule[:92].rstrip()+"..."
    supplied=str(analysis.get("input",{}).get("pdb_id") or "").upper()
    pdbid=header_id or (supplied if re.fullmatch(r"[0-9][A-Z0-9]{3}",supplied) else "")
    label=molecule or "Protein name not provided in original PDB"
    chain=str(analysis.get("input",{}).get("chain") or "ALL")
    return {
        "name":label,
        "accession":pdbid or "not verified",
        "accession_evidence":"original PDB HEADER" if header_id else ("requested PDB ID" if pdbid else "not available"),
        "name_evidence":"original PDB COMPND" if molecule else "not available",
        "source_label":source_name,
        "chain":chain,
        "source_sha256":str(analysis.get("input",{}).get("source_pdb_sha256") or ""),
        "model_version":str(analysis.get("version") or "unrecorded"),
        "export_timestamp_utc":datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    }
