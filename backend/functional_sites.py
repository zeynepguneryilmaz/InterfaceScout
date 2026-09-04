"""Functional-site annotation helpers for InterfaceScout 2.0.

SITE records from the deposited PDB are treated only as generic structural-site
annotations. They are never automatically labelled catalytic/active sites.
Users can additionally provide protected residues explicitly; those annotations
are carried into patch-overlap reporting without affecting compatibility scores.
"""
from __future__ import annotations

import re
import urllib.request
from typing import Any, Dict, List, Optional


def parse_residue_key(text: str) -> Optional[str]:
    """Normalize user residue notation to CHAIN:SEQ:ICODE.

    Accepted examples: A:42, A:42:, A:42:B.  Returns None for invalid input.
    """
    s = str(text or "").strip()
    m = re.fullmatch(r"([^:]):(-?\d+)(?::([^:]*))?", s)
    if not m:
        return None
    return f"{m.group(1)}:{int(m.group(2))}:{(m.group(3) or '').strip()}"


def normalize_protected_residues(values: Optional[List[str]]) -> List[str]:
    out: List[str] = []
    for value in values or []:
        key = parse_residue_key(value)
        if key and key not in out:
            out.append(key)
    return out


def _source_text(pdb_id: Optional[str], pdb_text: Optional[str]) -> str:
    if pdb_text:
        return pdb_text
    if pdb_id:
        try:
            with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id.strip().upper()}.pdb", timeout=30) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception:
            return ""
    return ""


def pdb_site_annotations(pdb_id: Optional[str], pdb_text: Optional[str]) -> List[Dict[str, Any]]:
    """Parse legacy PDB SITE records without assigning biological function."""
    text = _source_text(pdb_id, pdb_text)
    sites: Dict[str, Dict[str, Any]] = {}
    for line in text.splitlines():
        if not line.startswith("SITE"):
            continue
        site_id = line[11:14].strip() or "SITE"
        row = sites.setdefault(site_id, {"site_id": site_id, "residues": []})
        # SITE records contain up to four residue fields, 11 chars each,
        # starting at column 19 (1-based). Parse by fixed-width layout.
        for start in (18, 29, 40, 51):
            field = line[start:start+10]
            if not field.strip():
                continue
            res_name = field[0:3].strip()
            chain = field[4:5].strip()
            seq_txt = field[5:9].strip()
            icode = field[9:10].strip()
            if not res_name or not chain or not seq_txt:
                continue
            try:
                seq = int(seq_txt)
            except ValueError:
                continue
            item = {"res_name": res_name, "chain": chain, "res_seq": seq, "icode": icode}
            key = (res_name, chain, seq, icode)
            if key not in {(r["res_name"], r["chain"], r["res_seq"], r.get("icode", "")) for r in row["residues"]}:
                row["residues"].append(item)
    return list(sites.values())
