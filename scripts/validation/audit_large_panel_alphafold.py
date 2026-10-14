#!/usr/bin/env python3
"""Uniform AlphaFold DB audit for the frozen InterfaceScout 285-protein panel.

Purpose
-------
Resolve the frozen source names to reviewed human UniProtKB accessions and query
AlphaFold DB for the corresponding predicted monomer models. This script does
NOT use InterfaceScout scores, adsorption labels, or experimental PDB
availability to include/exclude proteins.

AlphaFold confidence is recorded during the later batch run as QC metadata. It is not used as a\nweight in InterfaceScout V1.
"""
from __future__ import annotations

import csv
import json
import re
import sys
import time
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "validation/large_panel_285_2026-10-13/source_285_master.csv"
OUT = ROOT / "validation/large_panel_285_2026-10-13/alphafold_structure_audit.csv"
SUMMARY = ROOT / "validation/large_panel_285_2026-10-13/alphafold_structure_audit_summary.json"

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb"
AF_API = "https://alphafold.ebi.ac.uk/api/prediction"

CURATED_SOURCE_ALIASES = {
    "Serum albumin": "P02768",
    "Insulin-like growth factor-binding protein complex acid labile chain": "P35858",
    "Hepatocyte growth factor activator": "Q04756",
    "78kDa glucose-regulated protein": "P11021",
    "Inactive cell surface hyaluronidase CEMIP2": "Q9UHN6",
    "Coiled coil domain containing protein 144B": "Q3MJ40",
}

session = requests.Session()
retry = Retry(
    total=5, connect=5, read=5, status=5,
    backoff_factor=0.6,
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=frozenset(["GET"]),
)
session.mount("https://", HTTPAdapter(max_retries=retry))
session.headers.update({"User-Agent": "InterfaceScout-AF-audit/2026-10-14"})


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())


def uniprot_names(rec: dict[str, Any]) -> list[str]:
    pd = rec.get("proteinDescription") or {}
    names: list[str] = []
    rn = ((pd.get("recommendedName") or {}).get("fullName") or {}).get("value")
    if rn:
        names.append(rn)
    for k in ("alternativeNames", "submissionNames"):
        for x in pd.get(k) or []:
            v = ((x or {}).get("fullName") or {}).get("value")
            if v:
                names.append(v)
    return names


def fetch_accession(accession: str) -> dict[str, Any]:
    r = session.get(f"{UNIPROT_ENTRY}/{accession}.json", timeout=45)
    r.raise_for_status()
    rec = r.json()
    if ((rec.get("organism") or {}).get("taxonId")) != 9606:
        raise RuntimeError(f"{accession} is not human")
    if "reviewed" not in (rec.get("entryType") or "").lower():
        raise RuntimeError(f"{accession} is not reviewed UniProtKB")
    return rec


def resolve_identity(source_name: str):
    if source_name in CURATED_SOURCE_ALIASES:
        rec = fetch_accession(CURATED_SOURCE_ALIASES[source_name])
        names = uniprot_names(rec)
        return (
            "RESOLVED_CURATED_ALIAS_REVIEWED_HUMAN",
            rec["primaryAccession"],
            rec.get("uniProtkbId", ""),
            names[0] if names else "",
            int((rec.get("sequence") or {}).get("length") or 0),
        )

    queries = [
        f'(protein_name:"{source_name.replace(chr(34), "")}") AND (organism_id:9606) AND (reviewed:true)',
        f'("{source_name.replace(chr(34), "")}") AND (organism_id:9606) AND (reviewed:true)',
    ]
    records = []
    for q in queries:
        r = session.get(
            UNIPROT_SEARCH,
            params={"query": q, "format": "json", "size": 25},
            timeout=45,
        )
        r.raise_for_status()
        records = r.json().get("results", [])
        if records:
            break

    target = norm(source_name)
    exact: dict[str, dict[str, Any]] = {}
    for rec in records:
        acc = rec.get("primaryAccession")
        if acc and any(norm(n) == target for n in uniprot_names(rec)):
            exact[acc] = rec

    if len(exact) == 1:
        rec = next(iter(exact.values()))
        names = uniprot_names(rec)
        return (
            "RESOLVED_EXACT_REVIEWED_HUMAN",
            rec["primaryAccession"],
            rec.get("uniProtkbId", ""),
            names[0] if names else "",
            int((rec.get("sequence") or {}).get("length") or 0),
        )
    if len(exact) > 1:
        return ("REVIEW_REQUIRED_MULTIPLE_EXACT", "", "", "", 0)
    return ("REVIEW_REQUIRED_NO_EXACT", "", "", "", 0)


def af_records(accession: str) -> list[dict[str, Any]]:
    r = session.get(f"{AF_API}/{accession}", timeout=60)
    if r.status_code == 404:
        return []
    r.raise_for_status()
    payload = r.json()
    if isinstance(payload, dict):
        payload = [payload]
    return payload if isinstance(payload, list) else []


def choose_primary(records: list[dict[str, Any]], uniprot_length: int):
    """Choose the model record with maximal sequence coverage, deterministic ties."""
    parsed = []
    for rec in records:
        s = int(rec.get("sequenceStart") or 1)
        e = int(rec.get("sequenceEnd") or uniprot_length or 0)
        span = max(0, e - s + 1)
        cov = span / uniprot_length if uniprot_length else 0.0
        parsed.append((cov, span, -s, str(rec.get("entryId") or ""), rec))
    parsed.sort(reverse=True)
    return parsed[0][-1] if parsed else None


def main():
    with PANEL.open(newline="", encoding="utf-8-sig") as fh:
        source = list(csv.DictReader(fh))
    if len(source) != 285:
        raise RuntimeError(f"Expected frozen n=285, got {len(source)}")

    fields = [
        "panel_id","protein_name","source_detected_count","source_detection_frequency",
        "identity_status","uniprot_accession","uniprot_id","uniprot_recommended_name",
        "uniprot_length","alphafold_status","alphafold_record_count","af_entry_id",
        "af_sequence_start","af_sequence_end","af_coverage_fraction","af_pdb_url",
        "af_cif_url","af_pae_doc_url","af_latest_version","af_is_reviewed",
        "af_is_reference_proteome","n_ca_confidence","mean_plddt","median_plddt",
        "frac_plddt_ge90","frac_plddt_ge70","frac_plddt_lt50",
        "confidence_qc_flag","primary_large_panel_status","notes",
    ]

    out = []
    for i, src in enumerate(source, 1):
        name = src["protein_name"].strip()
        print(f"[{i:03d}/285] {name}", flush=True)
        try:
            istatus, acc, uid, recname, ulen = resolve_identity(name)
        except Exception as exc:
            print(f"WARN identity {name}: {exc}", file=sys.stderr)
            istatus, acc, uid, recname, ulen = "LOOKUP_ERROR", "", "", "", 0

        row = {
            "panel_id": src["panel_id"],
            "protein_name": name,
            "source_detected_count": src["source_detected_count"],
            "source_detection_frequency": src["source_detection_frequency"],
            "identity_status": istatus,
            "uniprot_accession": acc,
            "uniprot_id": uid,
            "uniprot_recommended_name": recname,
            "uniprot_length": ulen or "",
            "alphafold_status": "NOT_QUERIED_UNRESOLVED_IDENTITY",
            "alphafold_record_count": 0,
            "af_entry_id": "",
            "af_sequence_start": "",
            "af_sequence_end": "",
            "af_coverage_fraction": "",
            "af_pdb_url": "",
            "af_cif_url": "",
            "af_pae_doc_url": "",
            "af_latest_version": "",
            "af_is_reviewed": "",
            "af_is_reference_proteome": "",
            "n_ca_confidence": "",
            "mean_plddt": "",
            "median_plddt": "",
            "frac_plddt_ge90": "",
            "frac_plddt_ge70": "",
            "frac_plddt_lt50": "",
            "confidence_qc_flag": "UNAVAILABLE",
            "primary_large_panel_status": "PENDING",
            "notes": "",
        }

        if acc:
            try:
                records = af_records(acc)
                row["alphafold_record_count"] = len(records)
                if not records:
                    row["alphafold_status"] = "NO_ALPHAFOLD_DB_MODEL"
                    row["primary_large_panel_status"] = "NO_MODEL"
                else:
                    primary = choose_primary(records, ulen)
                    s = int(primary.get("sequenceStart") or 1)
                    e = int(primary.get("sequenceEnd") or ulen or 0)
                    span = max(0, e - s + 1)
                    cov = span / ulen if ulen else 0.0
                    row.update({
                        "alphafold_status": "MODEL_AVAILABLE",
                        "af_entry_id": primary.get("entryId", ""),
                        "af_sequence_start": s,
                        "af_sequence_end": e,
                        "af_coverage_fraction": f"{cov:.6f}",
                        "af_pdb_url": primary.get("pdbUrl", ""),
                        "af_cif_url": primary.get("cifUrl", ""),
                        "af_pae_doc_url": primary.get("paeDocUrl", ""),
                        "af_latest_version": primary.get("latestVersion", ""),
                        "af_is_reviewed": primary.get("isReviewed", ""),
                        "af_is_reference_proteome": primary.get("isReferenceProteome", ""),
                        "primary_large_panel_status": "MODEL_AVAILABLE_FULL_OR_NEAR_FULL" if cov >= 0.95 else "MODEL_AVAILABLE_FRAGMENT_OR_PARTIAL",
                    })
                    row["confidence_qc_flag"] = "DEFERRED_TO_BATCH"
                    row["notes"] = "pLDDT/PAE QC deferred to batch run to avoid downloading each model twice."
            except Exception as exc:
                print(f"WARN AlphaFold {name} {acc}: {exc}", file=sys.stderr)
                row["alphafold_status"] = "ALPHAFOLD_LOOKUP_ERROR"
                row["primary_large_panel_status"] = "LOOKUP_ERROR"
                row["notes"] = str(exc)

        out.append(row)
        time.sleep(0.03)

    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(out)

    def counts(col):
        d = {}
        for r in out:
            x = r[col] or "EMPTY"
            d[x] = d.get(x, 0) + 1
        return dict(sorted(d.items()))

    summary = {
        "source_panel_n": 285,
        "identity_counts": counts("identity_status"),
        "alphafold_counts": counts("alphafold_status"),
        "primary_large_panel_status_counts": counts("primary_large_panel_status"),
        "confidence_qc_counts": counts("confidence_qc_flag"),
        "reference_condition": {"pH": 7.4, "ionic_mM": 150.0, "temp_K": 298.0},
        "confidence_policy": "AlphaFold confidence is descriptive QC metadata only; it is not an InterfaceScout V1 weight or fitted threshold.",
        "selection_policy": "No protein is removed on the basis of InterfaceScout output or PC-DB adsorption label.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("IDENTITY_COUNTS=" + json.dumps(summary["identity_counts"], sort_keys=True))
    print("ALPHAFOLD_COUNTS=" + json.dumps(summary["alphafold_counts"], sort_keys=True))
    print("PRIMARY_STATUS_COUNTS=" + json.dumps(summary["primary_large_panel_status_counts"], sort_keys=True))
    print("CONFIDENCE_QC_COUNTS=" + json.dumps(summary["confidence_qc_counts"], sort_keys=True))
    print(f"WROTE {OUT}")
    print(f"WROTE {SUMMARY}")


if __name__ == "__main__":
    main()
