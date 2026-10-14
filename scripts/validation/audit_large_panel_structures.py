#!/usr/bin/env python3
"""Conservative automated identity + experimental-PDB audit for the frozen 285-protein panel.

This script NEVER changes panel membership and NEVER scores InterfaceScout.
It resolves only unique exact-name reviewed human UniProtKB matches, then queries
RCSB PDB for experimental polymer entities linked to that accession. PDB hits are
candidates for manual construct/coverage review, not final structure selections.
"""
from __future__ import annotations

import csv
import json
import math
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
OUT = ROOT / "validation/large_panel_285_2026-10-13/structure_audit_auto.csv"

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core"

session = requests.Session()
retry = Retry(
    total=5, connect=5, read=5, status=5,
    backoff_factor=0.5,
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=frozenset(["GET", "POST"]),
)
session.mount("https://", HTTPAdapter(max_retries=retry))
session.headers.update({"User-Agent": "InterfaceScout-structure-audit/2026-10-13"})


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


def search_uniprot(source_name: str) -> tuple[str, str, str, str, int | None, list[str]]:
    # First use protein_name-scoped query; if no results, broaden but still require
    # an exact normalized protein name before auto-resolution.
    queries = [
        f'(protein_name:"{source_name.replace(chr(34), "")}") AND (organism_id:9606) AND (reviewed:true)',
        f'("{source_name.replace(chr(34), "")}") AND (organism_id:9606) AND (reviewed:true)',
    ]
    records: list[dict[str, Any]] = []
    for q in queries:
        r = session.get(UNIPROT, params={"query": q, "format": "json", "size": 25}, timeout=45)
        r.raise_for_status()
        records = r.json().get("results", [])
        if records:
            break

    target = norm(source_name)
    exact: dict[str, dict[str, Any]] = {}
    for rec in records:
        accession = rec.get("primaryAccession")
        if not accession:
            continue
        if any(norm(n) == target for n in uniprot_names(rec)):
            exact[accession] = rec

    if len(exact) == 1:
        accession, rec = next(iter(exact.items()))
        names = uniprot_names(rec)
        recommended = names[0] if names else ""
        length = ((rec.get("sequence") or {}).get("length"))
        return (
            "RESOLVED_EXACT_REVIEWED_HUMAN",
            accession,
            rec.get("uniProtkbId", ""),
            recommended,
            int(length) if isinstance(length, int) else None,
            names,
        )
    if len(exact) > 1:
        return (
            "REVIEW_REQUIRED_MULTIPLE_EXACT",
            "",
            "",
            "",
            None,
            sorted(exact.keys()),
        )
    candidates = [x.get("primaryAccession", "") for x in records if x.get("primaryAccession")]
    return (
        "REVIEW_REQUIRED_NO_EXACT",
        "",
        "",
        "",
        None,
        candidates[:10],
    )


def rcsb_entities(accession: str, rows: int = 25) -> list[str]:
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": rows},
            "sort": [{"sort_by": "entity_poly.rcsb_sample_sequence_length", "direction": "desc"}],
        },
    }
    r = session.post(RCSB_SEARCH, json=query, timeout=60)
    if r.status_code == 400:
        query["request_options"].pop("sort", None)
        r = session.post(RCSB_SEARCH, json=query, timeout=60)
    if r.status_code in (204, 404):
        return []
    r.raise_for_status()
    data = r.json()
    return [x.get("identifier") for x in data.get("result_set", []) if x.get("identifier")]


def coverage_from_alignment(entity: dict[str, Any], accession: str, uniprot_length: int | None) -> tuple[float | None, str]:
    aligns = entity.get("rcsb_polymer_entity_align") or []
    lengths: list[int] = []
    for a in aligns:
        ref = a.get("reference_database_accession")
        if ref and ref != accession:
            continue
        for reg in a.get("aligned_regions") or []:
            L = reg.get("length")
            if isinstance(L, int) and L > 0:
                lengths.append(L)
    if lengths and uniprot_length:
        # Alignments may occasionally overlap. Summed coverage is capped at 1.0;
        # manual QC remains mandatory.
        return min(sum(lengths) / uniprot_length, 1.0), "RCSB_UNIPROT_ALIGNMENT"

    seq = (entity.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can", "")
    seq_len = len(re.sub(r"[^A-Za-z]", "", seq))
    if seq_len and uniprot_length:
        return min(seq_len / uniprot_length, 1.0), "ENTITY_LENGTH_PROXY"
    return None, "UNAVAILABLE"


_entry_cache: dict[str, dict[str, Any]] = {}


def get_entry(entry_id: str) -> dict[str, Any]:
    if entry_id not in _entry_cache:
        r = session.get(f"{RCSB_DATA}/entry/{entry_id}", timeout=45)
        r.raise_for_status()
        _entry_cache[entry_id] = r.json()
    return _entry_cache[entry_id]


def inspect_entity(entity_id: str, accession: str, uniprot_length: int | None) -> dict[str, Any] | None:
    if "_" not in entity_id:
        return None
    entry_id, ent = entity_id.split("_", 1)
    r = session.get(f"{RCSB_DATA}/polymer_entity/{entry_id}/{ent}", timeout=45)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    pe = r.json()
    coverage, basis = coverage_from_alignment(pe, accession, uniprot_length)
    ids = pe.get("rcsb_polymer_entity_container_identifiers") or {}
    chains = ids.get("auth_asym_ids") or ids.get("asym_ids") or []
    poly = pe.get("entity_poly") or {}
    seq = poly.get("pdbx_seq_one_letter_code_can", "")
    entity_len = len(re.sub(r"[^A-Za-z]", "", seq))

    entry = get_entry(entry_id)
    methods = ";".join(x.get("method", "") for x in (entry.get("exptl") or []) if x.get("method"))
    resolutions = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
    resolution = min([float(x) for x in resolutions if isinstance(x, (int, float))], default=None)
    title = ((entry.get("struct") or {}).get("title") or "").replace("\n", " ").strip()
    mut = ((pe.get("rcsb_polymer_entity") or {}).get("rcsb_mutation_count"))
    mutation_count = int(mut) if isinstance(mut, int) else None

    return {
        "entity_id": entity_id,
        "entry_id": entry_id,
        "chains": ",".join(chains),
        "coverage": coverage,
        "coverage_basis": basis,
        "entity_length": entity_len,
        "method": methods,
        "resolution": resolution,
        "mutation_count": mutation_count,
        "title": title,
    }


def candidate_sort_key(x: dict[str, Any]):
    cov = x["coverage"] if x["coverage"] is not None else -1.0
    mut = x["mutation_count"] if x["mutation_count"] is not None else 999
    res = x["resolution"] if x["resolution"] is not None else 999.0
    return (-cov, mut, res, -x["entity_length"], x["entry_id"])


def pdb_audit(accession: str, length: int | None) -> tuple[str, dict[str, Any] | None, list[dict[str, Any]]]:
    ids = rcsb_entities(accession)
    hits: list[dict[str, Any]] = []
    for entity_id in ids:
        try:
            hit = inspect_entity(entity_id, accession, length)
            if hit:
                hits.append(hit)
        except Exception as e:
            print(f"WARN RCSB entity {entity_id}: {e}", file=sys.stderr)
        time.sleep(0.02)
    hits.sort(key=candidate_sort_key)
    if not hits:
        return "NO_EXPERIMENTAL_PDB_LINKED", None, []
    top = hits[0]
    cov = top["coverage"]
    if cov is None:
        status = "PDB_CANDIDATE_COVERAGE_UNKNOWN"
    elif cov >= 0.80:
        status = "FULL_OR_NEAR_FULL_PDB_CANDIDATE"
    elif cov >= 0.50:
        status = "PARTIAL_PDB_CANDIDATE"
    else:
        status = "DOMAIN_FRAGMENT_PDB_CANDIDATE"
    return status, top, hits[:5]


def compact_hit(h: dict[str, Any]) -> str:
    cov = "" if h["coverage"] is None else f'{h["coverage"]:.3f}'
    res = "" if h["resolution"] is None else f'{h["resolution"]:.2f}'
    mut = "" if h["mutation_count"] is None else str(h["mutation_count"])
    return "|".join([
        h["entity_id"], h["chains"], cov, h["coverage_basis"],
        h["method"], res, mut, str(h["entity_length"]),
    ])


def main() -> None:
    with PANEL.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if len(rows) != 285:
        raise RuntimeError(f"Expected frozen n=285, got {len(rows)}")

    fieldnames = [
        "panel_id", "protein_name", "source_detected_count", "source_detection_frequency",
        "identity_status", "uniprot_accession", "uniprot_id", "uniprot_recommended_name",
        "uniprot_length", "identity_candidates_or_names",
        "pdb_auto_status", "top_pdb_candidate", "top_entity", "top_chain",
        "top_coverage", "top_coverage_basis", "top_method", "top_resolution",
        "top_mutation_count", "top_entity_length", "top_title", "top5_candidates",
        "manual_structure_qc_required", "final_structure_selection_status",
    ]

    output: list[dict[str, Any]] = []
    for i, src in enumerate(rows, 1):
        name = src["protein_name"].strip()
        print(f"[{i:03d}/285] {name}", flush=True)
        try:
            istatus, acc, uid, recname, ulen, names = search_uniprot(name)
        except Exception as e:
            print(f"WARN UniProt {name}: {e}", file=sys.stderr)
            istatus, acc, uid, recname, ulen, names = "LOOKUP_ERROR", "", "", "", None, []

        result = {
            "panel_id": src["panel_id"],
            "protein_name": name,
            "source_detected_count": src["source_detected_count"],
            "source_detection_frequency": src["source_detection_frequency"],
            "identity_status": istatus,
            "uniprot_accession": acc,
            "uniprot_id": uid,
            "uniprot_recommended_name": recname,
            "uniprot_length": "" if ulen is None else ulen,
            "identity_candidates_or_names": ";".join(names),
            "pdb_auto_status": "NOT_QUERIED_UNRESOLVED_IDENTITY",
            "top_pdb_candidate": "",
            "top_entity": "",
            "top_chain": "",
            "top_coverage": "",
            "top_coverage_basis": "",
            "top_method": "",
            "top_resolution": "",
            "top_mutation_count": "",
            "top_entity_length": "",
            "top_title": "",
            "top5_candidates": "",
            "manual_structure_qc_required": "YES",
            "final_structure_selection_status": "PENDING_MANUAL_QC",
        }

        if acc:
            try:
                pstatus, top, hits = pdb_audit(acc, ulen)
                result["pdb_auto_status"] = pstatus
                if top:
                    result.update({
                        "top_pdb_candidate": top["entry_id"],
                        "top_entity": top["entity_id"],
                        "top_chain": top["chains"],
                        "top_coverage": "" if top["coverage"] is None else f'{top["coverage"]:.6f}',
                        "top_coverage_basis": top["coverage_basis"],
                        "top_method": top["method"],
                        "top_resolution": "" if top["resolution"] is None else f'{top["resolution"]:.3f}',
                        "top_mutation_count": "" if top["mutation_count"] is None else top["mutation_count"],
                        "top_entity_length": top["entity_length"],
                        "top_title": top["title"],
                        "top5_candidates": ";".join(compact_hit(h) for h in hits),
                    })
            except Exception as e:
                print(f"WARN PDB {name} {acc}: {e}", file=sys.stderr)
                result["pdb_auto_status"] = "PDB_LOOKUP_ERROR"
        output.append(result)
        time.sleep(0.03)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(output)

    counts: dict[str, int] = {}
    pdb_counts: dict[str, int] = {}
    for r in output:
        counts[r["identity_status"]] = counts.get(r["identity_status"], 0) + 1
        pdb_counts[r["pdb_auto_status"]] = pdb_counts.get(r["pdb_auto_status"], 0) + 1
    print("IDENTITY_COUNTS=" + json.dumps(counts, sort_keys=True))
    print("PDB_COUNTS=" + json.dumps(pdb_counts, sort_keys=True))
    print(f"WROTE {OUT}")


if __name__ == "__main__":
    main()
