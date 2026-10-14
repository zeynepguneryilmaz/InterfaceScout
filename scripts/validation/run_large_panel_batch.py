#!/usr/bin/env python3
"""Batch InterfaceScout V1 across the frozen 285-protein AlphaFold panel.

This is an applicability / structural-landscape analysis, not an adsorption
classifier. Panel membership, structure source, reference environment and
11-channel V1 model are frozen independently of InterfaceScout outcomes.

Outputs per shard:
- one compact CSV with one row per protein x chemistry channel;
- one status CSV with one row per source protein;
- one JSON file per successfully analyzed protein containing the complete
  InterfaceScout all-map response plus AlphaFold QC metadata.

AlphaFold pLDDT is recorded as QC metadata only and does not modify V1 scores.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import statistics
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[2]
BACKEND = ROOT / "backend"
sys.path.insert(0, str(BACKEND))

from interface_engine import analyze_all_maps  # noqa: E402
from model_settings import MODEL_VERSION  # noqa: E402

PANEL = ROOT / "validation/large_panel_285_2026-10-13/source_285_master.csv"
OUTROOT = ROOT / "validation/large_panel_285_2026-10-13/batch_outputs"

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb"
AF_API = "https://alphafold.ebi.ac.uk/api/prediction"

REFERENCE_PH = 7.4
REFERENCE_IONIC_MM = 150.0
REFERENCE_TEMP_K = 298.0

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
session.headers.update({"User-Agent": "InterfaceScout-large-panel-batch/2026-10-14"})


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())


def uniprot_names(rec: dict[str, Any]) -> list[str]:
    pd = rec.get("proteinDescription") or {}
    out: list[str] = []
    rn = ((pd.get("recommendedName") or {}).get("fullName") or {}).get("value")
    if rn:
        out.append(rn)
    for key in ("alternativeNames", "submissionNames"):
        for x in pd.get(key) or []:
            v = ((x or {}).get("fullName") or {}).get("value")
            if v:
                out.append(v)
    return out


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
            rec["primaryAccession"],
            rec.get("uniProtkbId", ""),
            names[0] if names else "",
            int((rec.get("sequence") or {}).get("length") or 0),
            "CURATED_ALIAS",
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
    exact = {}
    for rec in records:
        acc = rec.get("primaryAccession")
        if acc and any(norm(n) == target for n in uniprot_names(rec)):
            exact[acc] = rec
    if len(exact) != 1:
        raise RuntimeError(f"Could not uniquely resolve reviewed human UniProt for {source_name}")
    rec = next(iter(exact.values()))
    names = uniprot_names(rec)
    return (
        rec["primaryAccession"],
        rec.get("uniProtkbId", ""),
        names[0] if names else "",
        int((rec.get("sequence") or {}).get("length") or 0),
        "EXACT_NAME",
    )


def af_records(accession: str) -> list[dict[str, Any]]:
    r = session.get(f"{AF_API}/{accession}", timeout=60)
    if r.status_code == 404:
        return []
    r.raise_for_status()
    x = r.json()
    if isinstance(x, dict):
        x = [x]
    return x if isinstance(x, list) else []


def choose_primary(records: list[dict[str, Any]], uniprot_length: int):
    ranked = []
    for rec in records:
        start = int(rec.get("sequenceStart") or 1)
        end = int(rec.get("sequenceEnd") or uniprot_length or 0)
        span = max(0, end - start + 1)
        cov = span / uniprot_length if uniprot_length else 0.0
        ranked.append((cov, span, -start, str(rec.get("entryId") or ""), rec))
    ranked.sort(reverse=True)
    return ranked[0][-1] if ranked else None


def download_text(url: str) -> str:
    r = session.get(url, timeout=180)
    r.raise_for_status()
    return r.text


def plddt_qc(pdb_text: str) -> dict[str, Any]:
    vals = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except Exception:
                pass
    if not vals:
        return {
            "n_ca": 0, "mean": None, "median": None,
            "frac_ge90": None, "frac_ge70": None, "frac_lt50": None,
            "flag": "UNAVAILABLE",
        }
    n = len(vals)
    frac_lt50 = sum(v < 50 for v in vals) / n
    med = statistics.median(vals)
    if frac_lt50 >= 0.30:
        flag = "SUBSTANTIAL_LOW_CONFIDENCE"
    elif med < 70:
        flag = "LOW_MEDIAN_CONFIDENCE"
    else:
        flag = "STANDARD"
    return {
        "n_ca": n,
        "mean": statistics.fmean(vals),
        "median": med,
        "frac_ge90": sum(v >= 90 for v in vals) / n,
        "frac_ge70": sum(v >= 70 for v in vals) / n,
        "frac_lt50": frac_lt50,
        "flag": flag,
    }


def scalar(x):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.6f}"
    return x


def top_patch_summary(map_result: dict[str, Any]) -> dict[str, Any]:
    patches = map_result.get("patches") or []
    if not patches:
        return {
            "top_patch_rank": "", "top_patch_front": "", "top_patch_center": "",
            "top_patch_size": "", "top_patch_members": "",
            "top_patch_coherence": "", "top_patch_chemistry_support": "",
        }
    p = patches[0]
    members = p.get("members") or p.get("member_keys") or []
    if members and isinstance(members[0], dict):
        members = [m.get("key", "") for m in members]
    return {
        "top_patch_rank": 1,
        "top_patch_front": p.get("pareto_front", ""),
        "top_patch_center": p.get("center_key", p.get("center", "")),
        "top_patch_size": p.get("size", len(members)),
        "top_patch_members": ";".join(str(x) for x in members),
        "top_patch_coherence": p.get("patch_coherence", p.get("coherence", "")),
        "top_patch_chemistry_support": p.get("chemistry_support", ""),
    }


def selected_rows(all_rows, shard_index: int, shard_count: int):
    return [r for i, r in enumerate(all_rows) if i % shard_count == shard_index]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard-index", type=int, required=True)
    ap.add_argument("--shard-count", type=int, required=True)
    args = ap.parse_args()
    if args.shard_count < 1 or not (0 <= args.shard_index < args.shard_count):
        raise SystemExit("Invalid shard settings")

    with PANEL.open(newline="", encoding="utf-8-sig") as fh:
        all_rows = list(csv.DictReader(fh))
    if len(all_rows) != 285:
        raise RuntimeError(f"Expected frozen source n=285, got {len(all_rows)}")

    rows = selected_rows(all_rows, args.shard_index, args.shard_count)
    shard_dir = OUTROOT / f"shard_{args.shard_index:02d}_of_{args.shard_count:02d}"
    json_dir = shard_dir / "protein_json"
    json_dir.mkdir(parents=True, exist_ok=True)

    status_rows = []
    summary_rows = []

    for pos, src in enumerate(rows, 1):
        panel_id = src["panel_id"]
        name = src["protein_name"].strip()
        print(f"[shard {args.shard_index}/{args.shard_count} {pos}/{len(rows)}] {panel_id} {name}", flush=True)
        status = {
            "panel_id": panel_id,
            "protein_name": name,
            "source_detected_count": src["source_detected_count"],
            "source_detection_frequency": src["source_detection_frequency"],
            "status": "PENDING",
            "uniprot_accession": "",
            "uniprot_id": "",
            "uniprot_length": "",
            "identity_basis": "",
            "af_entry_id": "",
            "af_sequence_start": "",
            "af_sequence_end": "",
            "af_coverage_fraction": "",
            "plddt_mean": "",
            "plddt_median": "",
            "plddt_frac_lt50": "",
            "plddt_qc_flag": "",
            "n_surface_residues": "",
            "n_maps": "",
            "error": "",
        }
        try:
            acc, uid, recname, ulen, basis = resolve_identity(name)
            status.update({
                "uniprot_accession": acc, "uniprot_id": uid,
                "uniprot_length": ulen, "identity_basis": basis,
            })
            records = af_records(acc)
            if not records:
                status["status"] = "NO_ALPHAFOLD_MODEL"
                status_rows.append(status)
                continue
            rec = choose_primary(records, ulen)
            start = int(rec.get("sequenceStart") or 1)
            end = int(rec.get("sequenceEnd") or ulen or 0)
            cov = (max(0, end - start + 1) / ulen) if ulen else 0.0
            pdb_url = rec.get("pdbUrl", "")
            if not pdb_url:
                raise RuntimeError("AlphaFold record has no PDB URL")

            pdb_text = download_text(pdb_url)
            qc = plddt_qc(pdb_text)
            status.update({
                "af_entry_id": rec.get("entryId", ""),
                "af_sequence_start": start,
                "af_sequence_end": end,
                "af_coverage_fraction": f"{cov:.6f}",
                "plddt_mean": scalar(qc["mean"]),
                "plddt_median": scalar(qc["median"]),
                "plddt_frac_lt50": scalar(qc["frac_lt50"]),
                "plddt_qc_flag": qc["flag"],
            })

            result = analyze_all_maps(
                pH=REFERENCE_PH,
                ionic_mM=REFERENCE_IONIC_MM,
                temp_K=REFERENCE_TEMP_K,
                pdb_text=pdb_text,
                chain=None,
            )
            result["large_panel_metadata"] = {
                "panel_id": panel_id,
                "protein_name": name,
                "source_detected_count": int(src["source_detected_count"]),
                "source_detection_frequency": float(src["source_detection_frequency"]),
                "uniprot_accession": acc,
                "uniprot_id": uid,
                "uniprot_recommended_name": recname,
                "uniprot_length": ulen,
                "identity_basis": basis,
                "structure_source": "AlphaFold DB",
                "af_entry_id": rec.get("entryId", ""),
                "af_latest_version": rec.get("latestVersion", ""),
                "af_sequence_start": start,
                "af_sequence_end": end,
                "af_coverage_fraction": cov,
                "af_pdb_url": pdb_url,
                "plddt_qc": qc,
                "plddt_used_in_v1_scoring": False,
                "reference_condition": {
                    "pH": REFERENCE_PH,
                    "ionic_mM": REFERENCE_IONIC_MM,
                    "temp_K": REFERENCE_TEMP_K,
                },
                "adsorption_label_used_as_target": False,
            }

            out_json = json_dir / f"{panel_id}_{acc}.json"
            out_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

            status["status"] = "SUCCESS"
            status["n_surface_residues"] = result.get("diagnostics", {}).get("n_surface_residues", "")
            status["n_maps"] = result.get("n_maps", "")

            for map_key in result.get("map_order", []):
                m = result["maps"][map_key]
                tops = top_patch_summary(m)
                primary = m.get("primary_patches") or []
                props = [r.get("propensity") for r in (m.get("residues") or []) if isinstance(r.get("propensity"), (int, float))]
                summary_rows.append({
                    "panel_id": panel_id,
                    "protein_name": name,
                    "uniprot_accession": acc,
                    "af_entry_id": rec.get("entryId", ""),
                    "af_coverage_fraction": f"{cov:.6f}",
                    "plddt_qc_flag": qc["flag"],
                    "chemistry_channel": map_key,
                    "n_surface_residues": result.get("diagnostics", {}).get("n_surface_residues", ""),
                    "n_patches": m.get("n_patches", ""),
                    "n_pareto_primary_patches": m.get("n_pareto_primary_patches", ""),
                    "residue_propensity_mean": f"{statistics.fmean(props):.6f}" if props else "",
                    "residue_propensity_max": f"{max(props):.6f}" if props else "",
                    **{k: scalar(v) for k, v in tops.items()},
                })
        except Exception as exc:
            status["status"] = "ERROR"
            status["error"] = repr(exc)
            print(f"ERROR {panel_id} {name}: {exc}", file=sys.stderr)
        status_rows.append(status)

    status_fields = list(status_rows[0].keys()) if status_rows else []
    with (shard_dir / "protein_status.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=status_fields)
        w.writeheader()
        w.writerows(status_rows)

    summary_fields = [
        "panel_id","protein_name","uniprot_accession","af_entry_id",
        "af_coverage_fraction","plddt_qc_flag","chemistry_channel",
        "n_surface_residues","n_patches","n_pareto_primary_patches",
        "residue_propensity_mean","residue_propensity_max",
        "top_patch_rank","top_patch_front","top_patch_center","top_patch_size",
        "top_patch_members","top_patch_coherence","top_patch_chemistry_support",
    ]
    with (shard_dir / "channel_summary.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=summary_fields)
        w.writeheader()
        w.writerows(summary_rows)

    success = sum(r["status"] == "SUCCESS" for r in status_rows)
    errors = sum(r["status"] == "ERROR" for r in status_rows)
    no_model = sum(r["status"] == "NO_ALPHAFOLD_MODEL" for r in status_rows)
    meta = {
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "source_panel_n": 285,
        "assigned_n": len(rows),
        "success_n": success,
        "error_n": errors,
        "no_model_n": no_model,
        "model_version": MODEL_VERSION,
        "reference_condition": {
            "pH": REFERENCE_PH,
            "ionic_mM": REFERENCE_IONIC_MM,
            "temp_K": REFERENCE_TEMP_K,
        },
        "channels_per_success": 11,
    }
    (shard_dir / "shard_summary.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print("SHARD_SUMMARY=" + json.dumps(meta, sort_keys=True))

    # Do not fail the entire shard for per-protein errors. Errors remain explicit
    # in protein_status.csv and are auditable/retriable without selection bias.


if __name__ == "__main__":
    main()
