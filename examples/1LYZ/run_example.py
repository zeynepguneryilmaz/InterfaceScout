#!/usr/bin/env python3
"""Export a real 1LYZ analysis from an already running local InterfaceScout.

Uses the application's /analyze endpoint. No dependency installation, simulated
scores, example fallback values, or PyMOL requirement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import urllib.error
import urllib.request

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def post_analysis(base_url: str, request_data: dict) -> dict:
    req = urllib.request.Request(
        base_url.rstrip("/") + "/analyze",
        data=json.dumps(request_data).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=180) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"InterfaceScout HTTP {exc.code}: {detail[:1200]}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"Cannot reach InterfaceScout at {base_url}. Start the app first. "
            f"Underlying error: {exc.reason}"
        ) from exc


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--url", default="http://127.0.0.1:8000")
    ap.add_argument("--pdb", type=Path, help="Optional local PDB; avoids backend network download")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--out", type=Path, default=RESULTS)
    args = ap.parse_args()
    payload = {"chain": args.chain, "pH": args.ph, "ionic_mM": 150.0, "temp_K": 298.0}
    if args.pdb:
        if not args.pdb.is_file():
            ap.error(f"PDB file not found: {args.pdb}")
        pdb = args.pdb.read_text(encoding="utf-8", errors="replace")
        if not any(x.startswith("ATOM  ") for x in pdb.splitlines()):
            ap.error("The selected file contains no PDB ATOM records")
        payload["pdb_text"] = pdb
        source = str(args.pdb.resolve())
        original_sha256 = hashlib.sha256(pdb.encode("utf-8")).hexdigest()
    else:
        payload["pdb_id"] = "1LYZ"
        source = "RCSB PDB ID 1LYZ (retrieved by InterfaceScout backend)"
        original_sha256 = None

    result = post_analysis(args.url, payload)
    maps = result.get("maps")
    if not isinstance(maps, dict) or len(maps) != 11:
        raise RuntimeError(f"Expected 11 canonical maps, got {len(maps) if isinstance(maps, dict) else 'no map data'}")
    order = result.get("map_order") or list(maps.keys())
    if len(order) != 11 or set(order) != set(maps):
        raise RuntimeError("Inconsistent map_order and maps in analysis response")

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    (out / "1LYZ_analysis.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    residue_fields = [
        "map_key","map_label","residue_key","chain","res_seq","icode",
        "res_name","scrsa","local_score","propensity","multiscale_persistence",
    ]
    residue_rows = []
    patch_fields = [
        "map_key","map_label","patch_rank","pareto_front","center_key",
        "member_key","chain","res_seq","icode","res_name","is_chemistry_seed",
    ]
    patch_rows = []
    for key in order:
        m = maps[key]
        for r in m.get("residues", []):
            residue_rows.append({
                "map_key": key,"map_label": m.get("label",""),
                "residue_key": r.get("key",""),"chain": r.get("chain",""),
                "res_seq": r.get("res_seq",""),"icode": r.get("icode",""),
                "res_name": r.get("res_name",""),"scrsa": r.get("scrsa",""),
                "local_score": r.get("local_score",""),
                "propensity": r.get("propensity",""),
                "multiscale_persistence": r.get("multiscale_persistence",""),
            })
        for p in m.get("patches", []):
            seed_keys = set(p.get("seed_members") or [])
            rank = p.get("display_rank")
            for r in p.get("member_residues") or []:
                patch_rows.append({
                    "map_key": key,"map_label": m.get("label",""),
                    "patch_rank": rank,"pareto_front": p.get("pareto_front",""),
                    "center_key": p.get("center_key",""),
                    "member_key": r.get("key",""),"chain": r.get("chain",""),
                    "res_seq": r.get("res_seq",""),"icode": r.get("icode",""),
                    "res_name": r.get("res_name",""),
                    "is_chemistry_seed": r.get("key") in seed_keys,
                })
    write_csv(out / "1LYZ_residues.csv", residue_fields, residue_rows)
    write_csv(out / "1LYZ_patch_members.csv", patch_fields, patch_rows)
    manifest = {
        "structure": "1LYZ",
        "source": source,
        "local_input_sha256": original_sha256,
        "chain": args.chain,
        "pH": args.ph,
        "ionic_mM": 150.0,
        "temp_K": 298.0,
        "server_url": args.url,
        "model_version": result.get("version"),
        "public_version": result.get("public_version"),
        "source_pdb_sha256": result.get("input", {}).get("source_pdb_sha256"),
        "map_keys": order,
        "n_residue_rows": len(residue_rows),
        "n_patch_member_rows": len(patch_rows),
        "exports": [
            "1LYZ_analysis.json","1LYZ_residues.csv","1LYZ_patch_members.csv",
        ],
        "note": "Generated from a real local InterfaceScout /analyze response; not archived validation results.",
    }
    (out / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Saved 11 maps, {len(residue_rows)} residue rows and "
          f"{len(patch_rows)} patch-member rows to {out.resolve()}")
    print("For 3D screenshots, PDB/Excel/PDF files use Complete export in the app.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
