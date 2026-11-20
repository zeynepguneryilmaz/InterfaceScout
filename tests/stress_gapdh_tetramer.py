"""Actual experimental GAPDH stress: monomer vs homotetramer.

PDB 1U8F, human placental GAPDH, X-ray 1.75 A, 4 chains and
approximately 1332 modeled residues. The original structure is obtained
directly from RCSB. NAD/other non-protein molecules are removed by the
existing protein preparation logic, with no changes to the frozen maps.

This validates both physical multimer-dependent accessibility and
cross-platform large PDB export; it does not model dynamics or adsorption.
"""
from __future__ import annotations

import base64
import csv
import hashlib
from io import BytesIO, StringIO
import json
import platform
from pathlib import Path
import sys
from time import perf_counter
from urllib.request import Request, urlopen
import zipfile

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))
from interface_engine import analyze_all_maps
from export_bundle import MAP_NAMES, build_export_zip, _bfactor_pdb

PDB_ID = "1U8F"
PDB_URLS = [
    "https://files.rcsb.org/download/1U8F.pdb",
    "https://files.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb1u8f.ent.gz",
]


def fetch_pdb() -> bytes:
    last = None
    for url in PDB_URLS:
        try:
            with urlopen(Request(url, headers={"User-Agent": "InterfaceScout-QA/1.0"}), timeout=65) as response:
                data = response.read()
            if url.endswith(".gz"):
                import gzip
                data = gzip.decompress(data)
            if len(data) < 400_000 or b"ATOM  " not in data:
                raise ValueError("Downloaded GAPDH is empty or incomplete")
            return data
        except Exception as exc:
            last = exc
    raise RuntimeError(f"Unable to download experimentally deposited 1U8F PDB from wwPDB: {last}")


def determine_chains(pdb: str) -> tuple[list[str], dict[str, int]]:
    residues: dict[str, set[tuple[str, str]]] = {}
    for line in pdb.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM  ") or len(line) < 27:
            continue
        chain = line[21]
        residues.setdefault(chain, set()).add((line[22:26], line[26:27]))
    names = [x for x in residues if x.strip()]
    return names, {x: len(residues[x]) for x in names}


def png_uri(seed: int) -> str:
    im = Image.new("RGB", (700, 390), (221, 235, 248))
    for x in range(60, 650):
        for y in range(50, 350):
            if (x - 350)**2 / 265**2 + (y - 195)**2 / 135**2 < 1:
                im.putpixel((x, y), (50 + (seed * 9) % 70, 76 + y // 10, 160))
    stream = BytesIO()
    im.save(stream, format="PNG")
    return "data:image/png;base64," + base64.b64encode(stream.getvalue()).decode("ascii")


def main():
    raw = fetch_pdb()
    pdb = raw.decode("utf-8", errors="replace")
    chains, raw_counts = determine_chains(pdb)
    assert len(chains) == 4, f"GAPDH must contain 4 protein chains; got {chains}"
    assert sum(raw_counts.values()) >= 1250, f"Too few GAPDH tetramer residues: {raw_counts}"
    print(f"GAPDH PDB {PDB_ID}: {len(raw):,} bytes, 4 chains, {raw_counts} residues", flush=True)
    t0 = perf_counter()
    monomer = analyze_all_maps(pdb_text=pdb, chain=chains[0])
    print(f"Monomer: {len(monomer['residue_profiles'])} residues; "
          f"analysis {perf_counter() - t0:.1f}s", flush=True)
    t1 = perf_counter()
    tetramer = analyze_all_maps(pdb_text=pdb, chain=None)
    print(f"Tetramer: {len(tetramer['residue_profiles'])} residues; "
          f"analysis {perf_counter() - t1:.1f}s", flush=True)

    assert monomer["n_maps"] == tetramer["n_maps"] == 11
    assert len(monomer["residue_profiles"]) >= 300
    assert len(tetramer["residue_profiles"]) >= 1250
    assert len(set(r["chain"] for r in tetramer["residue_profiles"])) == 4
    assert monomer["input"]["source_pdb_sha256"] == tetramer["input"]["source_pdb_sha256"]
    atomically_exposed_mono = {
        r["key"]: float(r["sidechain_sasa"])
        for r in monomer["residue_profiles"]
    }
    tetramer_profile = {r["key"]: r for r in tetramer["residue_profiles"]}
    common = atomically_exposed_mono.keys() & tetramer_profile.keys()
    assert len(common) >= 300
    total_surface_occlusion = sum(
        atomically_exposed_mono[key] -
        float(tetramer_profile[key]["sidechain_sasa"])
        for key in common
    )
    assert total_surface_occlusion > 50, (
        "Tetrameric assembly should bury some surfaces relative to the isolated subunit; "
        f"measured side-chain SASA occlusion={total_surface_occlusion:.1f} A^2"
    )
    print("Inter-subunit surface occlusion: "
          f"{total_surface_occlusion:.1f} A^2 in tested chain", flush=True)

    started = perf_counter()
    zip_data = build_export_zip({
        "structures": [{
            "name": PDB_ID+"_tetramer",
            "pdb_text": pdb,
            "analysis": tetramer,
            "screenshots": {k: png_uri(i) for i, k in enumerate(MAP_NAMES)},
            "combined_maps": None,
        }],
        "comparison": [], "consensus": None,
    })
    elapsed_export = perf_counter() - started
    with zipfile.ZipFile(BytesIO(zip_data)) as package:
        names = package.namelist()
        assert package.testzip() is None
        prefix = f"structures/{PDB_ID}_tetramer/"
        assert sum(n.startswith(prefix+"pdf_maps/") and n.endswith(".pdf")
                   for n in names) == 11
        assert sum(n.startswith(prefix+"screenshots/") and n.endswith(".png")
                   for n in names) == 11
        assert sum(n.startswith(prefix+"pdb_bfactor/") and n.endswith(".pdb")
                   for n in names) == 44
        assert package.read(prefix+"full_report.pdf").count(b"/Subtype /Image") >= 11
        report = package.read(prefix+"csv/Map_Residue_Scores.csv").decode("utf-8-sig")
        rows = list(csv.DictReader(StringIO(report)))
        assert len(rows) == 11 * len(tetramer["residue_profiles"]), (
            len(rows), len(tetramer["residue_profiles"])
        )
        for i, key in enumerate(MAP_NAMES, 1):
            pdf = package.read(prefix+f"pdf_maps/{i:02d}_{key}.pdf")
            assert pdf.startswith(b"%PDF") and b"/Subtype /Image" in pdf
        bfactor = package.read(prefix+"pdb_bfactor/01_anionic_Bfactor.pdb").decode("utf-8")
        exported_chains, outcounts = determine_chains(bfactor)
        assert set(exported_chains) == set(chains), (exported_chains, chains)
        assert outcounts == raw_counts, (outcounts, raw_counts)
        assert "NOT EXPERIMENTAL" in bfactor

    print(
        "GAPDH TETRAMER STRESS PASS: "
        f"OS={platform.platform()}, "
        f"tetramer residues={len(tetramer['residue_profiles'])}, "
        f"11 full PDFs, 11 screenshots, 44 B-factor PDBs, "
        f"11*N CSV rows={len(rows)}, "
        f"ZIP bytes={len(zip_data):,}, export seconds={elapsed_export:.1f}",
        flush=True
    )


if __name__ == "__main__":
    main()
