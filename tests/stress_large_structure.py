"""Real larger-protein stress and full-data export on the three supported OSes.

Uses the 2.7 MB bundled PDB (1JNJ), not a small synthetic PDB. Every output
is checked for structural completeness. Prints timings for diagnosis, without
asserting arbitrary platform-dependent performance thresholds.
"""
from __future__ import annotations

import base64
import csv
from io import BytesIO, StringIO
import os
from pathlib import Path
import platform
import sys
from time import perf_counter
import zipfile

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from interface_engine import analyze_all_maps
from export_bundle import MAP_NAMES, build_export_zip

CANDIDATES = [
    ROOT / "examples" / "4F5S" / "InterfaceScout_4F5S_anionic_Bfactor.pdb",
    ROOT / "examples" / "1SUV" / "InterfaceScout_1SUV_hydrophobic_Bfactor.pdb",
    ROOT / "examples" / "1TGU" / "InterfaceScout_1TGU_hydrophobic_Bfactor.pdb",
    ROOT / "examples" / "3ECA" / "InterfaceScout_3ECA_cationic_Bfactor.pdb",
    ROOT / "examples" / "1JNJ" / "InterfaceScout_1JNJ_anionic_Bfactor.pdb",
]


def first_model_residue_count(pdb: str) -> int:
    rows = set()
    for line in pdb.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith("ATOM  "):
            rows.add((line[21:22], line[22:26], line[26:27]))
    return len(rows)


def choose_large_first_model() -> tuple[Path, str, int]:
    options = []
    for path in CANDIDATES:
        if path.is_file():
            pdb = path.read_text("utf-8")
            options.append((path, pdb, first_model_residue_count(pdb)))
    assert options, "No example PDB available for stress test"
    path, pdb, n = max(options, key=lambda p: p[2])
    assert n > 300, f"No stress fixture with >300 first-model residues; counts: " + (
        ", ".join(f"{p.name}: {c}" for p, _, c in options)
    )
    return path, pdb, n


def png_uri(seed: int) -> str:
    image = Image.new("RGB", (800, 430), (233, 240, 251))
    for x in range(110, 700):
        for y in range(80, 350):
            if ((x - 400) / 255)**2 + ((y - 215) / 126)**2 < 1.0:
                image.putpixel((x, y), ((35 + seed * 17) % 220, 85 + y // 15, 175))
    stream = BytesIO()
    image.save(stream, "PNG")
    return "data:image/png;base64," + base64.b64encode(stream.getvalue()).decode("ascii")


def main():
    fixture, pdb, model_residues = choose_large_first_model()
    print(f"Selected {fixture.name}, first model has {model_residues} residues "
          f"({len(pdb):,} file bytes)", flush=True)
    assert len(pdb) > 500_000, "Selected stress fixture is too small"

    started = perf_counter()
    analysis = analyze_all_maps(pdb_text=pdb, pH=7.4)
    duration_analyze = perf_counter() - started
    assert analysis["n_maps"] == 11
    n = len(analysis["residue_profiles"])
    assert n > 300, f"Unexpectedly few standard residues after preparation: {n}"
    assert len(analysis["map_order"]) == 11
    assert analysis["input"]["source_pdb_sha256"]
    print(f"Large PDB {len(pdb):,} bytes, {n} residues; 11 maps: {duration_analyze:.1f} sec",
          flush=True)

    started = perf_counter()
    blob = build_export_zip({
        "structures": [{
            "name": "large_protein",
            "analysis": analysis,
            "pdb_text": pdb,
            "screenshots": {
                name: png_uri(i + 1) for i, name in enumerate(MAP_NAMES)
            },
            "combined_maps": None,
        }],
        "comparison": [],
        "consensus": None,
        "consensus_screenshots": {},
    })
    duration_export = perf_counter() - started
    assert len(blob) > 100_000
    with zipfile.ZipFile(BytesIO(blob)) as z:
        assert z.testzip() is None
        names = z.namelist()
        base = "structures/large_protein/"
        assert sum(f.startswith(base + "pdf_maps/") and f.endswith(".pdf")
                   for f in names) == 11
        assert sum(f.startswith(base + "screenshots/") and f.endswith(".png")
                   for f in names) == 11
        assert sum(f.startswith(base + "pdb_bfactor/") and f.endswith(".pdb")
                   for f in names) == 44
        for i, name in enumerate(MAP_NAMES, start=1):
            report = z.read(base + f"pdf_maps/{i:02d}_{name}.pdf")
            assert report.startswith(b"%PDF")
            assert b"/Subtype /Image" in report
        full = z.read(base + "full_report.pdf")
        assert full.startswith(b"%PDF")
        assert full.count(b"/Subtype /Image") >= 11, \
            "Full PDF should include distinct 3D images for all 11 maps"
        text = z.read(base + "csv/Map_Residue_Scores.csv").decode("utf-8-sig")
        reader = csv.DictReader(StringIO(text))
        count = sum(1 for _ in reader)
        assert count == 11 * n, (count, n)
        assert z.read(base + "original.pdb").decode("utf-8") == pdb

    print(
        f"LARGE STRUCTURE STRESS PASS: {platform.platform()}, "
        f"{n} residues, {len(names)} output files, "
        f"{len(blob):,} ZIP bytes, "
        f"analyze={duration_analyze:.2f}s, "
        f"export={duration_export:.2f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
