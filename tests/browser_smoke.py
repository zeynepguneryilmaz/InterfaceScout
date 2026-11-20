"""Real-browser end-to-end smoke test (GitHub Actions ubuntu browser job only).

Loads a bundled PDB into the exact user interface, exercises 3Dmol, exports
the actual viewer screenshots, checks ZIP contents and a two-PDB consensus.
"""
from __future__ import annotations

from io import BytesIO
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from urllib.request import urlopen
import zipfile

from PIL import Image, ImageChops
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
PDB = ROOT / "examples" / "1CRN" / "InterfaceScout_1CRN_anionic_Bfactor.pdb"
SERVER = "http://127.0.0.1:8765"


def assert_export(download, label: str, multi: bool):
    path = download.path()
    if not path:
        raise AssertionError(label + ": browser has no downloaded ZIP")
    b = Path(path).read_bytes()
    # Optional machine-readable release evidence: archive exactly the ZIP that
    # passed the structural/content checks, never placeholder renderings.
    evidence = os.environ.get("INTERFACESCOUT_ARTIFACT_DIR")
    if evidence:
        dest = Path(evidence)
        dest.mkdir(parents=True, exist_ok=True)
        target = dest / ("multi_pdb_browser_export.zip" if multi else "single_pdb_browser_export.zip")
        target.write_bytes(b)
        record = {"workflow": label, "file": target.name,
                  "sha256": hashlib.sha256(b).hexdigest(),
                  "bytes": len(b), "passed_export_checks": False}
        (dest / (target.stem + ".json")).write_text(
            json.dumps(record, indent=2), encoding="utf-8")
    with zipfile.ZipFile(BytesIO(b)) as archive:
        names = archive.namelist()
        bases = ["structures/InterfaceScout_1CRN_anionic_Bfactor"]
        if multi:
            bases = ["structures/InterfaceScout_1CRN_anionic_Bfactor",
                     "structures/InterfaceScout_1CRN_anionic_Bfactor_2"]
        for base in bases:
            pdfs = [p for p in names if p.startswith(base + "/pdf_maps/") and p.endswith(".pdf")]
            screenshots = [p for p in names if p.startswith(base + "/screenshots/") and p.endswith(".png")]
            bfactor = [p for p in names if p.startswith(base + "/pdb_bfactor/") and p.endswith(".pdb")]
            if len(pdfs) != 11 or len(screenshots) != 11 or len(bfactor) != 44:
                raise AssertionError(
                    f"{label}: unexpected map outputs for {base}, "
                    f"PDF={len(pdfs)}, screenshot={len(screenshots)}, PDB={len(bfactor)}"
                )
            a = archive.read(base + "/screenshots/01_anionic.png")
            c = archive.read(base + "/screenshots/03_hydrophobic.png")
            with Image.open(BytesIO(a)) as im_a, Image.open(BytesIO(c)) as im_c:
                if im_a.width < 200 or im_a.height < 100:
                    raise AssertionError(label + ": PNG 3D viewport too small")
                if im_a.size != im_c.size:
                    raise AssertionError(label + ": inconsistent 3D viewport screenshots")
                rgb_a = im_a.convert("RGB")
                rgb_c = im_c.convert("RGB")
                if ImageChops.difference(rgb_a, rgb_c).getbbox() is None:
                    raise AssertionError(label + ": different maps exported identical screenshots")
            if base + "/combined/selected_combination.png" not in names:
                raise AssertionError(label + ": selected-map screenshot is missing")
            if base + "/csv/Map_Colocation.csv" not in names:
                raise AssertionError(label + ": map co-location CSV is missing")
            with zipfile.ZipFile(BytesIO(b)) as second:
                map_pdf = second.read(base + "/pdf_maps/01_anionic.pdf")
                full_pdf = second.read(base + "/full_report.pdf")
                if b"/Subtype /Image" not in map_pdf or b"/Subtype /Image" not in full_pdf:
                    raise AssertionError(label + ": PDF reports lack embedded 3D screenshots")
        if multi:
            if "consensus/consensus_report.pdf" not in names:
                raise AssertionError(label + ": consensus report absent")
            if "comparison/comparison_report.pdf" not in names:
                raise AssertionError(label + ": comparison report absent")
    if evidence:
        record["passed_export_checks"] = True
        record["n_archive_members"] = len(names)
        (dest / (target.stem + ".json")).write_text(
            json.dumps(record, indent=2), encoding="utf-8")
    print(f"{label}: PASS ({len(b):,} ZIP bytes; {len(names)} files)")


def main():
    if not PDB.exists():
        raise AssertionError("Bundled source PDB missing")
    proc = subprocess.Popen(
        [sys.executable, "-m", "uvicorn", "app:app", "--app-dir", str(ROOT / "backend"),
         "--host", "127.0.0.1", "--port", "8765"],
        cwd=str(ROOT), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT,
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
    )
    try:
        for attempt in range(120):
            if proc.poll() is not None:
                raise RuntimeError("InterfaceScout server failed to start")
            try:
                with urlopen(SERVER + "/health", timeout=1) as r:
                    if r.status == 200:
                        break
            except Exception:
                time.sleep(.5)
        else:
            raise RuntimeError("InterfaceScout server did not become healthy")

        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(
                headless=True,
                args=["--use-gl=angle", "--use-angle=swiftshader",
                      "--enable-unsafe-swiftshader", "--enable-webgl",
                      "--disable-dev-shm-usage"],
            )
            page = browser.new_page(viewport={"width": 1200, "height": 850},
                                    accept_downloads=True)
            page.on("pageerror", lambda error: print("BROWSER PAGE ERROR:", error))
            page.goto(SERVER, wait_until="domcontentloaded")
            page.wait_for_function("() => document.querySelector('#run') && !document.querySelector('#run').disabled")
            page.set_input_files("#file", str(PDB))
            page.wait_for_function("() => typeof pdbText === 'string' && pdbText.includes('ATOM')")
            page.locator("#run").click()
            page.wait_for_function(
                "() => typeof data !== 'undefined' && data !== null && data.n_maps === 11",
                timeout=120_000,
            )
            print("Single-PDB UI: 11 maps computed")
            # Exercise actual multi-chemistry intersection selection in the UI.
            page.evaluate("() => { if(selected.has('anionic')) toggleChem('anionic'); }")
            page.evaluate("() => toggleChem('hydrophobic')")
            page.evaluate("() => toggleChem('pi_carbon')")
            page.wait_for_function(
                "() => decisionCombo !== null && decisionCombo.selected_maps.length === 2",
                timeout=60_000,
            )
            print("Two-map intersection/co-location completed")
            try:
                with page.expect_download(timeout=180_000) as dinfo:
                    page.locator("#isFullExport").click()
                assert_export(dinfo.value, "Single-PDB full-export", multi=False)
            except Exception as e:
                print("UI error:", page.locator("#isError").inner_text())
                raise

            page.set_input_files("#isMultiFiles", [str(PDB), str(PDB)])
            page.wait_for_function("() => ensembleInputs.length === 2")
            page.locator("#isAnalyzeSet").click()
            page.wait_for_function(
                "() => ensembleResult !== null && ensembleResult.n_structures === 2 "
                "&& ensembleResult.consensus !== null",
                timeout=150_000,
            )
            print("Multi-PDB UI: both structures and matched consensus computed")
            try:
                with page.expect_download(timeout=240_000) as dinfo:
                    page.locator("#isFullExport").click()
                assert_export(dinfo.value, "Multi-PDB full-export", multi=True)
            except Exception:
                print("UI error:", page.locator("#isError").inner_text())
                raise
            browser.close()
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


if __name__ == "__main__":
    main()
