"""Real-browser end-to-end smoke test (GitHub Actions ubuntu browser job only).

Loads a bundled PDB into the exact user interface, exercises 3Dmol, exports
the actual viewer screenshots, checks ZIP contents and a two-PDB consensus.
"""
from __future__ import annotations

from io import BytesIO
import os
from pathlib import Path
import subprocess
import sys
import time
from urllib.request import urlopen
import zipfile

from PIL import Image, ImageChops
from openpyxl import load_workbook
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
PDB = ROOT / "examples" / "1CRN" / "InterfaceScout_1CRN_anionic_Bfactor.pdb"
SERVER = "http://127.0.0.1:8765"


def assert_export(download, label: str, multi: bool):
    path = download.path()
    if not path:
        raise AssertionError(label + ": browser has no downloaded ZIP")
    b = Path(path).read_bytes()
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
            # Selection changes presentation, never the full scientific export:
            # every Excel workbook must contain scores for all eleven maps.
            workbook = load_workbook(BytesIO(archive.read(base + "/complete_results.xlsx")), read_only=True, data_only=True)
            assert "Map_Residue_Scores" in workbook.sheetnames
            sheet = workbook["Map_Residue_Scores"]
            headers = [cell.value for cell in next(sheet.iter_rows(max_row=1))]
            map_col = headers.index("Map")
            exported_maps = {row[map_col] for row in sheet.iter_rows(min_row=2, values_only=True)}
            expected = {
                "anionic", "cationic", "hydrophobic", "pi_carbon", "hbond_donor",
                "hbond_acceptor", "oxide", "calcium_phosphate", "metal_coord",
                "soft_metal_sulfur", "phosphate",
            }
            assert exported_maps == expected, (label, exported_maps)
            workbook.close()
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
            # Presets are transparent user choices, not unexplained predictions.
            presets = page.locator("#researchGoal option").all_text_contents()
            assert all("Suggest" not in p for p in presets)
            assert not any("donor + acceptor" in p.lower() for p in presets)
            page.locator("#researchGoal").select_option("hydrophobic_aromatic")
            page.get_by_role("button", name="Apply preset").click()
            page.wait_for_function("() => selected.has('hydrophobic') && selected.has('pi_carbon') && selected.size === 2")
            page.wait_for_function("() => decisionCombo !== null && decisionCombo.selected_maps.length === 2")
            page.locator("#researchGoal").select_option("individual")
            page.get_by_role("button", name="Apply preset").click()
            page.wait_for_function("() => selected.size === 1")
            assert page.locator('.workspace-mode[aria-pressed="true"]').inner_text() == "Single protein"
            assert page.locator("#isFullExport").is_visible()
            assert not page.locator("#isMultiFiles").is_visible()
            assert page.locator("#run").is_visible()
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
                # Retry manually with a *real user click* on a persistent link.
                # Windows Chrome/Edge can block the programmatic click following
                # a long screenshot/PDF generation process.
                ready = page.locator("#isReadyDownload")
                assert ready.is_visible()
                assert ready.get_attribute("href").startswith("blob:")
                assert ready.get_attribute("download").endswith("_Complete_Results.zip")
                with page.expect_download(timeout=30_000) as manual_download:
                    ready.click()
                assert_export(manual_download.value, "Manual save retry", multi=False)
            except Exception as e:
                print("UI error:", page.locator("#isError").inner_text())
                raise

            # Compare workflow must export BOTH structures and all eleven maps,
            # independent of whichever chemistry chips are currently selected.
            page.locator('[data-workspace-target="compare"]').click()
            assert page.locator("#secondFile").is_visible()
            assert not page.locator("#isMultiFiles").is_visible()
            page.set_input_files("#secondFile", str(PDB))
            page.wait_for_function("() => typeof decisionSecondPdb === 'string'")
            page.locator("#compareDecisionButton").click()
            page.wait_for_function(
                "() => ensembleResult !== null && ensembleResult.n_structures === 2 "
                "&& decisionCompare !== null",
                timeout=150_000,
            )
            try:
                with page.expect_download(timeout=240_000) as dinfo:
                    page.locator("#isFullExport").click()
                assert_export(dinfo.value, "Compare-two full export", multi=True)
            except Exception:
                print("Compare UI error:", page.locator("#isError").inner_text())
                raise

            page.locator('[data-workspace-target="multi"]').click()
            assert page.locator("#isMultiFiles").is_visible()
            assert not page.locator("#secondFile").is_visible()
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
