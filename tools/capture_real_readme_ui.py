"""Capture *real* InterfaceScout UI (not a mockup) for the README.
Uses Playwright Chromium/WebGL and bundled two-state adenylate kinase PDBs.
Three screenshots are captured after successful 11-map calculations. No
invented rankings, values, or UI controls.
"""
from __future__ import annotations
import os,subprocess,sys,time
from pathlib import Path
from urllib.request import urlopen
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
PDB1=ROOT/"tests"/"fixtures"/"1AKE_A.pdb"
PDB2=ROOT/"tests"/"fixtures"/"4AKE_A.pdb"
OUT=ROOT/"docs"/"figures";OUT.mkdir(parents=True,exist_ok=True)
def shot(page,name):
    page.wait_for_timeout(850)
    page.screenshot(path=str(OUT/name),full_page=True,animations="disabled")
    print("REAL_UI_SCREENSHOT",name,flush=True)
def main():
    assert PDB1.is_file() and PDB2.is_file()
    proc=subprocess.Popen([sys.executable,"-m","uvicorn","app:app","--app-dir",str(ROOT/"backend"),
       "--host","127.0.0.1","--port","8765"],cwd=str(ROOT),stdout=subprocess.DEVNULL,
       stderr=subprocess.STDOUT,env={**os.environ,"PYTHONUNBUFFERED":"1"})
    try:
        for _ in range(120):
            if proc.poll() is not None:raise RuntimeError("Backend process stopped")
            try:
                with urlopen("http://127.0.0.1:8765/health",timeout=1) as r:
                    if r.status==200:break
            except Exception:time.sleep(.5)
        else:raise RuntimeError("Backend did not become healthy")
        with sync_playwright() as p:
            browser=p.chromium.launch(headless=True,args=[
             "--use-gl=angle","--use-angle=swiftshader","--enable-unsafe-swiftshader",
             "--enable-webgl","--disable-dev-shm-usage"])
            page=browser.new_page(viewport={"width":1540,"height":980},device_scale_factor=1)
            page.on("pageerror",lambda e:print("BROWSER_JS_ERROR",str(e),flush=True))
            page.goto("http://127.0.0.1:8765/",wait_until="domcontentloaded")
            page.wait_for_function("() => document.querySelector('#run') && !document.querySelector('#run').disabled")
            page.set_input_files("#file",str(PDB1))
            page.wait_for_function("() => typeof pdbText === 'string' && pdbText.includes('ATOM')")
            page.locator("#run").click()
            page.wait_for_function("() => typeof data !== 'undefined' && data && data.n_maps === 11",timeout=120000)
            page.evaluate("() => {openChem('hydrophobic');setRepresentation('surface')}")
            page.wait_for_timeout(1100)
            shot(page,"interfacescout-real-single-protein.png")
            # Compare precisely two structures from the same protein, with
            # matched residue positions; identical pH in both calculations.
            page.locator('[data-workspace-target="compare"]').click()
            page.set_input_files("#secondFile",str(PDB2))
            page.wait_for_function("() => typeof decisionSecondPdb === 'string'")
            page.locator("#compareDecisionButton").click()
            page.wait_for_function("() => typeof ensembleResult !== 'undefined' && ensembleResult && ensembleResult.n_structures === 2 && decisionCompare !== null",timeout=180000)
            page.wait_for_timeout(700)
            shot(page,"interfacescout-real-compare-two.png")
            # Actual Multi-PDB collection, two supplied related structures.
            page.locator('[data-workspace-target="multi"]').click()
            page.set_input_files("#isMultiFiles",[str(PDB1),str(PDB2)])
            page.wait_for_function("() => typeof ensembleInputs !== 'undefined' && ensembleInputs.length === 2")
            page.locator("#isAnalyzeSet").click()
            page.wait_for_function("() => ensembleResult && ensembleResult.n_structures === 2",timeout=180000)
            page.wait_for_timeout(700)
            shot(page,"interfacescout-real-multi-pdb.png")
            browser.close()
    finally:
        proc.terminate()
        try:proc.wait(timeout=5)
        except subprocess.TimeoutExpired:proc.kill()
if __name__=="__main__":main()
