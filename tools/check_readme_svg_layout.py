"""Render SVG figures in Chromium and fail on actual text-text overlap.

Measures SVGTextElement.getBBox() in browser coordinates, tests margins and
intersections. Also saves PNG proof images as CI artifacts for visual review.
"""
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
FIGURES=["interfacescout-workflows.svg","interfacescout-method.svg"]
def main():
    with sync_playwright() as playwright:
        browser=playwright.chromium.launch(headless=True,args=["--no-sandbox"])
        for name in FIGURES:
            source=(ROOT/"docs"/"figures"/name).resolve()
            page=browser.new_page(viewport={"width":1600,"height":1100},device_scale_factor=1)
            page.goto(source.as_uri())
            page.wait_for_selector("svg")
            measurements=page.evaluate("""() => {
              const root=document.querySelector('svg');
              const bb=root.viewBox.baseVal;
              const texts=Array.from(root.querySelectorAll('text')).map((e,i)=>{
                 const b=e.getBBox();return {i,text:e.textContent.trim(),
                    x:b.x,y:b.y,w:b.width,h:b.height}
              });
              return {view:{x:bb.x,y:bb.y,w:bb.width,h:bb.height},texts};
            }""")
            issues=[];texts=measurements["texts"];v=measurements["view"]
            for a in texts:
                if a["x"] < v["x"]-1 or a["y"]<v["y"]-1 or a["x"]+a["w"]>v["x"]+v["w"]+1 or a["y"]+a["h"]>v["y"]+v["h"]+1:
                    issues.append(f"Outside SVG: {a['text']} {a}")
            for i,a in enumerate(texts):
                for b in texts[i+1:]:
                    ix=min(a["x"]+a["w"],b["x"]+b["w"])-max(a["x"],b["x"])
                    iy=min(a["y"]+a["h"],b["y"]+b["h"])-max(a["y"],b["y"])
                    if ix>2 and iy>2:issues.append(f"Text overlap: '{a['text']}' / '{b['text']}' (intersection {ix:.1f}x{iy:.1f})")
            png=ROOT/"docs"/"figures"/name.replace(".svg","-layout-proof.png")
            page.locator("svg").screenshot(path=str(png))
            print(f"SVG_LAYOUT {name} text_nodes={len(texts)} overlaps={len(issues)} proof={png.name}",flush=True)
            for issue in issues:print("SVG_LAYOUT_ERROR",issue,flush=True)
            page.close()
            if issues:raise AssertionError(f"{name}: {len(issues)} layout issues")
        browser.close()
if __name__=="__main__": main()
