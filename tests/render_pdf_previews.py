"""Render and sanity-check the real 1CRN PDF previews produced by browser_smoke."""
from pathlib import Path
from PIL import ImageStat
import pypdfium2 as pdfium
from pypdf import PdfReader

DIRECTORY=Path(__file__).resolve().parents[1]/"dist"/"pdf_preview"
for filename in ("1CRN_hydrophobic_map.pdf", "1CRN_Complete_Surface_Atlas.pdf"):
    src=DIRECTORY/filename
    if not src.is_file():
        raise AssertionError(f"Real-browser PDF preview missing: {src}")
    doc=pdfium.PdfDocument(str(src))
    pages=len(doc)
    if filename.endswith("Complete_Surface_Atlas.pdf"):
        assert pages >= 13, (filename,pages)
    else:
        assert pages >= 2, (filename,pages)
    screenshot=doc[0].render(scale=1.7).to_pil().convert("RGB")
    assert screenshot.width>500 and screenshot.height>700
    assert max(ImageStat.Stat(screenshot).stddev)>15, "Unexpected blank/low-contrast PDF first page"
    screenshot.save(DIRECTORY/(src.stem+"_first_page.png"))
    if "Complete_Surface_Atlas" in filename:
        # Previous cover overflow made page 2 almost empty.
        second=doc[1].render(scale=1.4).to_pil().convert("RGB")
        width,height=second.size
        area=second.crop((int(width*.08),int(height*.19),
                          int(width*.93),int(height*.84)))
        nonwhite=sum(1 for pixel in area.getdata() if min(pixel)<230)
        coverage=nonwhite/(area.width*area.height)
        second.save(DIRECTORY/(src.stem+"_second_page.png"))
        second_text=doc[1].get_textpage().get_text_range()
        assert "Contents" in second_text, "Navigation table must appear on page 2"
        assert "Anionic surface" in second_text and "Hydrophobic surface" in second_text, (
            "TOC must contain the actual clickable map entries, not only its heading")
        assert coverage>0.012, f"Atlas contents page is unexpectedly blank: {coverage:.3%}"
        third_text=doc[2].get_textpage().get_text_range()
        assert "Anionic surface" in third_text, "Eleven-map gallery must begin on page 3"
        assert "Chem. support" in third_text, "Ranked summary must display chemistry support"
        cover=doc[0].get_textpage().get_text_range()
        assert "name" in cover.lower() or "surface" in cover.lower()
        assert "Version" in cover, "Atlas cover must identify the calculation version"
        # Read PDF outline entries via pypdf; Pdfium has changed bookmark
        # classes between versions, but native PDF /Outlines remain stable.
        def flatten(items):
            for element in items:
                if isinstance(element,list):
                    yield from flatten(element)
                else:
                    yield element.title
        outline=list(flatten(PdfReader(str(src)).outline))
        names=outline
        for title in ("Anionic surface","Hydrophobic surface","Complete model parameters"):
            assert any(title in name for name in names), f"PDF bookmark absent: {title}"
        print(f"Atlas second-page navigation coverage: {coverage:.1%}; {len(outline)} outline entries")
    print(f"{filename}: {pages} pages; first page rendered {screenshot.size}")
