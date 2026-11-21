"""Surface Atlas PDF regression: honest legends, real screenshot, full data."""
from __future__ import annotations

import re
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"backend"))
sys.path.insert(0,str(ROOT/"tests"))
from test_ensemble_exports import minimal_analysis,png_data_uri
from export_bundle import MAP_NAMES, RANK_COLORS, _image_bytes, map_pdf, report_pdf


class SurfaceAtlasPdfTests(unittest.TestCase):
    def test_rank_legend_matches_actual_viewer(self):
        ui=(ROOT/"frontend"/"index.html").read_text(encoding="utf-8")
        match=re.search(r"const PATCH_COLORS=\[([^]]+)\]",ui)
        self.assertIsNotNone(match)
        ui_colors=tuple(re.findall(r"#[0-9a-fA-F]{6}",match.group(1)))
        self.assertEqual(ui_colors, RANK_COLORS)

    def test_map_pdf_is_readable_and_contains_real_image(self):
        analysis=minimal_analysis()
        photo=_image_bytes(png_data_uri())
        pdf=map_pdf("fixture",analysis,"hydrophobic",photo)
        self.assertTrue(pdf.startswith(b"%PDF-"))
        self.assertGreater(len(pdf),1500)
        self.assertIn(b"/Subtype /Image",pdf)

    def test_opposite_view_adds_separate_page_without_losing_primary(self):
        analysis=minimal_analysis()
        photo=_image_bytes(png_data_uri())
        baseline=map_pdf("test",analysis,"hydrophobic",photo)
        with_second=map_pdf("test",analysis,"hydrophobic",photo,
                            alternate_png=photo)
        self.assertGreater(with_second.count(b"/Type /Page"),
                           baseline.count(b"/Type /Page"))
        self.assertIn(b"/Subtype /Image",with_second)

    def test_complete_11_map_atlas_keeps_all_sections(self):
        analysis=minimal_analysis()
        photo=_image_bytes(png_data_uri())
        pdf=report_pdf("fixture",analysis,{k:photo for k in MAP_NAMES},None)
        self.assertTrue(pdf.startswith(b"%PDF-"))
        self.assertGreater(len(pdf),5000)
        self.assertIn(b"/Subtype /Image",pdf)
        self.assertGreaterEqual(pdf.count(b"/Type /Page"),12)


if __name__=="__main__":
    unittest.main()
