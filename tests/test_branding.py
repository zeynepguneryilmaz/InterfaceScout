"""Regression checks for the single approved InterfaceScout brand mark."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


class BrandingTests(unittest.TestCase):
    def test_app_header_and_tab_icon_match_approved_master(self):
        logo = ROOT / "frontend" / "logo.png"
        favicon = ROOT / "frontend" / "favicon.png"
        self.assertTrue(logo.is_file())
        self.assertTrue(favicon.is_file())
        with Image.open(logo) as im:
            self.assertEqual(im.size, (256, 256))
            self.assertIn("A", im.getbands())
            self.assertEqual(im.getpixel((0, 0))[3], 0)
            self.assertEqual(im.getpixel((128, 128))[3], 255)
        with Image.open(favicon) as im:
            self.assertEqual(im.size, (48, 48))
            self.assertIn("A", im.getbands())
            self.assertEqual(im.getpixel((0, 0))[3], 0)

    def test_windows_icon_and_desktop_launcher_exist(self):
        icon = ROOT / "assets" / "InterfaceScout.ico"
        self.assertTrue(icon.is_file())
        with Image.open(icon) as im:
            self.assertEqual(im.format, "ICO")
            self.assertEqual(im.size, (256, 256))
        self.assertTrue((ROOT / "Create Desktop Shortcut.bat").is_file())
        self.assertTrue((ROOT / "launchers" / "create_windows_shortcut.ps1").is_file())

    def test_distribution_includes_only_platform_relevant_launchers(self):
        sys.path.insert(0, str(ROOT / "tools"))
        from build_release_archives import included_files

        names = {
            platform: {p.relative_to(ROOT).as_posix() for p in included_files(platform)}
            for platform in ("Windows", "macOS", "Linux")
        }
        for platform in names:
            self.assertTrue({"frontend/logo.png", "frontend/favicon.png"} <= names[platform])
        self.assertIn("assets/InterfaceScout.ico", names["Windows"])
        self.assertNotIn("assets/InterfaceScout.ico", names["macOS"])
        self.assertNotIn("assets/InterfaceScout.ico", names["Linux"])
        self.assertIn("Create Desktop Shortcut.bat", names["Windows"])
        self.assertIn("launchers/create_windows_shortcut.ps1", names["Windows"])
        self.assertNotIn("Create Desktop Shortcut.bat", names["Linux"])
        self.assertIn("Create Desktop App.command", names["macOS"])
        self.assertNotIn("Create Desktop App.command", names["Windows"])
        self.assertIn("Create Desktop Shortcut.sh", names["Linux"])
        self.assertNotIn("Create Desktop Shortcut.sh", names["macOS"])

    def test_frontend_uses_only_canonical_logo_and_favicon(self):
        html = (ROOT / "frontend" / "index.html").read_text(encoding="utf-8")
        self.assertIn('href="/favicon.png"', html)
        self.assertIn('src="/logo.png"', html)


if __name__ == "__main__":
    unittest.main()
