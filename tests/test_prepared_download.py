"""Prepared-download contract: native attachment, secure opaque token and ZIP persistence."""
from __future__ import annotations
from io import BytesIO
from pathlib import Path
from unittest.mock import patch
import sys
import unittest
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend"))
import app


class PreparedDownloadTests(unittest.TestCase):
    def test_server_prepares_direct_zip_attachment(self):
        b = BytesIO()
        with zipfile.ZipFile(b, "w") as z:
            z.writestr("README_EXPORT.txt", "All 11 chemistry maps preserved.")
        with patch.object(app, "build_export_zip", return_value=b.getvalue()):
            result = app.prepare_export_package(app.ExportRequest(
                structures=[{"name": "test", "pdb_text": "ATOM"}]))
        url = result["download_url"]
        self.assertRegex(url, r"^/export-package/download/[0-9a-f]{32}$")
        self.assertGreater(result["bytes"], 100)
        response = app.download_prepared_export(url.rsplit("/", 1)[-1])
        self.assertEqual(response.media_type, "application/zip")
        self.assertEqual(Path(response.path).read_bytes(), b.getvalue())
        self.assertEqual(response.headers.get("cache-control"), "no-store")
        self.assertEqual(response.filename, "InterfaceScout_Complete_Results.zip")

    def test_unknown_token_not_downloadable(self):
        from fastapi import HTTPException
        with self.assertRaises(HTTPException) as context:
            app.download_prepared_export("a" * 32)
        self.assertEqual(context.exception.status_code, 404)


if __name__ == "__main__":
    unittest.main()
