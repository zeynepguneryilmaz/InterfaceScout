"""Serve complete atlas at http://127.0.0.1:8765 (Python standard library only)."""
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
import os
os.chdir(Path(__file__).resolve().parent)
print("Open http://127.0.0.1:8765 in a browser",flush=True)
ThreadingHTTPServer(("127.0.0.1",8765),SimpleHTTPRequestHandler).serve_forever()
