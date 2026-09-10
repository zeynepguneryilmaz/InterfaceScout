#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$ROOT/backend"
if [ ! -f "$BACKEND/v2/api.py" ] || [ ! -f "$ROOT/frontend/index.html" ]; then
  echo "ERROR: InterfaceScout files are incomplete." >&2
  exit 1
fi
PY=""
for c in python3.12 python3.11 python3.10 python3; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c "import ssl" >/dev/null 2>&1; then PY="$c"; break; fi
done
if [ -z "$PY" ]; then echo "ERROR: Python 3.10-3.12 with SSL is required." >&2; exit 1; fi
cd "$BACKEND"
[ -d .venv ] || "$PY" -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
chmod +x "$ROOT/start.sh" "$ROOT/start.command" 2>/dev/null || true
if [ "$(uname -s)" = "Darwin" ]; then exec "$ROOT/start.command"; else exec "$ROOT/start.sh"; fi
