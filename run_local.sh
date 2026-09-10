#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$ROOT/backend"
if [ ! -f "$BACKEND/app.py" ] || [ ! -f "$ROOT/frontend/index.html" ]; then
  echo "ERROR: InterfaceScout files are incomplete." >&2
  exit 1
fi
PY=""
for c in python3.12 python3.11 python3.10 python3; do
  if command -v "$c" >/dev/null 2>&1 && "$c" - <<'PY' >/dev/null 2>&1
import sys, ssl
raise SystemExit(0 if (3,10) <= sys.version_info[:2] <= (3,12) else 1)
PY
  then PY="$c"; break; fi
done
if [ -z "$PY" ]; then
  echo "ERROR: Python 3.10, 3.11, or 3.12 with SSL support is required." >&2
  exit 1
fi
if [ "${1:-}" = "--check" ]; then
  "$PY" -c "import sys; print('InterfaceScout launcher check:', sys.version.split()[0])"
  exit 0
fi
cd "$BACKEND"
[ -d .venv ] || "$PY" -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
if [ "$(uname -s)" = "Darwin" ]; then exec bash "$ROOT/start.command"; else exec bash "$ROOT/start.sh"; fi
