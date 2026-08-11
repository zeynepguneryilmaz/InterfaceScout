#!/usr/bin/env bash
# ============================================================
# InterfaceScout - macOS daily launcher
# ============================================================

HERE="$(cd "$(dirname "$0")" && pwd)"
BACKEND="$HERE/backend"
FRONTEND="$HERE/frontend"

if [ ! -f "$BACKEND/main.py" ] || [ ! -f "$FRONTEND/index.html" ]; then
  osascript -e 'display dialog "InterfaceScout files are incomplete. Keep start.command next to backend/ and frontend/." buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1
  exit 1
fi

if [ ! -f "$BACKEND/.venv/bin/activate" ]; then
  osascript -e 'display dialog "First-time setup is required. Run run_local.sh once." buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1
  exit 1
fi

# If InterfaceScout is already running, only open the browser.
if curl -fsS "http://127.0.0.1:8000/health" >/dev/null 2>&1; then
  open "http://localhost:8000" >/dev/null 2>&1 || true
  exit 0
fi

# Do not kill an unrelated service using port 8000.
if command -v lsof >/dev/null 2>&1 && lsof -ti tcp:8000 >/dev/null 2>&1; then
  osascript -e 'display dialog "Port 8000 is already in use by another process." buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1
  exit 1
fi

cd "$BACKEND" || exit 1
# shellcheck disable=SC1091
source .venv/bin/activate

nohup python main.py >/tmp/interfacescout.log 2>&1 &
sleep 2

if curl -fsS "http://127.0.0.1:8000/health" >/dev/null 2>&1; then
  open "http://localhost:8000" >/dev/null 2>&1 || true
else
  osascript -e 'display dialog "InterfaceScout did not start. Check /tmp/interfacescout.log." buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1
fi

exit 0
