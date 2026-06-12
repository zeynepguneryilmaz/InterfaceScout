#!/usr/bin/env bash
# ============================================================
# InterfaceScout - Linux background launcher (no terminal window)
# Target of the InterfaceScout.desktop icon. Starts the backend
# detached and opens the browser at http://localhost:8000.
# Run run_local.sh once first for setup.
# This file sits in the InterfaceScout folder, next to backend/.
# ============================================================
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/backend" || exit 1

if [ ! -f ".venv/bin/activate" ]; then
  if command -v zenity >/dev/null 2>&1; then
    zenity --error --title="InterfaceScout" \
      --text="First-time setup needed.\nPlease run run_local.sh once, then use this icon." 2>/dev/null
  else
    echo "First-time setup needed. Run run_local.sh once." >&2
  fi
  exit 1
fi

# shellcheck disable=SC1091
source .venv/bin/activate

# free port 8000 if still bound from a previous run
if command -v lsof >/dev/null 2>&1; then
  PID=$(lsof -ti tcp:8000 2>/dev/null || true)
  [ -n "$PID" ] && kill -9 $PID 2>/dev/null || true
fi

# start backend detached; it opens the browser itself via webbrowser.open
nohup python main.py >/tmp/interfacescout.log 2>&1 &

sleep 2
if command -v xdg-open >/dev/null 2>&1; then
  xdg-open "http://localhost:8000" >/dev/null 2>&1 || true
fi
exit 0
