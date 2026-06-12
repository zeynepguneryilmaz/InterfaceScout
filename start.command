#!/usr/bin/env bash
# ============================================================
# InterfaceScout - macOS launcher (double-clickable in Finder)
# Starts the backend in the background (no Terminal stays open)
# and opens the browser at http://localhost:8000.
# Run run_local.sh once first for setup.
# This file sits in the InterfaceScout folder, next to backend/.
# ============================================================
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/backend" || exit 1

if [ ! -f ".venv/bin/activate" ]; then
  osascript -e 'display dialog "First-time setup needed.\n\nPlease run run_local.sh once (right-click → Open), then use this launcher." buttons {"OK"} default button "OK" with title "InterfaceScout"' >/dev/null 2>&1
  exit 1
fi

# shellcheck disable=SC1091
source .venv/bin/activate

# free port 8000 if still bound from a previous run
if command -v lsof >/dev/null 2>&1; then
  PID=$(lsof -ti tcp:8000 2>/dev/null || true)
  [ -n "$PID" ] && kill -9 $PID 2>/dev/null || true
fi

# start backend in the background (detached); it opens the browser itself
nohup python main.py >/tmp/interfacescout.log 2>&1 &

sleep 2
open "http://localhost:8000" 2>/dev/null || true

# close the Terminal window that Finder opened for this .command
osascript -e 'tell application "Terminal" to close (every window whose name contains "start.command")' >/dev/null 2>&1 &
exit 0
