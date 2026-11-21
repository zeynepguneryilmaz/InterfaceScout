#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
if [ "${1:-}" = "--check" ]; then bash -n "$0"; echo "Linux desktop shortcut check passed"; exit 0; fi
if [ "$(uname -s)" != "Linux" ]; then echo "This launcher is for Linux." >&2; exit 1; fi
if [ ! -f "$ROOT/frontend/logo.png" ] || [ ! -f "$ROOT/start.sh" ]; then
  echo "Extract the complete InterfaceScout Linux ZIP first." >&2; exit 1
fi
DESKTOP_DIR="${XDG_DESKTOP_DIR:-$HOME/Desktop}"
mkdir -p "$DESKTOP_DIR" "$HOME/.local/share/applications"
# Write a freedesktop launcher. Quote the installation path for spaces.
ENTRY="$HOME/.local/share/applications/interfacescout.desktop"
{
  printf '[Desktop Entry]\nVersion=1.0\nType=Application\nName=InterfaceScout\nComment=Physicochemical protein-surface mapping\n'
  printf 'Exec=/bin/bash "%s/start.sh"\n' "$ROOT"
  printf 'Path=%s\n' "$ROOT"
  printf 'Icon=%s/frontend/logo.png\n' "$ROOT"
  printf 'Terminal=false\nCategories=Science;Education;\nStartupNotify=true\n'
} > "$ENTRY"
chmod +x "$ENTRY"
cp "$ENTRY" "$DESKTOP_DIR/InterfaceScout.desktop"
chmod +x "$DESKTOP_DIR/InterfaceScout.desktop"
command -v gio >/dev/null 2>&1 && gio set "$DESKTOP_DIR/InterfaceScout.desktop" metadata::trusted true 2>/dev/null || true
echo "InterfaceScout desktop shortcut created: $DESKTOP_DIR/InterfaceScout.desktop"
echo "On some desktops, select 'Allow Launching' the first time."
