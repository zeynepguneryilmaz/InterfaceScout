#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
if [ "${1:-}" = "--check" ]; then bash -n "$0"; echo "macOS app shortcut check passed"; exit 0; fi
if [ "$(uname -s)" != "Darwin" ]; then echo "This launcher is for macOS." >&2; exit 1; fi
if [ ! -f "$ROOT/frontend/logo.png" ] || [ ! -f "$ROOT/start.command" ]; then
  echo "Extract the complete InterfaceScout macOS ZIP first." >&2; exit 1
fi
APP="$HOME/Desktop/InterfaceScout.app"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
ICONSET="$(mktemp -d)/AppIcon.iconset"
mkdir -p "$ICONSET"
for n in 16 32 128 256 512; do
  sips -s format png -z "$n" "$n" "$ROOT/frontend/logo.png" --out "$ICONSET/icon_${n}x${n}.png" >/dev/null
done
for n in 16 32 128 256 512; do
  if [ "$n" -le 512 ]; then
    sips -s format png -z "$((2*n))" "$((2*n))" "$ROOT/frontend/logo.png" --out "$ICONSET/icon_${n}x${n}@2x.png" >/dev/null
  fi
done
iconutil -c icns "$ICONSET" -o "$APP/Contents/Resources/AppIcon.icns"
cat > "$APP/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
<key>CFBundleName</key><string>InterfaceScout</string>
<key>CFBundleDisplayName</key><string>InterfaceScout</string>
<key>CFBundleIdentifier</key><string>org.interfacescout.local</string>
<key>CFBundleExecutable</key><string>InterfaceScout</string>
<key>CFBundlePackageType</key><string>APPL</string>
<key>CFBundleIconFile</key><string>AppIcon</string>
<key>CFBundleVersion</key><string>1.1.0</string>
</dict></plist>
PLIST
printf '#!/usr/bin/env bash\nexec /bin/bash %q\n' "$ROOT/start.command" > "$APP/Contents/MacOS/InterfaceScout"
chmod +x "$APP/Contents/MacOS/InterfaceScout"
echo "Desktop app created: $APP"
echo "Keep the extracted InterfaceScout folder in place; the app launcher points there."
