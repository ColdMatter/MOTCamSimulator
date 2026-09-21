#!/usr/bin/env bash
# Rebuild CaF_MOT_Camera_Comparison.nb, the thresholds CSV and verification/*.png
# from CaFCameraNoise.wl, running the verification suite. Exit code 1 if any test fails.
set -euo pipefail
cd "$(dirname "$0")"

find_wolframscript() {
  if command -v wolframscript >/dev/null 2>&1; then command -v wolframscript; return; fi
  for app in /Applications/Wolfram.app /Applications/Mathematica.app; do
    [ -x "$app/Contents/MacOS/wolframscript" ] && { echo "$app/Contents/MacOS/wolframscript"; return; }
  done
  for p in /usr/local/bin/wolframscript /opt/Wolfram/*/Executables/wolframscript /usr/local/Wolfram/*/*/Executables/wolframscript; do
    [ -x "$p" ] && { echo "$p"; return; }
  done
  return 1
}

WS="$(find_wolframscript)" || { echo "wolframscript not found: install Wolfram/Mathematica or add it to PATH" >&2; exit 2; }
echo "Using $WS"
exec "$WS" -file GenerateCaFComparisonNotebook.wls
