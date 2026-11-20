#!/usr/bin/env bash
# Render every docs/img/*.mmd to SVG with the system Chrome (no puppeteer download).
set -euo pipefail
cd "$(dirname "$0")"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
echo '{"executablePath":"/usr/bin/google-chrome","args":["--no-sandbox","--disable-gpu"]}' > "$T/p.json"
echo '{"theme":"base","themeVariables":{"fontSize":"20px"}}' > "$T/c.json"
for f in *.mmd; do
  PUPPETEER_SKIP_DOWNLOAD=1 npx -y @mermaid-js/mermaid-cli@11 -q -i "$f" -o "${f%.mmd}.svg" -p "$T/p.json" -c "$T/c.json" -b transparent
done
