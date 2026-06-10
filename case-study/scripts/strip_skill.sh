#!/usr/bin/env bash
# Remove every co-pilot asset from a checkout so the "no-skill" arms (A, B) run
# against a clean repo. Idempotent. Usage: strip_skill.sh <checkout_dir>
set -euo pipefail

DEST="${1:?usage: strip_skill.sh <checkout_dir>}"
cd "$DEST"

rm -rf \
  AGENTS.md CLAUDE.md GEMINI.md QUICKSTART.md \
  skills/ context/ installer/ \
  .claude/ .codex/ .cursor/ .gemini/ .agents/ .claude-plugin/ \
  .github/copilot-instructions.md \
  .mcp.json

echo "Stripped co-pilot assets from $DEST."
echo "Remaining agent-config files (should be none):"
ls -1d AGENTS.md CLAUDE.md skills context .mcp.json .cursor .claude 2>/dev/null || \
  echo "  (clean)"
