# Provider-agnostic skill installer

Author the co-pilot's rules, skills, and context **once** in this repo, then
project them into whichever AI coding agent the user runs. The MCP server is
already protocol-standard and works everywhere; this installer handles the
*non*-standard part — each agent's rules/skills/config files.

## Source of truth

| Asset | File(s) | Role |
|---|---|---|
| Persistent rules | `AGENTS.md` (repo root) | Canonical, provider-neutral instructions |
| Skills | `skills/<name>/SKILL.md` | On-demand task playbooks (Claude/Codex format) |
| Domain context | `context/*.md` | OpenProblems facts, data-format & pipeline contracts |

Edit those. **Do not** hand-edit the generated provider files — re-run the
installer instead.

## Usage

```bash
pip install -e .                       # provides the `spatialai-install` command

spatialai-install list                 # show targets and discovered skills
spatialai-install install -t all       # install for every agent, into ./
spatialai-install install -t claude -t cursor
spatialai-install install -t cursor --dest /path/to/another/project
```

`python -m installer ...` works identically if you'd rather not install the
console script.

## What each target gets

| Target | Generated files |
|---|---|
| `claude` (Claude Code) | `CLAUDE.md` (pointer → `AGENTS.md` via `@import`), `.claude/skills/<name>/SKILL.md`, `.mcp.json` |
| `codex` (OpenAI Codex) | `AGENTS.md` (read natively), `.codex/skills/<name>/SKILL.md`, `.codex/mcp-config.toml` snippet |
| `cursor` (Cursor) | `.cursor/rules/openproblems.mdc` (always-applied, embeds `AGENTS.md`), `.cursor/rules/<skill>.mdc` (description-gated), `.cursor/mcp.json` |
| `copilot` (GitHub Copilot) | `.github/copilot-instructions.md` (self-contained) |

Claude Code and Codex support importing/native-reading `AGENTS.md`, so they get a
thin pointer. Cursor and Copilot have no import mechanism, so the installer
embeds the full `AGENTS.md` content into their files to keep them self-contained.

## Design notes

- Stdlib only — no third-party dependency, so the installer runs anywhere Python
  3.8+ does. SKILL.md frontmatter is parsed without a YAML library.
- Adding a skill = add `skills/<name>/SKILL.md` with `name`/`description`
  frontmatter, then re-run the installer; every target picks it up automatically.
- Adding a target = add one `install_*` function and a `Target` entry in
  `installer/targets.py`.
