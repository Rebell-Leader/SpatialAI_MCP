# Case-study runner

Drives the experiment arms through a coding-agent harness and grades the output
automatically. One command per config; results land as JSON + a markdown table.

## What it does per arm

1. Copies the task repo into an isolated `workdir/<arm>` (skipping `.git`).
2. **Skill arm:** copies `AGENTS.md` + `skills/` + `context/` in, then runs
   `python -m installer install --target <agent_target>` to write the
   harness-native files (`GEMINI.md` + `.gemini/`, `CLAUDE.md` + `.claude/`,
   etc.) and `.mcp.json`.
   **No-skill arm:** runs `scripts/strip_skill.sh` to guarantee a clean repo.
3. Runs the harness headless on the frozen prompt (`task/prompt.txt`).
4. Finds the produced `config.vsh.yaml` for the component and grades it with
   `grade.py` (optionally `--run-viash` for build/test).

It measures the **zero-intervention** scenario: one prompt, no operator help.
The **human-validations** metric comes from *interactive* runs you do by hand —
this runner deliberately doesn't fake those.

## Quick start

```bash
cp case-study/runner/arms.example.json case-study/runner/arms.json
# edit arms.json: set task_repo to your task_ist_preprocessing clone, set models

# validate the plan without running anything:
python case-study/runner/run.py --config case-study/runner/arms.json --dry-run

# run one arm, or all:
python case-study/runner/run.py --config case-study/runner/arms.json --only C_gemini_skill
python case-study/runner/run.py --config case-study/runner/arms.json
```

Output: `case-study/runs/results.json` and `results.md`.

## Smoke-test offline first (no tokens, no setup)

Before spending real tokens, validate the whole prepare → run → grade → table
pipeline with a fake agent:

```bash
python case-study/runner/run.py --config case-study/runner/arms.example.json \
  --mock-harness skill-aware
```

`--mock-harness`:

- **`skill-aware`** (recommended) — the fake agent writes the *good* reference
  component when skill assets are present in the workdir, and the *bad* one
  otherwise. So a full run reproduces the expected pattern (skill arms 11/11
  PASS, no-skill arms 0/11 FAIL) — proving the install/strip, grading, and
  reporting all work end-to-end.
- **`good`** / **`bad`** — always produce that fixture, regardless of skill.

If the configured `task_repo` doesn't exist, mock mode synthesizes a minimal
stand-in repo automatically, so the smoke test needs nothing but this repo.
Swap `--mock-harness ...` for the real run once it looks right.

## Harness setup

The `command` templates in the config are substituted with `{prompt}` and
`{model}`. Defaults provided:

- **Gemini CLI** — `["gemini", "-y", "-p", "{prompt}"]`. Headless via `-p`; `-y`
  auto-approves tool calls. Reads `GEMINI.md` + `.gemini/settings.json` (which the
  `gemini` installer target writes). Auth via `gemini` login or `GEMINI_API_KEY`.
- **opencode** — `["opencode", "run", "{prompt}", "--model", "{model}"]`. Reads
  `AGENTS.md` natively (the skill arm copies it in). Use an OpenRouter model id
  like `openrouter/qwen/qwen3-coder`; set your OpenRouter key per opencode's docs.

Both must be installed and on `PATH` with credentials configured. If a harness
binary is missing, the arm is recorded as not-run rather than crashing the batch.

## Recommended arms (matches case-study/README.md)

| Arm | Harness | Model | Skill |
| --- | --- | --- | :---: |
| A | gemini | commercial Gemini | ❌ |
| C | gemini | **same** Gemini | ✅ |
| B | opencode | OpenRouter OSS | ❌ |
| C′ | opencode | **same** OSS | ✅ |

The decisive contrast is **same model, skill toggled** (A vs C, B vs C′): it
isolates the skill's effect. Run N≥3 repeats and grade blind. The runner records
the conformance score, build/test (with `--run-viash`), and the agent log per
arm; copy the numbers into `results/RESULTS.md`.

## Notes

- Stdlib only — no extra dependencies.
- `run_viash: true` in the config makes the grader also run `viash ns build` /
  `viash test` from each workdir (needs Viash + Docker on the runner).
- Component discovery matches a `config.vsh.yaml` whose `name:` is the component,
  falling back to a path/text match — so it finds the output wherever the agent
  placed it under `src/`.
