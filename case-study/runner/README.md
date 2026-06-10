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
# REQUIRED: edit arms.json -> set "task_repo" to your task_ist_preprocessing
# clone (and set the models). A real run errors clearly if it's left unset.

# validate the plan without running anything:
python case-study/runner/run.py --config case-study/runner/arms.json --dry-run

# run one arm, or all:
python case-study/runner/run.py --config case-study/runner/arms.json --only C_gemini_skill
python case-study/runner/run.py --config case-study/runner/arms.json
```

Output: `case-study/runs/results.json` and `results.md`.

> **Run it from anywhere.** Paths in the config (`copilot_repo`, `workdir`,
> `prompt_file`) resolve against the repo root, not your current directory — so
> `python runner/run.py ...` from inside `case-study/` works too. Only the
> `--config` path itself is taken relative to your shell's working directory.

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

## How a "skill arm" gets the skill

For skill arms the runner does two things: (1) installs the agent-native files
(via the installer `agent_target`) and copies `AGENTS.md` + `skills/` + `context/`
into the workdir, and (2) **injects the skill playbooks into the prompt**
(`inject_skills: true`). The injection matters because agents differ in whether
they auto-load skills — Gemini/Antigravity CLI won't open `SKILL.md` files on
their own. Injecting guarantees the skill *content* reaches the model, so the
experiment measures the value of the content rather than each agent's discovery
mechanism. Set `inject_skill_names` to a subset, or `inject_skills: false` to
test pure auto-load behavior.

## Repeats (beat single-run noise)

Set `repeats` in the config or pass `--repeats N` (N≥3 recommended). Each arm runs
N times in isolated `…__rep<i>` workdirs; the table reports the **median**, the
range, and the pass-rate. A single run is an anecdote — a 2-point swing is well
within sampling noise.

## Harness setup

The `command` templates in the config are substituted with `{prompt}` and
`{model}`. Defaults provided:

- **Antigravity CLI** — `["agy", "--yolo", "-p", "{prompt}"]`. Headless via `-p`;
  `--yolo` auto-approves tool calls. Reads `AGENTS.md` + `.agents/` (the
  `antigravity` installer target writes `.agents/mcp_config.json` + skills).
  Note: `agy -p` can emit no stdout when spawned as a subprocess, but it still
  writes files — which is what the grader reads, so grading is unaffected.
  (Gemini CLI is still available as the legacy `gemini` harness/target.)
- **opencode** — `["opencode", "run", "{prompt}", "--model", "{model}"]`. Reads
  `AGENTS.md` natively. Use an OpenRouter model id like
  `openrouter/qwen/qwen3-coder`; set your OpenRouter key per opencode's docs.

Both must be installed and on `PATH` with credentials configured. If a harness
binary is missing, the arm is recorded as not-run rather than crashing the batch.

## Recommended arms (matches case-study/README.md)

| Arm | Harness | Model | Skill |
| --- | --- | --- | :---: |
| A | antigravity | commercial Gemini | ❌ |
| C | antigravity | **same** Gemini | ✅ |
| B | opencode | OpenRouter OSS | ❌ |
| C′ | opencode | **same** OSS | ✅ |

The decisive contrast is **same model, skill toggled** (A vs C, B vs C′): it
isolates the skill's effect. Run N≥3 repeats and grade blind. Turn on
`run_viash: true` so `viash build`/`viash test` give the real correctness signal
alongside the static rubric. Copy the numbers into `results/RESULTS.md`.

> **Headroom matters.** If your commercial model already scores near 11/11
> without the skill (ceiling), the skill can't show a lift — focus the
> comparison on the weaker/OSS model, or choose a less common pipeline stage.

## Notes

- Stdlib only — no extra dependencies.
- `run_viash: true` in the config makes the grader also run `viash ns build` /
  `viash test` from each workdir (needs Viash + Docker on the runner).
- Component discovery matches a `config.vsh.yaml` whose `name:` is the component,
  falling back to a path/text match — so it finds the output wherever the agent
  placed it under `src/`.
