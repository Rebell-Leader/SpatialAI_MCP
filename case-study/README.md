# Case study: does the skill let a smaller model match a bigger one?

**Question.** Can an open-source coding model **with** the OpenProblems spatial
co-pilot skill match a frontier commercial model **without** it — on a real
`task_ist_preprocessing` task — while needing **fewer agent steps** and **fewer
human validations**?

This directory is a runnable experiment, not a claim. It ships the task, the
grading rubric, an automated grader, and reference outputs. You run the four
arms and fill in `results/RESULTS.md`; the grader scores them objectively.

## The task (one realistic pipeline component)

Implement a **log-normalization method** as a Viash component for the
**normalization stage** of
[`task_ist_preprocessing`](https://github.com/openproblems-bio/task_ist_preprocessing).

Why this task:

- **Real and representative.** Normalization is an actual stage in the iST
  pipeline; contributing a method is the canonical OpenProblems workflow.
- **Trivial algorithm, hard harness.** The math (normalize to counts-per-10k,
  then `log1p`) is a few lines. So success is almost entirely about *getting the
  OpenProblems/Viash conventions and the SpatialData contract right* — exactly
  what the skill encodes and a general agent gets wrong in predictable ways.
- **Objectively gradeable.** "Correct" is checkable: does it conform to the
  component API, use the right base image and data format, ship a test, and
  build? See the rubric.

The exact, frozen user prompt is in [`task/USER_PROMPT.md`](task/USER_PROMPT.md);
the shared environment/dataset every arm gets is in
[`task/SETUP.md`](task/SETUP.md).

## Experimental arms (2×2)

| Arm | Model | Skill installed? | Role |
| --- | --- | --- | --- |
| **A** | Commercial frontier (e.g. Claude Opus 4.x / GPT-class) | ❌ | Baseline — "what you get today" |
| **B** | Open-source via OpenRouter (e.g. Qwen3-Coder, DeepSeek-V3, Llama-3.x) | ❌ | Weak baseline |
| **C** | **Same open-source model as B** | ✅ | **Hero condition** |
| **D** | Commercial frontier (same as A) | ✅ | Ceiling |

**Primary hypothesis:** **C ≥ A** on conformance and success, at **fewer steps
and fewer human validations**. Secondary: D is the ceiling; C ≫ B isolates the
skill's contribution on a fixed model.

The model is held constant between B and C, so any difference is attributable to
the skill, not the model. The skill is held constant between A/D and B/C contrasts.

## What "with skill" vs "without skill" means (kept fair)

Both conditions get the **same** agent harness, the **same** repo checkout, the
**same** local tools (Viash/Docker), and the **same** max-turn budget. The only
difference:

- **Without skill:** the agent has the repo and a general web/search capability,
  but none of our assets. (This is the honest baseline — a competent agent can
  still read the repo.)
- **With skill:** the agent additionally has the installed skill + MCP server,
  i.e. `spatialai-install install --target <agent>` (or the plugin). That gives
  it `skills/viash-component-authoring`, the `context/` contracts, and the
  `validate_spatial_data` / `analyze_workflow_configuration` MCP tools.

To avoid leakage, "without skill" runs in a checkout with `AGENTS.md`, `skills/`,
`context/`, `CLAUDE.md`, `.cursor/`, `.github/copilot-instructions.md`, and
`.mcp.json` **removed** (see `scripts/strip_skill.sh`).

## Metrics (all objective)

| Metric | How measured | Why it matters |
| --- | --- | --- |
| **Conformance score** | `grade.py` static rubric (0–11) | Did it follow OpenProblems/Viash/SpatialData conventions? |
| **Build success** | `viash ns build` exit 0 | Does it actually compile? |
| **Test success** | `viash test` exit 0 | Does the shipped unit test pass? |
| **Task success (binary)** | conformance ≥ 9 **and** build **and** test pass | The headline pass/fail |
| **Agent steps** | count of agent turns (tool calls + messages) to first success | Efficiency / cost |
| **Human validations** | count of times the operator had to correct, redirect, or approve a risky step | The real UX cost to a biologist |
| **Tokens / wall-clock** (optional) | harness logs | Cost proxy |

"Human validations" is the metric biologists care about most: a plain agent that
*eventually* succeeds after five corrections is worse than one that succeeds in
two clean steps. We count an intervention whenever the operator types anything
other than "continue/approve" — e.g. "no, that's the old Viash syntax", "this
should read zarr, not h5ad".

## Predicted failure modes without the skill (the rubric items)

From comparing a real upstream component to typical agent output, a no-skill run
predictably:

1. emits **legacy Viash** (`functionality:` / `platforms:`) instead of the flat
   `engines:`/`runners:` config;
2. uses a **raw `python:3.11`** image instead of `openproblems/base_python:1` +
   setup partial;
3. **omits `__merge__`** of the stage's component API;
4. drops the **metadata block** (`label`/`summary`/`links`/`references`);
5. treats data as **AnnData/`.h5ad`** instead of **SpatialData/zarr**;
6. ships **no `viash test`**;
7. doesn't **validate the input** or distinguish raw vs. normalized.

Each is a rubric line. The skill addresses each directly, so the experiment
measures whether that knowledge transfer actually closes the model gap.

## How to run

```bash
# 0. Install the package and a Viash/Docker toolchain on the runner.
pip install -e ".[dev]"

# 1. For each arm, prepare a checkout:
#    - skill arms (C, D): spatialai-install install --target <agent>
#    - no-skill arms (A, B): bash case-study/scripts/strip_skill.sh <checkout>

# 2. Give the agent the frozen prompt (task/USER_PROMPT.md) in task/SETUP.md's
#    environment. Log every turn and every human intervention.

# 3. Grade the produced component directory:
python case-study/grade.py path/to/produced/component --json

# 4. Record build/test:
( cd <task_repo> && viash ns build && viash test path/to/config.vsh.yaml )

# 5. Fill results/RESULTS.md with the numbers from each arm.
```

The grader is demonstrated on two reference components in `examples/`:
`plain_agent_typical/` (what arm A/B tends to produce) scores low;
`skill_guided/` (what arm C/D should produce) scores high. Run
`python case-study/grade.py case-study/examples/skill_guided` to see it.

## Honesty notes

- The `examples/` are **reference fixtures** illustrating the rubric extremes and
  validating the grader — they are not model outputs and not results.
- `results/RESULTS.md` ships **empty** (a template). Do not report numbers you
  did not measure.
- The grader scores **static** conformance offline; **build/test** require a real
  Viash/Docker runner and are recorded separately.
- Confirm the exact component-API filename (`/src/api/comp_*`) against the live
  task repo before running — the grader matches the `comp_` prefix, not an exact
  name, to stay robust to upstream renames.
