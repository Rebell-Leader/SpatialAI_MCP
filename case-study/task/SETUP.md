# Shared environment (identical for every arm)

Keep everything below constant across arms A–D so the only variables are
**model** and **skill present/absent**.

## Runner

- OS with Docker available and running.
- [Viash](https://viash.io) installed (pin one version, e.g. `viash 0.9.x`) and
  on `PATH`.
- Python 3.10+.
- The agent harness of your choice, but the **same harness** for all arms, with
  the same tool permissions (file read/write, terminal, web).

## Task repository

- A fresh clone of `openproblems-bio/task_ist_preprocessing` at a **pinned
  commit** (record the SHA in `results/RESULTS.md`).
- `git submodule update --init --recursive` so `common/` is present.
- Confirm the normalization component API exists (e.g.
  `src/api/comp_*normalization*.yaml`) and note its exact path; the agent is
  expected to discover it.

## Skill present vs. absent

- **Skill arms (C, D):** from the co-pilot repo, run
  `spatialai-install install --target <agent> --dest <task_repo_checkout>`
  (or install the plugin). This adds the rules, the `viash-component-authoring`
  skill, the `context/` contracts, and the MCP server config.
- **No-skill arms (A, B):** run `bash case-study/scripts/strip_skill.sh
  <task_repo_checkout>` to guarantee none of our assets are present.

## Budgets and stopping

- **Max agent turns:** 12 (a turn = one model response, including its tool calls).
- **Operator script:** after the frozen prompt, only say `continue` to let the
  agent proceed. Every other message is logged as a human validation. Stop at
  task success, budget exhaustion, or give-up.

## Models (swap in your choices, record exact IDs)

- **Commercial frontier (A, D):** e.g. `claude-opus-4-x` or a GPT-class model.
- **Open-source via OpenRouter (B, C):** e.g. `qwen/qwen3-coder`,
  `deepseek/deepseek-v3`, or `meta-llama/llama-3.x-70b-instruct`. Use the **same**
  ID for B and C.

## Reproducibility

- Fixed temperature (e.g. 0) where the API allows it.
- Run **N ≥ 3 seeds/repeats** per arm; report median and range. Single runs are
  anecdotes.
- Blind the grading: grade produced components without knowing which arm made them.
