# AGENTS.md — OpenProblems Spatial Transcriptomics Co-pilot

> Canonical, provider-neutral instructions for any AI coding agent working in this
> repository (Claude Code, OpenAI Codex, Cursor, GitHub Copilot, and others).
> This is the single source of truth. Tool-specific files (`CLAUDE.md`,
> `.cursor/rules/`, `.github/copilot-instructions.md`) are **generated** from this
> file and the `skills/` directory by `installer/` — do not edit them by hand.

## What this project is

A co-pilot for computational biologists working on the
[OpenProblems](https://openproblems.bio) spatial transcriptomics benchmarks. It
ships two things:

1. **An MCP server** (`src/openproblems_mcp/`) — domain-aware validation and
   analysis tools for spatial data and bioinformatics workflows, exposed over the
   Model Context Protocol so any MCP-capable agent can call them.
2. **A provider-agnostic skill installer** (`skills/`, `context/`, `installer/`) —
   author rules/skills/context once here, then project them into whichever agent
   the user runs.

## Ground truth: which MCP tools actually exist

Only call tools that are registered. As of now the server is **read-only
analysis** — it validates and inspects, it does **not** execute pipelines.

**Available now:**
`health_check`, `list_tools_status`, `get_server_info`,
`validate_spatial_data`, `validate_multiple_spatial_files`,
`analyze_spatial_metadata`, `check_spatial_data_compatibility`,
`extract_bioinformatics_metadata`, `analyze_workflow_configuration`,
`assess_data_quality`, `analyze_workflow_dependencies`.

**Roadmap — NOT implemented, do not call:**
`run_nextflow_workflow`, `run_viash_component`, `build_docker_image`,
`run_nf_test`, `analyze_nextflow_log`, `create_spatial_component`,
`setup_spatial_environment`, and the OpenProblems build/benchmark/submission
tools. For these actions, drive the user's local `nextflow` / `viash` / `docker`
CLIs directly via the terminal instead, and say so plainly.

## Domain rules (spatial transcriptomics)

- **Validate before processing.** Run `validate_spatial_data` on any input before
  reasoning about it. Distinguish raw counts vs. normalized/log-transformed data —
  methods have explicit expectations.
- **Prefer SpatialData/zarr.** The current OpenProblems imaging-based spatial
  transcriptomics (iST) benchmark is built on the `SpatialData` (zarr) format;
  `.h5ad` AnnData is still used for some tables. See
  `context/data-formats.md`.
- **Know the pipeline stages.** The `task_ist_preprocessing` benchmark is a
  pipeline: preprocessing → segmentation → transcript assignment → cell-type
  annotation → expression correction → normalization → QC filtering → count
  aggregation → metrics. See `context/ist-preprocessing-pipeline.md`.
- **Components are Viash.** Methods are packaged as Viash components with a
  `config.vsh.yaml` and a script; the namespace is built with `viash ns build`.
  See `skills/viash-component-authoring/SKILL.md`.

## Coding standards

- **Python**: target 3.10+ (the `fastmcp` dependency requires ≥3.10), format with
  `black` (line length 88), lint with `ruff`, type-check with `mypy`. Match
  existing module style.
- **Viash components (OpenProblems conventions)**: use the modern flat config
  (top-level `name`, `arguments`, `resources`, `engines`, `runners`) — not the
  legacy `functionality:` / `platforms:` blocks. Conform to the task's component
  API via `__merge__: /src/api/comp_<stage>.yaml`. Build on the shared base
  images (e.g. `openproblems/base_python:1`) with `__merge__`-ed setup partials,
  not raw `python:*` images. Include the standard metadata (`label`, `summary`,
  `description`, `links`, `references`). In scripts, keep the
  `## VIASH START` / `## VIASH END` block, use `SpatialData` objects for I/O,
  handle coordinate systems explicitly, and fail loudly with a non-zero exit.
- **Reproducibility & verification**: pin package versions and base-image tags
  (no unpinned `latest`); add a `viash test` unit test (`test_resources`) per
  component; document parameters with biological rationale. See
  `skills/viash-component-authoring/SKILL.md` and
  `skills/openproblems-submission/SKILL.md`.

## Working in this repo

- Install dev environment: `pip install -e ".[dev]"` (add `[spatial]` for the
  heavy scientific stack: spatialdata, zarr, anndata, h5py).
- Run tests: `python -m pytest` (config lives in `pyproject.toml`; there is no
  separate `pytest.ini` or `setup.py`).
- Lint/format/type-check: `ruff check src/`, `black src/`, `mypy src/`.
- Health check the server: `openproblems-mcp check`.
- Run the server (normally launched by the agent host over stdio):
  `openproblems-mcp-server`.

## Skills

Task-specific playbooks live in `skills/<name>/SKILL.md` (shared Claude/Codex
format: YAML frontmatter + markdown). Load the one matching the task:

- `spatial-data-validation` — validating and inspecting spatial data files.
- `viash-component-authoring` — writing a new Viash method component.
- `nextflow-debugging` — diagnosing Nextflow pipeline failures.
- `openproblems-submission` — preparing a method for an OpenProblems benchmark.

## Conventions for changes

- Keep `README.md` and this file honest about what is actually implemented. If
  you add or remove a tool, update both (and the tool lists above).
- Do not advertise execution tools as available until they have code and tests.
- New skills go in `skills/` and are surfaced everywhere by re-running the
  installer; do not hand-edit the generated provider files.
