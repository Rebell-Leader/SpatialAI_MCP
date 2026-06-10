# Quickstart — a real win in ~5 minutes

This walks you from zero to a working spatial-data validation and a scaffolded
Viash component. Pick the install path for your agent, then follow the worked
example.

## 1. Install

```bash
git clone https://github.com/Rebell-Leader/SpatialAI_MCP.git
cd SpatialAI_MCP
pip install -e .
```

For **deep** spatial-data validation (reading SpatialData/zarr/AnnData
internals), also install the scientific stack:

```bash
pip install -e ".[spatial]"   # spatialdata, zarr, anndata, h5py
```

> Without `[spatial]`, validation still works but falls back to lightweight
> structural checks — which is the right default for terabyte-scale iST data.

## 2. Set up your agent

**Either** project the rules/skills/MCP config into your agent:

```bash
spatialai-install install --target claude    # or codex / cursor / copilot / all
```

**Or**, for Claude Code, install the whole thing as a plugin:

```text
/plugin marketplace add Rebell-Leader/SpatialAI_MCP
/plugin install openproblems-spatial@openproblems-spatial
```

The plugin bundles the four skills and the MCP server config. The MCP server
itself is the `openproblems-mcp-server` console script, so the `pip install`
from step 1 is still the prerequisite that makes the tools runnable.

## 3. Verify the server

```bash
openproblems-mcp check
```

You'll see which local tools (nextflow, viash, docker, …) are detected. Missing
tools are fine — the validation/analysis tools don't need them.

## 4. Worked example: validate, then scaffold

### a) Validate a spatial dataset

Point the validator at any SpatialData/zarr or AnnData file. In your agent, ask:

> "Validate `path/to/dataset.zarr` and tell me if it's ready for the
> segmentation stage."

The agent calls the `validate_spatial_data` MCP tool (and
`analyze_spatial_metadata`). You get back the detected format, structural
integrity, which SpatialData elements are present (`images/`, `labels/`,
`points/`, …), and whether the data looks like raw counts or normalized — the
distinction that silently breaks methods expecting one or the other.

Don't have data handy? Any zarr store works for a structure check, and the
[`task_ist_preprocessing`](https://github.com/openproblems-bio/task_ist_preprocessing)
test resources are a good source of real SpatialData objects.

### b) Scaffold a Viash component

Ask:

> "Scaffold a new transcript-assignment component for task_ist_preprocessing."

The `viash-component-authoring` skill activates and produces a `config.vsh.yaml`
that follows current OpenProblems conventions — flat config, `__merge__` of the
stage API, the shared `openproblems/base_python:1` base image with a setup
partial, the standard metadata block (`label`/`summary`/`links`/`references`),
and a nextflow runner with resource `directives` — plus a `script.py` with the
`## VIASH START` / `## VIASH END` block. Then build and test it with your local
Viash:

```bash
viash ns build
viash test config.vsh.yaml
```

## What just happened (and why it beats a bare agent)

- The agent did **not** guess the data format — it ran a real validator.
- It did **not** emit obsolete Viash syntax (`functionality:` / `platforms:`) or
  a raw `python:*` image — the skill pinned the modern, OpenProblems-aligned shape.
- It knew the iST pipeline stages and which elements each needs, from
  `context/ist-preprocessing-pipeline.md`.

## Where to go next

- `AGENTS.md` — the canonical rules every target is generated from.
- `skills/` — the task playbooks (validation, authoring, debugging, submission).
- `context/` — OpenProblems facts, the data-format contract, the pipeline map.
- `installer/README.md` — how the provider-agnostic projection works.
