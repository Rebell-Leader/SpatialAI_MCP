# Testing Guide

How to run and extend the test suite. Configuration lives in `pyproject.toml`
under `[tool.pytest.ini_options]` (there is no `pytest.ini` or `setup.py`); it
sets `pythonpath = ["src"]`, so no manual path juggling is needed.

## Setup

```bash
pip install -e ".[dev]"          # pytest, pytest-mock, ruff, black, mypy
pip install -e ".[spatial]"      # optional: spatialdata, zarr, anndata, h5py
```

The suite is designed to run **without** the heavy `[spatial]` stack — validators
degrade to lightweight structural checks and tests use mocking where needed.

## Run the tests

```bash
python -m pytest                 # everything
python -m pytest -q              # quiet
python -m pytest tests/test_installer.py        # one file
python -m pytest -k validation                  # by keyword
```

## What's covered

| File | Area |
| --- | --- |
| `tests/test_spatial_validation.py` | `SpatialDataValidator` (format detection, integrity, levels) |
| `tests/test_metadata_analysis.py` | `BioinformaticsMetadataExtractor` |
| `tests/test_spatial_tools.py` | MCP tool interface (`SpatialMCPTools`) |
| `tests/test_fastmcp_server.py` | FastMCP server wiring |
| `tests/test_installer.py` | Provider-agnostic installer (every target + CLI) |
| `tests/test_plugin_manifest.py` | Claude Code plugin marketplace manifests |
| `tests/test_case_study.py` | Case-study grader (rubric discrimination) |
| `tests/test_case_study_runner.py` | Case-study runner (dry-run, mock harness, grading) |
| `tests/conftest.py` | Shared fixtures |

## Lint, format, type-check

```bash
ruff check installer/ src/
black --check installer/
mypy src/
```

The MCP server (`src/`) carries some pre-existing line-length debt; CI lints it
**report-only**. New code under `installer/` and the case-study tooling is held
to a clean bar (CI fails on violations there).

## CI (GitHub Actions, `.github/workflows/ci.yml`)

On every push/PR:

1. **Tests** on Python 3.10, 3.11, 3.12.
2. **Generated-files-in-sync** — re-runs `spatialai-install install --target all`
   and fails if the committed per-agent files drift from `AGENTS.md` + `skills/`.
   So after editing those, regenerate and commit, or CI will flag it.
3. **Lint** — strict on new code, report-only on legacy `src/`.

Weekly (scheduled): an **upstream-drift** check comparing our `context/`
assumptions to the live `task_ist_preprocessing` repo (report-only).

## Smoke-testing the case study offline

The case-study runner has a no-tokens mock mode that exercises the full
prepare → run → grade → table pipeline:

```bash
python case-study/runner/run.py --config case-study/runner/arms.example.json \
  --mock-harness skill-aware
```

## Adding tests

1. Put new tests in `tests/test_*.py` with `test_*` functions.
2. Prefer lightweight checks; mock heavy I/O rather than requiring `[spatial]`.
3. If you change `AGENTS.md` or `skills/`, re-run the installer so the in-sync
   CI job stays green.
4. Run `python -m pytest` before committing.
