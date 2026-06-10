# Test suite

Unit and wiring tests for the MCP server, the provider-agnostic installer, and
the case-study tooling. For setup, how to run, lint, and CI, see the root
[`TESTING.md`](../TESTING.md).

## Files

| File | Area |
| --- | --- |
| `test_spatial_validation.py` | `SpatialDataValidator` — format detection, validation levels, integrity |
| `test_metadata_analysis.py` | `BioinformaticsMetadataExtractor` — Nextflow/Viash/spatial metadata |
| `test_spatial_tools.py` | `SpatialMCPTools` — MCP tool interface, result formatting, errors |
| `test_fastmcp_server.py` | FastMCP server wiring |
| `test_installer.py` | Installer — every target projection + CLI |
| `test_plugin_manifest.py` | Claude Code plugin marketplace manifests |
| `test_case_study.py` | Case-study grader (rubric discrimination) |
| `test_case_study_runner.py` | Case-study runner (dry-run, mock harness, grading) |
| `conftest.py` | Shared fixtures |

## Conventions

- Lightweight by design: the suite runs **without** the heavy `[spatial]` stack
  (spatialdata/zarr/anndata/h5py) by mocking external libraries, so it's fast and
  CI-friendly.
- Put new tests in `test_*.py` with `test_*` functions; mock expensive I/O.
- If you change `AGENTS.md` or `skills/`, re-run `spatialai-install` so the
  generated-files-in-sync CI job stays green.

Quick run:

```bash
python -m pytest -q
```
