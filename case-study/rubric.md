# Grading rubric

Each item is 1 point, scored automatically by [`grade.py`](grade.py). Maximum
**11** static points. **Task success** = static ≥ 9 **and** `viash ns build`
passes **and** `viash test` passes.

| # | Check (grade.py key) | Pass condition | Skill item that supplies it |
| - | --- | --- | --- |
| 1 | `modern_flat_config` | No legacy `functionality:` / `platforms:` keys | viash-component-authoring |
| 2 | `merges_comp_api` | `__merge__` references `/src/api/comp_*` | viash-component-authoring + pipeline context |
| 3 | `openproblems_base_image` | Docker engine uses `openproblems/base_*` | viash-component-authoring + AGENTS.md |
| 4 | `metadata_core` | `label` + `summary` + `description` present | viash-component-authoring |
| 5 | `metadata_links_references` | `links` + `references` present | viash-component-authoring |
| 6 | `pinned_versions` | No `:latest`; base image is tagged | AGENTS.md reproducibility rule |
| 7 | `nextflow_directives` | nextflow runner has `directives.label` | viash-component-authoring |
| 8 | `viash_param_block` | Script keeps `## VIASH START`/`## VIASH END` | viash-component-authoring |
| 9 | `spatialdata_io` | Script uses SpatialData/zarr I/O, not AnnData-only | data-formats context + spatial-data-validation |
| 10 | `fails_loudly` | Script exits non-zero / raises on error | AGENTS.md scripting rule |
| 11 | `has_unit_test` | Config declares a `test_resources` script | openproblems-submission + AGENTS.md |

## Runtime checks (need Viash/Docker; recorded separately)

| Check | Command | Pass condition |
| --- | --- | --- |
| Build | `viash ns build` | exit 0 |
| Test | `viash test config.vsh.yaml` | exit 0 |

`grade.py --run-viash <task_repo>` runs these when `viash` is on `PATH` and folds
the result into the report's `runtime` field.

## Reference scores (sanity check)

Running the grader on the shipped fixtures should give:

- `examples/plain_agent_typical/` → **0/11** (every convention violated)
- `examples/skill_guided/` → **11/11**

These bracket the rubric and prove the grader discriminates. They are *fixtures*,
not model outputs.
