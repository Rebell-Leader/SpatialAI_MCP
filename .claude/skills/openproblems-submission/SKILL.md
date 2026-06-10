---
name: openproblems-submission
description: Prepare a spatial transcriptomics method for contribution to an OpenProblems benchmark (e.g. task_ist_preprocessing). Use when the user wants to submit/contribute a method, open a PR to a task_* repo, or check submission readiness.
---

# OpenProblems submission

Goal: a Viash component that conforms to a task's API, builds cleanly, passes the
test profile, and is reproducible — ready for a PR.

## Checklist

1. **Repo setup.** Clone the target `task_*` repo and initialize the `common/`
   git submodule (`git submodule update --init --recursive`). Read the repo's
   `CONTRIBUTING` / `README` for current conventions.

2. **Conform to the component API.** Your `config.vsh.yaml` must `__merge__` the
   correct `comp_<stage>.yaml` interface and match its inputs/outputs and
   `info`/metadata fields. Run `analyze_workflow_configuration` on the config.

3. **Data contract.** Inputs/outputs are SpatialData (zarr) / AnnData as the stage
   requires. Validate test data with `validate_spatial_data` and confirm the
   elements the stage needs are present (`check_spatial_data_compatibility` for
   multi-file). State raw-vs-normalized expectations explicitly.

4. **Reproducibility.** Pin the engine image and package versions in the config.
   No unpinned `latest`. Document every parameter with a default and biological
   rationale.

5. **Build & test locally** (server doesn't run these — use the terminal):
   - `viash ns build`
   - `viash run config.vsh.yaml -- --input <test>.zarr --output tmp/out.zarr`
   - `nextflow run main.nf -profile test,docker`

6. **Dependency sanity.** Run `analyze_workflow_dependencies` across the workflow
   files to catch missing/conflicting containers or libraries.

7. **PR hygiene.** Include test data + expected output, document parameter
   choices, and confirm CI (GitHub Actions) builds and runs the component.

## Don't

- Don't claim the server can build/benchmark/validate-submission automatically —
  `build_openproblems_method`, `run_openproblems_benchmark`, and
  `validate_openproblems_submission` are roadmap. Use the validation/analysis
  tools that exist, and run the build/benchmark steps via local CLIs.
