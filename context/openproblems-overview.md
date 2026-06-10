# OpenProblems: context for agents

[OpenProblems](https://openproblems.bio) is a community-driven, extensible
benchmarking platform for single-cell and spatial genomics. Benchmarks are
organized as independent `task_*` repositories under the
[`openproblems-bio`](https://github.com/openproblems-bio) GitHub org, each
formalizing one open problem with standardized datasets, methods, control
methods, and metrics.

## Why this matters for the co-pilot

Computational biologists contributing to OpenProblems spend disproportionate
effort on *auxiliary* tooling — Viash packaging, Nextflow orchestration, Docker
environments, data-format wrangling — rather than the science. This co-pilot's
job is to absorb that auxiliary complexity: validate data, check workflow
configs, scaffold components, and surface the project's conventions so the
researcher stays focused on method development.

## Technology stack (shared across tasks)

- **Viash** — meta-framework that turns a script + a `config.vsh.yaml` into a
  standalone CLI, a Docker image, and a Nextflow module. Components are the unit
  of contribution.
- **Nextflow** — orchestrates components into benchmark pipelines (DSL2). Runs
  locally, on HPC, and on AWS Batch.
- **Docker** — per-component reproducible environments.
- **Data formats** — `SpatialData` (zarr) for imaging-based spatial data;
  `AnnData` (`.h5ad`) for cell-by-gene tables. See `data-formats.md`.

## Relevant spatial tasks

- **`task_ist_preprocessing`** — benchmarking preprocessing of imaging-based
  spatial transcriptomics (iST). The primary target of this co-pilot. See
  `ist-preprocessing-pipeline.md`.
- **`task_spatial_simulators`** — benchmarking spatial transcriptomics
  simulators.
- **`task_spatially_variable_genes`** — detecting genes whose expression varies
  across spatial regions.

## Contribution model (high level)

1. Fork/clone the relevant `task_*` repo; it pulls shared resources via a
   `common/` git submodule and defines components under `src/`.
2. Add a method as a Viash component (`config.vsh.yaml` + script), conforming to
   the task's component API (the `__merge__`-ed `comp_*` interface).
3. Build the namespace (`viash ns build`), run on the test profile, iterate.
4. Open a PR; CI builds and runs the component against test data.

> Keep this file updated as the upstream tasks evolve — it is the agent's mental
> model of OpenProblems. Verify specifics against the live repos before relying
> on exact paths or component names.
