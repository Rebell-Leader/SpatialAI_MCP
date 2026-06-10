# task_ist_preprocessing: pipeline contract

Reference: <https://github.com/openproblems-bio/task_ist_preprocessing>

Benchmarks preprocessing of **imaging-based spatial transcriptomics (iST)** data
(e.g. Xenium, MERFISH, CosMx). Input is raw imaging-derived data (images +
detected transcripts); the pipeline produces a clean cell-by-gene table suitable
for downstream analysis, and scores the preprocessing choices.

## Repository shape (verify against upstream before relying on paths)

- `src/` — component implementations (the processing stages below).
- `common/` — shared OpenProblems resources, pulled as a **git submodule**.
- `main.nf` — Nextflow workflow tying components together.
- `_viash.yaml` — project-level Viash config (namespace, package settings).
- `scripts/` — helper scripts; `.github/` — CI workflows.
- Languages: Python (majority), Shell, Nextflow, some R.

## Processing stages (the component pipeline)

1. **Data preprocessing** — ingest/standardize raw iST data into SpatialData.
2. **Segmentation** — define cell boundaries (e.g. from nuclear/membrane stains).
3. **Transcript assignment** — assign detected transcripts to segmented cells.
4. **Cell-type annotation** — label cells by type.
5. **Expression correction** — correct for technical artifacts (e.g. spillover).
6. **Normalization** — normalize the resulting expression matrix.
7. **QC filtering** — drop low-quality cells/genes.
8. **Cell volume calculation** — compute per-cell volume.
9. **Count aggregation** — aggregate transcripts into the cell-by-gene matrix.
10. **Metrics** — evaluate the preprocessed result against the benchmark metrics.

Each stage is one or more Viash components conforming to a `comp_*` API
interface (merged into each component config via `__merge__`). A method
contribution implements one stage's interface.

## What the co-pilot should do here

- Validate that inputs/outputs are SpatialData (zarr) with the elements a given
  stage expects (e.g. segmentation needs `images`; assignment needs `points` +
  `labels`).
- Check `config.vsh.yaml` against the stage's `comp_*` interface before building.
- Help scaffold a new component for one stage (see
  `skills/viash-component-authoring/SKILL.md`).
- When a Nextflow run fails, triage logs (see `skills/nextflow-debugging`).

> The stage list and repo layout reflect the upstream repository and will drift.
> Treat this as orientation, and confirm exact component names, the `comp_*`
> interface, and paths against the live repo before generating code.
