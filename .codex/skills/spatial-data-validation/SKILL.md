---
name: spatial-data-validation
description: Validate and inspect spatial transcriptomics data files (SpatialData/zarr/AnnData) before any analysis. Use when the user points at a spatial data file, reports a loading/format error, or asks whether data is suitable for a method.
---

# Spatial data validation

Validate first, reason second. Never assume a file's format or contents from its
extension.

## Steps

1. **Detect & validate format.** Call the MCP tool `validate_spatial_data` with
   the file path. Choose `validation_level`:
   - `basic` — extension/existence only.
   - `structure` (default) — format + structural integrity.
   - `integrity` — deeper content checks.
   - `domain` — spatial-biology-specific expectations.
   For several files, use `validate_multiple_spatial_files`.

2. **Read the metadata.** Call `analyze_spatial_metadata` to extract dimensions,
   spatial coordinates, and gene/feature info. Confirm `obsm['spatial']` (AnnData)
   or coordinate systems (SpatialData) are present when spatial analysis is
   intended.

3. **Check raw vs. processed.** Integer-valued matrices are likely raw counts;
   non-integer/log-scaled values are processed. State which the downstream method
   expects and flag a mismatch.

4. **Check compatibility** when combining files: `check_spatial_data_compatibility`
   (coordinate systems + gene overlap).

## Key facts (see `context/data-formats.md`)

- A `.zarr` store is SpatialData **only** if it has SpatialData markers or ≥2 of
  `images/ labels/ points/ shapes/ tables/`. Otherwise it's a plain zarr array.
- iST datasets can be terabytes. Default to lightweight structural checks; avoid
  loading whole objects. The heavy scientific stack (`pip install -e ".[spatial]"`)
  is optional — the validators degrade gracefully without it.

## Don't

- Don't claim a file is valid SpatialData without inspecting contents.
- Don't call execution tools (e.g. `run_nextflow_workflow`) — they aren't
  implemented. Drive local CLIs directly if execution is needed.
