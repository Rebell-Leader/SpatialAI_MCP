# Spatial data formats: contract for agents

The server's validators (`src/openproblems_mcp/spatial_validation.py`) recognize
three formats. This file is the reference for what "valid" means and which format
to expect where.

## SpatialData (zarr) — preferred for imaging-based spatial

A `SpatialData` object is persisted as a **zarr directory** (often named
`*.zarr` or `*.spatialdata`). It composes several element types:

- `images/` — raster image layers (e.g. DAPI, stains).
- `labels/` — segmentation masks (raster, integer-labeled).
- `points/` — transcript locations (per-molecule coordinates).
- `shapes/` — polygons/boundaries (e.g. cell shapes).
- `tables/` — `AnnData` tables annotating elements (cell-by-gene matrices).

Detection: a `.zarr` store is only classified as SpatialData after inspecting
its contents (SpatialData markers in `.zattrs`, or ≥2 of the expected element
directories). A plain zarr array without those markers is classified as `zarr`,
not `spatialdata` — see `_detect_format` / `_is_spatialdata_zarr`. **Do not
assume every `.zarr` is a SpatialData object.**

Coordinate systems are first-class in SpatialData. Elements are related through
named coordinate systems (e.g. `global`); validation should confirm they exist
and are consistent across elements being analyzed together.

## AnnData (.h5ad / .h5 / .hdf5)

An annotated cell-by-gene matrix: `X` (counts/expression), `obs` (cell
metadata), `var` (gene metadata), `obsm` (e.g. `spatial` coordinates), `layers`,
`uns`. Used for the table layer inside SpatialData and for some standalone
benchmark datasets.

Key checks: matrix shape consistency, presence of `obsm['spatial']` for spatial
data, raw-counts vs. normalized distinction (integer-valued `X` ⇒ likely raw).

## zarr (plain array store)

A generic zarr store (`.zarray` / `.zgroup` metadata) that is not a SpatialData
object. Validated structurally (presence of `.zarray`, readable chunks) without
spatial-domain assumptions.

## Raw vs. processed data

Many methods require **raw counts**; others require normalized/log-transformed
input. Always state which a method expects and check the input accordingly —
silently feeding normalized data to a method expecting raw counts is a common,
hard-to-spot error.

## Scale warning

iST datasets are large (often 10–100× scRNA-seq, up to terabytes per run).
Prefer lazy/structural inspection over loading whole objects into memory; the
server's lightweight checks work without the heavy `[spatial]` dependency stack
and should be the default for large files.
