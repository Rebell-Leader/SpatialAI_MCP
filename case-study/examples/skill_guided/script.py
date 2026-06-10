import sys

import numpy as np
import scanpy as sc
import spatialdata as sd

## VIASH START
par = {
    "input": "resources_test/task_ist_preprocessing/dataset.zarr",
    "output": "output.zarr",
    "n_target": 10000,
}
## VIASH END


def main() -> None:
    sdata = sd.read_zarr(par["input"])

    # The normalization stage operates on the cell-by-gene table.
    if "table" not in sdata.tables:
        raise ValueError(
            "Input SpatialData has no 'table' element to normalize. "
            "Run count aggregation first."
        )
    adata = sdata.tables["table"]

    # CP10k + log1p. Guard against the common raw-vs-normalized error: this method
    # expects raw counts.
    if not np.allclose(adata.X.data % 1, 0) if hasattr(adata.X, "data") else False:
        print("WARNING: input does not look like raw counts.", file=sys.stderr)

    sc.pp.normalize_total(adata, target_sum=par["n_target"])
    sc.pp.log1p(adata)

    sdata.tables["table"] = adata
    sdata.write(par["output"])
    print("Normalization complete.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # fail loudly with a non-zero exit
        print(f"Normalization failed: {exc}", file=sys.stderr)
        sys.exit(1)
