"""Minimal viash test: the normalized table should be log-scaled (non-integer,
bounded) rather than raw counts."""

import subprocess
import sys
from os import path

import numpy as np
import spatialdata as sd

input_path = meta["resources_dir"] + "/dataset.zarr"  # noqa: F821 (viash-injected)
output_path = "output.zarr"

subprocess.run(
    [meta["executable"], "--input", input_path, "--output", output_path],  # noqa: F821
    check=True,
)

assert path.exists(output_path), "no output written"
sdata = sd.read_zarr(output_path)
x = sdata.tables["table"].X
maximum = x.max()
assert maximum < 50, f"values look un-normalized (max={maximum})"
assert not np.allclose(
    x.data % 1 if hasattr(x, "data") else x % 1, 0
), "output still looks like raw integer counts"
print("test passed")
sys.exit(0)
