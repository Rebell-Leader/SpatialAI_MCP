#!/usr/bin/env python3
"""Debug script for format detection."""

import sys
from pathlib import Path

# Add src to Python path
project_root = Path(__file__).parent
src_path = project_root / "src"
sys.path.insert(0, str(src_path))

def debug_format_detection():
    """Debug the format detection logic."""
    from openproblems_mcp.spatial_validation import SpatialDataValidator

    validator = SpatialDataValidator()

    # Test the problematic case
    path = Path("test.spatialdata")
    print(f"Testing path: {path}")
    print(f"Path suffix: {path.suffix}")
    print(f"Path suffix lower: {path.suffix.lower()}")
    print(f"SPATIALDATA_EXTENSIONS: {validator.SPATIALDATA_EXTENSIONS}")
    print(f"ZARR_EXTENSIONS: {validator.ZARR_EXTENSIONS}")

    # Check conditions
    suffix = path.suffix.lower()
    print(f"suffix in SPATIALDATA_EXTENSIONS: {suffix in validator.SPATIALDATA_EXTENSIONS}")
    print(f"file_path.name.endswith('.spatialdata'): {path.name.endswith('.spatialdata')}")

    # Test detection
    result = validator._detect_format(path)
    print(f"Detection result: {result}")

    return result

if __name__ == "__main__":
    debug_format_detection()
