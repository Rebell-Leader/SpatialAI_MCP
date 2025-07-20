"""Setup script for development installation."""

from setuptools import setup, find_packages

setup(
    name="openproblems-spatial-mcp",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "fastmcp>=2.0.0",
        "pyyaml>=6.0",
        "click>=8.1.0",
        "psutil>=5.9.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-asyncio>=0.21.0",
            "pytest-mock>=3.10.0",
        ],
        "spatial": [
            "spatialdata",
            "zarr",
            "anndata",
            "h5py",
        ]
    }
)
