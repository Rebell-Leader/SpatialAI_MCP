"""Pytest configuration and fixtures for spatial transcriptomics tests."""

import pytest
import tempfile
import json
import yaml
from pathlib import Path
from typing import Dict, Any


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)


@pytest.fixture
def sample_nextflow_config():
    """Sample Nextflow configuration content."""
    return """
    // Nextflow configuration file

    params {
        input = 'data/*.fastq'
        output = 'results'
        genome = 'hg38'
    }

    process {
        executor = 'local'
        cpus = 2
        memory = '4.GB'
        container = 'biocontainers/fastqc:v0.11.9_cv8'
    }

    profiles {
        standard {
            process.executor = 'local'
        }

        docker {
            docker.enabled = true
            process.container = 'ubuntu:20.04'
        }

        cluster {
            process.executor = 'slurm'
            process.queue = 'compute'
        }
    }
    """


@pytest.fixture
def sample_nextflow_workflow():
    """Sample Nextflow workflow content."""
    return """
    #!/usr/bin/env nextflow

    params.input = 'data/*.fastq'
    params.output = 'results'

    process QUALITY_CONTROL {
        container 'biocontainers/fastqc:v0.11.9_cv8'
        publishDir params.output, mode: 'copy'

        input:
        path reads

        output:
        path "*.html"
        path "*.zip"

        script:
        '''
        fastqc !{reads}
        '''
    }

    process TRIMMING {
        container 'biocontainers/trimmomatic:v0.39_cv1'

        input:
        path reads

        output:
        path "*_trimmed.fastq"

        script:
        '''
        trimmomatic SE !{reads} !{reads.baseName}_trimmed.fastq TRAILING:20 MINLEN:50
        '''
    }

    workflow {
        input:
            reads_ch = Channel.fromPath(params.input)

        main:
            QUALITY_CONTROL(reads_ch)
            TRIMMING(reads_ch)

        emit:
            qc_reports = QUALITY_CONTROL.out
            trimmed_reads = TRIMMING.out
    }
    """


@pytest.fixture
def sample_viash_config():
    """Sample Viash configuration."""
    return {
        "name": "spatial_analysis_component",
        "description": "Component for spatial transcriptomics analysis",
        "version": "1.0.0",
        "authors": [
            {"name": "Test Author", "email": "test@example.com"}
        ],
        "functionality": {
            "arguments": [
                {
                    "name": "input",
                    "type": "file",
                    "direction": "input",
                    "description": "Input spatial data file",
                    "required": True
                },
                {
                    "name": "output",
                    "type": "file",
                    "direction": "output",
                    "description": "Output analysis results",
                    "required": True
                },
                {
                    "name": "method",
                    "type": "string",
                    "default": "leiden",
                    "description": "Clustering method to use"
                }
            ],
            "resources": [
                {
                    "type": "python_script",
                    "path": "script.py"
                }
            ]
        },
        "platforms": [
            {
                "type": "docker",
                "image": "python:3.9-slim",
                "setup": [
                    {"type": "python", "packages": ["scanpy", "spatialdata", "pandas"]}
                ]
            },
            {
                "type": "native"
            }
        ]
    }


@pytest.fixture
def sample_python_script():
    """Sample Python script content."""
    return """
    #!/usr/bin/env python3

    import pandas as pd
    import numpy as np
    import scanpy as sc
    import spatialdata as sd
    from pathlib import Path
    import argparse

    def load_spatial_data(input_path):
        '''Load spatial transcriptomics data.'''
        if input_path.suffix == '.h5ad':
            return sc.read_h5ad(input_path)
        elif input_path.suffix == '.zarr':
            return sd.read_zarr(input_path)
        else:
            raise ValueError(f"Unsupported format: {input_path.suffix}")

    def preprocess_data(adata):
        '''Preprocess spatial data.'''
        # Basic preprocessing
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)

        # Normalization
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        return adata

    def perform_clustering(adata, method='leiden'):
        '''Perform clustering analysis.'''
        sc.pp.highly_variable_genes(adata)
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)

        if method == 'leiden':
            sc.tl.leiden(adata)
        elif method == 'louvain':
            sc.tl.louvain(adata)

        return adata

    def main():
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', required=True)
        parser.add_argument('--output', required=True)
        parser.add_argument('--method', default='leiden')
        args = parser.parse_args()

        # Load and process data
        adata = load_spatial_data(Path(args.input))
        adata = preprocess_data(adata)
        adata = perform_clustering(adata, args.method)

        # Save results
        adata.write_h5ad(args.output)
        print(f"Analysis complete. Results saved to {args.output}")

    if __name__ == '__main__':
        main()
    """


@pytest.fixture
def create_test_files(temp_dir):
    """Create various test files in temporary directory."""
    def _create_files(file_configs: Dict[str, Any]):
        """
        Create test files based on configuration.

        Args:
            file_configs: Dict with filename as key and content/config as value
        """
        created_files = {}

        for filename, config in file_configs.items():
            file_path = temp_dir / filename

            if filename.endswith('.yaml') or filename.endswith('.yml'):
                with open(file_path, 'w') as f:
                    yaml.dump(config, f)
            elif filename.endswith('.json'):
                with open(file_path, 'w') as f:
                    json.dump(config, f, indent=2)
            elif isinstance(config, str):
                with open(file_path, 'w') as f:
                    f.write(config)
            elif config is None:
                # Create empty file
                file_path.touch()
            else:
                # Create directory
                file_path.mkdir(parents=True, exist_ok=True)

            created_files[filename] = file_path

        return created_files

    return _create_files


@pytest.fixture
def mock_spatial_data_structure(temp_dir):
    """Create a mock SpatialData directory structure."""
    spatialdata_dir = temp_dir / "test.spatialdata"
    spatialdata_dir.mkdir()

    # Create SpatialData components
    (spatialdata_dir / "images").mkdir()
    (spatialdata_dir / "labels").mkdir()
    (spatialdata_dir / "points").mkdir()
    (spatialdata_dir / "shapes").mkdir()
    (spatialdata_dir / "table").mkdir()

    # Create metadata files
    zattrs = spatialdata_dir / ".zattrs"
    with open(zattrs, 'w') as f:
        json.dump({
            "spatialdata": "0.1.0",
            "created_by": "test_suite"
        }, f)

    # Create table structure
    table_dir = spatialdata_dir / "table"
    table_zarray = table_dir / ".zarray"
    with open(table_zarray, 'w') as f:
        json.dump({
            "shape": [1000, 2000],
            "chunks": [100, 200],
            "dtype": "float32",
            "compressor": None
        }, f)

    return spatialdata_dir


@pytest.fixture
def mock_anndata_file(temp_dir):
    """Create a mock AnnData file."""
    h5ad_file = temp_dir / "test.h5ad"

    # Create a minimal file (won't be valid HDF5, but tests file existence)
    with open(h5ad_file, 'wb') as f:
        # Write minimal HDF5-like header
        f.write(b'\x89HDF\r\n\x1a\n')
        f.write(b'\x00' * 100)  # Padding

    return h5ad_file


# Test markers
pytest.mark.unit = pytest.mark.unit
pytest.mark.integration = pytest.mark.integration
pytest.mark.slow = pytest.mark.slow
