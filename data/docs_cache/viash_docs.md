# Viash Component Architecture Guide

## Overview
Viash enables building reusable, portable components across Docker, native, and Nextflow platforms.

## Component Structure

### Configuration File (config.vsh.yaml)
```yaml
name: "spatial_qc"
description: "Spatial transcriptomics quality control component"

argument_groups:
  - name: "Input/Output"
    arguments:
      - name: "--input"
        type: "file"
        description: "Input spatial data (h5ad format)"
        required: true
      - name: "--output"
        type: "file"
        direction: "output"
        description: "Output filtered data"
        required: true

  - name: "Parameters"
    arguments:
      - name: "--min_genes"
        type: "integer"
        description: "Minimum genes per cell"
        default: 200

resources:
  - type: "python_script"
    path: "script.py"

platforms:
  - type: "docker"
    image: "quay.io/biocontainers/scanpy:1.9.1--pyhd8ed1ab_0"
  - type: "nextflow"
```

### Script Implementation
```python
import argparse
import scanpy as sc
import json

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--min_genes', type=int, default=200)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)
sc.pp.filter_cells(adata, min_genes=args.min_genes)
adata.write(args.output)
```

## Development Workflow
```bash
# Build component
viash build config.vsh.yaml -p docker

# Test component
viash test config.vsh.yaml

# Build for Nextflow
viash build config.vsh.yaml -p nextflow -o target/nextflow/
```

## Best Practices
1. **Single Responsibility**: Each component should do one thing well
2. **Clear Interfaces**: Well-defined inputs and outputs
3. **Comprehensive Testing**: Unit tests for all functionality
4. **Documentation**: Clear descriptions and examples
