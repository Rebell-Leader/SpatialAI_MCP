# OpenProblems Framework Guide

## Overview
OpenProblems is a community effort to benchmark single-cell and spatial transcriptomics methods.

## Project Architecture

### Repository Structure
```
src/
├── tasks/                    # Benchmark tasks
│   ├── spatial_decomposition/
│   │   ├── methods/         # Benchmark methods
│   │   ├── metrics/         # Evaluation metrics
│   │   └── datasets/        # Task datasets
│   └── other_tasks/
├── common/                  # Shared components
└── workflows/              # Nextflow workflows
```

### Component Types

#### Dataset Components
Load benchmark datasets with standardized formats.

#### Method Components
Implement spatial analysis methods following OpenProblems standards.

#### Metric Components
Evaluate method performance with standardized metrics.

## Data Formats

### AnnData Structure
```python
import anndata as ad

# Spatial data structure
adata_spatial = ad.read_h5ad('spatial_data.h5ad')
# adata_spatial.X: expression matrix
# adata_spatial.obs: spot metadata
# adata_spatial.var: gene metadata
# adata_spatial.obsm['spatial']: spatial coordinates

# Reference single-cell data
adata_reference = ad.read_h5ad('reference_data.h5ad')
# adata_reference.obs['cell_type']: cell type annotations
```

### Standard Metadata Fields
- **Cell types**: obs['cell_type']
- **Spatial coordinates**: obsm['spatial']
- **Batch information**: obs['batch']

## Best Practices
- Follow OpenProblems naming conventions
- Use standard data formats (AnnData h5ad)
- Include comprehensive documentation
- Ensure reproducibility across platforms
