#!/usr/bin/env python3
"""
Documentation Generator for OpenProblems MCP Server

Generates comprehensive, curated documentation for:
- Nextflow best practices and DSL2 patterns
- Viash component architecture and workflows
- OpenProblems project structure and guidelines
- Docker optimization for bioinformatics
- Spatial transcriptomics pipeline templates

This provides structured knowledge that complements Continue.dev's
real-time documentation access.
"""

import asyncio
import json
from pathlib import Path
from typing import Dict

class DocumentationGenerator:
    def __init__(self, cache_dir: str = "data/docs_cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    async def generate_all_documentation(self) -> Dict[str, str]:
        """Generate comprehensive curated documentation."""
        print("📚 Generating curated documentation for OpenProblems MCP Server...")

        documentation = {
            "nextflow": await self._generate_nextflow_docs(),
            "viash": await self._generate_viash_docs(),
            "openproblems": await self._generate_openproblems_docs(),
            "docker": await self._generate_docker_docs(),
            "spatial_templates": await self._generate_spatial_templates()
        }

        # Save to cache
        print("🔄 Saving documentation to cache...")
        await self._save_documentation_cache(documentation)

        return documentation

    async def _generate_nextflow_docs(self) -> str:
        """Generate comprehensive Nextflow DSL2 documentation and best practices."""
        return """# Nextflow DSL2 Best Practices Guide

## Overview
Nextflow enables scalable and reproducible scientific workflows using software containers.

## Essential DSL2 Patterns

### Basic Pipeline Structure
```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Pipeline parameters
params.input = './data/*.fastq'
params.output_dir = './results'

// Import modules
include { QUALITY_CONTROL } from './modules/qc.nf'
include { ALIGNMENT } from './modules/align.nf'

// Main workflow
workflow {
    // Create input channel
    input_ch = Channel.fromPath(params.input)

    // Execute processes
    QUALITY_CONTROL(input_ch)
    ALIGNMENT(QUALITY_CONTROL.out.trimmed)
}
```

### Process Definition Best Practices
```nextflow
process SPATIAL_ANALYSIS {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/scanpy:1.9.1--pyhd8ed1ab_0'
    publishDir "${params.output_dir}/spatial_analysis", mode: 'copy'

    input:
    tuple val(sample_id), path(spatial_data)

    output:
    tuple val(sample_id), path("${sample_id}_analyzed.h5ad"), emit: analyzed
    path "${sample_id}_metrics.json", emit: metrics

    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    import json

    # Load and analyze spatial data
    adata = sc.read_h5ad('${spatial_data}')

    # Spatial analysis workflow
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Save results
    adata.write('${sample_id}_analyzed.h5ad')

    # Generate metrics
    metrics = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'sample_id': '${sample_id}'
    }

    with open('${sample_id}_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    """
}
```

## Resource Management
```nextflow
// nextflow.config
process {
    withLabel: 'process_low' {
        cpus = 2
        memory = '4.GB'
        time = '1.h'
    }
    withLabel: 'process_medium' {
        cpus = 4
        memory = '8.GB'
        time = '2.h'
    }
    withLabel: 'process_high' {
        cpus = 8
        memory = '16.GB'
        time = '4.h'
    }
    withLabel: 'process_spatial' {
        cpus = 6
        memory = '12.GB'
        time = '3.h'
    }
}

docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}
```

## Error Handling and Retry Strategies
```nextflow
process ROBUST_PROCESS {
    errorStrategy 'retry'
    maxRetries 3

    script:
    '''
    # Process implementation with error handling
    set -euo pipefail

    # Your analysis code here
    '''
}
```

## Channel Operations for Spatial Data
```nextflow
// Pair spatial data with metadata
Channel.fromPath('*.h5ad')
    .map { file ->
        def sample_id = file.baseName
        return [sample_id, file]
    }
    .set { spatial_data_ch }

// Combine with reference data
spatial_data_ch
    .combine(Channel.fromPath(params.reference_data))
    .set { analysis_input_ch }
```

## Debugging and Monitoring
```bash
# Run with comprehensive logging
nextflow run pipeline.nf -with-trace -with-report -with-timeline -with-dag

# Resume interrupted runs
nextflow run pipeline.nf -resume

# Check specific work directory
ls work/a1/b2c3d4*/
```

## Common Issues and Solutions
1. **Out of Memory**: Increase memory allocation or use dynamic resources
2. **File Not Found**: Check file paths and ensure proper input staging
3. **Container Issues**: Verify container accessibility and user permissions
4. **Process Hanging**: Check resource requirements and time limits
"""

    async def _generate_viash_docs(self) -> str:
        """Generate comprehensive Viash component documentation."""
        return """# Viash Component Architecture Guide

## Overview
Viash enables building reusable, portable components that work across Docker, native, and Nextflow platforms.

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
        example: "spatial_data.h5ad"
      - name: "--output"
        type: "file"
        direction: "output"
        description: "Output filtered data"
        required: true
        example: "filtered_spatial.h5ad"
      - name: "--metrics_output"
        type: "file"
        direction: "output"
        description: "QC metrics JSON file"
        required: true

  - name: "Parameters"
    arguments:
      - name: "--min_genes"
        type: "integer"
        description: "Minimum genes per cell"
        default: 200
      - name: "--min_cells"
        type: "integer"
        description: "Minimum cells per gene"
        default: 3

resources:
  - type: "python_script"
    path: "script.py"

platforms:
  - type: "docker"
    image: "quay.io/biocontainers/scanpy:1.9.1--pyhd8ed1ab_0"
    setup:
      - type: "python"
        packages: ["anndata>=0.8.0", "pandas>=1.5.0"]
  - type: "nextflow"
    directives:
      label: ["process_medium"]
```

### Script Implementation
```python
# script.py
import argparse
import scanpy as sc
import pandas as pd
import json

# Parse arguments
parser = argparse.ArgumentParser(description='Spatial QC component')
parser.add_argument('--input', required=True, help='Input spatial data')
parser.add_argument('--output', required=True, help='Output filtered data')
parser.add_argument('--metrics_output', required=True, help='Metrics output')
parser.add_argument('--min_genes', type=int, default=200, help='Min genes per cell')
parser.add_argument('--min_cells', type=int, default=3, help='Min cells per gene')

args = parser.parse_args()

# Load spatial data
adata = sc.read_h5ad(args.input)

# Quality control
n_cells_before = adata.n_obs
n_genes_before = adata.n_vars

# Filter cells and genes
sc.pp.filter_cells(adata, min_genes=args.min_genes)
sc.pp.filter_genes(adata, min_cells=args.min_cells)

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# Save results
adata.write(args.output)

# Generate metrics
metrics = {
    'n_cells_before': int(n_cells_before),
    'n_cells_after': int(adata.n_obs),
    'n_genes_before': int(n_genes_before),
    'n_genes_after': int(adata.n_vars),
    'median_genes_per_cell': float(adata.obs['n_genes_by_counts'].median()),
    'median_counts_per_cell': float(adata.obs['total_counts'].median())
}

with open(args.metrics_output, 'w') as f:
    json.dump(metrics, f, indent=2)
```

## Development Workflow
```bash
# Build component for Docker
viash build config.vsh.yaml -p docker -o spatial_qc_docker

# Test component
viash test config.vsh.yaml

# Build for Nextflow
viash build config.vsh.yaml -p nextflow -o target/nextflow/

# Build all components in namespace
viash ns build --parallel
```

## Integration Patterns

### With Nextflow
```nextflow
// Include built Viash component
include { SPATIAL_QC } from './target/nextflow/spatial_qc/main.nf'

workflow {
    input_ch = Channel.fromPath(params.input)
    SPATIAL_QC(input_ch)
}
```

### Component Testing
```yaml
# Add to config.vsh.yaml
test_resources:
  - type: "python_script"
    path: "test_component.py"
  - path: "test_data.h5ad"
    dest: "test_data.h5ad"

tests:
  - name: "basic_test"
    script: "test_component.py"
    expect:
      - type: "file"
        name: "output.h5ad"
```

## Best Practices
1. **Single Responsibility**: Each component should do one thing well
2. **Clear Interfaces**: Well-defined inputs, outputs, and parameters
3. **Comprehensive Testing**: Unit tests for all functionality
4. **Documentation**: Clear descriptions, examples, and parameter explanations
5. **Version Control**: Use semantic versioning for component releases
"""

    async def _generate_openproblems_docs(self) -> str:
        """Generate OpenProblems project documentation."""
        return """# OpenProblems Framework Guide

## Overview
OpenProblems is a community effort to benchmark single-cell and spatial transcriptomics analysis methods.

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
│   ├── datasets/           # Common dataset loaders
│   └── metrics/            # Shared metrics
└── workflows/              # Nextflow workflows
```

### Component Types

#### Dataset Components
```yaml
name: "openproblems_spatial_dataset"
description: "Load spatial transcriptomics benchmark dataset"

argument_groups:
  - name: "Output"
    arguments:
      - name: "--output_spatial"
        type: "file"
        direction: "output"
        description: "Spatial expression matrix (h5ad)"
      - name: "--output_reference"
        type: "file"
        direction: "output"
        description: "Reference single-cell data (h5ad)"
      - name: "--output_solution"
        type: "file"
        direction: "output"
        description: "Ground truth solution (h5ad)"

platforms:
  - type: "docker"
    image: "openproblems/base_python:1.0.0"
  - type: "nextflow"
```

#### Method Components
```yaml
name: "spatial_decomposition_method"
description: "Spatial cell type decomposition method"

argument_groups:
  - name: "Input"
    arguments:
      - name: "--input_spatial"
        type: "file"
        description: "Spatial expression data"
        required: true
      - name: "--input_reference"
        type: "file"
        description: "Reference single-cell data"
        required: true

  - name: "Output"
    arguments:
      - name: "--output_proportions"
        type: "file"
        direction: "output"
        description: "Cell type proportions per spot"
        required: true
```

#### Metric Components
```yaml
name: "spatial_decomposition_metric"
description: "Evaluate spatial decomposition accuracy"

argument_groups:
  - name: "Input"
    arguments:
      - name: "--input_proportions"
        type: "file"
        description: "Predicted proportions"
      - name: "--input_solution"
        type: "file"
        description: "Ground truth proportions"

  - name: "Output"
    arguments:
      - name: "--output_scores"
        type: "file"
        direction: "output"
        description: "Evaluation scores"
```

## Data Formats

### AnnData Structure
```python
import anndata as ad

# Spatial data structure
adata_spatial = ad.read_h5ad('spatial_data.h5ad')
# adata_spatial.X: expression matrix
# adata_spatial.obs: spot metadata (including spatial coordinates)
# adata_spatial.var: gene metadata
# adata_spatial.obsm['spatial']: spatial coordinates

# Reference single-cell data
adata_reference = ad.read_h5ad('reference_data.h5ad')
# adata_reference.obs['cell_type']: cell type annotations
```

### Standard Metadata Fields
- **Cell types**: `obs['cell_type']`
- **Spatial coordinates**: `obsm['spatial']`
- **Batch information**: `obs['batch']`
- **Dataset information**: `uns['dataset_id']`

## Development Guidelines

### Component Implementation
```python
# Standard imports for OpenProblems
import anndata as ad
import pandas as pd
import numpy as np
from scipy import sparse

def main(input_spatial, input_reference, output_proportions):
    # Load data
    adata_spatial = ad.read_h5ad(input_spatial)
    adata_reference = ad.read_h5ad(input_reference)

    # Get common genes
    common_genes = adata_spatial.var_names.intersection(adata_reference.var_names)
    adata_spatial = adata_spatial[:, common_genes]
    adata_reference = adata_reference[:, common_genes]

    # Method implementation here
    # ...

    # Create output proportions matrix
    cell_types = adata_reference.obs['cell_type'].unique()
    proportions = pd.DataFrame(
        data=predicted_proportions,  # Your method output
        index=adata_spatial.obs_names,
        columns=cell_types
    )

    # Save as AnnData
    adata_out = ad.AnnData(
        X=proportions.values,
        obs=adata_spatial.obs,
        var=pd.DataFrame(index=cell_types)
    )
    adata_out.write(output_proportions)
```

### Testing Framework
```bash
# Test individual component
viash test src/tasks/spatial_decomposition/methods/method_name/config.vsh.yaml

# Run full benchmark pipeline
nextflow run . \\
  --input datasets/spatial_dataset.h5ad \\
  --output results/ \\
  --publish_dir_mode copy

# Evaluate results
python scripts/evaluate_benchmark.py --results results/
```

## Contribution Workflow
1. **Fork repository** from GitHub
2. **Create feature branch** for your method/metric
3. **Implement component** following templates
4. **Add comprehensive tests** and documentation
5. **Submit pull request** with benchmark results
6. **Participate in review** process with community

## Best Practices
- Follow OpenProblems naming conventions
- Use standard data formats (AnnData h5ad)
- Include comprehensive documentation
- Provide example data and expected outputs
- Ensure reproducibility across platforms
"""

    async def _generate_docker_docs(self) -> str:
        """Generate Docker best practices for bioinformatics."""
        return """# Docker Best Practices for Bioinformatics

## Multi-stage Builds for Spatial Analysis

### Optimized Python + R Environment
```dockerfile
# Build stage - compile dependencies
FROM python:3.9-slim as builder
WORKDIR /build

# Install build dependencies
RUN apt-get update && apt-get install -y \\
    build-essential \\
    gcc \\
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt

# Production stage - minimal runtime
FROM python:3.9-slim
WORKDIR /app

# Copy only installed packages
COPY --from=builder /root/.local /root/.local

# Install R and system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \\
    r-base \\
    procps \\
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('Seurat', 'SingleCellExperiment'), repos='https://cloud.r-project.org')"

# Create non-root user for security
RUN groupadd -g 1000 biouser && useradd -u 1000 -g biouser biouser
USER biouser
```

### Bioinformatics-Specific Patterns

#### Scanpy + Spatial Analysis Stack
```dockerfile
FROM python:3.9-slim

# System dependencies for spatial analysis
RUN apt-get update && apt-get install -y --no-install-recommends \\
    libhdf5-dev \\
    libffi-dev \\
    libblas-dev \\
    liblapack-dev \\
    gfortran \\
    && rm -rf /var/lib/apt/lists/*

# Python spatial transcriptomics stack
RUN pip install --no-cache-dir \\
    scanpy>=1.9.0 \\
    squidpy>=1.2.0 \\
    anndata>=0.8.0 \\
    pandas>=1.5.0 \\
    numpy>=1.21.0 \\
    scipy>=1.9.0 \\
    matplotlib>=3.5.0 \\
    seaborn>=0.11.0

WORKDIR /app
```

#### Conda-based Environment
```dockerfile
FROM continuumio/miniconda3:latest

# Copy environment specification
COPY environment.yml /tmp/environment.yml

# Create conda environment
RUN conda env create -f /tmp/environment.yml && \\
    conda clean -afy

# Activate environment in shell
SHELL ["conda", "run", "-n", "spatial-env", "/bin/bash", "-c"]

# Set environment as default
ENV PATH /opt/conda/envs/spatial-env/bin:$PATH
```

#### OpenProblems Compatible Container
```dockerfile
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \\
    procps \\
    curl \\
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics Python stack
RUN pip install --no-cache-dir \\
    anndata>=0.8.0 \\
    scanpy>=1.9.0 \\
    pandas>=1.5.0 \\
    numpy>=1.21.0 \\
    scipy>=1.9.0 \\
    scikit-learn>=1.1.0

# Create non-root user (required for Nextflow)
RUN groupadd -g 1000 nextflow && \\
    useradd -u 1000 -g nextflow -s /bin/bash nextflow

USER nextflow
WORKDIR /app

# Set Python entrypoint
ENTRYPOINT ["python"]
```

## Security and Performance Best Practices

### Dockerfile Optimization
```dockerfile
# Use specific versions for reproducibility
FROM python:3.9.7-slim

# Combine RUN commands to reduce layers
RUN apt-get update && apt-get install -y --no-install-recommends \\
    package1 \\
    package2 \\
    && rm -rf /var/lib/apt/lists/* \\
    && pip install --no-cache-dir package3

# Use .dockerignore to reduce build context
# Add to .dockerignore:
# .git
# __pycache__
# *.pyc
# .pytest_cache
# work/
# results/
```

### Resource Management
```dockerfile
# Add health check for long-running containers
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \\
    CMD python -c "import scanpy; print('healthy')" || exit 1

# Use init system for proper signal handling
RUN apt-get update && apt-get install -y --no-install-recommends tini
ENTRYPOINT ["tini", "--"]
CMD ["python", "analysis.py"]
```

### Memory and Storage Optimization
```dockerfile
# Use multi-stage builds to reduce final image size
FROM python:3.9-slim as deps
RUN pip install large-package

FROM python:3.9-slim as runtime
COPY --from=deps /usr/local/lib/python3.9/site-packages /usr/local/lib/python3.9/site-packages

# For large datasets, use volume mounts
VOLUME ["/data", "/results"]
```

## Container Usage Examples

### Local Development
```bash
# Build spatial analysis container
docker build -t spatial-analysis:latest .

# Run with volume mounts for data
docker run -v $(pwd)/data:/data -v $(pwd)/results:/results \\
    spatial-analysis:latest script.py --input /data/spatial.h5ad
```

### Nextflow Integration
```nextflow
process SPATIAL_ANALYSIS {
    container 'spatial-analysis:latest'

    input:
    path spatial_data

    output:
    path "analysis_results.h5ad"

    script:
    """
    python /app/spatial_analysis.py \\
        --input ${spatial_data} \\
        --output analysis_results.h5ad
    """
}
```

### Production Considerations
- Pin all software versions for reproducibility
- Use official base images when possible
- Minimize attack surface with minimal base images
- Implement proper logging and monitoring
- Use health checks for service containers
- Set appropriate resource limits in orchestration
"""

    async def _generate_spatial_templates(self) -> str:
        """Generate spatial transcriptomics workflow templates."""
        return """# Spatial Transcriptomics Pipeline Templates

## 1. Complete Quality Control Workflow

```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Pipeline parameters
params.input_pattern = "*.h5ad"
params.output_dir = "./results"
params.min_genes_per_cell = 200
params.min_cells_per_gene = 3
params.max_pct_mt = 20

process SPATIAL_QC {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/scanpy:1.9.1--pyhd8ed1ab_0'
    publishDir "${params.output_dir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(spatial_data)

    output:
    tuple val(sample_id), path("${sample_id}_qc.h5ad"), emit: filtered_data
    path "${sample_id}_qc_metrics.json", emit: metrics
    path "${sample_id}_qc_plots.pdf", emit: plots

    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    import pandas as pd
    import json
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    # Configure scanpy
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=80, facecolor='white')

    # Load spatial data
    adata = sc.read_h5ad('${spatial_data}')

    # Store original counts
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # Generate QC plots
    with PdfPages('${sample_id}_qc_plots.pdf') as pdf:
        # Basic statistics
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Total counts per cell
        sc.pl.violin(adata, ['total_counts'], jitter=0.4, ax=axes[0,0])
        axes[0,0].set_title('Total counts per cell')

        # Number of genes per cell
        sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4, ax=axes[0,1])
        axes[0,1].set_title('Number of genes per cell')

        # Mitochondrial gene percentage
        sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, ax=axes[1,0])
        axes[1,0].set_title('Mitochondrial gene %')

        # Ribosomal gene percentage
        sc.pl.violin(adata, ['pct_counts_ribo'], jitter=0.4, ax=axes[1,1])
        axes[1,1].set_title('Ribosomal gene %')

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

        # Spatial plots if coordinates available
        if 'spatial' in adata.obsm:
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))

            sc.pl.spatial(adata, color='total_counts', ax=axes[0,0], show=False)
            axes[0,0].set_title('Total counts')

            sc.pl.spatial(adata, color='n_genes_by_counts', ax=axes[0,1], show=False)
            axes[0,1].set_title('Number of genes')

            sc.pl.spatial(adata, color='pct_counts_mt', ax=axes[1,0], show=False)
            axes[1,0].set_title('Mitochondrial %')

            sc.pl.spatial(adata, color='pct_counts_ribo', ax=axes[1,1], show=False)
            axes[1,1].set_title('Ribosomal %')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

    # Apply filters
    sc.pp.filter_cells(adata, min_genes=${params.min_genes_per_cell})
    sc.pp.filter_genes(adata, min_cells=${params.min_cells_per_gene})

    # Filter by mitochondrial percentage
    adata = adata[adata.obs.pct_counts_mt < ${params.max_pct_mt}].copy()

    # Save filtered data
    adata.write('${sample_id}_qc.h5ad')

    # Generate summary metrics
    metrics = {
        'sample_id': '${sample_id}',
        'n_cells_before': int(n_cells_before),
        'n_cells_after': int(adata.n_obs),
        'n_genes_before': int(n_genes_before),
        'n_genes_after': int(adata.n_vars),
        'cells_filtered': int(n_cells_before - adata.n_obs),
        'genes_filtered': int(n_genes_before - adata.n_vars),
        'median_genes_per_cell': float(adata.obs['n_genes_by_counts'].median()),
        'median_counts_per_cell': float(adata.obs['total_counts'].median()),
        'median_mt_percent': float(adata.obs['pct_counts_mt'].median())
    }

    with open('${sample_id}_qc_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    """
}

workflow SPATIAL_QC_WORKFLOW {
    take:
    spatial_files_ch

    main:
    // Execute QC for each sample
    SPATIAL_QC(spatial_files_ch)

    emit:
    filtered_data = SPATIAL_QC.out.filtered_data
    metrics = SPATIAL_QC.out.metrics
    plots = SPATIAL_QC.out.plots
}

workflow {
    // Create input channel from file pattern
    input_ch = Channel.fromPath(params.input_pattern)
        .map { file ->
            def sample_id = file.baseName.replaceAll(/\\.h5ad$/, '')
            return [sample_id, file]
        }

    // Run QC workflow
    SPATIAL_QC_WORKFLOW(input_ch)

    // Collect metrics for summary report
    SPATIAL_QC_WORKFLOW.out.metrics
        .collectFile(name: 'qc_summary.json', storeDir: params.output_dir)
}
```

## 2. Spatial Cell Type Decomposition Pipeline

```nextflow
process SPATIAL_DECOMPOSITION {
    tag "$sample_id"
    label 'process_high'
    container 'openproblems/spatial-decomposition:latest'
    publishDir "${params.output_dir}/decomposition", mode: 'copy'

    input:
    tuple val(sample_id), path(spatial_data), path(reference_data)

    output:
    tuple val(sample_id), path("${sample_id}_decomposition.h5ad"), emit: results
    path "${sample_id}_proportions.csv", emit: proportions
    path "${sample_id}_decomp_metrics.json", emit: metrics

    script:
    """
    #!/usr/bin/env python
    import anndata as ad
    import pandas as pd
    import numpy as np
    import scanpy as sc
    from scipy.spatial.distance import pdist, squareform
    import json

    # Load data
    adata_spatial = ad.read_h5ad('${spatial_data}')
    adata_reference = ad.read_h5ad('${reference_data}')

    print(f"Spatial data: {adata_spatial.shape}")
    print(f"Reference data: {adata_reference.shape}")

    # Find common genes
    common_genes = adata_spatial.var_names.intersection(adata_reference.var_names)
    print(f"Common genes: {len(common_genes)}")

    adata_spatial = adata_spatial[:, common_genes].copy()
    adata_reference = adata_reference[:, common_genes].copy()

    # Get cell types from reference
    cell_types = adata_reference.obs['cell_type'].unique()
    print(f"Cell types: {cell_types}")

    # Placeholder decomposition (replace with actual method)
    # In practice, use methods like Cell2location, SpatialDWLS, etc.
    n_spots = adata_spatial.n_obs
    n_cell_types = len(cell_types)

    # Generate random proportions (replace with actual algorithm)
    np.random.seed(42)
    proportions_matrix = np.random.dirichlet(np.ones(n_cell_types), size=n_spots)

    # Create proportions DataFrame
    proportions_df = pd.DataFrame(
        proportions_matrix,
        columns=cell_types,
        index=adata_spatial.obs_names
    )

    # Add spatial coordinates if available
    if 'spatial' in adata_spatial.obsm:
        coords = adata_spatial.obsm['spatial']
        proportions_df['x_coord'] = coords[:, 0]
        proportions_df['y_coord'] = coords[:, 1]

    # Save proportions
    proportions_df.to_csv('${sample_id}_proportions.csv')

    # Add proportions to spatial data
    for cell_type in cell_types:
        adata_spatial.obs[f'prop_{cell_type}'] = proportions_df[cell_type].values

    # Calculate spatial autocorrelation if coordinates available
    spatial_metrics = {}
    if 'spatial' in adata_spatial.obsm:
        coords = adata_spatial.obsm['spatial']

        # Calculate pairwise distances
        distances = squareform(pdist(coords))

        # Simple spatial autocorrelation for each cell type
        for cell_type in cell_types:
            props = proportions_df[cell_type].values
            # Simplified Moran's I calculation
            n = len(props)
            mean_prop = np.mean(props)

            # Weight matrix (inverse distance, with cutoff)
            W = 1.0 / (distances + 1e-10)
            W[distances > np.percentile(distances, 10)] = 0  # Keep only close neighbors
            W = W / W.sum(axis=1, keepdims=True)  # Normalize

            # Moran's I
            numerator = np.sum(W * np.outer(props - mean_prop, props - mean_prop))
            denominator = np.sum((props - mean_prop) ** 2)

            if denominator > 0:
                morans_i = (n / np.sum(W)) * (numerator / denominator)
                spatial_metrics[f'morans_i_{cell_type}'] = float(morans_i)

    # Save results
    adata_spatial.write('${sample_id}_decomposition.h5ad')

    # Generate metrics
    metrics = {
        'sample_id': '${sample_id}',
        'n_spots': int(adata_spatial.n_obs),
        'n_genes': int(adata_spatial.n_vars),
        'n_cell_types': int(len(cell_types)),
        'cell_types': list(cell_types),
        'mean_entropy': float(np.mean(-np.sum(proportions_matrix * np.log(proportions_matrix + 1e-10), axis=1))),
        **spatial_metrics
    }

    with open('${sample_id}_decomp_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    """
}

workflow SPATIAL_DECOMPOSITION_WORKFLOW {
    take:
    spatial_ch
    reference_ch

    main:
    // Combine spatial data with reference
    input_ch = spatial_ch.combine(reference_ch)

    // Run decomposition
    SPATIAL_DECOMPOSITION(input_ch)

    emit:
    results = SPATIAL_DECOMPOSITION.out.results
    proportions = SPATIAL_DECOMPOSITION.out.proportions
    metrics = SPATIAL_DECOMPOSITION.out.metrics
}
```

## 3. Comprehensive Spatial Analysis Configuration

```nextflow
// nextflow.config
params {
    // Input/Output
    input_dir = './data'
    output_dir = './results'
    reference_data = './reference/reference_atlas.h5ad'

    // QC parameters
    min_genes_per_cell = 200
    min_cells_per_gene = 3
    max_pct_mt = 20

    // Analysis parameters
    n_top_genes = 2000
    resolution = 0.5

    // Visualization
    generate_plots = true
    plot_format = 'pdf'
}

// Process resource allocation
process {
    withLabel: 'process_low' {
        cpus = 2
        memory = '4.GB'
        time = '1.h'
    }

    withLabel: 'process_medium' {
        cpus = 4
        memory = '8.GB'
        time = '2.h'
    }

    withLabel: 'process_high' {
        cpus = 8
        memory = '16.GB'
        time = '4.h'
    }

    withLabel: 'process_spatial' {
        cpus = 6
        memory = '12.GB'
        time = '3.h'
    }
}

// Execution profiles
profiles {
    standard {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    cluster {
        process.executor = 'slurm'
        process.queue = 'compute'
        singularity.enabled = true
    }

    test {
        params.input_dir = './test_data'
        params.output_dir = './test_results'
    }
}

// Resource monitoring
trace {
    enabled = true
    file = "${params.output_dir}/trace.txt"
}

report {
    enabled = true
    file = "${params.output_dir}/report.html"
}

timeline {
    enabled = true
    file = "${params.output_dir}/timeline.html"
}

dag {
    enabled = true
    file = "${params.output_dir}/dag.svg"
}
```

## 4. Integration with OpenProblems Benchmarking

```nextflow
// OpenProblems-compatible spatial workflow
include { LOAD_DATASET } from './modules/openproblems/datasets.nf'
include { RUN_METHOD } from './modules/openproblems/methods.nf'
include { CALCULATE_METRICS } from './modules/openproblems/metrics.nf'

workflow OPENPROBLEMS_SPATIAL_BENCHMARK {
    // Load benchmark datasets
    LOAD_DATASET()

    // Run multiple methods
    methods_ch = Channel.from(['cell2location', 'rctd', 'spatialdecon'])

    methods_ch
        .combine(LOAD_DATASET.out.spatial)
        .combine(LOAD_DATASET.out.reference)
        .set { method_input_ch }

    RUN_METHOD(method_input_ch)

    // Calculate evaluation metrics
    RUN_METHOD.out.results
        .combine(LOAD_DATASET.out.solution)
        .set { metrics_input_ch }

    CALCULATE_METRICS(metrics_input_ch)

    // Aggregate results
    CALCULATE_METRICS.out.scores
        .collectFile(name: 'benchmark_results.csv', storeDir: params.output_dir)
}
```

This comprehensive set of templates provides:

1. **Production-ready QC pipeline** with comprehensive filtering and reporting
2. **Spatial decomposition workflow** with built-in evaluation metrics
3. **Flexible configuration** for different computing environments
4. **OpenProblems integration** for standardized benchmarking
5. **Comprehensive monitoring** and resource tracking
"""

    async def _save_documentation_cache(self, documentation: Dict[str, str]):
        """Save documentation to cache files."""
        for source, content in documentation.items():
            cache_file = self.cache_dir / f"{source}_docs.md"
            with open(cache_file, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"   💾 Cached {source} documentation ({len(content):,} chars)")

    async def load_cached_documentation(self) -> Dict[str, str]:
        """Load documentation from cache if available."""
        documentation = {}

        for source in ["nextflow", "viash", "openproblems", "docker", "spatial_templates"]:
            cache_file = self.cache_dir / f"{source}_docs.md"
            if cache_file.exists():
                with open(cache_file, 'r', encoding='utf-8') as f:
                    documentation[source] = f.read()

        return documentation

async def main():
    """Main function to generate and cache documentation."""
    print("📚 OpenProblems Documentation Generator")
    print("=" * 50)

    generator = DocumentationGenerator()

    print("🔄 Generating curated documentation...")
    documentation = await generator.generate_all_documentation()

    print(f"\n📊 Documentation generation complete!")
    total_chars = 0
    for source, content in documentation.items():
        chars = len(content)
        total_chars += chars
        print(f"   ✅ {source}: {chars:,} characters")

    print(f"\n🎉 Total: {total_chars:,} characters of documentation cached!")
    print("   💾 Documentation saved to: data/docs_cache/")
    print("   🔗 Now available via MCP Resources in your server")

    return documentation

if __name__ == "__main__":
    asyncio.run(main())
