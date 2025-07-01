#!/usr/bin/env python3
"""
Simple Documentation Generator for OpenProblems MCP Server

Generates curated documentation for:
- Nextflow best practices
- Viash components
- OpenProblems guidelines
- Docker patterns
- Spatial workflow templates
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
            "nextflow": self._generate_nextflow_docs(),
            "viash": self._generate_viash_docs(),
            "openproblems": self._generate_openproblems_docs(),
            "docker": self._generate_docker_docs(),
            "spatial_templates": self._generate_spatial_templates()
        }

        # Save to cache
        print("🔄 Saving documentation to cache...")
        await self._save_documentation_cache(documentation)

        return documentation

    def _generate_nextflow_docs(self) -> str:
        """Generate Nextflow documentation as a JSON string."""
        nextflow_docs = {
            "overview": "Nextflow is a workflow framework for bioinformatics pipelines",
            "best_practices": {
                "dsl_version": "Use DSL2 for all new workflows",
                "resource_management": "Specify memory and CPU requirements for each process",
                "error_handling": "Implement retry strategies and error handling",
                "containerization": "Use Docker/Singularity containers for reproducibility",
            },
            "common_patterns": {
                "input_channels": "Use Channel.fromPath() for file inputs",
                "output_publishing": "Use publishDir directive for results",
                "conditional_execution": "Use when clause for conditional processes",
            },
            "troubleshooting": {
                "oom_errors": "Increase memory allocation or implement dynamic resource allocation",
                "missing_files": "Check file paths and ensure proper input staging",
                "container_issues": "Verify container availability and permissions",
            },
            "code_examples": {
                "basic_pipeline": "#!/usr/bin/env nextflow\nnextflow.enable.dsl=2\n...",
                "process_definition": "process SPATIAL_ANALYSIS { ... }"
            }
        }
        return json.dumps(nextflow_docs, indent=2)

    def _generate_viash_docs(self) -> str:
        """Generate Viash documentation as a JSON string."""
        viash_docs = {
            "overview": "Viash is a meta-framework for building reusable workflow modules",
            "component_structure": {
                "config_file": "YAML configuration defining component metadata",
                "script": "Core functionality implementation",
                "platforms": "Target platforms (docker, native, nextflow)",
            },
            "best_practices": {
                "modularity": "Keep components focused on single tasks",
                "documentation": "Provide clear descriptions and examples",
                "testing": "Include unit tests for all components",
                "versioning": "Use semantic versioning for component releases",
            },
            "common_commands": {
                "build": "viash build config.vsh.yaml",
                "run": "viash run config.vsh.yaml",
                "test": "viash test config.vsh.yaml",
                "ns_build": "viash ns build",
            },
            "code_examples": {
                "config_file": "name: 'spatial_qc'\ndescription: '...'\n...",
                "script_implementation": "import argparse\n..."
            }
        }
        return json.dumps(viash_docs, indent=2)

    def _generate_openproblems_docs(self) -> str:
        """Generate OpenProblems documentation."""
        return """# OpenProblems Framework Guide

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
"""

    def _generate_docker_docs(self) -> str:
        """Generate Docker documentation."""
        docker_docs = {
            "overview": "Docker best practices for bioinformatics workflows",
            "dockerfile_optimization": {
                "multi_stage_builds": "Use multi-stage builds to reduce image size",
                "base_images": "Use minimal base images like python:3.9-slim",
                "layer_caching": "Combine RUN commands to reduce layers",
                "user_security": "Create non-root users for security"
            },
            "bioinformatics_specific": {
                "dependencies": "Install common bio packages: scanpy, anndata, pandas",
                "resources": "Set appropriate memory and CPU limits",
                "nextflow_compatibility": "Ensure containers work with Nextflow",
                "health_checks": "Include health checks for services"
            },
            "common_patterns": {
                "python_bio": "FROM python:3.9-slim + bio packages",
                "nextflow_user": "Create user with uid 1000 for Nextflow",
                "apt_cleanup": "Remove apt cache to reduce size"
            }
        }
        return json.dumps(docker_docs, indent=2)

    def _generate_spatial_templates(self) -> str:
        """Generate spatial workflow templates."""
        spatial_templates = {
            "basic_preprocessing": {
                "name": "Basic Spatial Preprocessing",
                "description": "Quality control and basic preprocessing for spatial transcriptomics data",
                "inputs": ["spatial_data.h5ad"],
                "outputs": ["filtered_data.h5ad", "qc_metrics.json"],
                "workflow": "nextflow spatial_qc.nf",
                "parameters": {
                    "min_genes_per_cell": 200,
                    "min_cells_per_gene": 3
                }
            },
            "spatially_variable_genes": {
                "name": "Spatially Variable Gene Detection",
                "description": "Identify genes with spatial expression patterns",
                "inputs": ["spatial_data.h5ad"],
                "outputs": ["svg_results.h5ad", "spatial_features.csv"],
                "workflow": "nextflow svg_detection.nf",
                "parameters": {
                    "n_top_genes": 2000,
                    "spatial_key": "spatial"
                }
            },
            "label_transfer": {
                "name": "Cell Type Label Transfer",
                "description": "Transfer cell type labels from reference to spatial data",
                "inputs": ["spatial_data.h5ad", "reference_data.h5ad"],
                "outputs": ["labeled_spatial.h5ad", "transfer_scores.csv"],
                "workflow": "nextflow label_transfer.nf",
                "parameters": {
                    "reference_key": "cell_type",
                    "confidence_threshold": 0.5
                }
            }
        }
        return json.dumps(spatial_templates, indent=2)

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
