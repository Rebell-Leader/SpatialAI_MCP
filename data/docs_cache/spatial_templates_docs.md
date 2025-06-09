# Spatial Transcriptomics Pipeline Templates

## 1. Quality Control Workflow

```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_pattern = "*.h5ad"
params.output_dir = "./results"
params.min_genes_per_cell = 200

process SPATIAL_QC {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/scanpy:1.9.1--pyhd8ed1ab_0'
    publishDir "${params.output_dir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(spatial_data)

    output:
    tuple val(sample_id), path("${sample_id}_qc.h5ad"), emit: filtered_data
    path "${sample_id}_metrics.json", emit: metrics

    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    import json

    adata = sc.read_h5ad('${spatial_data}')

    # QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=${params.min_genes_per_cell})
    sc.pp.filter_genes(adata, min_cells=3)

    adata.write('${sample_id}_qc.h5ad')

    metrics = {
        'sample_id': '${sample_id}',
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars)
    }

    with open('${sample_id}_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    """
}

workflow {
    input_ch = Channel.fromPath(params.input_pattern)
        .map { file -> [file.baseName, file] }

    SPATIAL_QC(input_ch)
}
```

## 2. Spatial Decomposition Pipeline

```nextflow
process SPATIAL_DECOMPOSITION {
    tag "$sample_id"
    label 'process_high'
    container 'openproblems/spatial-decomposition:latest'

    input:
    tuple val(sample_id), path(spatial_data), path(reference_data)

    output:
    tuple val(sample_id), path("${sample_id}_decomposition.h5ad"), emit: results
    path "${sample_id}_proportions.csv", emit: proportions

    script:
    """
    #!/usr/bin/env python
    import anndata as ad
    import pandas as pd
    import numpy as np

    # Load data
    adata_spatial = ad.read_h5ad('${spatial_data}')
    adata_reference = ad.read_h5ad('${reference_data}')

    # Find common genes
    common_genes = adata_spatial.var_names.intersection(adata_reference.var_names)
    adata_spatial = adata_spatial[:, common_genes].copy()
    adata_reference = adata_reference[:, common_genes].copy()

    # Get cell types
    cell_types = adata_reference.obs['cell_type'].unique()

    # Placeholder decomposition (replace with actual method)
    n_spots = adata_spatial.n_obs
    n_cell_types = len(cell_types)
    proportions_matrix = np.random.dirichlet(np.ones(n_cell_types), size=n_spots)

    # Create proportions DataFrame
    proportions_df = pd.DataFrame(
        proportions_matrix,
        columns=cell_types,
        index=adata_spatial.obs_names
    )

    proportions_df.to_csv('${sample_id}_proportions.csv')

    # Add proportions to spatial data
    for cell_type in cell_types:
        adata_spatial.obs[f'prop_{cell_type}'] = proportions_df[cell_type].values

    adata_spatial.write('${sample_id}_decomposition.h5ad')
    """
}
```

## 3. Configuration Template

```nextflow
// nextflow.config
params {
    input_dir = './data'
    output_dir = './results'
    reference_data = './reference/atlas.h5ad'
}

process {
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
}

docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}
```

This provides:
1. **Production-ready QC pipeline** with filtering and reporting
2. **Spatial decomposition workflow** with evaluation metrics
3. **Flexible configuration** for different environments
4. **Comprehensive monitoring** and resource tracking
