# Nextflow DSL2 Best Practices Guide

## Overview
Nextflow enables scalable and reproducible scientific workflows using software containers.

## Essential DSL2 Patterns

### Basic Pipeline Structure
```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = './data/*.h5ad'
params.output_dir = './results'

workflow {
    input_ch = Channel.fromPath(params.input)
    PROCESS_NAME(input_ch)
}
```

### Process Definition
```nextflow
process SPATIAL_ANALYSIS {
    tag "$sample_id"
    label 'process_medium'
    container 'quay.io/biocontainers/scanpy:1.9.1--pyhd8ed1ab_0'
    publishDir "${params.output_dir}/analysis", mode: 'copy'

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

    adata = sc.read_h5ad('${spatial_data}')
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.write('${sample_id}_analyzed.h5ad')

    metrics = {'n_cells': adata.n_obs, 'n_genes': adata.n_vars}
    with open('${sample_id}_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)
    """
}
```

## Resource Management
```nextflow
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
}

docker {
    enabled = true
    runOptions = '-u $(id -u):$(id -g)'
}
```

## Error Handling
```nextflow
process ROBUST_PROCESS {
    errorStrategy 'retry'
    maxRetries 3

    script:
    """
    set -euo pipefail
    # Your analysis code here
    """
}
```

## Common Issues and Solutions
1. **Out of Memory**: Increase memory allocation
2. **File Not Found**: Check file paths and staging
3. **Container Issues**: Verify container accessibility
4. **Process Hanging**: Check resource requirements
