{
  "basic_preprocessing": {
    "name": "Basic Spatial Preprocessing",
    "description": "Quality control and basic preprocessing for spatial transcriptomics data",
    "inputs": [
      "spatial_data.h5ad"
    ],
    "outputs": [
      "filtered_data.h5ad",
      "qc_metrics.json"
    ],
    "workflow": "nextflow spatial_qc.nf",
    "parameters": {
      "min_genes_per_cell": 200,
      "min_cells_per_gene": 3
    }
  },
  "spatially_variable_genes": {
    "name": "Spatially Variable Gene Detection",
    "description": "Identify genes with spatial expression patterns",
    "inputs": [
      "spatial_data.h5ad"
    ],
    "outputs": [
      "svg_results.h5ad",
      "spatial_features.csv"
    ],
    "workflow": "nextflow svg_detection.nf",
    "parameters": {
      "n_top_genes": 2000,
      "spatial_key": "spatial"
    }
  },
  "label_transfer": {
    "name": "Cell Type Label Transfer",
    "description": "Transfer cell type labels from reference to spatial data",
    "inputs": [
      "spatial_data.h5ad",
      "reference_data.h5ad"
    ],
    "outputs": [
      "labeled_spatial.h5ad",
      "transfer_scores.csv"
    ],
    "workflow": "nextflow label_transfer.nf",
    "parameters": {
      "reference_key": "cell_type",
      "confidence_threshold": 0.5
    }
  }
}