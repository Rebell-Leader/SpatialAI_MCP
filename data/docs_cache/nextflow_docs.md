{
  "overview": "Nextflow is a workflow framework for bioinformatics pipelines",
  "best_practices": {
    "dsl_version": "Use DSL2 for all new workflows",
    "resource_management": "Specify memory and CPU requirements for each process",
    "error_handling": "Implement retry strategies and error handling",
    "containerization": "Use Docker/Singularity containers for reproducibility"
  },
  "common_patterns": {
    "input_channels": "Use Channel.fromPath() for file inputs",
    "output_publishing": "Use publishDir directive for results",
    "conditional_execution": "Use when clause for conditional processes"
  },
  "troubleshooting": {
    "oom_errors": "Increase memory allocation or implement dynamic resource allocation",
    "missing_files": "Check file paths and ensure proper input staging",
    "container_issues": "Verify container availability and permissions"
  },
  "code_examples": {
    "basic_pipeline": "#!/usr/bin/env nextflow\nnextflow.enable.dsl=2\n...",
    "process_definition": "process SPATIAL_ANALYSIS { ... }"
  }
}