{
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