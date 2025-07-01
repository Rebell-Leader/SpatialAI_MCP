{
  "overview": "Viash is a meta-framework for building reusable workflow modules",
  "component_structure": {
    "config_file": "YAML configuration defining component metadata",
    "script": "Core functionality implementation",
    "platforms": "Target platforms (docker, native, nextflow)"
  },
  "best_practices": {
    "modularity": "Keep components focused on single tasks",
    "documentation": "Provide clear descriptions and examples",
    "testing": "Include unit tests for all components",
    "versioning": "Use semantic versioning for component releases"
  },
  "common_commands": {
    "build": "viash build config.vsh.yaml",
    "run": "viash run config.vsh.yaml",
    "test": "viash test config.vsh.yaml",
    "ns_build": "viash ns build"
  },
  "code_examples": {
    "config_file": "name: 'spatial_qc'\ndescription: '...'\n...",
    "script_implementation": "import argparse\n..."
  }
}