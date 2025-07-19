# OpenProblems Spatial Transcriptomics Local MCP Server

A local Model Context Protocol (MCP) server that enables AI agents like Continue.dev to seamlessly work with spatial transcriptomics workflows, Nextflow pipelines, Viash components, and OpenProblems benchmarks directly in your development environment.

## Quick Start

### Installation

```bash
# Install from PyPI (when published)
pip install openproblems-spatial-mcp

# Or install from source
git clone https://github.com/openproblems-bio/SpatialAI_MCP.git
cd SpatialAI_MCP
pip install -e .
```

### Continue.dev Configuration

Add this to your Continue.dev configuration (`~/.continue/config.json`):

```json
{
  "experimental": {
    "modelContextProtocolServers": [
      {
        "name": "openproblems-spatial",
        "transport": {
          "type": "stdio",
          "command": "openproblems-mcp-server",
          "cwd": "."
        }
      }
    ]
  }
}
```

### Usage

Once configured, the Continue.dev agent can:

- **File Operations**: Read, write, and manage your project files
- **Nextflow Execution**: Run and manage Nextflow pipelines
- **Viash Components**: Build and execute Viash components
- **Docker Management**: Build and manage Docker images
- **OpenProblems Integration**: Work with OpenProblems repositories and benchmarks
- **Spatial Data Processing**: Handle spatial transcriptomics data formats

## Architecture

The local MCP server runs as a background process in your development environment with direct access to:

- Your project files and directories
- Your local tool installations (Nextflow, Viash, Docker, Git)
- Your system resources and permissions
- Your development workflow and preferences

## Features

### Core Capabilities

- **Direct File Access**: No file system barriers or path translation
- **Local Tool Execution**: Uses your existing tool installations
- **Native Performance**: No containerization overhead
- **Simple Setup**: Single pip install with minimal configuration

### Bioinformatics Tools

- **Nextflow**: Execute pipelines with your local installation
- **Viash**: Build and run components using your local Viash
- **Docker**: Build images using your local Docker daemon
- **Git**: Repository operations with your local Git

### Spatial Transcriptomics

- **Data Validation**: Check SpatialData, zarr, and AnnData formats
- **Method Development**: Generate Viash components for spatial methods
- **Workflow Optimization**: Analyze and optimize spatial workflows
- **OpenProblems Integration**: Seamless benchmarking workflow

## Requirements

### System Requirements

- Python 3.8 or higher
- Operating System: Linux, macOS, or Windows

### Optional Tools

The server will detect and use these tools if available:

- **Nextflow**: For pipeline execution
- **Viash**: For component building and execution
- **Docker**: For container operations
- **Git**: For repository operations
- **Conda/Mamba**: For environment management

## Configuration

### Default Configuration

The server works out-of-the-box with sensible defaults. Optional configuration can be placed in:

- `~/.openproblems-mcp/config.yaml` (user-wide)
- `.openproblems-mcp.yaml` (project-specific)

### Example Configuration

```yaml
server:
  log_level: INFO
  max_concurrent_executions: 3
  default_timeout_seconds: 3600

tools:
  nextflow_executable: nextflow
  viash_executable: viash
  docker_executable: docker
  git_executable: git

workspace:
  allowed_extensions: [".py", ".R", ".nf", ".yaml", ".yml", ".json", ".txt", ".md"]
  # No blocked_paths in local mode - respects user's file permissions
```

## Development

### Setting up Development Environment

```bash
git clone https://github.com/openproblems-bio/SpatialAI_MCP.git
cd SpatialAI_MCP

# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest

# Format code
black src/

# Type checking
mypy src/
```

### Project Structure

```
src/openproblems_mcp/
├── __init__.py          # Package initialization
├── main.py              # Main entry point
├── server.py            # MCP server core
├── filesystem.py        # File system operations
├── execution.py         # Tool execution engine
├── state.py             # State management
├── config.py            # Configuration management
├── tools/               # MCP tool implementations
│   ├── __init__.py
│   ├── file_tools.py    # File operation tools
│   ├── nextflow_tools.py # Nextflow tools
│   ├── viash_tools.py   # Viash tools
│   ├── docker_tools.py  # Docker tools
│   └── spatial_tools.py # Spatial transcriptomics tools
└── templates/           # Component templates
    ├── viash_spatial.yaml
    └── nextflow_spatial.nf
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **GitHub Issues**: [Report bugs and request features](https://github.com/openproblems-bio/SpatialAI_MCP/issues)
- **Documentation**: [Full documentation](https://github.com/openproblems-bio/SpatialAI_MCP/docs)
- **OpenProblems**: [Learn about OpenProblems](https://openproblems.bio)

## Acknowledgments

This project is part of the OpenProblems initiative to standardize benchmarking in single-cell and spatial biology.
