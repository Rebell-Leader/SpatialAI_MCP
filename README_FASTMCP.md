# OpenProblems Spatial MCP Server (FastMCP)

A production-ready Model Context Protocol server for OpenProblems spatial transcriptomics workflows, built with [FastMCP 2.0](https://gofastmcp.com/).

## Features

- **Production-ready**: Built with FastMCP 2.0 for optimal performance and reliability
- **Minimal dependencies**: Clean, focused implementation without unnecessary bloat
- **Tool detection**: Automatic detection of Nextflow, Viash, Docker, Git, and Python
- **Health monitoring**: Built-in health checks and status reporting
- **Easy configuration**: YAML-based configuration with environment variable overrides
- **Continue.dev integration**: Designed to work seamlessly with Continue.dev in VSCode

## Quick Start

### Installation

```bash
pip install openproblems-spatial-mcp
```

### Basic Usage

1. **Start the server**:
   ```bash
   openproblems-mcp-server
   ```

2. **Check health status**:
   ```bash
   openproblems-mcp check
   ```

3. **Initialize configuration**:
   ```bash
   openproblems-mcp init
   ```

### Continue.dev Integration

Add to your Continue.dev configuration (`~/.continue/config.json`):

```json
{
  "mcpServers": {
    "openproblems-spatial": {
      "command": "openproblems-mcp-server",
      "args": [],
      "env": {
        "OPENPROBLEMS_MCP_LOG_LEVEL": "INFO"
      }
    }
  }
}
```

## Available Tools

- `health_check`: Check server and tool status
- `list_tools_status`: List bioinformatics tool availability
- `get_server_info`: Get server configuration details

## Available Resources

- `config://server`: Server configuration (JSON)
- `status://tools`: Tool detection status (JSON)
- `status://health`: Overall health status (JSON)

## Configuration

The server uses sensible defaults but can be customized via:

1. **Configuration file**: `~/.openproblems-mcp/config.yaml`
2. **Environment variables**: `OPENPROBLEMS_MCP_*`
3. **Command line**: `openproblems-mcp-server [config_path]`

### Environment Variables

- `OPENPROBLEMS_MCP_WORKSPACE_ROOT`: Workspace directory
- `OPENPROBLEMS_MCP_LOG_LEVEL`: Logging level (DEBUG, INFO, WARNING, ERROR)
- `OPENPROBLEMS_MCP_MAX_CONCURRENT`: Max concurrent executions
- `OPENPROBLEMS_MCP_TIMEOUT`: Default timeout in seconds
- `OPENPROBLEMS_MCP_NEXTFLOW_EXECUTABLE`: Nextflow executable path
- `OPENPROBLEMS_MCP_VIASH_EXECUTABLE`: Viash executable path
- `OPENPROBLEMS_MCP_DOCKER_EXECUTABLE`: Docker executable path

## Architecture

This server implements the core infrastructure for Task 1 of the production MCP server specification:

- ✅ Clean Python package structure
- ✅ FastMCP-based async architecture
- ✅ Comprehensive logging and configuration
- ✅ Local tool detection and validation
- ✅ Pip-installable with CLI commands

## Development

```bash
# Install in development mode
pip install -e .

# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Format code
black src/
ruff check src/
```

## Why FastMCP?

FastMCP 2.0 provides:
- **Minimal boilerplate**: Focus on functionality, not protocol details
- **Production-ready**: Built for real-world deployment scenarios
- **High performance**: Optimized for MCP protocol efficiency
- **Pythonic**: Clean, intuitive API design
- **Comprehensive**: Full MCP ecosystem support

This implementation is significantly cleaner and more maintainable than raw MCP protocol implementations.
