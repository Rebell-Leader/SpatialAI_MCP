# OpenProblems Spatial Transcriptomics MCP Server

A production-ready Model Context Protocol (MCP) server for OpenProblems spatial transcriptomics workflows, built with [FastMCP 2.0](https://gofastmcp.com/) and designed to work seamlessly with Continue.dev in VSCode.

## 🎯 Purpose

This MCP server provides **bioinformatics-specific tools** that complement Continue.dev's built-in capabilities, focusing on:

- **Spatial transcriptomics domain expertise** (data validation, method development)
- **Bioinformatics tool execution** (Nextflow, Viash, Docker)
- **OpenProblems ecosystem integration** (benchmarking, method validation)
- **Workflow state management** (long-running pipeline tracking)

> **Note**: This server does NOT duplicate Continue.dev's built-in file operations, terminal commands, or git functionality. Instead, it provides specialized bioinformatics capabilities that work alongside Continue.dev's existing tools.

## 🚀 Quick Start

### Installation

```bash
pip install openproblems-spatial-mcp
```

### Continue.dev Configuration

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

### Basic Usage

1. **Check server health**:
   ```bash
   openproblems-mcp check
   ```

2. **Start the server** (usually automatic via Continue.dev):
   ```bash
   openproblems-mcp-server
   ```

3. **Initialize configuration**:
   ```bash
   openproblems-mcp init
   ```

## 🏗️ Architecture

The system operates as invisible local infrastructure in your development environment:

```mermaid
graph TB
    subgraph "Development Environment"
        subgraph "IDE"
            VSCode[VSCode IDE]
            Continue[Continue.dev Extension]
            VSCode --> Continue
        end

        subgraph "MCP Server (Local Process)"
            MCP[FastMCP Server]
            Tools[Bioinformatics Tools]
            State[Workflow State]
            MCP --> Tools
            MCP --> State
        end

        subgraph "Local Tools"
            Nextflow[Nextflow]
            Viash[Viash]
            Docker[Docker]
            Git[Git]
        end

        subgraph "User Files"
            Workspace[Project Files]
            Data[Spatial Data]
            Configs[Configurations]
        end
    end

    Continue -.->|MCP Protocol| MCP
    Tools --> Nextflow
    Tools --> Viash
    Tools --> Docker
    State -.-> Workspace
    MCP -.-> Data

    style MCP fill:#e1f5fe
    style Tools fill:#f3e5f5
    style Continue fill:#e8f5e8
```

## 🔧 Available Tools

### Core Infrastructure
- `health_check`: Check server and tool status
- `list_tools_status`: List bioinformatics tool availability
- `get_server_info`: Get server configuration details

### Bioinformatics Execution (Planned)
- `run_nextflow_workflow`: Execute Nextflow pipelines locally
- `build_viash_component`: Build Viash components
- `execute_viash_component`: Run Viash components
- `build_docker_image`: Build Docker images

### Spatial Transcriptomics (Planned)
- `validate_spatial_data`: Validate SpatialData/zarr/AnnData formats
- `analyze_spatial_metadata`: Extract spatial data metadata
- `setup_spatial_environment`: Generate conda/pip environments
- `create_spatial_component`: Generate Viash components for spatial methods

### OpenProblems Integration (Planned)
- `analyze_openproblems_repo`: Analyze repository structure
- `build_openproblems_method`: Build methods with OpenProblems framework
- `run_openproblems_benchmark`: Execute benchmarks
- `validate_openproblems_submission`: Validate submissions

### Workflow Management (Planned)
- `get_execution_status`: Check workflow status
- `cancel_execution`: Cancel running workflows
- `get_execution_history`: View execution history

## 📊 Available Resources

- `config://server`: Server configuration (JSON)
- `status://tools`: Tool detection status (JSON)
- `status://health`: OverON)health status (JSON)

## ⚙️ Configuration

### Default Configuration

The server works out-of-the-box with sensible defaults. Optional configuration:

- `~/.openproblems-mcp/config.yaml` (user-wide)
- `.openproblems-mcp.yaml` (project-specific)

### Example Configuration

```yaml
server:
  log_level: INFO
  max_concurrent_executions: 3
  default_timeout_seconds: 3600
  workspace_root: "."

tools:
  nextflow_executable: nextflow
  viash_executable: viash
  docker_executable: docker
  git_executable: git
  python_executable: python
```

### Environment Variables

- `OPENPROBLEMS_MCP_WORKSPACE_ROOT`: Workspace directory
- `OPENPROBLEMS_MCP_LOG_LEVEL`: Logging level (DEBUG, INFO, WARNING, ERROR)
- `OPENPROBLEMS_MCP_MAX_CONCURRENT`: Max concurrent executions
- `OPENPROBLEMS_MCP_TIMEOUT`: Default timeout in seconds
- `OPENPROBLEMS_MCP_NEXTFLOW_EXECUTABLE`: Nextflow executable path
- `OPENPROBLEMS_MCP_VIASH_EXECUTABLE`: Viash executable path
- `OPENPROBLEMS_MCP_DOCKER_EXECUTABLE`: Docker executable path

## 🔄 Integration with Continue.dev

This server **complements** Continue.dev's built-in tools:

### Continue.dev Built-ins (We DON'T duplicate)
- File operations (`read_file`, `create_new_file`)
- Search operations (`exact_search`, `file_glob_search`)
- Terminal commands (`run_terminal_command`)
- Git operations (`view_diff`)
- Directory operations (`view_subdirectory`)

### Our Specialized Tools (Unique value)
- Bioinformatics tool execution and management
- Spatial transcriptomics domain expertise
- OpenProblems ecosystem integration
- Long-running workflow state management
- Domain-specific error analysis and remediation

## 🧬 Use Cases

### Spatial Transcriptomics Method Development
```
Continue.dev Agent: "I need to develop a spatial clustering method"
↓
1. Agent uses Continue.dev's file tools to read existing code
2. Agent calls our validate_spatial_data to check test data
3. Agent calls our create_spatial_component to generate Viash component
4. Agent uses Continue.dev's terminal to run tests
5. Agent calls our build_openproblems_method to prepare submission
```

### Nextflow Pipeline Development
```
Continue.dev Agent: "Let's optimize this Nextflow pipeline"
↓
1. Agent uses Continue.dev's file tools to read pipeline code
2. Agent calls our run_nextflow_workflow to test execution
3. Agent calls our get_execution_status to monitor progress
4. Agent uses Continue.dev's diff tools to review changes
5. Agent calls our validate_openproblems_submission to check compliance
```

## 🛠️ Development

### Local Development Setup

```bash
git clone https://github.com/openproblems-bio/SpatialAI_MCP.git
cd SpatialAI_MCP

# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest

# Format code
black src/
ruff check src/
```

### Project Structure

```
src/openproblems_mcp/
├── __init__.py          # Package initialization
├── main.py              # Main entry point
├── server.py            # FastMCP server core
├── cli.py               # CLI commands
├── config.py            # Configuration management
├── tool_detection.py    # Tool detection logic
└── exceptions.py        # Error handling
```

## 📋 Requirements

### System Requirements
- Python 3.8 or higher
- Operating System: Linux, macOS, or Windows

### Optional Tools (Auto-detected)
- **Nextflow**: For pipeline execution
- **Viash**: For component building and execution
- **Docker**: For container operations
- **Git**: For repository operations
- **Python/Conda**: For environment management

## 🚦 Current Status

### ✅ Completed (Task 1)
- ✅ Clean Python package structure
- ✅ FastMCP-based server architecture
- ✅ Comprehensive logging and configuration
- ✅ Local tool detection and validation
- ✅ Pip-installable with CLI commands
- ✅ Health monitoring and status reporting

### 🚧 In Development
- 🚧 Bioinformatics tool execution (Task 2)
- 🚧 Spatial data validation (Task 3)
- 🚧 OpenProblems integration (Task 4)
- 🚧 Workflow state management (Task 5)

## 🤝 Why FastMCP?

FastMCP 2.0 provides:
- **Minimal boilerplate**: Focus on functionality, not protocol details
- **Production-ready**: Built for real-world deployment scenarios
- **High performance**: Optimized for MCP protocol efficiency
- **Pythonic**: Clean, intuitive API design
- **Comprehensive**: Full MCP ecosystem support

This implementation is significantly cleaner and more maintainable than raw MCP protocol implementations.

## 📚 Documentation

- **Installation Guide**: See Quick Start section above
- **API Reference**: Coming soon
- **Continue.dev Integration**: See Integration section above
- **Troubleshooting**: Use `openproblems-mcp check` for diagnostics

## 🤝 Contributing

We welcome contributions! Please:

1. Focus on bioinformatics-specific functionality
2. Avoid duplicating Continue.dev built-in tools
3. Follow the FastMCP patterns established in the codebase
4. Add tests for new functionality

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **OpenProblems Initiative**: For standardizing benchmarking in spatial biology
- **FastMCP**: For providing an excellent MCP framework
- **Continue.dev**: For creating a powerful AI coding assistant platform

## 📞 Support

- **GitHub Issues**: [Report bugs and request features](https://github.com/openproblems-bio/SpatialAI_MCP/issues)
- **Health Check**: Run `openproblems-mcp check` for diagnostics
- **OpenProblems**: [Learn about OpenProblems](https://openproblems.bio)
