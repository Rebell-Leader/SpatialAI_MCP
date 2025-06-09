# Setup Guide - OpenProblems Spatial Transcriptomics MCP Server

This guide will help you set up and run the OpenProblems Spatial Transcriptomics MCP Server.

## Prerequisites

### System Requirements

- **Python**: 3.8 or higher
- **Operating System**: Linux, macOS, or Windows (with WSL2 recommended)
- **Memory**: Minimum 4GB RAM (8GB+ recommended for processing large datasets)
- **Storage**: 10GB+ free space for data and temporary files

### Required Tools

The MCP server integrates with these bioinformatics tools:

- **[Nextflow](https://www.nextflow.io/)**: Workflow orchestration
- **[Viash](https://viash.io/)**: Component framework
- **[Docker](https://www.docker.com/)**: Containerization
- **Java**: 11 or higher (required for Nextflow)

## Installation

### Option 1: Local Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/openproblems-bio/SpatialAI_MCP.git
   cd SpatialAI_MCP
   ```

2. **Create a Python virtual environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install the package**:
   ```bash
   pip install -e .
   ```

4. **Install external tools**:

   **Nextflow**:
   ```bash
   curl -s https://get.nextflow.io | bash
   sudo mv nextflow /usr/local/bin/
   ```

   **Viash**:
   ```bash
   curl -fsSL get.viash.io | bash -s -- --bin /usr/local/bin
   ```

   **Docker**: Follow the [official Docker installation guide](https://docs.docker.com/get-docker/)

### Option 2: Docker Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/openproblems-bio/SpatialAI_MCP.git
   cd SpatialAI_MCP
   ```

2. **Build the Docker image**:
   ```bash
   docker build -f docker/Dockerfile -t openproblems-spatial-mcp .
   ```

3. **Run with Docker Compose**:
   ```bash
   cd docker
   docker-compose up -d
   ```

### Option 3: Development Setup

For contributors and developers:

1. **Clone and install in development mode**:
   ```bash
   git clone https://github.com/openproblems-bio/SpatialAI_MCP.git
   cd SpatialAI_MCP
   pip install -e ".[dev]"
   ```

2. **Install pre-commit hooks**:
   ```bash
   pre-commit install
   ```

3. **Run tests**:
   ```bash
   pytest tests/
   ```

## Configuration

### Basic Configuration

The server uses `config/server_config.yaml` for configuration. Key settings:

```yaml
server:
  name: "OpenProblems-SpatialAI-MCP"
  transport:
    primary: "stdio"
    http_port: 8000

paths:
  data_dir: "./data"
  work_dir: "./work"
  logs_dir: "./logs"

tools:
  nextflow:
    default_profile: "docker"
  viash:
    default_engine: "docker"
```

### Environment Variables

You can override configuration with environment variables:

```bash
export MCP_SERVER_NAME="Custom-MCP-Server"
export MCP_DATA_DIR="/custom/data/path"
export MCP_LOG_LEVEL="DEBUG"
```

### Directory Structure

Create the required directories:

```bash
mkdir -p data work logs cache
chmod 755 data work logs cache
```

## Running the Server

### Method 1: Direct Python Execution

```bash
# Start the server
python -m mcp_server.main

# Or use the installed command
openproblems-mcp
```

### Method 2: Docker

```bash
# Run the container
docker run -it --rm \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/work:/app/work \
  -v $(pwd)/logs:/app/logs \
  -v /var/run/docker.sock:/var/run/docker.sock \
  openproblems-spatial-mcp
```

### Method 3: Docker Compose

```bash
cd docker
docker-compose up
```

## Testing the Installation

### Run the Test Suite

```bash
pytest tests/ -v
```

### Use the Example Client

```bash
python examples/simple_client.py
```

### Manual Testing

1. **Start the server** (in one terminal):
   ```bash
   python -m mcp_server.main
   ```

2. **Test with MCP client** (in another terminal):
   ```python
   import asyncio
   from mcp import ClientSession, StdioServerParameters
   from mcp.client.stdio import stdio_client

   async def test_connection():
       server_params = StdioServerParameters(
           command="python",
           args=["-m", "mcp_server.main"],
       )

       async with stdio_client(server_params) as (read, write):
           async with ClientSession(read, write) as session:
               await session.initialize()

               # Test echo
               result = await session.call_tool("echo_test", {"message": "Hello!"})
               print(f"Echo result: {result}")

               # List resources
               resources = await session.list_resources()
               print(f"Available resources: {len(resources)}")

   asyncio.run(test_connection())
   ```

## Troubleshooting

### Common Issues

1. **Import errors**:
   - Ensure the package is installed: `pip install -e .`
   - Check Python path: `python -c "import mcp_server; print('OK')"`

2. **Tool not found errors**:
   - Install missing tools (Nextflow, Viash, Docker)
   - Check PATH: `which nextflow`, `which viash`, `which docker`

3. **Permission errors**:
   - Ensure Docker daemon is running: `docker version`
   - Check directory permissions: `ls -la data/ work/ logs/`

4. **Port conflicts** (HTTP transport):
   - Change port in config: `transport.http_port: 8001`
   - Check port usage: `netstat -tulpn | grep 8000`

### Debug Mode

Enable debug logging:

```bash
export MCP_LOG_LEVEL=DEBUG
python -m mcp_server.main
```

### Log Files

Check server logs:

```bash
tail -f logs/mcp_server.log
```

### Health Check

Test server health:

```bash
# For Docker containers
docker exec openproblems-spatial-mcp python -c "import mcp; print('MCP SDK available')"

# For local installation
python -c "import mcp_server.main; print('Server module available')"
```

## Next Steps

1. **Read the [API Documentation](API.md)** to understand available tools and resources
2. **Explore [Examples](../examples/)** to see practical usage patterns
3. **Check the [Integration Guide](INTEGRATION.md)** for AI agent setup
4. **Review [Best Practices](BEST_PRACTICES.md)** for optimal usage

## Support

- **Issues**: [GitHub Issues](https://github.com/openproblems-bio/SpatialAI_MCP/issues)
- **Documentation**: [Project Docs](https://github.com/openproblems-bio/SpatialAI_MCP/docs)
- **Community**: [OpenProblems Discussions](https://github.com/openproblems-bio/openproblems/discussions)

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) for development guidelines and contribution instructions.
