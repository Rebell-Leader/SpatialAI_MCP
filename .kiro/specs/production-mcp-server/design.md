# Design Document

## Overview

This design document outlines the architecture for a production-ready, local OpenProblems Spatial Transcriptomics MCP Server that operates as invisible infrastructure supporting the three-part development environment: IDE (VSCode) + Continue.dev Extension + Local MCP Server.

The system is designed to be a robust, efficient, and maintainable local service that provides computational biology tools through the Model Context Protocol, enabling AI agents to seamlessly execute real bioinformatics workflows directly in the user's development environment.

## Architecture

### High-Level System Architecture

```mermaid
graph TB
    subgraph "User's Development Environment"
        subgraph "IDE"
            IDE[VSCode IDE]
            CONT[Continue.dev Extension]
            IDE --> CONT
        end

        subgraph "Local MCP Server Process"
            MCP[MCP Server Core]
            EXEC[Execution Engine]
            STATE[State Manager]
            FS[File System Handler]

            MCP --> EXEC
            MCP --> STATE
            MCP --> FS
        end

        subgraph "User's Local Tools"
            NF[Nextflow]
            VIASH[Viash]
            DOCKER[Docker]
            GIT[Git]
            PYTHON[Python/Conda]
        end

        subgraph "User's Project Files"
            WORKSPACE[Project Directory]
            LOGS[Local Logs]
            STATE_DB[Local State]
            CONFIG[User Config]
        end
    end

    CONT -.->|stdio/MCP Protocol| MCP
    EXEC --> NF
    EXEC --> VIASH
    EXEC --> DOCKER
    EXEC --> GIT
    EXEC --> PYTHON

    FS -.->|Direct Access| WORKSPACE
    STATE -.->|Direct Access| STATE_DB
    EXEC -.->|Direct Access| LOGS
    MCP -.->|Direct Access| CONFIG

    style MCP fill:#e1f5fe
    style EXEC fill:#f3e5f5
    style STATE fill:#e8f5e8
    style FS fill:#fff3e0
```

### Component Architecture

#### 1. MCP Server Core
- **Protocol Handler**: Implements MCP specification for tool and resource management
- **Request Router**: Routes incoming tool calls to appropriate execution handlers
- **Response Formatter**: Structures responses for Continue.dev agent consumption
- **Health Monitor**: Tracks system health and dependency status

#### 2. Execution Engine
- **Tool Executors**: Specialized handlers for each bioinformatics tool
- **Process Manager**: Manages concurrent executions with proper resource allocation
- **Timeout Handler**: Implements graceful timeouts and cleanup procedures
- **Result Collector**: Aggregates execution results and error information

#### 3. State Manager
- **Workflow Tracker**: Maintains state of long-running workflows
- **Progress Monitor**: Tracks execution progress and provides updates
- **History Manager**: Maintains execution history and logs
- **Configuration Store**: Manages persistent configuration and settings

#### 4. File System Handler
- **Path Resolver**: Handles secure path resolution and validation
- **Stream Processor**: Implements efficient file operations for large datasets
- **Permission Manager**: Enforces file system security and access controls
- **Workspace Manager**: Manages project workspace organization

## Components and Interfaces

### MCP Server Core Interface

```python
class MCPServer:
    """Main MCP server implementing the protocol specification"""

    async def handle_list_tools() -> List[Tool]
    async def handle_call_tool(name: str, arguments: Dict) -> List[TextContent]
    async def handle_list_resources() -> List[Resource]
    async def handle_read_resource(uri: str) -> str
    async def health_check() -> HealthStatus
```

### Execution Engine Interface

```python
class ExecutionEngine:
    """Manages execution of bioinformatics tools"""

    async def execute_nextflow(workflow_path: str, params: Dict) -> ExecutionResult
    async def build_viash_component(component_path: str) -> BuildResult
    async def build_docker_image(dockerfile_path: str, tag: str) -> BuildResult
    async def execute_viash_component(component: str, params: Dict) -> ExecutionResult
    async def clone_repository(repo_url: str, target_dir: str) -> CloneResult
    async def get_execution_status(execution_id: str) -> ExecutionStatus
    async def cancel_execution(execution_id: str) -> CancelResult
```

### State Manager Interface

```python
class StateManager:
    """Manages workflow state and execution history"""

    async def create_execution(tool: str, params: Dict) -> str  # Returns execution_id
    async def update_execution_status(execution_id: str, status: ExecutionStatus)
    async def get_execution_history(limit: int = 50) -> List[ExecutionRecord]
    async def cleanup_old_executions(retention_days: int = 7)
    async def get_active_executions() -> List[ExecutionRecord]
```

### File System Handler Interface

```python
class FileSystemHandler:
    """Handles secure file system operations"""

    async def read_file(path: str, encoding: str = 'utf-8') -> str
    async def write_file(path: str, content: str, encoding: str = 'utf-8')
    async def list_directory(path: str, recursive: bool = False) -> List[FileInfo]
    async def create_directory(path: str, parents: bool = True)
    async def delete_path(path: str, recursive: bool = False)
    async def get_file_info(path: str) -> FileInfo
    async def stream_large_file(path: str) -> AsyncIterator[bytes]
```

## Data Models

### Execution Models

```python
@dataclass
class ExecutionResult:
    execution_id: str
    status: ExecutionStatus
    return_code: int
    stdout: str
    stderr: str
    execution_time: float
    resource_usage: ResourceUsage
    artifacts: List[str]
    error_analysis: Optional[ErrorAnalysis]

@dataclass
class ExecutionStatus:
    state: ExecutionState  # PENDING, RUNNING, COMPLETED, FAILED, CANCELLED
    progress: float  # 0.0 to 1.0
    current_step: str
    estimated_completion: Optional[datetime]

@dataclass
class ResourceUsage:
    max_memory_mb: float
    cpu_time_seconds: float
    disk_io_mb: float
    network_io_mb: float

@dataclass
class ErrorAnalysis:
    error_type: str
    root_cause: str
    suggested_fixes: List[str]
    relevant_logs: List[str]
    documentation_links: List[str]
```

### File System Models

```python
@dataclass
class FileInfo:
    path: str
    name: str
    size: int
    modified_time: datetime
    is_directory: bool
    permissions: str
    mime_type: Optional[str]

@dataclass
class WorkspaceInfo:
    root_path: str
    total_size: int
    file_count: int
    directory_count: int
    last_modified: datetime
```

### Configuration Models

```python
@dataclass
class ServerConfig:
    workspace_root: str
    max_concurrent_executions: int
    default_timeout_seconds: int
    log_retention_days: int
    max_memory_per_execution_mb: int
    allowed_file_extensions: List[str]
    blocked_paths: List[str]

@dataclass
class ToolConfig:
    nextflow_executable: str
    viash_executable: str
    docker_executable: str
    git_executable: str
    default_nextflow_profile: str
    default_container_registry: str
```

## Error Handling

### Error Classification System

```python
class ErrorType(Enum):
    VALIDATION_ERROR = "validation_error"
    EXECUTION_ERROR = "execution_error"
    RESOURCE_ERROR = "resource_error"
    PERMISSION_ERROR = "permission_error"
    TIMEOUT_ERROR = "timeout_error"
    DEPENDENCY_ERROR = "dependency_error"
    CONFIGURATION_ERROR = "configuration_error"

@dataclass
class StructuredError:
    error_type: ErrorType
    error_code: str
    message: str
    context: Dict[str, Any]
    suggested_fixes: List[str]
    documentation_url: Optional[str]
    retry_possible: bool
```

### Error Recovery Strategies

1. **Automatic Retry**: For transient failures (network issues, temporary resource constraints)
2. **Graceful Degradation**: Fallback to alternative approaches when primary methods fail
3. **Resource Cleanup**: Automatic cleanup of partial executions and temporary files
4. **State Recovery**: Ability to resume interrupted workflows where possible

### Common Error Patterns and Solutions

```python
ERROR_PATTERNS = {
    "nextflow_oom": {
        "pattern": r"exit status 137|OutOfMemoryError",
        "analysis": "Nextflow process ran out of memory",
        "fixes": [
            "Increase memory allocation in nextflow.config",
            "Use process.memory directive for specific processes",
            "Consider data chunking for large datasets"
        ]
    },
    "viash_build_fail": {
        "pattern": r"viash build.*failed|docker build.*failed",
        "analysis": "Viash component build failed",
        "fixes": [
            "Check Dockerfile syntax and dependencies",
            "Verify base image availability",
            "Check network connectivity for package downloads"
        ]
    },
    "docker_permission": {
        "pattern": r"permission denied.*docker|Cannot connect to the Docker daemon",
        "analysis": "Docker permission or daemon issues",
        "fixes": [
            "Ensure Docker daemon is running",
            "Add user to docker group",
            "Check Docker socket permissions"
        ]
    }
}
```

## Testing Strategy

### Unit Testing
- **Component Isolation**: Test each component independently with mocked dependencies
- **Error Simulation**: Test error handling paths with simulated failures
- **Resource Management**: Test resource allocation and cleanup procedures
- **Configuration Validation**: Test configuration loading and validation

### Integration Testing
- **Tool Integration**: Test actual execution of Nextflow, Viash, and Docker
- **File System Operations**: Test file operations with various scenarios
- **State Management**: Test workflow state tracking and persistence
- **MCP Protocol**: Test MCP protocol compliance and message handling

### Performance Testing
- **Concurrent Execution**: Test multiple simultaneous tool executions
- **Large File Handling**: Test operations with large spatial transcriptomics datasets
- **Memory Management**: Test memory usage under various workloads
- **Timeout Handling**: Test graceful handling of long-running operations

### End-to-End Testing
- **Complete Workflows**: Test full spatial transcriptomics method development workflows
- **Error Recovery**: Test recovery from various failure scenarios
- **State Persistence**: Test state preservation across server restarts
- **Continue.dev Integration**: Test integration with actual Continue.dev agent

## Deployment Architecture

### Local Installation Structure

```bash
# Simple pip installation
pip install openproblems-spatial-mcp

# Package structure
openproblems-spatial-mcp/
├── src/
│   └── openproblems_mcp/
│       ├── __init__.py
│       ├── main.py              # Entry point
│       ├── server.py            # MCP server core
│       ├── execution.py         # Tool execution
│       ├── filesystem.py        # File operations
│       ├── state.py             # State management
│       └── config.py            # Configuration
├── config/
│   └── default_config.yaml      # Default configuration
└── setup.py                     # Package setup
```

### Continue.dev Configuration

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

### Local Configuration

```yaml
# ~/.openproblems-mcp/config.yaml
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
  # Uses current working directory by default
  allowed_extensions: [".py", ".R", ".nf", ".yaml", ".yml", ".json", ".txt", ".md"]
  blocked_paths: []  # No restrictions in local mode
```

### Security Considerations

1. **Local Process Security**:
   - Run with user's permissions (no privilege escalation)
   - Respect user's file system permissions
   - No network services or exposed ports

2. **File System Security**:
   - Operate within user's project directory by default
   - Respect user's file permissions and ownership
   - Optional path validation for enhanced security

3. **Resource Management**:
   - Configurable memory limits per execution
   - Respect user's system resource limits
   - Graceful handling of resource constraints

4. **Tool Execution Security**:
   - Execute tools with user's permissions
   - Use user's existing tool installations
   - No privilege escalation or system modification

## Performance Optimization

### Resource Management
- **Memory Pooling**: Reuse memory allocations for similar operations
- **Process Pooling**: Maintain pool of worker processes for tool execution
- **Caching**: Cache frequently accessed files and computation results
- **Streaming**: Use streaming for large file operations

### Concurrency Optimization
- **Async Operations**: Use asyncio for I/O-bound operations
- **Thread Pools**: Use thread pools for CPU-bound operations
- **Queue Management**: Implement fair queuing for concurrent executions
- **Resource Scheduling**: Intelligent scheduling based on resource requirements

### Monitoring and Metrics
- **Execution Metrics**: Track execution times, success rates, resource usage
- **System Metrics**: Monitor CPU, memory, disk usage
- **Error Metrics**: Track error rates and types
- **Performance Alerts**: Alert on performance degradation

This design provides a robust foundation for a production-ready MCP server that operates as invisible infrastructure while providing powerful bioinformatics capabilities to AI agents through the Continue.dev extension.
