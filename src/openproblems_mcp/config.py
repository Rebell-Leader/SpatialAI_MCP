"""Configuration management for the OpenProblems MCP Server."""

import os
import yaml
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
from pathlib import Path

from .exceptions import ConfigurationError


@dataclass
class ToolConfig:
    """Configuration for external bioinformatics tools."""
    nextflow_executable: str = "nextflow"
    viash_executable: str = "viash"
    docker_executable: str = "docker"
    git_executable: str = "git"
    python_executable: str = "python"
    default_nextflow_profile: str = "standard"
    default_container_registry: str = "docker.io"

    def validate(self) -> None:
        """Validate tool configuration."""
        # This will be implemented in later tasks
        pass


@dataclass
class ServerConfig:
    """Main server configuration."""
    workspace_root: str = "."
    max_concurrent_executions: int = 3
    default_timeout_seconds: int = 3600
    log_retention_days: int = 7
    max_memory_per_execution_mb: int = 4096
    allowed_file_extensions: List[str] = field(default_factory=lambda: [
        ".py", ".R", ".nf", ".yaml", ".yml", ".json", ".txt", ".md",
        ".csv", ".tsv", ".h5", ".h5ad", ".zarr"
    ])
    blocked_paths: List[str] = field(default_factory=list)
    log_level: str = "INFO"

    def validate(self) -> None:
        """Validate server configuration."""
        if self.max_concurrent_executions < 1:
            raise ConfigurationError("max_concurrent_executions must be at least 1")

        if self.default_timeout_seconds < 1:
            raise ConfigurationError("default_timeout_seconds must be at least 1")

        if self.max_memory_per_execution_mb < 128:
            raise ConfigurationError("max_memory_per_execution_mb must be at least 128")

        if self.log_level not in ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]:
            raise ConfigurationError(f"Invalid log_level: {self.log_level}")


@dataclass
class Config:
    """Complete configuration for the MCP server."""
    server: ServerConfig = field(default_factory=ServerConfig)
    tools: ToolConfig = field(default_factory=ToolConfig)

    def validate(self) -> None:
        """Validate complete configuration."""
        self.server.validate()
        self.tools.validate()


class ConfigManager:
    """Manages configuration loading and validation."""

    DEFAULT_CONFIG_PATHS = [
        "~/.openproblems-mcp/config.yaml",
        ".openproblems-mcp.yaml",
        "config/default_config.yaml"
    ]

    def __init__(self, config_path: Optional[str] = None):
        self.config_path = config_path
        self._config: Optional[Config] = None

    def load_config(self) -> Config:
        """Load configuration from file or environment variables."""
        if self._config is not None:
            return self._config

        config_data = {}

        # Try to load from file
        config_file = self._find_config_file()
        if config_file:
            try:
                with open(config_file, 'r', encoding='utf-8') as f:
                    config_data = yaml.safe_load(f) or {}
            except Exception as e:
                raise ConfigurationError(f"Failed to load config from {config_file}: {e}")

        # Override with environment variables
        config_data = self._apply_env_overrides(config_data)

        # Create config objects
        server_config = ServerConfig(**config_data.get('server', {}))
        tools_config = ToolConfig(**config_data.get('tools', {}))

        self._config = Config(server=server_config, tools=tools_config)
        self._config.validate()

        return self._config

    def _find_config_file(self) -> Optional[str]:
        """Find the first available configuration file."""
        if self.config_path:
            if os.path.exists(self.config_path):
                return self.config_path
            else:
                raise ConfigurationError(f"Specified config file not found: {self.config_path}")

        for path in self.DEFAULT_CONFIG_PATHS:
            expanded_path = os.path.expanduser(path)
            if os.path.exists(expanded_path):
                return expanded_path

        return None

    def _apply_env_overrides(self, config_data: Dict[str, Any]) -> Dict[str, Any]:
        """Apply environment variable overrides to configuration."""
        # Server configuration overrides
        server_config = config_data.setdefault('server', {})

        if 'OPENPROBLEMS_MCP_WORKSPACE_ROOT' in os.environ:
            server_config['workspace_root'] = os.environ['OPENPROBLEMS_MCP_WORKSPACE_ROOT']

        if 'OPENPROBLEMS_MCP_MAX_CONCURRENT' in os.environ:
            try:
                server_config['max_concurrent_executions'] = int(os.environ['OPENPROBLEMS_MCP_MAX_CONCURRENT'])
            except ValueError:
                raise ConfigurationError("OPENPROBLEMS_MCP_MAX_CONCURRENT must be an integer")

        if 'OPENPROBLEMS_MCP_TIMEOUT' in os.environ:
            try:
                server_config['default_timeout_seconds'] = int(os.environ['OPENPROBLEMS_MCP_TIMEOUT'])
            except ValueError:
                raise ConfigurationError("OPENPROBLEMS_MCP_TIMEOUT must be an integer")

        if 'OPENPROBLEMS_MCP_LOG_LEVEL' in os.environ:
            server_config['log_level'] = os.environ['OPENPROBLEMS_MCP_LOG_LEVEL']

        # Tools configuration overrides
        tools_config = config_data.setdefault('tools', {})

        for tool in ['nextflow', 'viash', 'docker', 'git', 'python']:
            env_var = f'OPENPROBLEMS_MCP_{tool.upper()}_EXECUTABLE'
            if env_var in os.environ:
                tools_config[f'{tool}_executable'] = os.environ[env_var]

        return config_data

    def create_default_config_file(self, path: str) -> None:
        """Create a default configuration file."""
        default_config = {
            'server': {
                'log_level': 'INFO',
                'max_concurrent_executions': 3,
                'default_timeout_seconds': 3600,
                'workspace_root': '.',
            },
            'tools': {
                'nextflow_executable': 'nextflow',
                'viash_executable': 'viash',
                'docker_executable': 'docker',
                'git_executable': 'git',
                'python_executable': 'python',
            }
        }

        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w', encoding='utf-8') as f:
            yaml.dump(default_config, f, default_flow_style=False, indent=2)


def setup_logging(config: Config) -> None:
    """Set up logging based on configuration."""
    log_level = getattr(logging, config.server.log_level.upper())

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
        ]
    )

    # Set specific logger levels
    logging.getLogger('openproblems_mcp').setLevel(log_level)
    logging.getLogger('fastmcp').setLevel(logging.WARNING)  # Reduce FastMCP library noise
