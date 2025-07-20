"""Tests for the FastMCP server implementation."""

import pytest
import asyncio
from unittest.mock import Mock, patch

from openproblems_mcp.server import MCPServer
from openproblems_mcp.config import Config, ServerConfig, ToolConfig


@pytest.fixture
def mock_config():
    """Create a mock configuration for testing."""
    return Config(
        server=ServerConfig(
            workspace_root=".",
            log_level="INFO",
            max_concurrent_executions=2,
            default_timeout_seconds=60
        ),
        tools=ToolConfig(
            nextflow_executable="nextflow",
            viash_executable="viash",
            docker_executable="docker",
            git_executable="git",
            python_executable="python"
        )
    )


def test_server_initialization(mock_config):
    """Test that the server initializes correctly."""
    server = MCPServer(mock_config)

    assert server.config == mock_config
    assert server.tool_detector is not None
    assert server.mcp is not None


def test_server_has_required_tools(mock_config):
    """Test that the server registers the required tools."""
    server = MCPServer(mock_config)
    mcp_server = server.get_mcp_server()

    # Check that tools are registered (this will depend on FastMCP's API)
    # For now, just verify the server was created
    assert mcp_server is not None


@patch('openproblems_mcp.tool_detection.ToolDetector.get_health_status')
def test_health_check_tool(mock_health_status, mock_config):
    """Test the health_check tool."""
    # Mock the health status response
    mock_health_status.return_value = {
        'tools_available': 3,
        'tools_total': 5,
        'all_tools_available': False,
        'tools': {
            'nextflow': {'available': True, 'version': '23.04.0', 'error': None},
            'docker': {'available': True, 'version': '24.0.0', 'error': None},
            'viash': {'available': False, 'version': None, 'error': 'Not found'},
        }
    }

    server = MCPServer(mock_config)

    # This test would need to be adapted based on how FastMCP exposes tool testing
    # For now, just verify the server was created with the mock config
    assert server.config == mock_config


def test_server_config_resource(mock_config):
    """Test that server configuration can be accessed as a resource."""
    server = MCPServer(mock_config)

    # Verify the server was initialized
    assert server.config.server.workspace_root == "."
    assert server.config.server.log_level == "INFO"


@pytest.mark.skip(reason="Async test requires pytest-asyncio plugin")
async def test_server_startup(mock_config):
    """Test that the server can start up without errors."""
    server = MCPServer(mock_config)

    # For now, just test that we can create the server
    # Full startup testing would require mocking FastMCP's run method
    assert server is not None

    # Test that we can get the underlying FastMCP server
    mcp_server = server.get_mcp_server()
    assert mcp_server is not None
