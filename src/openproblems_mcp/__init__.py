"""OpenProblems Spatial Transcriptomics MCP Server.

A production-ready Model Context Protocol server for spatial transcriptomics
workflows, designed to work seamlessly with Continue.dev and VSCode.
"""

__version__ = "0.1.0"
__author__ = "OpenProblems MCP Contributors"
__email__ = "info@openproblems.bio"

from .server import MCPServer
from .config import Config, ServerConfig, ToolConfig, ConfigManager
from .exceptions import MCPServerError, ExecutionError, ValidationError
from .tool_detection import ToolDetector, ToolInfo

__all__ = [
    "MCPServer",
    "Config",
    "ServerConfig",
    "ToolConfig",
    "ConfigManager",
    "ToolDetector",
    "ToolInfo",
    "MCPServerError",
    "ExecutionError",
    "ValidationError",
]
