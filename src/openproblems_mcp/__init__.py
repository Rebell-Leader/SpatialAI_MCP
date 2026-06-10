"""OpenProblems Spatial Transcriptomics MCP Server.

A Model Context Protocol server providing domain-aware validation and analysis
tools for spatial transcriptomics workflows. Provider-neutral: usable from any
MCP-capable agent (Claude Code, Codex, Cursor, Copilot, Gemini CLI, ...).
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
