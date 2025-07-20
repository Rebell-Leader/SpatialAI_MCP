"""Core FastMCP server implementation for OpenProblems spatial transcriptomics workflows."""

import logging
from typing import Any, Dict, List, Optional
import json

from fastmcp import FastMCP

from .config import Config, ConfigManager, setup_logging
from .tool_detection import ToolDetector
from .exceptions import MCPServerError, DependencyError


logger = logging.getLogger(__name__)


class MCPServer:
    """Main FastMCP server for OpenProblems spatial transcriptomics workflows."""

    def __init__(self, config: Optional[Config] = None):
        self.config = config or ConfigManager().load_config()
        self.tool_detector = ToolDetector(self.config.tools)

        # Initialize FastMCP server
        self.mcp = FastMCP("OpenProblems Spatial MCP")
        self._setup_tools()
        self._setup_resources()

    def _setup_tools(self) -> None:
        """Set up MCP tools using FastMCP decorators."""

        @self.mcp.tool()
        def health_check() -> str:
            """Check the health status of the MCP server and its dependencies."""
            try:
                health_status = self.tool_detector.get_health_status()

                # Add server-specific health information
                health_status.update({
                    "server_status": "healthy",
                    "config_valid": True,
                    "workspace_accessible": True
                })

          ync def atus_text = "✅ Server Health Check\n\n"
                status_text += f"Tools Available: {health_status['tools_available']}/{health_status['tools_total']}\n"
                status_text += f"All Tools Available: {'Yes' if health_status['all_tools_available'] else 'No'}\n\n"

                status_text += "Tool Details:\n"
                for tool_name, tool_info in health_status['tools'].items():
                    status = "✅" if tool_info['available'] else "❌"
                    version = f" (v{tool_info['version']})" if tool_info['version'] else ""
                    status_text += f"{status} {tool_name}{version}\n"
                    if tool_info['error']:
                        status_text += f"   Error: {tool_info['error']}\n"

                return status_text

            except Exception as e:
                logger.error(f"Health check failed: {e}")
                return f"❌ Health check failed: {e}"

        @self.mcp.tool()
        def list_tools_status() -> str:
            """List the status of all detected bioinformatics tools."""
            try:
                tools = self.tool_detector.detect_all_tools()

                status_text = "🔧 Bioinformatics Tools Status\n\n"

                for tool_name, tool_info in tools.items():
                    status = "✅ Available" if tool_info.available else "❌ Not Available"
                    status_text += f"**{tool_name.title()}**: {status}\n"
                    status_text += f"  Executable: {tool_info.executable}\n"

                    if tool_info.version:
                        status_text += f"  Version: {tool_info.version}\n"

                    if tool_info.error_message:
                        status_text += f"  Error: {tool_info.error_message}\n"

                    status_text += "\n"

                # Add installation suggestions for missing tools
                missing_tools = self.tool_detector.get_missing_tools()
                if missing_tools:
                    status_text += "📋 Installation Suggestions:\n\n"
                    suggestions = self.tool_detector._get_installation_suggestions([t.name for t in missing_tools])
                    for suggestion in suggestions:
                        status_text += f"• {suggestion}\n"

                return status_text

            except Exception as e:
                logger.error(f"Tools status check failed: {e}")
                return f"❌ Tools status check failed: {e}"

        @self.mcp.tool()
        def get_server_info() -> str:
            """Get information about the MCP server configuration and status."""
            try:
                info_text = "ℹ️ OpenProblems MCP Server Information\n\n"
                info_text += f"**Version**: 0.1.0\n"
                info_text += f"**Workspace Root**: {self.config.server.workspace_root}\n"
                info_text += f"**Max Concurrent Executions**: {self.config.server.max_concurrent_executions}\n"
                info_text += f"**Default Timeout**: {self.config.server.default_timeout_seconds}s\n"
                info_text += f"**Log Level**: {self.config.server.log_level}\n"
                info_text += f"**Memory Limit per Execution**: {self.config.server.max_memory_per_execution_mb}MB\n\n"

                info_text += "**Allowed File Extensions**:\n"
                for ext in self.config.server.allowed_file_extensions:
                    info_text += f"  • {ext}\n"

                if self.config.server.blocked_paths:
                    info_text += "\n**Blocked Paths**:\n"
                    for path in self.config.server.blocked_paths:
                        info_text += f"  • {path}\n"

                return info_text

            except Exception as e:
                logger.error(f"Server info retrieval failed: {e}")
                return f"❌ Server info retrieval failed: {e}"

    def _setup_resources(self) -> None:
        """Set up MCP resources using FastMCP decorators."""

        @self.mcp.resource("config://server")
        def server_config() -> str:
            """Current server configuration."""
            config_dict = {
                "server": {
                    "workspace_root": self.config.server.workspace_root,
                    "max_concurrent_executions": self.config.server.max_concurrent_executions,
                    "default_timeout_seconds": self.config.server.default_timeout_seconds,
                    "log_level": self.config.server.log_level,
                    "max_memory_per_execution_mb": self.config.server.max_memory_per_execution_mb,
                    "allowed_file_extensions": self.config.server.allowed_file_extensions,
                    "blocked_paths": self.config.server.blocked_paths
                },
                "tools": {
                    "nextflow_executable": self.config.tools.nextflow_executable,
                    "viash_executable": self.config.tools.viash_executable,
                    "docker_executable": self.config.tools.docker_executable,
                    "git_executable": self.config.tools.git_executable,
                    "python_executable": self.config.tools.python_executable,
                    "default_nextflow_profile": self.config.tools.default_nextflow_profile,
                    "default_container_registry": self.config.tools.default_container_registry
                }
            }
            return json.dumps(config_dict, indent=2)

        @self.mcp.resource("status://tools")
        def tools_status() -> str:
            """Status of detected bioinformatics tools."""
            return json.dumps(self.tool_detector.get_health_status(), indent=2)

        @self.mcp.resource("status://health")
        def health_status() -> str:
            """Overall health status of the server."""
            health_status = self.tool_detector.get_health_status()
            health_status.update({
                "server_status": "healthy",
                "config_valid": True,
                "workspace_accessible": True
            })
            return json.dumps(health_status, indent=2)

    async def run(self) -> None:
        """Run the FastMCP server."""
        logger.info("Starting OpenProblems MCP Server with FastMCP")

        # Perform initial health check
        try:
            health_status = self.tool_detector.get_health_status()
            logger.info(f"Initial health check: {health_status['tools_available']}/{health_status['tools_total']} tools available")

            missing_tools = self.tool_detector.get_missing_tools()
            if missing_tools:
                logger.warning(f"Missing tools: {[t.name for t in missing_tools]}")
        except Exception as e:
            logger.error(f"Initial health check failed: {e}")

        # Run the FastMCP server
        await self.mcp.run()

    def get_mcp_server(self) -> FastMCP:
        """Get the underlying FastMCP server instance."""
        return self.mcp
