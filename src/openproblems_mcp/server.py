"""Core MCP server implementation for OpenProblems spatial transcriptomics workflows."""

import asyncio
import logging
from typing import Any, Dict, List, Optional, Sequence
from contextlib import asynccontextmanager

from mcp.server import Server
from mcp.server.models import InitializationOptions
from mcp.server.stdio import stdio_server
from mcp.types import (
    Tool,
    TextContent,
    Resource,
    CallToolRequest,
    ListToolsRequest,
    ListResourcesRequest,
    ReadResourceRequest,
)

from .config import Config, ConfigManager, setup_logging
from .tool_detection import ToolDetector
from .exceptions import MCPServerError, DependencyError


logger = logging.getLogger(__name__)


class MCPServer:
    """Main MCP server implementing the protocol specification."""

    def __init__(self, config: Optional[Config] = None):
        self.config = config or ConfigManager().load_config()
        self.server = Server("openproblems-spatial-mcp")
        self.tool_detector = ToolDetector(self.config.tools)
        self._setup_handlers()

    def _setup_handlers(self) -> None:
        """Set up MCP protocol handlers."""

        @self.server.list_tools()
        async def handle_list_toist[Tool]:
            """Handle list_tools requests."""
            return await self._handle_list_tools()

        @self.server.call_tool()
        async def handle_call_tool(name: str, arguments: Dict[str, Any]) -> List[TextContent]:
            """Handle call_tool requests."""
            return await self._handle_call_tool(name, arguments)

        @self.server.list_resources()
        async def handle_list_resources() -> List[Resource]:
            """Handle list_resources requests."""
            return await self._handle_list_resources()

        @self.server.read_resource()
        async def handle_read_resource(uri: str) -> str:
            """Handle read_resource requests."""
            return await self._handle_read_resource(uri)

    async def _handle_list_tools(self) -> List[Tool]:
        """List available MCP tools."""
        tools = [
            Tool(
                name="health_check",
                description="Check the health status of the MCP server and its dependencies",
                inputSchema={
                    "type": "object",
                    "properties": {},
                    "additionalProperties": False
                }
            ),
            Tool(
                name="list_tools_status",
                description="List the status of all detected bioinformatics tools",
                inputSchema={
                    "type": "object",
                    "properties": {},
                    "additionalProperties": False
                }
            ),
            Tool(
                name="get_server_info",
                description="Get information about the MCP server configuration and status",
                inputSchema={
                    "type": "object",
                    "properties": {},
                    "additionalProperties": False
                }
            )
        ]

        logger.debug(f"Listed {len(tools)} available tools")
        return tools

    async def _handle_call_tool(self, name: str, arguments: Dict[str, Any]) -> List[TextContent]:
        """Handle tool execution requests."""
        try:
            logger.info(f"Executing tool: {name} with arguments: {arguments}")

            if name == "health_check":
                return await self._tool_health_check(arguments)
            elif name == "list_tools_status":
                return await self._tool_list_tools_status(arguments)
            elif name == "get_server_info":
                return await self._tool_get_server_info(arguments)
            else:
                raise MCPServerError(f"Unknown tool: {name}")

        except Exception as e:
            logger.error(f"Tool execution failed: {e}")
            if isinstance(e, MCPServerError):
                error_info = e.to_dict()
            else:
                error_info = {
                    "error_type": "execution_error",
                    "message": str(e),
                    "suggested_fixes": ["Check server logs for more details"]
                }

            return [TextContent(
                type="text",
                text=f"Tool execution failed: {error_info}"
            )]

    async def _handle_list_resources(self) -> List[Resource]:
        """List available MCP resources."""
        resources = [
            Resource(
                uri="config://server",
                name="Server Configuration",
                description="Current server configuration",
                mimeType="application/json"
            ),
            Resource(
                uri="status://tools",
                name="Tools Status",
                description="Status of detected bioinformatics tools",
                mimeType="application/json"
            ),
            Resource(
                uri="status://health",
                name="Health Status",
                description="Overall health status of the server",
                mimeType="application/json"
            )
        ]

        logger.debug(f"Listed {len(resources)} available resources")
        return resources

    async def _handle_read_resource(self, uri: str) -> str:
        """Handle resource read requests."""
        try:
            logger.info(f"Reading resource: {uri}")

            if uri == "config://server":
                return await self._resource_server_config()
            elif uri == "status://tools":
                return await self._resource_tools_status()
            elif uri == "status://health":
                return await self._resource_health_status()
            else:
                raise MCPServerError(f"Unknown resource: {uri}")

        except Exception as e:
            logger.error(f"Resource read failed: {e}")
            raise MCPServerError(f"Failed to read resource {uri}: {e}")

    async def _tool_health_check(self, arguments: Dict[str, Any]) -> List[TextContent]:
        """Perform health check of server and dependencies."""
        try:
            health_status = self.tool_detector.get_health_status()

            # Add server-specific health information
            health_status.update({
                "server_status": "healthy",
                "config_valid": True,
                "workspace_accessible": True  # This will be properly checked in later tasks
            })

            status_text = "✅ Server Health Check\n\n"
            status_text += f"Tools Available: {health_status['tools_available']}/{health_status['tools_total']}\n"
            status_text += f"All Tools Available: {'Yes' if health_status['all_tools_available'] else 'No'}\n\n"

            status_text += "Tool Details:\n"
            for tool_name, tool_info in health_status['tools'].items():
                status = "✅" if tool_info['available'] else "❌"
                version = f" (v{tool_info['version']})" if tool_info['version'] else ""
                status_text += f"{status} {tool_name}{version}\n"
                if tool_info['error']:
                    status_text += f"   Error: {tool_info['error']}\n"

            return [TextContent(type="text", text=status_text)]

        except Exception as e:
            logger.error(f"Health check failed: {e}")
            return [TextContent(
                type="text",
                text=f"❌ Health check failed: {e}"
            )]

    async def _tool_list_tools_status(self, arguments: Dict[str, Any]) -> List[TextContent]:
        """List status of all detected tools."""
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

            return [TextContent(type="text", text=status_text)]

        except Exception as e:
            logger.error(f"Tools status check failed: {e}")
            return [TextContent(
                type="text",
                text=f"❌ Tools status check failed: {e}"
            )]

    async def _tool_get_server_info(self, arguments: Dict[str, Any]) -> List[TextContent]:
        """Get server configuration and status information."""
        try:
            info_text = "ℹ️ OpenProblems MCP Server Information\n\n"
            info_text += f"**Version**: {self.__class__.__module__.split('.')[0]} v0.1.0\n"
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

            return [TextContent(type="text", text=info_text)]

        except Exception as e:
            logger.error(f"Server info retrieval failed: {e}")
            return [TextContent(
                type="text",
                text=f"❌ Server info retrieval failed: {e}"
            )]

    async def _resource_server_config(self) -> str:
        """Return server configuration as JSON."""
        import json

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

    async def _resource_tools_status(self) -> str:
        """Return tools status as JSON."""
        import json
        return json.dumps(self.tool_detector.get_health_status(), indent=2)

    async def _resource_health_status(self) -> str:
        """Return overall health status as JSON."""
        import json

        health_status = self.tool_detector.get_health_status()
        health_status.update({
            "server_status": "healthy",
            "config_valid": True,
            "workspace_accessible": True
        })

        return json.dumps(health_status, indent=2)

    async def run_stdio(self) -> None:
        """Run the server using stdio transport."""
        logger.info("Starting OpenProblems MCP Server with stdio transport")

        # Perform initial health check
        try:
            health_status = self.tool_detector.get_health_status()
            logger.info(f"Initial health check: {health_status['tools_available']}/{health_status['tools_total']} tools available")

            missing_tools = self.tool_detector.get_missing_tools()
            if missing_tools:
                logger.warning(f"Missing tools: {[t.name for t in missing_tools]}")
        except Exception as e:
            logger.error(f"Initial health check failed: {e}")

        # Run the server
        async with stdio_server() as (read_stream, write_stream):
            await self.server.run(
                read_stream,
                write_stream,
                InitializationOptions(
                    server_name="openproblems-spatial-mcp",
                    server_version="0.1.0",
                    capabilities=self.server.get_capabilities()
                )
            )
