#!/usr/bin/env python3
"""
OpenProblems Spatial Transcriptomics MCP Server

A Model Context Protocol server that provides AI agents with standardized access
to Nextflow pipelines, Viash components, and spatial transcriptomics workflows
within the OpenProblems project.
"""

import asyncio
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from mcp.server import Server
from mcp.server.models import InitializationOptions
from mcp.types import (
    GetPromptResult,
    Prompt,
    PromptArgument,
    PromptMessage,
    Resource,
    TextContent,
    Tool,
)
import mcp.server.stdio

# Import the new utils
from ..utils.ai_assistant import LLMAssistant
from ..utils.doc_indexer import EnhancedDocumentationIndexer

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize the MCP server
server = Server("OpenProblems-SpatialAI-MCP")

# Server configuration
SERVER_VERSION = "0.1.0"
SERVER_NAME = "OpenProblems Spatial Transcriptomics MCP"

# Initialize new utilities
llm_assistant = LLMAssistant()
doc_indexer = EnhancedDocumentationIndexer()


@server.list_resources()
async def handle_list_resources() -> List[Resource]:
    """List available resources for spatial transcriptomics workflows."""
    return [
        Resource(
            uri="server://status",
            name="Server Status",
            description="Current status and configuration of the MCP server",
            mimeType="application/json",
        ),
        Resource(
            uri="documentation://nextflow",
            name="Nextflow Documentation",
            description="Comprehensive documentation for Nextflow workflows and best practices",
            mimeType="application/json",
        ),
        Resource(
            uri="documentation://viash",
            name="Viash Documentation",
            description="Documentation for Viash components and configuration",
            mimeType="application/json",
        ),
        Resource(
            uri="documentation://docker",
            name="Docker Documentation",
            description="Docker best practices and optimization guidelines",
            mimeType="application/json",
        ),
        Resource(
            uri="templates://spatial-workflows",
            name="Spatial Transcriptomics Pipeline Templates",
            description="Curated Nextflow pipeline templates for spatial transcriptomics analysis",
            mimeType="application/json",
        ),
    ]


@server.read_resource()
async def handle_read_resource(uri: str) -> str:
    """Read and return resource content based on URI."""
    logger.info(f"Reading resource: {uri}")

    if uri == "server://status":
        status = {
            "server_name": SERVER_NAME,
            "version": SERVER_VERSION,
            "status": "running",
            "capabilities": {
                "nextflow_execution": True,
                "viash_components": True,
                "docker_builds": True,
                "automated_testing": True,
                "ai_code_analysis": True,
                "ai_log_debugging": True,
                "semantic_doc_search": True,
            },
            "supported_formats": ["h5ad", "json", "yaml", "nf", "vsh.yaml"],
            "documentation_available": True,
        }
        return json.dumps(status, indent=2)

    # For other resources, we can return a summary or a link to the Gradio UI
    # where the full functionality is available.
    return json.dumps({"message": f"Resource {uri} is available. Use the corresponding tool to interact with it."}, indent=2)


@server.list_tools()
async def handle_list_tools() -> List[Tool]:
    """List available tools for spatial transcriptomics workflows."""
    return [
        Tool(
            name="index_documentation",
            description="Scrapes and indexes documentation from official sources for semantic search.",
            inputSchema={}
        ),
        Tool(
            name="search_documentation",
            description="Performs a semantic search on the indexed documentation.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {"type": "string", "description": "The natural language query to search for."}
                },
                "required": ["query"]
            }
        ),
        Tool(
            name="analyze_code",
            description="Analyzes spatial transcriptomics code using an AI assistant with documentation context.",
            inputSchema={
                "type": "object",
                "properties": {
                    "code": {"type": "string", "description": "The Python code to analyze."},
                    "context": {"type": "string", "description": "A brief description of what the code should accomplish."}
                },
                "required": ["code", "context"]
            }
        ),
        Tool(
            name="debug_logs",
            description="Analyzes execution logs with an AI assistant to identify root causes and suggest fixes.",
            inputSchema={
                "type": "object",
                "properties": {
                    "logs": {"type": "string", "description": "The execution logs to debug."}
                },
                "required": ["logs"]
            }
        ),
        Tool(
            name="check_environment",
            description="Check if required tools and dependencies are installed",
            inputSchema={
                "type": "object",
                "properties": {
                    "tools": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of tools to check (nextflow, viash, docker, java, etc.)",
                        "default": ["nextflow", "viash", "docker", "java"]
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="validate_nextflow_config",
            description="Validate Nextflow configuration and pipeline syntax",
            inputSchema={
                "type": "object",
                "properties": {
                    "config_path": {
                        "type": "string",
                        "description": "Path to nextflow.config file"
                    },
                    "pipeline_path": {
                        "type": "string",
                        "description": "Path to main.nf or pipeline file"
                    }
                },
                "required": ["pipeline_path"]
            }
        ),
        Tool(
            name="analyze_nextflow_log",
            description="Analyze Nextflow execution logs for errors and troubleshooting",
            inputSchema={
                "type": "object",
                "properties": {
                    "log_file_path": {
                        "type": "string",
                        "description": "Path to the .nextflow.log file"
                    }
                },
                "required": ["log_file_path"]
            }
        ),
    ]


@server.call_tool()
async def handle_call_tool(name: str, arguments: Dict[str, Any]) -> List[TextContent]:
    """Handle tool execution requests."""
    logger.info(f"Executing tool: {name} with arguments: {arguments}")

    if name == "index_documentation":
        result = doc_indexer.index_documentation()
        return [TextContent(type="text", text=result)]

    elif name == "search_documentation":
        query = arguments.get("query", "")
        results = doc_indexer.search_documentation(query)
        return [TextContent(type="text", text=json.dumps(results, indent=2))]

    elif name == "analyze_code":
        code = arguments.get("code", "")
        context = arguments.get("context", "")
        docs_context = doc_indexer.search_documentation(f"spatial transcriptomics {context} {code[:150]}", limit=3)
        analysis = llm_assistant.analyze_with_context(code, "code", docs_context)
        return [TextContent(type="text", text=analysis)]

    elif name == "debug_logs":
        logs = arguments.get("logs", "")
        docs_context = doc_indexer.search_documentation(f"nextflow error debugging {logs[:200]}", limit=3)
        analysis = llm_assistant.analyze_with_context(logs, "logs", docs_context)
        return [TextContent(type="text", text=analysis)]

    elif name == "check_environment":
        return await _check_environment(arguments)

    elif name == "validate_nextflow_config":
        return await _validate_nextflow_config(arguments)

    elif name == "analyze_nextflow_log":
        return await _analyze_nextflow_log(arguments)

    else:
        raise ValueError(f"Unknown tool: {name}")


async def _check_environment(arguments: Dict[str, Any]) -> List[TextContent]:
    """Check if required tools and dependencies are installed."""
    tools = arguments.get("tools", ["nextflow", "viash", "docker", "java"])

    environment_status = {
        "overall_status": "ready",
        "tools": {},
        "recommendations": []
    }

    try:
        for tool in tools:
            tool_status = {"available": False, "version": None, "path": None}

            try:
                if tool == "nextflow":
                    result = subprocess.run(["nextflow", "-version"], capture_output=True, text=True, timeout=10)
                    if result.returncode == 0:
                        tool_status["available"] = True
                        tool_status["version"] = result.stdout.strip()
                        tool_status["path"] = subprocess.run(["which", "nextflow"], capture_output=True, text=True).stdout.strip()

                elif tool == "viash":
                    result = subprocess.run(["viash", "--version"], capture_output=True, text=True, timeout=10)
                    if result.returncode == 0:
                        tool_status["available"] = True
                        tool_status["version"] = result.stdout.strip()
                        tool_status["path"] = subprocess.run(["which", "viash"], capture_output=True, text=True).stdout.strip()

                elif tool == "docker":
                    result = subprocess.run(["docker", "--version"], capture_output=True, text=True, timeout=10)
                    if result.returncode == 0:
                        tool_status["available"] = True
                        tool_status["version"] = result.stdout.strip()
                        tool_status["path"] = subprocess.run(["which", "docker"], capture_output=True, text=True).stdout.strip()

                elif tool == "java":
                    result = subprocess.run(["java", "-version"], capture_output=True, text=True, timeout=10)
                    if result.returncode == 0:
                        tool_status["available"] = True
                        tool_status["version"] = result.stderr.strip()  # Java outputs version to stderr
                        tool_status["path"] = subprocess.run(["which", "java"], capture_output=True, text=True).stdout.strip()

                else:
                    # Generic tool check
                    result = subprocess.run([tool, "--version"], capture_output=True, text=True, timeout=10)
                    if result.returncode == 0:
                        tool_status["available"] = True
                        tool_status["version"] = result.stdout.strip()
                        tool_status["path"] = subprocess.run(["which", tool], capture_output=True, text=True).stdout.strip()

            except (subprocess.TimeoutExpired, FileNotFoundError):
                tool_status["available"] = False

            environment_status["tools"][tool] = tool_status

            # Add recommendations for missing tools
            if not tool_status["available"]:
                environment_status["overall_status"] = "incomplete"
                if tool == "nextflow":
                    environment_status["recommendations"].append("Install Nextflow: curl -s https://get.nextflow.io | bash")
                elif tool == "viash":
                    environment_status["recommendations"].append("Install Viash: curl -fsSL get.viash.io | bash")
                elif tool == "docker":
                    environment_status["recommendations"].append("Install Docker: https://docs.docker.com/get-docker/")
                elif tool == "java":
                    environment_status["recommendations"].append("Install Java: sudo apt install openjdk-17-jre-headless")

        return [TextContent(type="text", text=json.dumps(environment_status, indent=2))]

    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "status": "error",
                "error": f"Failed to check environment: {str(e)}"
            }, indent=2)
        )]


async def _validate_nextflow_config(arguments: Dict[str, Any]) -> List[TextContent]:
    """Validate Nextflow configuration and pipeline syntax."""
    pipeline_path = arguments["pipeline_path"]
    config_path = arguments.get("config_path")

    validation_results = {
        "pipeline_path": pipeline_path,
        "config_path": config_path,
        "issues": [],
        "warnings": [],
        "status": "valid"
    }

    try:
        # Check if pipeline file exists
        pipeline_file = Path(pipeline_path)
        if not pipeline_file.exists():
            validation_results["issues"].append(f"Pipeline file not found: {pipeline_path}")
            validation_results["status"] = "invalid"
            return [TextContent(type="text", text=json.dumps(validation_results, indent=2))]

        # Read and check pipeline content
        with open(pipeline_file, 'r') as f:
            pipeline_content = f.read()

        # Basic Nextflow syntax checks
        if 'nextflow.enable.dsl=2' not in pipeline_content and 'nextflow { dsl = 2 }' not in pipeline_content:
            validation_results["warnings"].append("DSL2 not explicitly enabled - recommend adding 'nextflow.enable.dsl=2'")

        if 'process ' not in pipeline_content and 'workflow ' not in pipeline_content:
            validation_results["issues"].append("No process or workflow blocks found in pipeline")
            validation_results["status"] = "invalid"

        # Check for common issues
        if 'publishDir' in pipeline_content and 'output:' not in pipeline_content:
            validation_results["warnings"].append("publishDir found but no output block - this may cause issues")

        # Check config file if provided
        if config_path:
            config_file = Path(config_path)
            if not config_file.exists():
                validation_results["warnings"].append(f"Config file not found: {config_path}")
            else:
                with open(config_file, 'r') as f:
                    config_content = f.read()

                # Basic config validation
                if 'process ' in config_content:
                    validation_results["warnings"].append("Config looks good - process configuration found")

        # Try to run nextflow validation if available
        try:
            result = subprocess.run(
                ["nextflow", "config", pipeline_path],
                capture_output=True, text=True, timeout=30
            )
            if result.returncode != 0:
                validation_results["issues"].append(f"Nextflow config validation failed: {result.stderr}")
                validation_results["status"] = "invalid"
        except (subprocess.TimeoutExpired, FileNotFoundError):
            validation_results["warnings"].append("Nextflow not available - performed basic syntax check only")

        return [TextContent(type="text", text=json.dumps(validation_results, indent=2))]

    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "status": "error",
                "error": f"Failed to validate Nextflow configuration: {str(e)}"
            }, indent=2)
        )]


async def _analyze_nextflow_log(arguments: Dict[str, Any]) -> List[TextContent]:
    """Analyze Nextflow execution logs for errors and troubleshooting."""
    log_file_path = arguments["log_file_path"]

    try:
        log_path = Path(log_file_path)
        if not log_path.exists():
            return [TextContent(
                type="text",
                text=json.dumps({
                    "status": "error",
                    "error": f"Log file not found: {log_file_path}"
                }, indent=2)
            )]

        # Read and analyze the log file
        with open(log_path, 'r') as f:
            log_content = f.read()

        analysis = {
            "log_file": str(log_path),
            "file_size": log_path.stat().st_size,
            "issues_found": [],
            "suggestions": [],
        }

        # Common error patterns and their solutions
        error_patterns = {
            "exit status 137": {
                "issue": "Out of memory (OOM) error",
                "suggestion": "Increase memory allocation for the process or implement dynamic resource allocation"
            },
            "error exit status (137)": {
                "issue": "Out of memory (OOM) error",
                "suggestion": "Increase memory allocation for the process or implement dynamic resource allocation"
            },
            "Command exit status:\n  137": {
                "issue": "Out of memory (OOM) error",
                "suggestion": "Increase memory allocation for the process or implement dynamic resource allocation"
            },
            "exit status 1": {
                "issue": "General execution error",
                "suggestion": "Check process logs for specific error details"
            },
            "command not found": {
                "issue": "Missing command or tool",
                "suggestion": "Ensure required tools are installed in the container or environment"
            },
            "No such file or directory": {
                "issue": "Missing input file",
                "suggestion": "Verify input file paths and ensure proper file staging"
            },
            "Permission denied": {
                "issue": "File permission error",
                "suggestion": "Check file permissions and container user settings"
            },
        }

        # Analyze log content for known patterns
        for pattern, info in error_patterns.items():
            if pattern.lower() in log_content.lower():
                analysis["issues_found"].append({
                    "pattern": pattern,
                    "issue": info["issue"],
                    "suggestion": info["suggestion"]
                })

        # Additional check for exit status 137 pattern variations
        if "137" in log_content:
            # Check if it's an OOM-related 137 exit code
            if any(phrase in log_content.lower() for phrase in ["exit status", "exit code", "terminated"]):
                # Only add if an OOM error is not already found
                is_oom_found = any("out of memory" in issue["issue"].lower() for issue in analysis["issues_found"])
                if not is_oom_found:
                    analysis["issues_found"].append({
                        "pattern": "exit status 137",
                        "issue": "Out of memory (OOM) error",
                        "suggestion": "Increase memory allocation for the process or implement dynamic resource allocation"
                    })

        # Extract execution statistics if available
        if "Execution completed" in log_content:
            analysis["execution_status"] = "completed"
        elif "Execution cancelled" in log_content:
            analysis["execution_status"] = "cancelled"
        elif "Execution failed" in log_content:
            analysis["execution_status"] = "failed"
        else:
            analysis["execution_status"] = "unknown"

        return [TextContent(
            type="text",
            text=json.dumps(analysis, indent=2)
        )]

    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "status": "error",
                "error": f"Failed to analyze log file: {str(e)}"
            }, indent=2)
        )]


async def main():
    """Main entry point for the MCP server."""
    logger.info(f"Starting {SERVER_NAME} v{SERVER_VERSION}")

    async with mcp.server.stdio.stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            InitializationOptions(
                server_name=SERVER_NAME,
                server_version=SERVER_VERSION,
                capabilities={
                    "resources": {},
                    "tools": {},
                    "prompts": {},
                    "logging": {}
                },
            ),
        )


if __name__ == "__main__":
    asyncio.run(main())