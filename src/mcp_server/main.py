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
from .documentation_generator_simple import DocumentationGenerator

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

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize the MCP server
server = Server("OpenProblems-SpatialAI-MCP")

# Server configuration
SERVER_VERSION = "0.1.0"
SERVER_NAME = "OpenProblems Spatial Transcriptomics MCP"

# Initialize documentation generator
doc_generator = DocumentationGenerator()


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
                "log_analysis": True,
            },
            "supported_formats": ["h5ad", "json", "yaml", "nf", "vsh.yaml"],
            "documentation_available": True,
        }
        return json.dumps(status, indent=2)

    elif uri == "documentation://nextflow":
        # Try to load cached documentation first
        cached_docs = await doc_generator.load_cached_documentation()
        if "nextflow" in cached_docs:
            return cached_docs["nextflow"]
        else:
            # Fallback to basic documentation
            nextflow_docs = {
                "overview": "Nextflow is a workflow framework for bioinformatics pipelines",
                "status": "Real documentation not yet cached - run 'python -m mcp_server.documentation_scraper' to download",
                "best_practices": {
                    "dsl_version": "Use DSL2 for all new workflows",
                    "resource_management": "Specify memory and CPU requirements for each process",
                    "error_handling": "Implement retry strategies and error handling",
                    "containerization": "Use Docker/Singularity containers for reproducibility",
                },
                "common_patterns": {
                    "input_channels": "Use Channel.fromPath() for file inputs",
                    "output_publishing": "Use publishDir directive for results",
                    "conditional_execution": "Use when clause for conditional processes",
                },
                "troubleshooting": {
                    "oom_errors": "Increase memory allocation or implement dynamic resource allocation",
                    "missing_files": "Check file paths and ensure proper input staging",
                    "container_issues": "Verify container availability and permissions",
                },
            }
            return json.dumps(nextflow_docs, indent=2)

    elif uri == "documentation://viash":
        # Try to load cached documentation first
        cached_docs = await doc_generator.load_cached_documentation()
        if "viash" in cached_docs:
            return cached_docs["viash"]
        else:
            # Fallback to basic documentation
            viash_docs = {
                "overview": "Viash is a meta-framework for building reusable workflow modules",
                "status": "Real documentation not yet cached - run 'python -m mcp_server.documentation_scraper' to download",
                "component_structure": {
                    "config_file": "YAML configuration defining component metadata",
                    "script": "Core functionality implementation",
                    "platforms": "Target platforms (docker, native, nextflow)",
                },
                "best_practices": {
                    "modularity": "Keep components focused on single tasks",
                    "documentation": "Provide clear descriptions and examples",
                    "testing": "Include unit tests for all components",
                    "versioning": "Use semantic versioning for component releases",
                },
                "common_commands": {
                    "build": "viash build config.vsh.yaml",
                    "run": "viash run config.vsh.yaml",
                    "test": "viash test config.vsh.yaml",
                    "ns_build": "viash ns build",
                },
            }
            return json.dumps(viash_docs, indent=2)

    elif uri == "documentation://docker":
        # Try to load cached documentation first
        cached_docs = await doc_generator.load_cached_documentation()
        if "docker" in cached_docs:
            return cached_docs["docker"]
        else:
            # Return generated Docker best practices
            return await doc_generator._generate_docker_docs()

    elif uri == "templates://spatial-workflows":
        # Try to load cached documentation first
        cached_docs = await doc_generator.load_cached_documentation()
        if "spatial_templates" in cached_docs:
            return cached_docs["spatial_templates"]
        else:
            # Return generated spatial workflow templates
            return await doc_generator._generate_spatial_templates()

    else:
        raise ValueError(f"Unknown resource URI: {uri}")


@server.list_tools()
async def handle_list_tools() -> List[Tool]:
    """List available tools for spatial transcriptomics workflows."""
    return [
        Tool(
            name="echo_test",
            description="Simple echo test to verify MCP communication",
            inputSchema={
                "type": "object",
                "properties": {
                    "message": {
                        "type": "string",
                        "description": "Message to echo back"
                    }
                },
                "required": ["message"]
            }
        ),
        Tool(
            name="list_available_tools",
            description="List all available MCP tools and their descriptions",
            inputSchema={
                "type": "object",
                "properties": {},
            }
        ),
        Tool(
            name="run_nextflow_workflow",
            description="Execute a Nextflow pipeline from OpenProblems repositories",
            inputSchema={
                "type": "object",
                "properties": {
                    "workflow_name": {
                        "type": "string",
                        "description": "Name of the Nextflow workflow (e.g., main.nf)"
                    },
                    "github_repo_url": {
                        "type": "string",
                        "description": "GitHub URL of the repository containing the workflow"
                    },
                    "profile": {
                        "type": "string",
                        "description": "Nextflow profile to use (e.g., docker, test)",
                        "default": "docker"
                    },
                    "params": {
                        "type": "object",
                        "description": "Key-value pairs for pipeline parameters",
                        "default": {}
                    },
                    "config_file": {
                        "type": "string",
                        "description": "Path to custom Nextflow configuration file"
                    }
                },
                "required": ["workflow_name", "github_repo_url"]
            }
        ),
        Tool(
            name="run_viash_component",
            description="Execute a Viash component with specified parameters",
            inputSchema={
                "type": "object",
                "properties": {
                    "component_name": {
                        "type": "string",
                        "description": "Name of the Viash component"
                    },
                    "component_config_path": {
                        "type": "string",
                        "description": "Path to the Viash config file (.vsh.yaml)"
                    },
                    "engine": {
                        "type": "string",
                        "description": "Execution engine (native, docker)",
                        "default": "docker"
                    },
                    "args": {
                        "type": "object",
                        "description": "Component-specific arguments",
                        "default": {}
                    }
                },
                "required": ["component_name", "component_config_path"]
            }
        ),
        Tool(
            name="build_docker_image",
            description="Build a Docker image from a Dockerfile",
            inputSchema={
                "type": "object",
                "properties": {
                    "dockerfile_path": {
                        "type": "string",
                        "description": "Path to the Dockerfile"
                    },
                    "image_tag": {
                        "type": "string",
                        "description": "Tag for the Docker image"
                    },
                    "context_path": {
                        "type": "string",
                        "description": "Build context directory",
                        "default": "."
                    }
                },
                "required": ["dockerfile_path", "image_tag"]
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
        Tool(
            name="read_file",
            description="Read contents of a file for analysis or editing",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the file to read"
                    }
                },
                "required": ["file_path"]
            }
        ),
        Tool(
            name="write_file",
            description="Write or create a file with specified content",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the file to write"
                    },
                    "content": {
                        "type": "string",
                        "description": "Content to write to the file"
                    }
                },
                "required": ["file_path", "content"]
            }
        ),
        Tool(
            name="list_directory",
            description="List contents of a directory",
            inputSchema={
                "type": "object",
                "properties": {
                    "directory_path": {
                        "type": "string",
                        "description": "Path to the directory to list"
                    },
                    "include_hidden": {
                        "type": "boolean",
                        "description": "Include hidden files and directories",
                        "default": False
                    }
                },
                "required": ["directory_path"]
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
    ]


@server.call_tool()
async def handle_call_tool(name: str, arguments: Dict[str, Any]) -> List[TextContent]:
    """Handle tool execution requests."""
    logger.info(f"Executing tool: {name} with arguments: {arguments}")

    if name == "echo_test":
        message = arguments.get("message", "")
        return [TextContent(type="text", text=f"Echo: {message}")]

    elif name == "list_available_tools":
        tools = await handle_list_tools()
        tool_list = []
        for tool in tools:
            tool_list.append({
                "name": tool.name,
                "description": tool.description,
                "required_params": tool.inputSchema.get("required", [])
            })
        return [TextContent(
            type="text",
            text=json.dumps(tool_list, indent=2)
        )]

    elif name == "run_nextflow_workflow":
        return await _execute_nextflow_workflow(arguments)

    elif name == "run_viash_component":
        return await _execute_viash_component(arguments)

    elif name == "build_docker_image":
        return await _build_docker_image(arguments)

    elif name == "analyze_nextflow_log":
        return await _analyze_nextflow_log(arguments)

    elif name == "read_file":
        return await _read_file(arguments)

    elif name == "write_file":
        return await _write_file(arguments)

    elif name == "list_directory":
        return await _list_directory(arguments)

    elif name == "validate_nextflow_config":
        return await _validate_nextflow_config(arguments)

    elif name == "check_environment":
        return await _check_environment(arguments)

    else:
        raise ValueError(f"Unknown tool: {name}")


async def _execute_nextflow_workflow(arguments: Dict[str, Any]) -> List[TextContent]:
    """Execute a Nextflow workflow."""
    workflow_name = arguments["workflow_name"]
    github_repo_url = arguments["github_repo_url"]
    profile = arguments.get("profile", "docker")
    params = arguments.get("params", {})
    config_file = arguments.get("config_file")

    # Build the command
    cmd = ["nextflow", "run", f"{github_repo_url}/{workflow_name}"]

    if profile:
        cmd.extend(["-profile", profile])

    if config_file:
        cmd.extend(["-c", config_file])

    # Add parameters
    for key, value in params.items():
        cmd.append(f"--{key}")
        cmd.append(str(value))

    try:
        # Execute the command
        logger.info(f"Executing command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )

        execution_result = {
            "command": " ".join(cmd),
            "exit_code": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "status": "completed" if result.returncode == 0 else "failed"
        }

        return [TextContent(
            type="text",
            text=json.dumps(execution_result, indent=2)
        )]

    except subprocess.TimeoutExpired:
        return [TextContent(
            type="text",
            text=json.dumps({
                "command": " ".join(cmd),
                "status": "timeout",
                "error": "Workflow execution timed out after 1 hour"
            }, indent=2)
        )]
    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "command": " ".join(cmd),
                "status": "error",
                "error": str(e)
            }, indent=2)
        )]


async def _execute_viash_component(arguments: Dict[str, Any]) -> List[TextContent]:
    """Execute a Viash component."""
    component_name = arguments["component_name"]
    component_config_path = arguments["component_config_path"]
    engine = arguments.get("engine", "docker")
    args = arguments.get("args", {})

    # Build the command
    cmd = ["viash", "run", component_config_path, "-p", engine]

    # Add component arguments
    if args:
        cmd.append("--")
        for key, value in args.items():
            cmd.append(f"--{key}")
            cmd.append(str(value))

    try:
        logger.info(f"Executing Viash component: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800  # 30 minutes timeout
        )

        execution_result = {
            "component": component_name,
            "command": " ".join(cmd),
            "exit_code": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "status": "completed" if result.returncode == 0 else "failed"
        }

        return [TextContent(
            type="text",
            text=json.dumps(execution_result, indent=2)
        )]

    except subprocess.TimeoutExpired:
        return [TextContent(
            type="text",
            text=json.dumps({
                "component": component_name,
                "command": " ".join(cmd),
                "status": "timeout",
                "error": "Component execution timed out after 30 minutes"
            }, indent=2)
        )]
    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "component": component_name,
                "command": " ".join(cmd),
                "status": "error",
                "error": str(e)
            }, indent=2)
        )]


async def _build_docker_image(arguments: Dict[str, Any]) -> List[TextContent]:
    """Build a Docker image."""
    dockerfile_path = arguments["dockerfile_path"]
    image_tag = arguments["image_tag"]
    context_path = arguments.get("context_path", ".")

    cmd = ["docker", "build", "-t", image_tag, "-f", dockerfile_path, context_path]

    try:
        logger.info(f"Building Docker image: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800  # 30 minutes timeout
        )

        build_result = {
            "image_tag": image_tag,
            "command": " ".join(cmd),
            "exit_code": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "status": "completed" if result.returncode == 0 else "failed"
        }

        return [TextContent(
            type="text",
            text=json.dumps(build_result, indent=2)
        )]

    except subprocess.TimeoutExpired:
        return [TextContent(
            type="text",
            text=json.dumps({
                "image_tag": image_tag,
                "command": " ".join(cmd),
                "status": "timeout",
                "error": "Docker build timed out after 30 minutes"
            }, indent=2)
        )]
    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "image_tag": image_tag,
                "command": " ".join(cmd),
                "status": "error",
                "error": str(e)
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


async def _read_file(arguments: Dict[str, Any]) -> List[TextContent]:
    """Read contents of a file for analysis or editing."""
    file_path = arguments["file_path"]

    try:
        with open(file_path, 'r') as f:
            file_content = f.read()
        return [TextContent(type="text", text=file_content)]
    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "status": "error",
                "error": f"Failed to read file: {str(e)}"
            }, indent=2)
        )]


async def _write_file(arguments: Dict[str, Any]) -> List[TextContent]:
    """Write or create a file with specified content."""
    file_path = arguments["file_path"]
    content = arguments["content"]

    try:
        with open(file_path, 'w') as f:
            f.write(content)
        return [TextContent(type="text", text="File written successfully")]
    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "status": "error",
                "error": f"Failed to write file: {str(e)}"
            }, indent=2)
        )]


async def _list_directory(arguments: Dict[str, Any]) -> List[TextContent]:
    """List contents of a directory."""
    directory_path = arguments["directory_path"]
    include_hidden = arguments.get("include_hidden", False)

    try:
        entries = []
        for entry in Path(directory_path).iterdir():
            if include_hidden or not entry.name.startswith('.'):
                entries.append({
                    "name": entry.name,
                    "is_directory": entry.is_dir(),
                    "size": entry.stat().st_size
                })
        return [TextContent(
            type="text",
            text=json.dumps(entries, indent=2)
        )]
    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({
                "status": "error",
                "error": f"Failed to list directory: {str(e)}"
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
