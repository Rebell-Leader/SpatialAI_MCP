#!/usr/bin/env python3
"""
Gradio Web Interface for OpenProblems MCP Server Tools

This module provides a visual web interface for testing and using our MCP tools
while maintaining the full MCP server functionality in parallel.
"""

import gradio as gr
import asyncio
import json
from typing import Any, Dict, List, Optional
from pathlib import Path

# Import our existing MCP server tools
from .main import (
    handle_call_tool,
    handle_list_tools,
    handle_read_resource,
    handle_list_resources
)


class OpenProblemsMCPInterface:
    """Gradio interface wrapper for OpenProblems MCP Server tools."""

    def __init__(self):
        self.tools = None
        self.resources = None
        # Run initialization in a new event loop
        asyncio.run(self.initialize())

    async def initialize(self):
        """Initialize tools and resources."""
        self.tools = await handle_list_tools()
        self.resources = await handle_list_resources()

    async def check_environment(self, tools_to_check: str = "nextflow,viash,docker,java") -> str:
        """
        Check if required bioinformatics tools are installed and available.

        Args:
            tools_to_check (str): Comma-separated list of tools to check

        Returns:
            str: Environment check results in JSON format
        """
        tools_list = [tool.strip() for tool in tools_to_check.split(",")]

        try:
            result = await handle_call_tool("check_environment", {
                "tools": tools_list
            })
            return result[0].text
        except Exception as e:
            return f"Error: {str(e)}"

    async def validate_nextflow_config(self, pipeline_path: str, config_path: str = "") -> str:
        """
        Validate Nextflow pipeline syntax and configuration.

        Args:
            pipeline_path (str): Path to the Nextflow pipeline file (.nf)
            config_path (str): Optional path to nextflow.config file

        Returns:
            str: Validation results in JSON format
        """
        args = {"pipeline_path": pipeline_path}
        if config_path:
            args["config_path"] = config_path

        try:
            result = await handle_call_tool("validate_nextflow_config", args)
            return result[0].text
        except Exception as e:
            return f"Error: {str(e)}"

    async def run_nextflow_workflow(
        self,
        workflow_name: str,
        github_repo_url: str,
        profile: str = "docker",
        params_json: str = "{}"
    ) -> str:
        """
        Execute a Nextflow workflow from OpenProblems repositories.

        Args:
            workflow_name (str): Name of the workflow (e.g., main.nf)
            github_repo_url (str): GitHub repository URL
            profile (str): Nextflow profile to use
            params_json (str): Pipeline parameters as JSON string

        Returns:
            str: Execution results in JSON format
        """
        try:
            params = json.loads(params_json) if params_json.strip() else {}
            result = await handle_call_tool("run_nextflow_workflow", {
                "workflow_name": workflow_name,
                "github_repo_url": github_repo_url,
                "profile": profile,
                "params": params
            })
            return result[0].text
        except Exception as e:
            return f"Error: {str(e)}"

    async def analyze_nextflow_log(self, log_file_path: str) -> str:
        """
        Analyze Nextflow execution logs for errors and troubleshooting insights.

        Args:
            log_file_path (str): Path to the .nextflow.log file

        Returns:
            str: Log analysis results in JSON format
        """
        try:
            result = await handle_call_tool("analyze_nextflow_log", {
                "log_file_path": log_file_path
            })
            return result[0].text
        except Exception as e:
            return f"Error: {str(e)}"

    async def read_file(self, file_path: str) -> str:
        """
        Read and display file contents for analysis.

        Args:
            file_path (str): Path to the file to read

        Returns:
            str: File contents or error message
        """
        try:
            result = await handle_call_tool("read_file", {
                "file_path": file_path
            })
            return result[0].text
        except Exception as e:
            return f"Error: {str(e)}"

    async def write_file(self, file_path: str, content: str) -> str:
        """
        Write content to a file.

        Args:
            file_path (str): Path where to write the file
            content (str): Content to write

        Returns:
            str: Success message or error
        """
        try:
            result = await handle_call_tool("write_file", {
                "file_path": file_path,
                "content": content
            })
            return result[0].text
        except Exception as e:
            return f"Error: {str(e)}"

    async def list_directory(self, directory_path: str, include_hidden: bool = False) -> str:
        """
        List contents of a directory.

        Args:
            directory_path (str): Path to the directory
            include_hidden (bool): Whether to include hidden files

        Returns:
            str: Directory listing in JSON format
        """
        try:
            result = await handle_call_tool("list_directory", {
                "directory_path": directory_path,
                "include_hidden": include_hidden
            })
            return result[0].text
        except Exception as e:
            return f"Error: {str(e)}"

    async def get_documentation(self, doc_type: str) -> str:
        """
        Get documentation resources.

        Args:
            doc_type (str): Type of documentation (nextflow, viash, docker, spatial-workflows)

        Returns:
            str: Documentation content
        """
        uri_mapping = {
            "nextflow": "documentation://nextflow",
            "viash": "documentation://viash",
            "docker": "documentation://docker",
            "spatial-workflows": "templates://spatial-workflows",
            "server-status": "server://status"
        }

        uri = uri_mapping.get(doc_type)
        if not uri:
            return f"Invalid documentation type. Available: {list(uri_mapping.keys())}"

        try:
            result = await handle_read_resource(uri)
            return result
        except Exception as e:
            return f"Error: {str(e)}"


def create_gradio_interface():
    """Create the Gradio interface for OpenProblems MCP Server."""

    mcp_interface = OpenProblemsMCPInterface()

    with gr.Blocks(
        title="OpenProblems Spatial Transcriptomics MCP Server",
        theme=gr.themes.Soft(),
        css="""
        .gradio-container { max-width: 1200px; margin: auto; }
        .tool-section { border: 1px solid #e0e0e0; border-radius: 8px; padding: 20px; margin: 10px 0; }
        """
    ) as demo:

        gr.Markdown("""
        # 🧬 OpenProblems Spatial Transcriptomics MCP Server

        **Visual interface for testing MCP tools and accessing documentation resources.**

        This interface provides access to the same tools available through the MCP protocol,
        allowing you to test functionality before integrating with AI agents like Continue.dev.
        """)

        with gr.Tabs():

            # Environment Tools Tab
            with gr.Tab("🔧 Environment & Validation"):
                gr.Markdown("### Environment Validation")
                with gr.Row():
                    tools_input = gr.Textbox(
                        value="nextflow,viash,docker,java",
                        label="Tools to Check",
                        placeholder="Comma-separated list of tools"
                    )
                    check_btn = gr.Button("Check Environment", variant="primary")

                env_output = gr.JSON(label="Environment Check Results")
                check_btn.click(mcp_interface.check_environment, tools_input, env_output)

                gr.Markdown("### Nextflow Configuration Validation")
                with gr.Row():
                    pipeline_path = gr.Textbox(label="Pipeline Path", placeholder="path/to/main.nf")
                    config_path = gr.Textbox(label="Config Path (optional)", placeholder="path/to/nextflow.config")

                validate_btn = gr.Button("Validate Configuration", variant="primary")
                validate_output = gr.JSON(label="Validation Results")
                validate_btn.click(
                    mcp_interface.validate_nextflow_config,
                    [pipeline_path, config_path],
                    validate_output
                )

            # Workflow Execution Tab
            with gr.Tab("⚡ Workflow Execution"):
                gr.Markdown("### Execute Nextflow Workflow")
                with gr.Row():
                    workflow_name = gr.Textbox(
                        label="Workflow Name",
                        value="main.nf",
                        placeholder="main.nf"
                    )
                    repo_url = gr.Textbox(
                        label="GitHub Repository URL",
                        placeholder="https://github.com/openproblems-bio/task_spatial_decomposition"
                    )

                with gr.Row():
                    profile = gr.Dropdown(
                        choices=["docker", "singularity", "conda", "test"],
                        value="docker",
                        label="Profile"
                    )
                    params_json = gr.Textbox(
                        label="Parameters (JSON)",
                        value='{"input": "data.h5ad", "output": "results/"}',
                        placeholder='{"key": "value"}'
                    )

                run_btn = gr.Button("Run Workflow", variant="primary")
                workflow_output = gr.JSON(label="Workflow Execution Results")
                run_btn.click(
                    mcp_interface.run_nextflow_workflow,
                    [workflow_name, repo_url, profile, params_json],
                    workflow_output
                )

            # File Management Tab
            with gr.Tab("📁 File Management"):
                with gr.Row():
                    with gr.Column():
                        gr.Markdown("### List Directory")
                        dir_path = gr.Textbox(label="Directory Path", value=".")
                        include_hidden = gr.Checkbox(label="Include Hidden Files")
                        list_btn = gr.Button("List Directory")
                        list_output = gr.JSON(label="Directory Contents")
                        list_btn.click(
                            mcp_interface.list_directory,
                            [dir_path, include_hidden],
                            list_output
                        )

                    with gr.Column():
                        gr.Markdown("### Read File")
                        read_path = gr.Textbox(label="File Path", placeholder="path/to/file.txt")
                        read_btn = gr.Button("Read File")
                        read_output = gr.Textbox(label="File Contents", lines=10)
                        read_btn.click(mcp_interface.read_file, read_path, read_output)

                gr.Markdown("### Write File")
                with gr.Row():
                    write_path = gr.Textbox(label="File Path", placeholder="path/to/new_file.txt")
                    write_content = gr.Textbox(label="Content", lines=5, placeholder="File content here...")

                write_btn = gr.Button("Write File", variant="primary")
                write_output = gr.Textbox(label="Write Result")
                write_btn.click(
                    mcp_interface.write_file,
                    [write_path, write_content],
                    write_output
                )

            # Log Analysis Tab
            with gr.Tab("🔍 Log Analysis"):
                gr.Markdown("### Nextflow Log Analysis")
                log_path = gr.Textbox(
                    label="Log File Path",
                    placeholder="path/to/.nextflow.log",
                    value="work/.nextflow.log"
                )
                analyze_btn = gr.Button("Analyze Log", variant="primary")
                log_output = gr.JSON(label="Log Analysis Results")
                analyze_btn.click(mcp_interface.analyze_nextflow_log, log_path, log_output)

            # Documentation Tab
            with gr.Tab("📚 Documentation & Resources"):
                gr.Markdown("### Access MCP Resources")
                doc_type = gr.Dropdown(
                    choices=["nextflow", "viash", "docker", "spatial-workflows", "server-status"],
                    value="nextflow",
                    label="Documentation Type"
                )
                doc_btn = gr.Button("Get Documentation", variant="primary")
                doc_output = gr.Textbox(label="Documentation Content", lines=20)
                doc_btn.click(mcp_interface.get_documentation, doc_type, doc_output)

        gr.Markdown("""
        ---
        ### 🤖 AI Agent Integration

        To use these tools with AI agents like Continue.dev, add this to your `~/.continue/config.json`:

        ```json
        {
          "experimental": {
            "modelContextProtocolServers": [
              {
                "name": "openproblems-spatial",
                "transport": {
                  "type": "stdio",
                  "command": "python",
                  "args": ["-m", "mcp_server.main"],
                  "cwd": "/path/to/your/SpatialAI_MCP"
                }
              }
            ]
          }
        }
        ```

        **📖 Documentation**: [Setup Guide](docs/CONTINUE_DEV_SETUP.md) | [Agent Rules](docs/AGENT_RULES.md)
        """)

    return demo


def launch_gradio_interface(share: bool = False, server_port: int = 7860):
    """Launch the Gradio interface."""
    demo = create_gradio_interface()

    print("🚀 Starting OpenProblems MCP Server Gradio Interface...")
    print(f"📱 Web Interface: http://localhost:{server_port}")
    print("🤖 MCP Server: Use 'python -m mcp_server.main' for AI agents")

    demo.launch(
        share=share,
        server_port=server_port,
        server_name="0.0.0.0",
        show_error=True,
        # Note: Not setting mcp_server=True to avoid conflicts with our main MCP server
    )


if __name__ == "__main__":
    launch_gradio_interface()
