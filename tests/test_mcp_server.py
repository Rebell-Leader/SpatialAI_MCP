#!/usr/bin/env python3
"""
Test suite for the OpenProblems Spatial Transcriptomics MCP Server.
"""

import asyncio
import json
import pytest
from unittest.mock import AsyncMock, MagicMock, patch

# Import the server components
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent / "src"))

from mcp_server.main import (
    handle_list_resources,
    handle_read_resource,
    handle_list_tools,
    handle_call_tool,
)


class TestMCPServer:
    """Test cases for the MCP server functionality."""

    @pytest.mark.asyncio
    async def test_list_resources(self):
        """Test that resources are properly listed."""
        resources = await handle_list_resources()

        assert len(resources) == 5
        resource_uris = [r.uri for r in resources]

        expected_uris = [
            "server://status",
            "documentation://nextflow",
            "documentation://viash",
            "documentation://docker",
            "templates://spatial-workflows"
        ]

        for uri in expected_uris:
            assert uri in resource_uris

    @pytest.mark.asyncio
    async def test_read_server_status_resource(self):
        """Test reading the server status resource."""
        status_content = await handle_read_resource("server://status")
        status_data = json.loads(status_content)

        assert status_data["server_name"] == "OpenProblems Spatial Transcriptomics MCP"
        assert status_data["version"] == "0.1.0"
        assert status_data["status"] == "running"
        assert "capabilities" in status_data
        assert status_data["capabilities"]["nextflow_execution"] is True

    @pytest.mark.asyncio
    async def test_read_documentation_resources(self):
        """Test reading documentation resources."""
        # Test Nextflow documentation
        nextflow_docs = await handle_read_resource("documentation://nextflow")
        nextflow_data = json.loads(nextflow_docs)
        assert "best_practices" in nextflow_data
        assert "dsl_version" in nextflow_data["best_practices"]

        # Test Viash documentation
        viash_docs = await handle_read_resource("documentation://viash")
        viash_data = json.loads(viash_docs)
        assert "component_structure" in viash_data
        assert "best_practices" in viash_data

        # Test Docker documentation
        docker_docs = await handle_read_resource("documentation://docker")
        docker_data = json.loads(docker_docs)
        assert "dockerfile_optimization" in docker_data
        assert "bioinformatics_specific" in docker_data

    @pytest.mark.asyncio
    async def test_read_templates_resource(self):
        """Test reading pipeline templates resource."""
        templates_content = await handle_read_resource("templates://spatial-workflows")
        templates_data = json.loads(templates_content)

        expected_templates = [
            "basic_preprocessing",
            "spatially_variable_genes",
            "label_transfer"
        ]

        for template in expected_templates:
            assert template in templates_data
            assert "name" in templates_data[template]
            assert "description" in templates_data[template]
            assert "inputs" in templates_data[template]
            assert "outputs" in templates_data[template]

    @pytest.mark.asyncio
    async def test_invalid_resource_uri(self):
        """Test handling of invalid resource URIs."""
        with pytest.raises(ValueError, match="Unknown resource URI"):
            await handle_read_resource("invalid://resource")

    @pytest.mark.asyncio
    async def test_list_tools(self):
        """Test that tools are properly listed."""
        tools = await handle_list_tools()

        expected_tools = [
            "echo_test",
            "list_available_tools",
            "run_nextflow_workflow",
            "run_viash_component",
            "build_docker_image",
            "analyze_nextflow_log"
        ]

        tool_names = [t.name for t in tools]

        for tool_name in expected_tools:
            assert tool_name in tool_names

        # Check that tools have proper schemas
        for tool in tools:
            assert hasattr(tool, 'inputSchema')
            assert 'type' in tool.inputSchema
            assert tool.inputSchema['type'] == 'object'

    @pytest.mark.asyncio
    async def test_echo_test_tool(self):
        """Test the echo test tool."""
        result = await handle_call_tool("echo_test", {"message": "Hello MCP!"})

        assert len(result) == 1
        assert result[0].type == "text"
        assert result[0].text == "Echo: Hello MCP!"

    @pytest.mark.asyncio
    async def test_list_available_tools_tool(self):
        """Test the list available tools tool."""
        result = await handle_call_tool("list_available_tools", {})

        assert len(result) == 1
        assert result[0].type == "text"

        tools_data = json.loads(result[0].text)
        assert isinstance(tools_data, list)
        assert len(tools_data) >= 6  # We have at least 6 tools

        # Check structure of tool entries
        for tool in tools_data:
            assert "name" in tool
            assert "description" in tool
            assert "required_params" in tool

    @pytest.mark.asyncio
    async def test_invalid_tool_name(self):
        """Test handling of invalid tool names."""
        with pytest.raises(ValueError, match="Unknown tool"):
            await handle_call_tool("invalid_tool", {})

    @pytest.mark.asyncio
    @patch('mcp_server.main.subprocess.run')
    async def test_nextflow_workflow_execution(self, mock_subprocess):
        """Test Nextflow workflow execution tool."""
        # Mock successful subprocess execution
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Nextflow execution completed successfully"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result

        arguments = {
            "workflow_name": "main.nf",
            "github_repo_url": "https://github.com/openproblems-bio/test-workflow",
            "profile": "docker",
            "params": {"input": "test.h5ad", "output": "results/"}
        }

        result = await handle_call_tool("run_nextflow_workflow", arguments)

        assert len(result) == 1
        assert result[0].type == "text"

        execution_data = json.loads(result[0].text)
        assert execution_data["status"] == "completed"
        assert execution_data["exit_code"] == 0

    @pytest.mark.asyncio
    @patch('mcp_server.main.subprocess.run')
    async def test_viash_component_execution(self, mock_subprocess):
        """Test Viash component execution tool."""
        # Mock successful subprocess execution
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Viash component executed successfully"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result

        arguments = {
            "component_name": "test_component",
            "component_config_path": "config.vsh.yaml",
            "engine": "docker",
            "args": {"input": "test.h5ad", "output": "result.h5ad"}
        }

        result = await handle_call_tool("run_viash_component", arguments)

        assert len(result) == 1
        assert result[0].type == "text"

        execution_data = json.loads(result[0].text)
        assert execution_data["status"] == "completed"
        assert execution_data["exit_code"] == 0
        assert execution_data["component"] == "test_component"

    @pytest.mark.asyncio
    @patch('mcp_server.main.subprocess.run')
    async def test_docker_image_build(self, mock_subprocess):
        """Test Docker image building tool."""
        # Mock successful subprocess execution
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Successfully built docker image"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result

        arguments = {
            "dockerfile_path": "Dockerfile",
            "image_tag": "openproblems/test:latest",
            "context_path": "."
        }

        result = await handle_call_tool("build_docker_image", arguments)

        assert len(result) == 1
        assert result[0].type == "text"

        build_data = json.loads(result[0].text)
        assert build_data["status"] == "completed"
        assert build_data["exit_code"] == 0
        assert build_data["image_tag"] == "openproblems/test:latest"

    @pytest.mark.asyncio
    @patch('mcp_server.main.Path')
    async def test_nextflow_log_analysis(self, mock_path):
        """Test Nextflow log analysis tool."""
        # Mock log file content
        mock_log_content = """
        N E X T F L O W  ~  version 23.04.0
        Launching `main.nf` [abc123] DSL2 - revision: def456

        executor >  local (4)
        [12/abc123] process > PROCESS_1 [100%] 2 of 2 ✓
        [34/def456] process > PROCESS_2 [100%] 2 of 2, failed: 1, retries: 1 ✗

        ERROR ~ Error executing process > 'PROCESS_2'

        Caused by:
          Process `PROCESS_2` terminated with an error exit status (137)

        Command executed:
          python script.py --input data.h5ad --output result.h5ad

        Command exit status:
          137

        Execution failed
        """

        # Mock file operations
        mock_log_path = MagicMock()
        mock_log_path.exists.return_value = True
        mock_log_path.stat.return_value.st_size = len(mock_log_content)
        mock_path.return_value = mock_log_path

        # Mock file reading
        with patch('builtins.open', mock_open(read_data=mock_log_content)):
            arguments = {"log_file_path": "/path/to/.nextflow.log"}
            result = await handle_call_tool("analyze_nextflow_log", arguments)

            assert len(result) == 1
            assert result[0].type == "text"

            analysis_data = json.loads(result[0].text)
            assert "issues_found" in analysis_data
            assert "execution_status" in analysis_data
            assert analysis_data["execution_status"] == "failed"

            # Check that OOM error was detected
            issues = analysis_data["issues_found"]
            oom_issue = next((issue for issue in issues if "exit status 137" in issue["pattern"]), None)
            assert oom_issue is not None
            assert "Out of memory" in oom_issue["issue"]


def mock_open(read_data):
    """Mock file opening for testing."""
    from unittest.mock import mock_open as mock_open_builtin
    return mock_open_builtin(read_data=read_data)


if __name__ == "__main__":
    pytest.main([__file__])
