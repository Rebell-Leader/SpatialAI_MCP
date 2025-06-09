#!/usr/bin/env python3
"""
Simple MCP Client Example for OpenProblems Spatial Transcriptomics

This example demonstrates how to connect to and interact with the
OpenProblems Spatial Transcriptomics MCP Server.
"""

import asyncio
import json
import subprocess
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client


async def demo_mcp_interaction():
    """Demonstrate basic interactions with the MCP server."""

    print("🚀 Starting OpenProblems Spatial Transcriptomics MCP Client Demo")
    print("=" * 60)

    # Configure server parameters
    server_params = StdioServerParameters(
        command="python",
        args=["-m", "mcp_server.main"],
        env=None,
    )

    try:
        # Connect to the MCP server
        async with stdio_client(server_params) as (read, write):
            async with ClientSession(read, write) as session:
                print("✅ Connected to MCP server")

                # Initialize the session
                await session.initialize()
                print("✅ Session initialized")

                # List available resources
                print("\n📚 Available Resources:")
                print("-" * 30)
                resources = await session.list_resources()
                for resource in resources:
                    print(f"  • {resource.name}: {resource.description}")

                # List available tools
                print("\n🛠️  Available Tools:")
                print("-" * 30)
                tools = await session.list_tools()
                for tool in tools:
                    print(f"  • {tool.name}: {tool.description}")

                # Test echo tool
                print("\n🔄 Testing Echo Tool:")
                print("-" * 30)
                echo_result = await session.call_tool(
                    "echo_test",
                    arguments={"message": "Hello from MCP client!"}
                )
                print(f"Echo response: {echo_result}")

                # Read server status
                print("\n📊 Server Status:")
                print("-" * 30)
                try:
                    status_content = await session.read_resource("server://status")
                    status_data = json.loads(status_content)
                    print(f"Server Name: {status_data['server_name']}")
                    print(f"Version: {status_data['version']}")
                    print(f"Status: {status_data['status']}")
                    print("Capabilities:")
                    for capability, enabled in status_data['capabilities'].items():
                        status_icon = "✅" if enabled else "❌"
                        print(f"  {status_icon} {capability}")
                except Exception as e:
                    print(f"Error reading server status: {e}")

                # Read documentation examples
                print("\n📖 Sample Documentation:")
                print("-" * 30)
                try:
                    nextflow_docs = await session.read_resource("documentation://nextflow")
                    docs_data = json.loads(nextflow_docs)
                    print("Nextflow Best Practices:")
                    for practice, description in docs_data['best_practices'].items():
                        print(f"  • {practice}: {description}")
                except Exception as e:
                    print(f"Error reading documentation: {e}")

                # List available tools using the MCP tool
                print("\n🔍 Detailed Tool Information:")
                print("-" * 30)
                try:
                    tools_result = await session.call_tool("list_available_tools", arguments={})
                    tools_data = json.loads(tools_result)
                    for tool in tools_data:
                        print(f"  • {tool['name']}")
                        print(f"    Description: {tool['description']}")
                        required_params = tool.get('required_params', [])
                        if required_params:
                            print(f"    Required params: {', '.join(required_params)}")
                        print()
                except Exception as e:
                    print(f"Error listing tools: {e}")

                # Read pipeline templates
                print("\n🧬 Spatial Transcriptomics Pipeline Templates:")
                print("-" * 30)
                try:
                    templates_content = await session.read_resource("templates://spatial-workflows")
                    templates_data = json.loads(templates_content)
                    for template_id, template_info in templates_data.items():
                        print(f"  • {template_info['name']}")
                        print(f"    Description: {template_info['description']}")
                        print(f"    Inputs: {', '.join(template_info['inputs'])}")
                        print(f"    Outputs: {', '.join(template_info['outputs'])}")
                        print()
                except Exception as e:
                    print(f"Error reading templates: {e}")

        print("✅ Demo completed successfully!")

    except Exception as e:
        print(f"❌ Error during demo: {e}")
        return False

    return True


async def demo_workflow_execution():
    """Demonstrate workflow execution capabilities (if tools are available)."""

    print("\n🧪 Workflow Execution Demo")
    print("=" * 60)

    # Check if required tools are available
    required_tools = ["nextflow", "docker"]
    missing_tools = []

    for tool in required_tools:
        try:
            result = subprocess.run([tool, "--version"],
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                print(f"✅ {tool} is available")
            else:
                missing_tools.append(tool)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_tools.append(tool)

    if missing_tools:
        print(f"⚠️  Missing tools: {', '.join(missing_tools)}")
        print("   Workflow execution demo skipped")
        return

    # Configure server parameters
    server_params = StdioServerParameters(
        command="python",
        args=["-m", "mcp_server.main"],
        env=None,
    )

    try:
        async with stdio_client(server_params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()

                # Example: Analyze a mock Nextflow log
                print("\n📋 Testing Log Analysis:")
                print("-" * 30)

                # Create a mock log file for testing
                mock_log_path = Path("/tmp/test_nextflow.log")
                mock_log_content = """
N E X T F L O W  ~  version 23.04.0
Launching `main.nf` [abc123] DSL2 - revision: def456

executor >  local (2)
[12/abc123] process > PROCESS_1 [100%] 1 of 1 ✓
[34/def456] process > PROCESS_2 [  0%] 0 of 1, failed: 1

ERROR ~ Error executing process > 'PROCESS_2'
Caused by:
  Process `PROCESS_2` terminated with an error exit status (137)

Command executed:
  python script.py --input data.h5ad

Command exit status:
  137

Execution failed
"""

                try:
                    with open(mock_log_path, 'w') as f:
                        f.write(mock_log_content)

                    # Analyze the log using MCP
                    log_analysis = await session.call_tool(
                        "analyze_nextflow_log",
                        arguments={"log_file_path": str(mock_log_path)}
                    )

                    analysis_data = json.loads(log_analysis)
                    print(f"Log analysis completed:")
                    print(f"  File size: {analysis_data['file_size']} bytes")
                    print(f"  Execution status: {analysis_data['execution_status']}")

                    if analysis_data['issues_found']:
                        print("  Issues found:")
                        for issue in analysis_data['issues_found']:
                            print(f"    • {issue['issue']}: {issue['suggestion']}")

                    # Clean up
                    mock_log_path.unlink(missing_ok=True)

                except Exception as e:
                    print(f"Error in log analysis demo: {e}")
                    mock_log_path.unlink(missing_ok=True)

    except Exception as e:
        print(f"❌ Error during workflow demo: {e}")


async def main():
    """Main function to run the demo."""

    print("🧬 OpenProblems Spatial Transcriptomics MCP Client")
    print("   Model Context Protocol Demo")
    print("   Version 0.1.0")
    print()

    # Run basic interaction demo
    success = await demo_mcp_interaction()

    if success:
        # Run workflow execution demo
        await demo_workflow_execution()

    print("\n" + "=" * 60)
    print("Demo completed! 🎉")
    print("\nTo use this MCP server with AI agents:")
    print("1. Start the server: python -m mcp_server.main")
    print("2. Configure your AI agent to connect via stdio transport")
    print("3. Use the available tools and resources for spatial transcriptomics workflows")


if __name__ == "__main__":
    # Check if the server module is available
    try:
        import mcp_server.main
    except ImportError:
        print("❌ MCP server module not found. Make sure you're in the project directory")
        print("   and have installed the package: pip install -e .")
        sys.exit(1)

    # Run the demo
    asyncio.run(main())
