#!/usr/bin/env python3
"""
Command-line interface for the OpenProblems Spatial Transcriptomics MCP Server.
"""

import asyncio
import click
import logging
import sys
from pathlib import Path

from .main import main as run_server


@click.group()
@click.version_option(version="0.1.0")
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
@click.option("--config", "-c", type=click.Path(exists=True), help="Configuration file path")
def cli(verbose, config):
    """OpenProblems Spatial Transcriptomics MCP Server CLI."""
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if config:
        # TODO: Load configuration from file
        click.echo(f"Using configuration from: {config}")


@cli.command()
@click.option("--host", default="localhost", help="Host to bind to (HTTP transport)")
@click.option("--port", default=8000, help="Port to bind to (HTTP transport)")
@click.option("--transport", default="stdio", type=click.Choice(["stdio", "http"]),
              help="Transport method")
def serve(host, port, transport):
    """Start the MCP server."""
    click.echo("🚀 Starting OpenProblems Spatial Transcriptomics MCP Server")
    click.echo(f"   Transport: {transport}")

    if transport == "http":
        click.echo(f"   Host: {host}")
        click.echo(f"   Port: {port}")
        click.echo("   Note: HTTP transport is not yet implemented")
        sys.exit(1)

    try:
        asyncio.run(run_server())
    except KeyboardInterrupt:
        click.echo("\n👋 Server stopped")
    except Exception as e:
        click.echo(f"❌ Server error: {e}", err=True)
        sys.exit(1)


@cli.command()
def test():
    """Run the test suite."""
    import subprocess

    click.echo("🧪 Running test suite...")

    try:
        result = subprocess.run(["pytest", "tests/", "-v"],
                              capture_output=True, text=True)

        click.echo(result.stdout)
        if result.stderr:
            click.echo(result.stderr, err=True)

        if result.returncode == 0:
            click.echo("✅ All tests passed!")
        else:
            click.echo("❌ Some tests failed")
            sys.exit(1)

    except FileNotFoundError:
        click.echo("❌ pytest not found. Install with: pip install pytest", err=True)
        sys.exit(1)


@cli.command()
def demo():
    """Run the interactive demo client."""
    click.echo("🎬 Starting MCP client demo...")

    try:
        import subprocess
        result = subprocess.run([sys.executable, "examples/simple_client.py"])
        sys.exit(result.returncode)
    except Exception as e:
        click.echo(f"❌ Demo error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--check-tools", is_flag=True, help="Check if external tools are available")
@click.option("--check-deps", is_flag=True, help="Check Python dependencies")
def doctor(check_tools, check_deps):
    """Diagnose installation and configuration issues."""
    click.echo("🔍 OpenProblems MCP Server Health Check")
    click.echo("=" * 50)

    all_good = True

    # Check Python imports
    click.echo("\n📦 Python Dependencies:")
    dependencies = [
        ("mcp", "MCP Python SDK"),
        ("yaml", "PyYAML"),
        ("docker", "Docker Python client"),
        ("pandas", "Pandas"),
        ("numpy", "NumPy"),
    ]

    for module, description in dependencies:
        try:
            __import__(module)
            click.echo(f"  ✅ {description}")
        except ImportError:
            click.echo(f"  ❌ {description} - not installed")
            all_good = False

    # Check external tools
    if check_tools:
        click.echo("\n🛠️  External Tools:")
        tools = [
            ("nextflow", "Nextflow workflow engine"),
            ("viash", "Viash component framework"),
            ("docker", "Docker containerization"),
            ("java", "Java runtime (required for Nextflow)"),
        ]

        import subprocess
        for tool, description in tools:
            try:
                result = subprocess.run([tool, "--version"],
                                      capture_output=True, timeout=10)
                if result.returncode == 0:
                    click.echo(f"  ✅ {description}")
                else:
                    click.echo(f"  ❌ {description} - not working properly")
                    all_good = False
            except (subprocess.TimeoutExpired, FileNotFoundError):
                click.echo(f"  ❌ {description} - not found")
                all_good = False

    # Check directories
    click.echo("\n📁 Directory Structure:")
    directories = ["data", "work", "logs", "cache"]

    for directory in directories:
        path = Path(directory)
        if path.exists():
            if path.is_dir():
                click.echo(f"  ✅ {directory}/ - exists")
            else:
                click.echo(f"  ❌ {directory} - exists but not a directory")
                all_good = False
        else:
            click.echo(f"  ⚠️  {directory}/ - missing (will be created)")
            try:
                path.mkdir(exist_ok=True)
                click.echo(f"     Created {directory}/")
            except Exception as e:
                click.echo(f"     Failed to create: {e}")
                all_good = False

    # Check server module
    click.echo("\n🖥️  Server Module:")
    try:
        from . import main
        click.echo("  ✅ MCP server module - importable")

        # Test basic functionality
        import asyncio
        async def test_handlers():
            try:
                resources = await main.handle_list_resources()
                tools = await main.handle_list_tools()
                click.echo(f"  ✅ Server handlers - working ({len(resources)} resources, {len(tools)} tools)")
            except Exception as e:
                click.echo(f"  ❌ Server handlers - error: {e}")
                return False
            return True

        handler_ok = asyncio.run(test_handlers())
        all_good = all_good and handler_ok

    except ImportError as e:
        click.echo(f"  ❌ MCP server module - import error: {e}")
        all_good = False

    # Summary
    click.echo("\n" + "=" * 50)
    if all_good:
        click.echo("✅ All checks passed! Your setup is ready.")
    else:
        click.echo("❌ Some issues found. Please fix them before running the server.")
        click.echo("\nFor help, see: docs/SETUP.md")
        sys.exit(1)


@cli.command()
def download_docs():
    """Download and cache documentation from OpenProblems, Nextflow, and Viash."""
    click.echo("📚 Downloading documentation from OpenProblems, Nextflow, and Viash...")

    async def download():
        from .documentation_generator_simple import DocumentationGenerator

        try:
            generator = DocumentationGenerator()
            documentation = await generator.generate_all_documentation()

            click.echo("\n📊 Documentation download complete!")
            total_chars = 0
            for source, content in documentation.items():
                chars = len(content)
                total_chars += chars
                click.echo(f"   ✅ {source}: {chars:,} characters")

            click.echo(f"\n🎉 Total: {total_chars:,} characters of documentation cached!")
            click.echo("   Documentation is now available in your MCP server resources.")

        except Exception as e:
            click.echo(f"❌ Failed to download documentation: {e}")
            sys.exit(1)

    asyncio.run(download())


@cli.command()
@click.argument("tool_name")
@click.argument("arguments", nargs=-1)
def tool(tool_name, arguments):
    """Execute a specific MCP tool directly."""
    click.echo(f"🔧 Executing tool: {tool_name}")

    # Parse arguments (simple key=value format)
    tool_args = {}
    for arg in arguments:
        if "=" in arg:
            key, value = arg.split("=", 1)
            tool_args[key] = value
        else:
            click.echo(f"❌ Invalid argument format: {arg}")
            click.echo("   Use: key=value format")
            sys.exit(1)

    click.echo(f"   Arguments: {tool_args}")

    async def run_tool():
        from .main import handle_call_tool
        try:
            result = await handle_call_tool(tool_name, tool_args)
            click.echo("\n📄 Result:")
            for item in result:
                click.echo(item.text)
        except Exception as e:
            click.echo(f"❌ Tool execution failed: {e}", err=True)
            sys.exit(1)

    asyncio.run(run_tool())


@cli.command()
def info():
    """Show server information and available tools/resources."""
    click.echo("📋 OpenProblems Spatial Transcriptomics MCP Server")
    click.echo("   Version: 0.1.0")
    click.echo("   Protocol: Model Context Protocol (MCP)")
    click.echo("   Purpose: Spatial transcriptomics workflow automation")

    async def show_info():
        from .main import handle_list_resources, handle_list_tools

        try:
            resources = await handle_list_resources()
            tools = await handle_list_tools()

            click.echo(f"\n📚 Available Resources ({len(resources)}):")
            for resource in resources:
                click.echo(f"  • {resource.name}")
                click.echo(f"    URI: {resource.uri}")
                click.echo(f"    Description: {resource.description}")
                click.echo()

            click.echo(f"🛠️  Available Tools ({len(tools)}):")
            for tool in tools:
                click.echo(f"  • {tool.name}")
                click.echo(f"    Description: {tool.description}")
                required = tool.inputSchema.get("required", [])
                if required:
                    click.echo(f"    Required parameters: {', '.join(required)}")
                click.echo()

        except Exception as e:
            click.echo(f"❌ Error getting server info: {e}", err=True)

    asyncio.run(show_info())


def main():
    """Main CLI entry point."""
    cli()


if __name__ == "__main__":
    main()
