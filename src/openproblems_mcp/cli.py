"""Command-line interface for the OpenProblems MCP Server."""

import click
import asyncio
import sys
import os
from pathlib import Path
from typing import Optional

from .config import ConfigManager, setup_logging
from .server import MCPServer
from .tool_detection import ToolDetector
from .exceptions import ConfigurationError, DependencyError


@click.group()
@click.version_option(version="0.1.0", prog_name="openproblems-mcp")
def cli():
    """OpenProblems Spatial Transcriptomics MCP Server CLI."""
    pass


@cli.command()
@click.option(
    "--config", "-c",
    type=click.Path(exists=True),
    help="Path to configuration file"
)
def serve(config: Optional[str]):
    """Start the MCP server."""
    from .main import main_async

    try:
        asyncio.run(main_async(config))
    except KeyboardInterrupt:
        click.echo("Server stopped.")


@cli.command()
@click.option(
    "--config", "-c",
    type=click.Path(exists=True),
    help="Path to configuration file"
)
def check(config: Optional[str]):
    """Check server health and tool availability."""
    try:
        # Load configuration
        config_manager = ConfigManager(config)
        server_config = config_manager.load_config()

        # Set up basic logging
        setup_logging(server_config)

        # Check tools
        tool_detector = ToolDetector(server_config.tools)
        health_status = tool_detector.get_health_status()

        click.echo("🔧 OpenProblems MCP Server Health Check")
        click.echo("=" * 50)

        # Overall status
        if health_status['all_tools_available']:
            click.echo("✅ All tools are available and ready!")
        else:
            click.echo(f"⚠️  {health_status['tools_available']}/{health_status['tools_total']} tools available")

        click.echo()

        # Tool details
        for tool_name, tool_info in health_status['tools'].items():
            status = "✅" if tool_info['available'] else "❌"
            version = f" (v{tool_info['version']})" if tool_info['version'] else ""
            click.echo(f"{status} {tool_name.title()}{version}")

            if tool_info['error']:
                click.echo(f"   Error: {tool_info['error']}")

        # Installation suggestions
        missing_tools = tool_detector.get_missing_tools()
        if missing_tools:
            click.echo()
            click.echo("📋 Installation Suggestions:")
            suggestions = tool_detector._get_installation_suggestions([t.name for t in missing_tools])
            for suggestion in suggestions:
                click.echo(f"  • {suggestion}")

        # Configuration info
        click.echo()
        click.echo("⚙️  Configuration:")
        click.echo(f"  Workspace: {server_config.server.workspace_root}")
        click.echo(f"  Max concurrent: {server_config.server.max_concurrent_executions}")
        click.echo(f"  Timeout: {server_config.server.default_timeout_seconds}s")
        click.echo(f"  Log level: {server_config.server.log_level}")

        # Exit with error code if tools are missing
        if not health_status['all_tools_available']:
            sys.exit(1)

    except ConfigurationError as e:
        click.echo(f"❌ Configuration error: {e}", err=True)
        if hasattr(e, 'suggested_fixes') and e.suggested_fixes:
            click.echo("Suggested fixes:", err=True)
            for fix in e.suggested_fixes:
                click.echo(f"  • {fix}", err=True)
        sys.exit(1)

    except Exception as e:
        click.echo(f"❌ Health check failed: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--output", "-o",
    type=click.Path(),
    default="~/.openproblems-mcp/config.yaml",
    help="Output path for configuration file"
)
@click.option(
    "--force", "-f",
    is_flag=True,
    help="Overwrite existing configuration file"
)
def init(output: str, force: bool):
    """Initialize a default configuration file."""
    output_path = Path(output).expanduser()

    if output_path.exists() and not force:
        click.echo(f"❌ Configuration file already exists: {output_path}")
        click.echo("Use --force to overwrite")
        sys.exit(1)

    try:
        config_manager = ConfigManager()
        config_manager.create_default_config_file(str(output_path))

        click.echo(f"✅ Created default configuration file: {output_path}")
        click.echo()
        click.echo("You can now:")
        click.echo(f"  1. Edit the configuration: {output_path}")
        click.echo("  2. Check server health: openproblems-mcp check")
        click.echo("  3. Start the server: openproblems-mcp serve")

    except Exception as e:
        click.echo(f"❌ Failed to create configuration file: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option(
    "--config", "-c",
    type=click.Path(exists=True),
    help="Path to configuration file"
)
def config_info(config: Optional[str]):
    """Show current configuration."""
    try:
        config_manager = ConfigManager(config)
        server_config = config_manager.load_config()

        click.echo("⚙️  OpenProblems MCP Server Configuration")
        click.echo("=" * 50)

        click.echo("Server Settings:")
        click.echo(f"  Workspace Root: {server_config.server.workspace_root}")
        click.echo(f"  Max Concurrent Executions: {server_config.server.max_concurrent_executions}")
        click.echo(f"  Default Timeout: {server_config.server.default_timeout_seconds}s")
        click.echo(f"  Log Level: {server_config.server.log_level}")
        click.echo(f"  Memory Limit per Execution: {server_config.server.max_memory_per_execution_mb}MB")
        click.echo(f"  Log Retention: {server_config.server.log_retention_days} days")

        click.echo()
        click.echo("Tool Executables:")
        click.echo(f"  Nextflow: {server_config.tools.nextflow_executable}")
        click.echo(f"  Viash: {server_config.tools.viash_executable}")
        click.echo(f"  Docker: {server_config.tools.docker_executable}")
        click.echo(f"  Git: {server_config.tools.git_executable}")
        click.echo(f"  Python: {server_config.tools.python_executable}")

        click.echo()
        click.echo("File System:")
        click.echo(f"  Allowed Extensions: {', '.join(server_config.server.allowed_file_extensions)}")
        if server_config.server.blocked_paths:
            click.echo(f"  Blocked Paths: {', '.join(server_config.server.blocked_paths)}")
        else:
            click.echo("  Blocked Paths: None")

    except ConfigurationError as e:
        click.echo(f"❌ Configuration error: {e}", err=True)
        sys.exit(1)

    except Exception as e:
        click.echo(f"❌ Failed to load configuration: {e}", err=True)
        sys.exit(1)


def main():
    """Main CLI entry point."""
    cli()


if __name__ == "__main__":
    main()
