"""Main entry point for the OpenProblems MCP Server."""

import asyncio
import logging
import sys
from typing import Optional

from .config import ConfigManager, setup_logging
from .server import MCPServer
from .exceptions import ConfigurationError, DependencyError


logger = logging.getLogger(__name__)


async def main_async(config_path: Optional[str] = None) -> None:
    """Main async entry point for the MCP server."""
    try:
        # Load configuration
        config_manager = ConfigManager(config_path)
        config = config_manager.load_config()

        # Set up logging
        setup_logging(config)

        logger.info("Starting OpenProblems Spatial Transcriptomics MCP Server")
        logger.info(f"Workspace root: {config.server.workspace_root}")
        logger.info(f"Log level: {config.server.log_level}")

        # Create and run server
        server = MCPServer(config)
        await server.run()

    except ConfigurationError as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        if hasattr(e, 'suggested_fixes') and e.suggested_fixes:
            print("Suggested fixes:", file=sys.stderr)
            for fix in e.suggested_fixes:
                print(f"  • {fix}", file=sys.stderr)
        sys.exit(1)

    except DependencyError as e:
        print(f"Dependency error: {e}", file=sys.stderr)
        if hasattr(e, 'suggested_fixes') and e.suggested_fixes:
            print("Installation suggestions:", file=sys.stderr)
            for fix in e.suggested_fixes:
                print(f"  • {fix}", file=sys.stderr)
        sys.exit(1)

    except KeyboardInterrupt:
        logger.info("Server shutdown requested")
        sys.exit(0)

    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


def main() -> None:
    """Main entry point for the MCP server."""
    # Handle command line arguments if needed
    config_path = None
    if len(sys.argv) > 1:
        if sys.argv[1] in ['-h', '--help']:
            print("OpenProblems Spatial Transcriptomics MCP Server")
            print("Usage: openproblems-mcp-server [config_path]")
            print("")
            print("Arguments:")
            print("  config_path    Optional path to configuration file")
            print("")
            print("Environment Variables:")
            print("  OPENPROBLEMS_MCP_WORKSPACE_ROOT    Workspace root directory")
            print("  OPENPROBLEMS_MCP_MAX_CONCURRENT    Max concurrent executions")
            print("  OPENPROBLEMS_MCP_TIMEOUT           Default timeout in seconds")
            print("  OPENPROBLEMS_MCP_LOG_LEVEL         Log level (DEBUG, INFO, WARNING, ERROR)")
            print("  OPENPROBLEMS_MCP_*_EXECUTABLE      Tool executable paths")
            sys.exit(0)
        else:
            config_path = sys.argv[1]

    # Run the async main function
    try:
        asyncio.run(main_async(config_path))
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
