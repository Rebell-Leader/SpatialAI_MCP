#!/usr/bin/env python3
"""Simple integration test for the FastMCP server."""

import asyncio
import sys
import logging
from pathlib import Path

# Add src to path for testing
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from openproblems_mcp.server import MCPServer
from openproblems_mcp.config import ConfigManager


async def test_server_creation():
    """Test that we can create and initialize the server."""
    print("🧪 Testing server creation...")

    try:
        # Load default configuration
        config_manager = ConfigManager()
        config = config_manager.load_config()

        print(f"✅ Configuration loaded successfully")
        print(f"   Workspace: {config.server.workspace_root}")
        print(f"   Log level: {config.server.log_level}")

        # Create server
        server = MCPServer(config)
        print(f"✅ Server created successfully")

        # Test tool detection
        health_status = server.tool_detector.get_health_status()
        print(f"✅ Tool detection completed")
        print(f"   Available tools: {health_status['tools_available']}/{health_status['tools_total']}")

        for tool_name, tool_info in health_status['tools'].items():
            status = "✅" if tool_info['available'] else "❌"
            version = f" (v{tool_info['version']})" if tool_info['version'] else ""
            print(f"   {status} {tool_name}{version}")

        # Get FastMCP server instance
        mcp_server = server.get_mcp_server()
        print(f"✅ FastMCP server instance obtained")

        print("\n🎉 All tests passed! Server is ready for production use.")
        return True

    except Exception as e:
        print(f"❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


async def main():
    """Run integration tests."""
    print("🚀 OpenProblems MCP Server Integration Test")
    print("=" * 50)

    success = await test_server_creation()

    if success:
        print("\n✅ Integration test completed successfully!")
        sys.exit(0)
    else:
        print("\n❌ Integration test failed!")
        sys.exit(1)


if __name__ == "__main__":
    # Set up basic logging
    logging.basicConfig(level=logging.INFO)

    # Run the test
    asyncio.run(main())
