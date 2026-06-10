#!/usr/bin/env python3
import asyncio
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def main():
    print("🚀 Starting Minimal MCP Client")
    server_params = StdioServerParameters(
        command="python",
        args=["-m", "mcp_server.main"],
    )
    try:
        async with stdio_client(server_params) as (read, write):
            async with ClientSession(read, write) as session:
                print("✅ Connected to MCP server")
                await session.initialize()
                print("✅ Session initialized")
                
                print("⏳ Listing resources...")
                resources = await session.list_resources()
                print("✅ list_resources() call completed.")
                print("Raw response:")
                print(resources)

    except Exception as e:
        print(f"❌ Error during minimal client execution: {e}")

if __name__ == "__main__":
    asyncio.run(main())
