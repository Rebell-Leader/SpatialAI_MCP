import asyncio
import json
import logging
from aiohttp import web
from mcp.server import Server
from mcp.server.models import InitializationOptions
from mcp.server.stdio import StdioStream
from mcp.types import TextContent

logger = logging.getLogger(__name__)

class HttpServer:
    def __init__(self, mcp_server: Server, server_name: str, server_version: str):
        self.mcp_server = mcp_server
        self.server_name = server_name
        self.server_version = server_version
        self.app = web.Application()
        self.app.router.add_post('/mcp', self.handle_mcp_request)
        self.app.router.add_get('/health', self.handle_health_check)

    async def handle_mcp_request(self, request: web.Request) -> web.Response:
        try:
            data = await request.json()
            logger.info(f"Received MCP request: {data}")

            # Create a mock StdioStream to process the request
            # This is a simplified approach; a more robust solution might involve
            # a dedicated message queue or direct function calls if the MCP server
            # can expose its internal methods.
            mock_stdin = asyncio.StreamReader()
            mock_stdin.feed_data(json.dumps(data).encode('utf-8') + b'\n')
            mock_stdin.feed_eof()

            mock_stdout_writer = StdioStream()
            mock_stdout_reader = asyncio.StreamReader()
            mock_stdout_writer.attach_reader(mock_stdout_reader)

            # Run the MCP server's request handling logic
            # This assumes the MCP server's run method can be called with mock streams
            # and will process a single request/response cycle.
            # In a real-world scenario, you might need to adapt the MCP server
            # to expose its internal request handling as a direct function call.
            await self.mcp_server.run(
                mock_stdin,
                mock_stdout_writer,
                InitializationOptions(
                    server_name=self.server_name,
                    server_version=self.server_version,
                    capabilities={
                        "resources": {},
                        "tools": {},
                        "prompts": {},
                        "logging": {}
                    },
                ),
                single_request=True # Assuming a mode for single request processing
            )

            response_data = await mock_stdout_reader.readuntil(b'\n')
            response_json = json.loads(response_data.decode('utf-8'))
            logger.info(f"Sending MCP response: {response_json}")
            return web.json_response(response_json)

        except json.JSONDecodeError:
            return web.json_response({"error": "Invalid JSON"}, status=400)
        except Exception as e:
            logger.exception("Error handling MCP request")
            return web.json_response({"error": str(e)}, status=500)

    async def handle_health_check(self, request: web.Request) -> web.Response:
        return web.json_response({"status": "ok", "server": self.server_name, "version": self.server_version})

    async def start(self, host: str = '0.0.0.0', port: int = 8080):
        logger.info(f"Starting HTTP server on http://{host}:{port}")
        runner = web.AppRunner(self.app)
        await runner.setup()
        site = web.TCPSite(runner, host, port)
        await site.start()
        logger.info("HTTP server started.")
        # Keep the server running indefinitely
        await asyncio.Event().wait()

if __name__ == '__main__':
    # Example usage (for testing purposes)
    # In a real application, the MCP server instance would be passed from main.py
    mcp_server_instance = Server("TestMCP")
    http_server = HttpServer(mcp_server_instance, "Test MCP HTTP Server", "1.0.0")
    asyncio.run(http_server.start())
