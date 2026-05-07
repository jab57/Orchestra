"""
Orchestra MCP Client

Manages subprocess lifecycle and MCP stdio transport for CASCADE and RegNetAgents
child servers. Each MCPClient instance is a persistent connection to one child server.

Usage:
    async with make_cascade_client() as cascade:
        result = await cascade.call_tool("get_gene_metadata", {"gene": "TP53", "cell_type": "epithelial_cell"})
"""

import json
import os
import sys
from contextlib import AsyncExitStack
from datetime import timedelta
from typing import Any

from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client


# Per-tool timeouts (seconds) matching the plan's latency budget
TIMEOUT_PERTURBATION = 60.0   # comprehensive_perturbation_analysis
TIMEOUT_PPI = 15.0            # get_protein_interactions
TIMEOUT_NETWORK = 60.0        # comprehensive_gene_analysis, pathway_focused_analysis
TIMEOUT_DEFAULT = 30.0


class MCPToolError(Exception):
    """Raised when a child server tool call returns isError=True."""


class MCPClient:
    """
    Persistent MCP connection to one child server subprocess via stdio transport.

    Spawns the child server as a subprocess on __aenter__, calls initialize(),
    and tears down cleanly on __aexit__. Keep one instance open per workflow run.
    """

    def __init__(
        self,
        server_name: str,
        command: str,
        args: list[str],
        cwd: str | None = None,
        env: dict[str, str] | None = None,
    ):
        self.server_name = server_name
        self._params = StdioServerParameters(
            command=command,
            args=args,
            cwd=cwd,
            env=env,
        )
        self._exit_stack: AsyncExitStack | None = None
        self._session: ClientSession | None = None

    async def __aenter__(self) -> "MCPClient":
        self._exit_stack = AsyncExitStack()
        await self._exit_stack.__aenter__()
        read, write = await self._exit_stack.enter_async_context(stdio_client(self._params))
        self._session = await self._exit_stack.enter_async_context(ClientSession(read, write))
        await self._session.initialize()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self._exit_stack is not None:
            result = await self._exit_stack.__aexit__(exc_type, exc_val, exc_tb)
            self._exit_stack = None
            self._session = None
            return result

    async def call_tool(
        self,
        name: str,
        arguments: dict[str, Any],
        timeout_seconds: float = TIMEOUT_DEFAULT,
    ) -> dict[str, Any]:
        """
        Call a tool on the child server and return its result as a dict.

        Raises MCPToolError if the server returns isError=True.
        Raises TimeoutError (from anyio) if the call exceeds timeout_seconds.
        """
        if self._session is None:
            raise RuntimeError(
                f"MCPClient '{self.server_name}' not connected — use 'async with'"
            )

        result = await self._session.call_tool(
            name,
            arguments,
            read_timeout_seconds=timedelta(seconds=timeout_seconds),
        )

        if result.isError:
            raise MCPToolError(
                f"{self.server_name}.{name} returned error: {_extract_text(result.content)}"
            )

        # Prefer structured content if available; fall back to JSON-parsing text
        # structuredContent added in mcp >1.9; use getattr for compatibility
        structured = getattr(result, "structuredContent", None)
        if structured:
            return structured

        text = _extract_text(result.content)
        try:
            return json.loads(text)
        except (json.JSONDecodeError, TypeError):
            return {"text": text}

    async def health_check(self) -> list[str]:
        """
        Verify the connection is live by listing available tools.
        Returns tool names. Raises if the connection is unhealthy.
        """
        if self._session is None:
            raise RuntimeError(
                f"MCPClient '{self.server_name}' not connected — use 'async with'"
            )
        result = await self._session.list_tools()
        return [t.name for t in result.tools]


def _extract_text(content: list) -> str:
    return "\n".join(c.text for c in content if hasattr(c, "text"))


# ---------------------------------------------------------------------------
# Pre-configured clients for the two child servers
# ---------------------------------------------------------------------------

def _venv_python(server_dir: str) -> str:
    """Return the venv Python for a child server, falling back to sys.executable."""
    candidate = os.path.join(server_dir, "env", "Scripts", "python.exe")
    return candidate if os.path.isfile(candidate) else sys.executable


def _child_env(extra: dict[str, str] | None = None) -> dict[str, str]:
    """Build subprocess env: inherit current environment, optionally propagate SSL bypass."""
    env = dict(os.environ)
    if os.environ.get("ORCHESTRA_SSL_NO_VERIFY") == "1":
        env["ORCHESTRA_SSL_NO_VERIFY"] = "1"
    if extra:
        env.update(extra)
    return env


def make_cascade_client(cwd: str = r"c:\Dev\CASCADE") -> MCPClient:
    """Returns a pre-configured MCPClient for the CASCADE child server."""
    return MCPClient(
        server_name="CASCADE",
        command=_venv_python(cwd),
        args=["cascade_langgraph_mcp_server.py"],
        cwd=cwd,
        env=_child_env(),
    )


def make_regnetagents_client(cwd: str = r"c:\Dev\RegNetAgents") -> MCPClient:
    """Returns a pre-configured MCPClient for the RegNetAgents child server."""
    return MCPClient(
        server_name="RegNetAgents",
        command=_venv_python(cwd),
        args=["regnetagents_langgraph_mcp_server.py"],
        cwd=cwd,
    )
