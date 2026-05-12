"""
Unit tests for MCPClient — mock MCP server, verify client round-trip.

No running child servers required. The MCP stdio transport and ClientSession are
mocked end-to-end, verifying that MCPClient correctly:
  - initializes the session on __aenter__ and clears it on __aexit__
  - returns structuredContent directly when present
  - falls back to JSON-parsing text content when structuredContent is absent
  - wraps plain text in {"text": ...} when the text is not valid JSON
  - raises MCPToolError when the server returns isError=True
  - raises RuntimeError when call_tool / health_check are called outside async with
"""

import json
import pytest
from contextlib import ExitStack
from unittest.mock import AsyncMock, MagicMock, patch

from mcp_client import MCPClient, MCPToolError, make_cascade_client, make_regnetagents_client


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class _AsyncCM:
    """Minimal async context manager returning a fixed value on __aenter__."""
    def __init__(self, value):
        self._value = value

    async def __aenter__(self):
        return self._value

    async def __aexit__(self, *args):
        return False


def _make_result(is_error: bool = False, structured=None, text: str | None = None):
    """Build a mock MCP tool-call result object."""
    result = MagicMock()
    result.isError = is_error
    # Explicitly set so getattr(result, "structuredContent", None) returns the right value.
    result.structuredContent = structured
    if text is not None:
        item = MagicMock()
        item.text = text
        result.content = [item]
    else:
        result.content = []
    return result


def _make_session(*, tool_result=None, tool_names: list[str] | None = None):
    """Return a mock ClientSession with configurable call_tool / list_tools."""
    session = MagicMock()
    session.initialize = AsyncMock()

    if tool_result is not None:
        session.call_tool = AsyncMock(return_value=tool_result)

    if tool_names is not None:
        tool_mocks = []
        for name in tool_names:
            t = MagicMock()
            t.name = name
            tool_mocks.append(t)
        list_result = MagicMock()
        list_result.tools = tool_mocks
        session.list_tools = AsyncMock(return_value=list_result)

    return session


def _enter_patches(stack: ExitStack, session):
    """Enter stdio_client and ClientSession patches via an ExitStack."""
    mock_read, mock_write = MagicMock(), MagicMock()
    stack.enter_context(
        patch("mcp_client.stdio_client", return_value=_AsyncCM((mock_read, mock_write)))
    )
    stack.enter_context(
        patch("mcp_client.ClientSession", return_value=_AsyncCM(session))
    )


# ---------------------------------------------------------------------------
# Connection lifecycle
# ---------------------------------------------------------------------------

class TestMCPClientLifecycle:
    async def test_enter_calls_initialize(self):
        session = _make_session()
        with ExitStack() as stack:
            _enter_patches(stack, session)
            async with MCPClient("Test", "python", ["server.py"]):
                session.initialize.assert_called_once()

    async def test_exit_clears_session_and_stack(self):
        session = _make_session()
        with ExitStack() as stack:
            _enter_patches(stack, session)
            client = MCPClient("Test", "python", ["server.py"])
            async with client:
                pass
        assert client._session is None
        assert client._exit_stack is None

    async def test_call_tool_outside_context_raises_runtime_error(self):
        client = MCPClient("Test", "python", ["server.py"])
        with pytest.raises(RuntimeError, match="not connected"):
            await client.call_tool("some_tool", {})

    async def test_health_check_outside_context_raises_runtime_error(self):
        client = MCPClient("Test", "python", ["server.py"])
        with pytest.raises(RuntimeError, match="not connected"):
            await client.health_check()


# ---------------------------------------------------------------------------
# call_tool — result routing (structuredContent → JSON → plain text)
# ---------------------------------------------------------------------------

class TestCallToolResultRouting:
    async def test_structured_content_returned_directly(self):
        expected = {"gene": "TP53", "gene_type": "master_regulator", "num_targets": 200}
        session = _make_session(tool_result=_make_result(structured=expected))
        with ExitStack() as stack:
            _enter_patches(stack, session)
            async with MCPClient("Test", "python", ["server.py"]) as client:
                result = await client.call_tool("get_gene_metadata", {"gene": "TP53"})
        assert result == expected

    async def test_json_text_fallback_when_no_structured_content(self):
        payload = {"interactions": [{"partner": "CTNNB1", "combined_score": 0.999}]}
        session = _make_session(
            tool_result=_make_result(structured=None, text=json.dumps(payload))
        )
        with ExitStack() as stack:
            _enter_patches(stack, session)
            async with MCPClient("Test", "python", ["server.py"]) as client:
                result = await client.call_tool("get_protein_interactions", {"gene": "APC"})
        assert result == payload

    async def test_plain_text_wrapped_in_dict(self):
        session = _make_session(
            tool_result=_make_result(structured=None, text="plain text response")
        )
        with ExitStack() as stack:
            _enter_patches(stack, session)
            async with MCPClient("Test", "python", ["server.py"]) as client:
                result = await client.call_tool("some_tool", {})
        assert result == {"text": "plain text response"}

    async def test_error_result_raises_mcp_tool_error(self):
        session = _make_session(
            tool_result=_make_result(is_error=True, text="tool execution failed")
        )
        with ExitStack() as stack:
            _enter_patches(stack, session)
            async with MCPClient("Test", "python", ["server.py"]) as client:
                with pytest.raises(MCPToolError):
                    await client.call_tool("failing_tool", {})

    async def test_error_message_contains_server_name(self):
        session = _make_session(
            tool_result=_make_result(is_error=True, text="bad args")
        )
        with ExitStack() as stack:
            _enter_patches(stack, session)
            async with MCPClient("CASCADE", "python", ["server.py"]) as client:
                with pytest.raises(MCPToolError, match="CASCADE"):
                    await client.call_tool("get_gene_metadata", {})


# ---------------------------------------------------------------------------
# health_check
# ---------------------------------------------------------------------------

class TestHealthCheck:
    async def test_returns_tool_names(self):
        names = ["get_gene_metadata", "comprehensive_perturbation_analysis", "get_protein_interactions"]
        session = _make_session(tool_names=names)
        with ExitStack() as stack:
            _enter_patches(stack, session)
            async with MCPClient("CASCADE", "python", ["server.py"]) as client:
                tools = await client.health_check()
        assert set(names) == set(tools)

    async def test_returns_list_of_strings(self):
        session = _make_session(tool_names=["tool_a", "tool_b"])
        with ExitStack() as stack:
            _enter_patches(stack, session)
            async with MCPClient("Test", "python", ["server.py"]) as client:
                tools = await client.health_check()
        assert isinstance(tools, list)
        assert all(isinstance(t, str) for t in tools)


# ---------------------------------------------------------------------------
# Factory functions
# ---------------------------------------------------------------------------

class TestFactoryFunctions:
    def test_make_cascade_client_returns_mcp_client(self):
        assert isinstance(make_cascade_client(), MCPClient)

    def test_make_cascade_client_server_name(self):
        assert make_cascade_client().server_name == "CASCADE"

    def test_make_cascade_client_correct_script(self):
        client = make_cascade_client()
        assert "cascade_langgraph_mcp_server.py" in client._params.args

    def test_make_regnetagents_client_returns_mcp_client(self):
        assert isinstance(make_regnetagents_client(), MCPClient)

    def test_make_regnetagents_client_server_name(self):
        assert make_regnetagents_client().server_name == "RegNetAgents"

    def test_make_regnetagents_client_correct_script(self):
        client = make_regnetagents_client()
        assert "regnetagents_langgraph_mcp_server.py" in client._params.args
