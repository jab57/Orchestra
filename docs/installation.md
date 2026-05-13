# Installation

## Prerequisites

- Python 3.10 or later
- [RegNetAgents](https://github.com/jab57/RegNetAgents) installed and working in its own directory
- [CASCADE](https://github.com/jab57/CASCADE) installed and working in its own directory
- Both child servers must have their own virtual environments (`env/`) with all dependencies installed

Orchestra spawns RegNetAgents and CASCADE as subprocesses via the MCP stdio transport. It uses each child server's own `env/Scripts/python.exe` (Windows) or `env/bin/python` (Unix), so the child servers' Python environments are independent of Orchestra's.

## Install Orchestra

```bash
git clone https://github.com/jab57/Orchestra.git
cd Orchestra
python -m venv env
# Windows:
env\Scripts\activate
# macOS/Linux:
source env/bin/activate
pip install -r requirements.txt
```

## Configure Environment

```bash
cp .env.example .env
```

Open `.env` and set the paths to your child server installations:

```env
REGNETAGENTS_SERVER_PATH=C:\Dev\RegNetAgents\regnetagents_langgraph_mcp_server.py
CASCADE_SERVER_PATH=C:\Dev\CASCADE\cascade_langgraph_mcp_server.py
```

The server path variables are informational — Orchestra derives the working directory from them. The actual paths used internally are hardcoded defaults (`c:\Dev\RegNetAgents` and `c:\Dev\CASCADE`). If your installations are in non-default locations, edit `mcp_client.py`:

```python
def make_cascade_client(cwd: str = r"c:\Dev\CASCADE") -> MCPClient:
def make_regnetagents_client(cwd: str = r"c:\Dev\RegNetAgents") -> MCPClient:
```

## Verify Installation

Run the unit tests (no child servers required):

```bash
pytest tests/
```

All 108 unit tests should pass. The 17 integration tests are skipped by default.

To run integration tests (requires running RegNetAgents and CASCADE):

```bash
set ORCHESTRA_INTEGRATION_TESTS=1   # Windows
# or: export ORCHESTRA_INTEGRATION_TESTS=1  (macOS/Linux)
pytest tests/
```

## Run a Validation Case (Optional)

Verify end-to-end with a live run:

```bash
python run_validation.py apc
```

Expected output: an APC→CTNNB1→Wnt signaling analysis report saved to `outputs/`.

## Connect to Claude Desktop (Optional)

To expose Orchestra as an MCP server to Claude Desktop, add it to your Claude Desktop configuration file:

**macOS:** `~/Library/Application Support/Claude/claude_desktop_config.json`  
**Windows:** `%APPDATA%\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "orchestra": {
      "command": "C:\\Dev\\Orchestra\\env\\Scripts\\python.exe",
      "args": ["C:\\Dev\\Orchestra\\orchestra_mcp_server.py"],
      "env": {
        "PYTHONPATH": "C:\\Dev\\Orchestra"
      }
    }
  }
}
```

Adjust paths for your installation. After restarting Claude Desktop, the three Orchestra tools (`causal_chain_analysis`, `validate_therapeutic_targets`, `effector_analysis`) will appear in the tools list.

## Troubleshooting

### SSL verification errors

On networks with corporate SSL inspection, connections to external APIs (STRING, GitHub) may fail with SSL certificate errors. Add this to your `.env`:

```env
ORCHESTRA_SSL_NO_VERIFY=1
```

This propagates the bypass to child server subprocesses. Note: only use this on internal/trusted networks.

### Child server not found

If Orchestra cannot start a child server, verify:
1. The directory exists: `c:\Dev\RegNetAgents` and `c:\Dev\CASCADE`
2. Each has a virtual environment: `env\Scripts\python.exe` (Windows)
3. Each child server runs independently: `python regnetagents_langgraph_mcp_server.py` should start without errors

### NumPy version conflict

Orchestra uses each child server's own Python environment (`env/Scripts/python.exe`) to avoid NumPy 2.x vs 1.x conflicts between projects. If you see import errors, confirm that `env/Scripts/python.exe` exists in both child server directories.

### Subprocess timeout

Default timeouts: 60s for perturbation analysis, 15s for PPI, 60s for network analysis. If you see `TimeoutError` on slow hardware, increase them in `.env`:

```env
CASCADE_TIMEOUT=60
REGNETAGENTS_TIMEOUT=60
```
