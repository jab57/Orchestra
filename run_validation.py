"""
Issue #6 Biological Validation Runner

Runs the three canonical JOSS validation cases against live child servers
and writes results to outputs/validation_<case>_<timestamp>.txt

Usage:
    python run_validation.py                  # run all three cases
    python run_validation.py apc              # run only APC case
    python run_validation.py tp53
    python run_validation.py brd4
"""

import asyncio
import io
import sys
import json
from datetime import datetime
from pathlib import Path

# Force UTF-8 stdout so unicode chars (arrows, checkmarks) don't crash on cp1252 terminals
if hasattr(sys.stdout, "buffer"):
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# Load .env before importing workflow (sets ORCHESTRA_SSL_NO_VERIFY if present)
from dotenv import load_dotenv
load_dotenv(Path(__file__).parent / ".env")

from orchestra_langgraph_workflow import OrchestraWorkflow

CASES = {
    "apc": {
        "gene": "APC",
        "cell_type": "epithelial_cell",
        "analysis_type": "causal_chain",
        "description": "APC / epithelial_cell — effector path (strict-necessity demonstration)",
    },
    "tp53": {
        "gene": "TP53",
        "cell_type": "epithelial_cell",
        "analysis_type": "causal_chain",
        "description": "TP53 / epithelial_cell — TF path (sanity check)",
    },
    "brd4": {
        "gene": "MYC",
        "cell_type": "cd4_t_cells",
        "analysis_type": "therapeutic_validation",
        "description": "MYC / cd4_t_cells — validate_therapeutic_targets (BRD4→MYC complementarity)",
    },
}

OUTPUT_DIR = Path(__file__).parent / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)


async def run_case(_name: str, cfg: dict) -> dict:
    print(f"\n{'='*60}")
    print(f"Running: {cfg['description']}")
    print(f"{'='*60}")
    workflow = OrchestraWorkflow()
    result = await workflow.run_analysis(
        gene=cfg["gene"],
        cell_type=cfg["cell_type"],
        analysis_type=cfg["analysis_type"],
    )
    return result


def save_result(name: str, cfg: dict, result: dict) -> Path:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = OUTPUT_DIR / f"validation_{name}_{ts}.txt"

    lines = [
        f"# Orchestra Issue #6 Validation — {name.upper()}",
        f"# {cfg['description']}",
        f"# Run: {datetime.now().isoformat()}",
        "",
        "## Final Report",
        result.get("final_report", "(no final_report in result)"),
        "",
        "## Synthesis (raw)",
        json.dumps(result.get("synthesis", {}), indent=2, default=str),
        "",
        "## PPI Interactions (raw)",
        json.dumps(result.get("ppi_interactions", {}), indent=2, default=str),
        "",
        "## Network Analysis (raw, truncated)",
        json.dumps(result.get("network_analysis", {}), indent=2, default=str)[:2000],
        "",
        "## Perturbation Result (raw, truncated)",
        json.dumps(result.get("perturbation_result", {}), indent=2, default=str)[:2000],
        "",
        "## Validated Targets (raw)",
        json.dumps(result.get("validated_targets", []), indent=2, default=str),
        "",
        "## Errors",
        json.dumps(result.get("errors", {}), indent=2, default=str),
    ]

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nSaved: {out_path}")
    return out_path


async def main():
    requested = sys.argv[1].lower() if len(sys.argv) > 1 else "all"
    cases_to_run = (
        list(CASES.items())
        if requested == "all"
        else [(requested, CASES[requested])]
        if requested in CASES
        else []
    )

    if not cases_to_run:
        print(f"Unknown case '{requested}'. Choose from: {', '.join(CASES)} or 'all'")
        sys.exit(1)

    for name, cfg in cases_to_run:
        try:
            result = await run_case(name, cfg)
            save_result(name, cfg, result)
            report = result.get("final_report", "")
            print("\n--- REPORT PREVIEW (first 80 lines) ---")
            for line in report.splitlines()[:80]:
                print(line.encode("utf-8", errors="replace").decode("utf-8"))
        except Exception as e:
            print(f"\nERROR running {name}: {e}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    asyncio.run(main())
