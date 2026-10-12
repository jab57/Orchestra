"""
Provenance helper for the corroboration data-collection scripts: records which
RegNetAgents / CASCADE / Orchestra versions produced a result file, so a reader can
tell later exactly what was run.
"""
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPOS = {
    "orchestra": Path(__file__).resolve().parent.parent,
    "regnetagents": Path(r"c:\Dev\RegNetAgents"),
    "cascade": Path(r"c:\Dev\CASCADE"),
}


def _git(repo, *args):
    try:
        out = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, timeout=30)
        return out.stdout.strip() or None
    except Exception:
        return None


def collect():
    info: dict = {"collected_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    for name, repo in REPOS.items():
        info[name] = {
            "describe": _git(repo, "describe", "--tags", "--always", "--dirty"),
            "commit": _git(repo, "rev-parse", "HEAD"),
        }
    return info


def write_manifest(path, script, panels, extra=None):
    manifest = {"script": script, "panels": panels, **collect(), **(extra or {})}
    Path(path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
