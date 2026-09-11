"""Step 05's sidecar shipped ``git_head: null`` on every 2026-09-10 reference table: its helper called
``git -C <dir>``, which the cluster's system git does not support (exit 129, empty stdout, swallowed). The
helper must return the 40-hex commit of the checkout that holds the running package.
"""
import importlib.util
import re
from pathlib import Path

STEP05 = Path(__file__).resolve().parents[3] / "src/magestic/pipelines/bc1_donor_bc0/snakemake/05_assign_sublibrary.py"


def _load_step05():
    spec = importlib.util.spec_from_file_location("bdb_step05_under_test", STEP05)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_step05_git_state_carries_the_commit_sha():
    g = _load_step05()._git_state()
    assert g["git_head"] is not None, "sidecar would ship git_head: null"
    assert re.fullmatch(r"[0-9a-f]{40}", g["git_head"]), g
    assert isinstance(g["git_tracked_dirty"], bool)
    assert g["git_repo"] and Path(g["git_repo"]).is_dir()
