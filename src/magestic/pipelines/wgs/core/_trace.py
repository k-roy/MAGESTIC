"""
magestic.pipelines.wgs.core._trace — the single import site for the pinned TRACE.

Every TRACE symbol used by `magestic.pipelines.wgs` is reached through this module,
for three reasons:

1.  **One-directional, pinned dependency** (ARCHITECTURE_POLICY §11). TRACE owns the
    algorithms and never imports magestic. Importing it in exactly one place makes the
    dependency auditable — `grep -rn "from . import _trace" core/` lists every consumer.
2.  **`import magestic.pipelines.wgs.core` must work without the extra installed.** The
    import is deferred to call time, so ruff, mypy and pytest collection do not require
    `magestic[wgs]`; a caller that actually needs an algorithm gets a loud ImportError
    naming the install command instead of a `ModuleNotFoundError` from three frames down.
3.  **Which TRACE ran is a recorded fact, not an assumption** — see `trace_provenance()`.
    More than one TRACE working tree exists on this program's machines and they are not
    byte-identical, so "it imported" never establishes which algorithm executed.

Never invoke TRACE as a SUBPROCESS. Its per-module `main()` entry points exist and are
reachable as runnable modules, but shelling out to one would re-create the inverted
dependency through the back door (ARCHITECTURE_POLICY §11). Call the library functions.

(The forbidden invocation string is deliberately not spelled out here: the repo guard
`test_wgs_never_shells_out_to_a_trace_cli` greps this package's source text for it, and
a literal example inside a prohibition would trip the guard it documents — a keyword
grep answers "is this word present", never "is this capability present".)
"""

from __future__ import annotations

import importlib
import hashlib
import os
import subprocess
from pathlib import Path
from types import ModuleType

__all__ = [
    "bc1_consensus",
    "donor_outcome",
    "reconcile",
    "sv_plasmid",
    "trace_provenance",
]

_EXTRA_HINT = (
    "TRACE (trace-crispr) is required by magestic.pipelines.wgs but is not importable. "
    "The WGS module is orchestration over a PINNED TRACE (ARCHITECTURE_POLICY §11); "
    "install it with:  pip install 'magestic[wgs]'"
)


def _module(name: str) -> ModuleType:
    """Import one `trace_crispr.wgs` submodule, or raise with the install hint."""
    try:
        return importlib.import_module(f"trace_crispr.wgs.{name}")
    except ImportError as exc:  # pragma: no cover - environment-dependent
        raise ImportError(f"{_EXTRA_HINT}  (importing trace_crispr.wgs.{name})") from exc


def bc1_consensus() -> ModuleType:
    """`trace_crispr.wgs.bc1_consensus` — Stage 1 bc1 calling + bc1->oligo bridge."""
    return _module("bc1_consensus")


def donor_outcome() -> ModuleType:
    """`trace_crispr.wgs.donor_outcome` — Stage 2 donor-window outcome pre-calling."""
    return _module("donor_outcome")


def reconcile() -> ModuleType:
    """`trace_crispr.wgs.reconcile` — Stage 4 per-cell lineage x edit-VAF reconciliation."""
    return _module("reconcile")


def sv_plasmid() -> ModuleType:
    """`trace_crispr.wgs.sv_plasmid` — Stage 3 plasmid integration + donor-read cleaning."""
    return _module("sv_plasmid")


def trace_provenance() -> dict[str, str]:
    """Record which TRACE actually loaded, for the run's `PROVENANCE.json`.

    Returns
    -------
    dict
        `trace_version` (`trace_crispr.__version__`), `trace_path` (the imported
        package file), and `trace_wgs_api_names` (the length of
        `trace_crispr.wgs.__all__`, as a string).

    Notes
    -----
    `trace_wgs_api_names` is the cheapest discriminator between TRACE working trees:
    the pinned commit exports 60 names, an adjacent development tree exports 63. An
    import that succeeds proves only that *a* TRACE is present — never which one — so
    this triple belongs in every provenance sidecar the module writes.

    `trace_crispr.__version__` is known to disagree with the packaging metadata
    (`0.6.1` in-package vs `0.6.3` in `pyproject.toml`/dist-info). Both are recorded as
    found; reconciling them is a TRACE-side fix, not a magestic one.
    """
    trace = importlib.import_module("trace_crispr")
    wgs = importlib.import_module("trace_crispr.wgs")
    package = Path(trace.__file__).resolve().parent
    source_hash = hashlib.sha256()
    for source in sorted((package/'wgs').rglob('*.py')):
        source_hash.update(str(source.relative_to(package)).encode())
        source_hash.update(source.read_bytes())
    def git(*args):
        try:
            return subprocess.check_output(['git', *args], cwd=package,
                                           stderr=subprocess.DEVNULL, text=True).strip()
        except (OSError, subprocess.CalledProcessError):
            return ''
    commit = git('rev-parse', 'HEAD')
    staged = os.environ.get('TRACE_BUILD_SHA', '').strip()
    if not commit and not staged:
        for root in (package, package.parent):
            identity_file = root/'BUILD_SHA'
            if identity_file.is_file():
                staged = identity_file.read_text().strip()
                break
    return {
        "trace_version": str(getattr(trace, "__version__", "unknown")),
        "trace_path": str(getattr(trace, "__file__", "unknown")),
        "trace_wgs_api_names": str(len(getattr(wgs, "__all__", ()))),
        'trace_git_commit': commit or staged or 'unresolved',
        'trace_commit_source': 'git' if commit else 'staged' if staged else 'unresolved',
        'trace_worktree_dirty': str(bool(git('status', '--porcelain'))) if commit else 'unknown',
        'trace_wgs_source_sha256': source_hash.hexdigest(),
    }
