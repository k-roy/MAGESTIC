"""
magestic.pipelines.wgs.qc.provenance — the §8 PROVENANCE sidecars, with a real git sha.

ARCHITECTURE_POLICY §8: "Every data output carries a PROVENANCE.json." This module is a
thin adapter over `magestic.utils.provenance` (`StepProvenance`, `write_run_provenance`),
adding the one thing the shared helper cannot know: **which commit the code came from
when the code is not running out of a git checkout.**

🔴 THE NULL-SHA TRAP, and why this file exists rather than a direct `StepProvenance` call.

`utils.provenance.git_state()` shells out to `git rev-parse` with `cwd=` set to the
running package's directory. Two independent ways that yields an EMPTY sha, both of
which apply to this module's actual execution environment:

1.  **The code is staged, not cloned.** The cluster run copies `src/` to Oak and runs it
    from there. A staged tree is not a git checkout, so `git rev-parse --show-toplevel`
    fails and `git_state()` returns `{"repo": None, "commit": None}` — correctly, but the
    sidecar then carries a null sha and acceptance criterion P6 fails on a technicality
    while the run itself was fine.
2.  **Sherlock's git predates `git -C`.** A helper written with `git -C <dir>` exits 129
    with EMPTY STDOUT there, which a `returncode != 0` guard turns into `None` — a silent
    null sha rather than an error. `utils.provenance` uses `cwd=` and so is not affected,
    but the same trap is one refactor away and it has already cost this program a wrong
    record once.

So the sha is resolved from an explicit, ordered chain and **verified non-empty before
anything is written** (fail condition F7 says to check the sha in the actual sidecar, not
that the run exited 0):

    1. `$MAGESTIC_BUILD_SHA`            — what the cluster driver exports
    2. `BUILD_SHA` beside the package    — what the staging step writes into the tree
    3. `utils.provenance.git_state()`    — a real checkout (the local/dev case)

`require_sha=True` (the default for `write_step_provenance`) raises when all three come
back empty, so a run that cannot record its own provenance fails loudly rather than
shipping an unattributable table.

Author: Kevin R. Roy
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from magestic.utils.provenance import (
    RUN_SIDECAR_NAME,
    SIDECAR_SUFFIX,
    StepProvenance,
    git_state,
)
from magestic.utils.provenance import (
    write_run_provenance as _write_run_provenance,
)

__all__ = [
    "BUILD_SHA_ENV",
    "BUILD_SHA_FILE",
    "PIPELINE",
    "RUN_SIDECAR_NAME",
    "SIDECAR_SUFFIX",
    "ProvenanceError",
    "build_sha",
    "code_state",
    "preflight_provenance",
    "write_run_provenance",
    "write_step_provenance",
]

#: Pipeline name recorded in every sidecar.
PIPELINE = "wgs"

#: Environment variable the cluster driver exports with the committed sha.
BUILD_SHA_ENV = "MAGESTIC_BUILD_SHA"

#: File the staging step writes beside the `magestic` package, holding the sha.
BUILD_SHA_FILE = "BUILD_SHA"


class ProvenanceError(RuntimeError):
    """Provenance could not be recorded truthfully (typically: no resolvable git sha)."""


def _sha_from_file() -> str:
    """The sha from a `BUILD_SHA` file beside the installed `magestic` package."""
    try:
        import magestic

        pkg = Path(magestic.__file__).resolve().parent
    except Exception:  # pragma: no cover - magestic is importable by construction here
        return ""
    # Look beside the package and one and two levels up (`src/` and the repo root),
    # which covers the staged layout and an installed-into-site-packages layout.
    for candidate in (pkg, pkg.parent, pkg.parent.parent):
        path = candidate / BUILD_SHA_FILE
        try:
            if path.is_file():
                return path.read_text().strip()
        except OSError:
            continue
    return ""


def build_sha() -> str:
    """Resolve the code sha, empty string if every source is empty.

    Order: `$MAGESTIC_BUILD_SHA`, then a `BUILD_SHA` file beside the package, then
    `git_state()["commit"]` for a genuine checkout. Never guesses and never returns a
    placeholder — an empty result is reported as empty so F7 can fire.
    """
    env = os.environ.get(BUILD_SHA_ENV, "").strip()
    if env:
        return env
    from_file = _sha_from_file()
    if from_file:
        return from_file
    return (git_state().get("commit") or "").strip()


def _verify_against_checkout(sha: str, state: Mapping[str, Any]) -> tuple[bool, str]:
    """`(verified, reason)` for a sha taken from `$MAGESTIC_BUILD_SHA` or `BUILD_SHA`.

    planning/284 X-PROV-ENV-SHA: every plate-5031 receipt named `9d654b8`, the sha the
    driver exported, while candidate trees with uncommitted repairs ran. A sha is
    verified only when the package sits in a clean checkout whose HEAD IS that sha.
    A checkout at a DIFFERENT commit raises (`preflight_provenance`), because the
    receipt would name code that did not run.
    """
    head = (state.get("commit") or "").strip()
    if not head:
        return False, ("no git checkout holds the package, so the exported sha cannot be "
                       "checked; magestic_wgs_source_sha256 identifies the code that ran")
    if not (head.startswith(sha) or sha.startswith(head)):
        raise ProvenanceError(
            f"the exported code sha {sha} is not the checkout's HEAD {head} "
            f"({state.get('repo')}). A receipt would name code that did not run "
            f"(planning/284 X-PROV-ENV-SHA). Export `git rev-parse HEAD` of the tree that runs, "
            f"or unset ${BUILD_SHA_ENV}.")
    if state.get("dirty"):
        return False, "the checkout at that sha has uncommitted changes (dirty_files)"
    return True, "clean checkout at that sha"


def preflight_provenance(step: str) -> str:
    """Resolve and check the code sha BEFORE a step computes or publishes anything.

    planning/284 X-S1-PARTIAL-PUBLISH: Stage 1 published a 96-row `bc1_calls.tsv` and
    only then raised `ProvenanceError`, leaving a complete-looking table with no
    receipt. Every `run_step_N` calls this first, so an unresolvable or mismatched sha
    fails before any output exists. Returns the sha.
    """
    sha = build_sha()
    if not sha:
        raise ProvenanceError(
            f"{step}: no code sha could be resolved (checked ${BUILD_SHA_ENV}, a "
            f"{BUILD_SHA_FILE} file beside the package, and git_state()), so its outputs "
            "could not carry a receipt. Refused before computing (planning/284 "
            "X-S1-PARTIAL-PUBLISH). Export the sha in the run driver.")
    code_state()  # raises when an exported sha contradicts the checkout
    return sha


def code_state() -> dict[str, Any]:
    """`git_state()` with the resolved sha merged in, plus how the sha was found.

    Both are recorded: `commit` is what P6/F7 read, `commit_source` is what makes a
    later reader able to tell a checkout sha from a staged build sha.
    """
    import hashlib
    state = dict(git_state())
    # A commit alone does not identify dirty or staged source trees. Hash the
    # actual imported algorithm/orchestration source in a path-independent order.
    directory = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    for source in sorted(directory.rglob("*.py")):
        digest.update(str(source.relative_to(directory.parent)).encode())
        digest.update(source.read_bytes())
    state["magestic_wgs_source_sha256"] = digest.hexdigest()
    sha = build_sha()
    if sha:
        source = (
            BUILD_SHA_ENV
            if os.environ.get(BUILD_SHA_ENV, "").strip()
            else (BUILD_SHA_FILE if _sha_from_file() else "git_state")
        )
        if source == "git_state":
            verified, reason = not state.get("dirty"), (
                "clean checkout" if not state.get("dirty") else
                "the checkout has uncommitted changes (dirty_files)")
        else:
            verified, reason = _verify_against_checkout(sha, state)
        state["checkout_commit"] = state.get("commit")
        state["commit"] = sha
        state["commit_source"] = source
        state["commit_verified"] = verified
        state["commit_verification"] = reason
    else:
        state["commit_source"] = "unresolved"
    return state


def _trace_state() -> dict[str, str]:
    """Which TRACE actually loaded — recorded, never assumed.

    More than one TRACE working tree exists on this program's machines and they are not
    byte-identical, so "it imported" establishes nothing about which algorithm ran.
    Failure to import is recorded rather than raised: provenance must still be writable
    for a step that needed no TRACE call.
    """
    try:
        from ..core._trace import trace_provenance

        return trace_provenance()
    except Exception as exc:  # pragma: no cover - environment-dependent
        return {"trace_version": "unavailable", "trace_error": str(exc)}


def write_step_provenance(
    step: str,
    output: str | os.PathLike,
    *,
    inputs: Mapping[str, str | os.PathLike] | None = None,
    parameters: Mapping[str, Any] | None = None,
    counts: Mapping[str, Any] | None = None,
    require_sha: bool = True,
    input_md5: Mapping[str, str] | None = None,
) -> Path:
    """Write `<output>.provenance.json` beside a step's output and return its path.

    Parameters
    ----------
    step : str
        Step name, e.g. `"step_1_bc1_call"`.
    output : path
        The data file this sidecar describes. The sidecar lands beside it with
        `SIDECAR_SUFFIX` appended.
    inputs, parameters, counts : mapping, optional
        Recorded verbatim. `inputs` values are paths and are md5'd by the shared helper;
        `counts` is where a step's row counts belong, so a later reader can check the
        emitted denominator without re-reading the table.
    require_sha : bool, default True
        Raise `ProvenanceError` when no sha can be resolved. Default True because a
        sidecar with a null sha passes a naive "the file exists" check and fails F7.
    input_md5 : mapping, optional
        `{role: raw md5}` already computed by the caller for the CURRENT bytes of that
        input, so a large file shared by many receipts is hashed once. Roles absent
        here are hashed as usual.

    Raises
    ------
    ProvenanceError
        No resolvable git sha (and `require_sha`), or the output does not exist — a
        sidecar for a file that was never written is a false provenance claim.
    """
    out = Path(output)
    if not out.exists():
        raise ProvenanceError(
            f"{step}: refusing to write provenance for {out}, which does not exist. "
            "A sidecar beside a missing output is a false claim (§7 no silent success)."
        )
    sha = build_sha()
    if require_sha and not sha:
        raise ProvenanceError(
            f"{step}: no code sha could be resolved (checked ${BUILD_SHA_ENV}, a "
            f"{BUILD_SHA_FILE} file beside the package, and git_state()). A sidecar with "
            "an empty sha fails acceptance criterion P6 / fail condition F7, so this is "
            "raised rather than written. Export the sha in the run driver."
        )

    prov = StepProvenance(PIPELINE, step, out)
    for role, path in (inputs or {}).items():
        prov.add_input(path, role=role, md5=(input_md5 or {}).get(role))
    if parameters:
        prov.set_parameters(**dict(parameters))
    if counts:
        prov.set_counts(**dict(counts))
    prov.add_output(out)

    # Overwrite the helper's own code block with the resolved one, so the staged-tree
    # case records the real commit instead of `None`, and record which TRACE ran.
    record = prov.build()
    record["code"] = code_state()
    record["trace"] = _trace_state()

    from magestic.utils.provenance import write_json_atomic

    target = Path(str(out) + SIDECAR_SUFFIX)
    write_json_atomic(record, target)
    return target


def write_run_provenance(
    directory: str | os.PathLike,
    *,
    run_name: str,
    step_sidecars: Iterable[str | os.PathLike],
    description: str = "",
) -> Path:
    """Write the run-level `PROVENANCE.json` indexing every step sidecar.

    Acceptance criterion P6 names `PROVENANCE.json` specifically while the dual-path
    parity test globs `*.provenance.json`, so BOTH are emitted and both carry the
    resolved sha.
    """
    target = _write_run_provenance(
        directory,
        pipeline=PIPELINE,
        run_name=run_name,
        step_sidecars=[str(s) for s in step_sidecars],
        description=description,
        extra={"code": code_state(), "trace": _trace_state()},
    )
    return Path(target)
