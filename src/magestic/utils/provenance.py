"""
Step-level provenance sidecars (planning/policies/PIPELINE_CORRECTNESS_POLICY.md §2).

Every step that writes a data file writes ``<output>.provenance.json`` beside it, and every
multi-step run writes a run-level ``PROVENANCE.json``. Minimum contents, all captured AT RUN
TIME by the code that produced the output (never reconstructed afterwards):

    inputs      every input path WITH its md5 (a path is not a version)
    code        repo, branch, commit sha, dirty flag (+ the dirty file list)
    command     the exact argv
    parameters  every threshold and flag that affects the result, by name and value
    environment host, python, key library versions, SLURM context
    started / finished
    counts      rows in, rows out, anything set aside
    decisions   any human ruling the run depends on, quoted, with date and source

This module is dependency-free (stdlib only) so any step can import it.

Author: Kevin R. Roy
"""

from __future__ import annotations

import hashlib
import importlib
import json
import os
import platform
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

SCHEMA = "magestic-step-provenance/1"
SIDECAR_SUFFIX = ".provenance.json"
RUN_SIDECAR_NAME = "PROVENANCE.json"

# The ruling this policy implements; quoted verbatim in every sidecar's `decisions` so a
# reader of the JSON alone knows why the file was rebuilt the way it was.
POLICY_DECISION = {
    "date": "2026-09-10",
    "who": "Kevin Roy",
    "source": "planning/policies/PIPELINE_CORRECTNESS_POLICY.md (SGTC workspace; mirrored to Oak projects/QTL/policies/)",
    "quote": ("Please do this all the correct way — correct our code and rerun rather than hodge-podge patching. "
              "This goes for B, C, and the table versions. We need provenance and sidecars for every step."),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def utc_of(ts: float) -> str:
    return datetime.fromtimestamp(ts, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def md5_file(path: str | os.PathLike, chunk: int = 1 << 22) -> str:
    """Raw md5 of the bytes as stored (md5_kind = 'raw')."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def file_record(path: str | os.PathLike, *, role: Optional[str] = None, md5: Optional[str] = None,
                compute_md5: bool = True, **extra: Any) -> Dict[str, Any]:
    """One ``inputs[]``/``outputs[]`` row: path, bytes, mtime, md5 (raw) and any extra fields."""
    p = Path(path)
    st = p.stat()
    rec: Dict[str, Any] = {"path": str(p), "bytes": st.st_size, "mtime_utc": utc_of(st.st_mtime)}
    if role:
        rec["role"] = role
    if md5 is None and compute_md5:
        md5 = md5_file(p)
    if md5 is not None:
        rec["md5"] = md5
        rec["md5_kind"] = "raw"
    rec.update(extra)
    return rec


def _run(cmd: List[str], cwd: Optional[str] = None, *, strip: bool = True) -> Optional[str]:
    try:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=30)
    except Exception:
        return None
    if r.returncode != 0:
        return None
    return r.stdout.strip() if strip else r.stdout


def git_state(repo_hint: Optional[str | os.PathLike] = None) -> Dict[str, Any]:
    """Commit sha, branch and dirty flag of the repo that holds the running ``magestic`` package."""
    if repo_hint is None:
        try:
            import magestic  # noqa: WPS433
            repo_hint = Path(magestic.__file__).resolve().parent
        except Exception:
            repo_hint = Path(__file__).resolve().parent
    top = _run(["git", "rev-parse", "--show-toplevel"], cwd=str(repo_hint))
    if not top:
        return {"repo": None, "commit": None, "dirty": None, "note": "not a git checkout"}
    # NUL records preserve spaces/newlines and avoid Git's quoted-path syntax.
    # In -z mode a rename has destination first, then a separate source record.
    raw_status = _run(["git", "status", "--porcelain=v1", "-z", "--untracked-files=all"],
                      cwd=top, strip=False)
    dirty_files = []
    records = iter((raw_status or "").split("\0"))
    for record in records:
        if not record:
            continue
        dirty_files.append(record[3:])
        if "R" in record[:2] or "C" in record[:2]:
            next(records, None)
    return {
        "repo": top,
        "branch": _run(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=top),
        "commit": _run(["git", "rev-parse", "HEAD"], cwd=top),
        "dirty": bool(dirty_files) if raw_status is not None else None,
        "dirty_files": dirty_files,
    }


def _pkg_version(name: str) -> Optional[str]:
    try:
        return getattr(importlib.import_module(name), "__version__", "unknown")
    except Exception:
        return None


def environment(packages: Iterable[str] = ("magestic", "pysam", "pandas", "numpy", "rapidfuzz")) -> Dict[str, Any]:
    slurm = {k[6:].lower(): v for k, v in os.environ.items()
             if k in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_JOB_NAME",
                      "SLURM_JOB_PARTITION", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_JOB_NODELIST")}
    return {
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "python_executable": sys.executable,
        "packages": {p: _pkg_version(p) for p in packages},
        "slurm": slurm or None,
        "MAGESTIC_COMMON_DIR": os.environ.get("MAGESTIC_COMMON_DIR"),
    }


class StepProvenance:
    """Collect a step's provenance while it runs, then write ``<output>.provenance.json``.

    Usage::

        prov = StepProvenance(pipeline="guide_donor_bc0", step="02_combine", output=out_path)
        prov.add_input(fastq, role="collapsed_fastq")           # md5 computed unless given
        prov.set_parameters(**parse_parameters())
        ... do the work ...
        prov.set_counts(reads_in=n, rows_out=m)
        prov.add_output(out_path, rows=m)                        # md5 of the output, computed
        prov.write()

    ``decisions`` always carries the 2026-09-10 policy ruling; add more with ``add_decision``.
    """

    def __init__(self, pipeline: str, step: str, output: str | os.PathLike, *,
                 description: str = "", decisions: Optional[List[Dict[str, Any]]] = None,
                 argv: Optional[List[str]] = None):
        self.pipeline = pipeline
        self.step = step
        self.output = str(output)
        self.description = description
        self.started = utc_now()
        self.argv = list(sys.argv if argv is None else argv)
        self.inputs: List[Dict[str, Any]] = []
        self.outputs: List[Dict[str, Any]] = []
        self.parameters: Dict[str, Any] = {}
        self.counts: Dict[str, Any] = {}
        self.decisions: List[Dict[str, Any]] = [dict(POLICY_DECISION)] + list(decisions or [])
        self.notes: List[str] = []
        self.code = git_state()
        self.env = environment()

    # -- builders --------------------------------------------------------------------
    def add_input(self, path: str | os.PathLike, *, role: Optional[str] = None,
                  md5: Optional[str] = None, compute_md5: bool = True, **extra: Any) -> "StepProvenance":
        self.inputs.append(file_record(path, role=role, md5=md5, compute_md5=compute_md5, **extra))
        return self

    def add_output(self, path: str | os.PathLike, *, rows: Optional[int] = None,
                   compute_md5: bool = True, **extra: Any) -> "StepProvenance":
        rec = file_record(path, compute_md5=compute_md5, **extra)
        if rows is not None:
            rec["rows"] = int(rows)
        self.outputs.append(rec)
        return self

    def set_parameters(self, **kw: Any) -> "StepProvenance":
        self.parameters.update(kw)
        return self

    def set_counts(self, **kw: Any) -> "StepProvenance":
        self.counts.update(kw)
        return self

    def add_decision(self, *, date: str, who: str, source: str, quote: str) -> "StepProvenance":
        self.decisions.append({"date": date, "who": who, "source": source, "quote": quote})
        return self

    def note(self, text: str) -> "StepProvenance":
        self.notes.append(text)
        return self

    # -- output ----------------------------------------------------------------------
    def build(self) -> Dict[str, Any]:
        return {
            "schema": SCHEMA,
            "pipeline": self.pipeline,
            "step": self.step,
            "description": self.description,
            "output": self.output,
            "started": self.started,
            "finished": utc_now(),
            "command": {"argv": self.argv, "cwd": os.getcwd()},
            "code": self.code,
            "environment": self.env,
            "parameters": self.parameters,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "counts": self.counts,
            "decisions": self.decisions,
            "notes": self.notes,
        }

    def sidecar_path(self) -> str:
        return self.output + SIDECAR_SUFFIX

    def write(self, path: Optional[str | os.PathLike] = None) -> str:
        target = str(path) if path else self.sidecar_path()
        write_json_atomic(self.build(), target)
        return target


def write_json_atomic(obj: Any, target: str | os.PathLike) -> None:
    import tempfile
    target = str(target)
    os.makedirs(os.path.dirname(os.path.abspath(target)), exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(target)+'.', suffix='.tmp',
                               dir=os.path.dirname(os.path.abspath(target)))
    try:
        with os.fdopen(fd, 'w') as fh:
            json.dump(obj, fh, indent=1, default=str)
            fh.write("\n")
        os.replace(tmp, target)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def read_sidecar(path: str | os.PathLike) -> Dict[str, Any]:
    with open(path) as fh:
        return json.load(fh)


def write_run_provenance(directory: str | os.PathLike, *, pipeline: str, run_name: str,
                         step_sidecars: Iterable[str | os.PathLike], description: str = "",
                         extra: Optional[Dict[str, Any]] = None) -> str:
    """Run-level ``PROVENANCE.json``: an index of every step sidecar (path + md5 of the
    sidecar + its output's md5/rows) plus the aggregate code/environment context."""
    steps = []
    for sc in step_sidecars:
        sc = str(sc)
        d = read_sidecar(sc)
        steps.append({
            "sidecar": sc, "sidecar_md5": md5_file(sc),
            "step": d.get("step"), "output": d.get("output"),
            "outputs": d.get("outputs"), "counts": d.get("counts"),
            "code_commit": (d.get("code") or {}).get("commit"),
            "code_dirty": (d.get("code") or {}).get("dirty"),
            "started": d.get("started"), "finished": d.get("finished"),
        })
    obj = {
        "schema": SCHEMA.replace("step", "run"),
        "pipeline": pipeline,
        "run_name": run_name,
        "description": description,
        "written": utc_now(),
        "command": {"argv": list(sys.argv), "cwd": os.getcwd()},
        "code": git_state(),
        "environment": environment(),
        "decisions": [dict(POLICY_DECISION)],
        "steps": steps,
    }
    if extra:
        obj.update(extra)
    target = os.path.join(str(directory), RUN_SIDECAR_NAME)
    write_json_atomic(obj, target)
    return target
