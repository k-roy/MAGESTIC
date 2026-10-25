"""
magestic.pipelines.wgs.pipeline._common — private helpers shared by the step orchestrators.

No analysis logic: these are the orchestration invariants every `step_*.py` has to
hold, factored out once so they cannot drift apart.

1. **Atomic output.** Every artifact is written to `<path>.tmp` and `os.replace`d
   into place, so a file that EXISTS is a file that is COMPLETE. Under a preempted
   and requeued SLURM array task (`ARCHITECTURE_POLICY` §6) a plain `open(path,"w")`
   can leave a truncated TSV that a naive `exists()` skip-check happily accepts —
   idempotent, but idempotently wrong. `core/emit.py` owns the schema of what goes
   into the file; this module owns the guarantee that it lands whole.
2. **Skip-checks read the output.** Combined with (1), `exists()` is a sound
   completion signal, so nothing here writes a sentinel.
3. **Fan-in walks the MANIFEST, never a glob.** Cohort row order is therefore
   manifest order under both the DIRECT path (A) and the SNAKEMAKE path (B); a
   `sorted(dir.iterdir())` fan-in would make §4 byte-identity depend on directory
   listing order. A missing shard raises; it never yields a shorter table.
4. **Sample identity comes from the manifest**, not from a BAM basename. TRACE's
   `bam_suffix` default matches none of this program's BAMs, and an empty glob
   there writes a header-only TSV and returns 0 (fail condition F1). Taking
   `sample` and `bam` from one manifest row makes that unreachable.

Author: Kevin R. Roy
"""

from __future__ import annotations

import csv
import errno
import hashlib
import json
from dataclasses import asdict, is_dataclass
from functools import lru_cache
import itertools
import logging
import os
import shutil
import time
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path

from ..config import WGSConfig
from ..core.manifest import WGSSample, load_manifest

logger = logging.getLogger(__name__)

#: Per-process counter making every `atomic()` temp file name unique, so two writers of
#: the same destination never share a temp path (see `atomic`).
_TMP_COUNTER = itertools.count()


class GateError(RuntimeError):
    """A stage that cannot run on this data because a documented GATE is open (§8).

    Raised — never swallowed — by `step_2_donor_precall` (GATEs 1/2/3) and by any
    step whose input is itself gated. `run_pipeline.run()` records it as a
    `"GATED: ..."` string in the artifacts dict, never as a path, so a reader can
    never mistake a gate for "ran, and found nothing".
    """


class StepIncompleteError(RuntimeError):
    """An output was shorter than its manifest denominator, or a shard was missing.

    The no-silent-success guard (`ARCHITECTURE_POLICY` §7): a short table written
    with exit code 0 is the most dangerous outcome this module can produce.
    """


# --------------------------------------------------------------------------- #
# Paths — every one derived from config, none hardcoded (§1)                    #
# --------------------------------------------------------------------------- #


def processed_root(config: WGSConfig) -> Path:
    """Scratch-aware equivalent of `config.paths.processed`.

    Mirrors `WGSConfig.get_working_dir` exactly (same `use_scratch` / `scratch_dir`
    rule) so a cohort table at the processed ROOT and a shard in a stage
    subdirectory always land on the same filesystem.
    """
    if config.use_scratch and config.scratch_dir is not None:
        return Path(config.scratch_dir)
    return Path(config.paths.processed)


def shard_dir(config: WGSConfig, subdir: str) -> Path:
    """Directory holding one stage's per-sample shards."""
    return config.get_working_dir(subdir)


def cohort_table(config: WGSConfig, name: str) -> Path:
    """A cohort-level table at the processed root, e.g. `bc1_calls.tsv`."""
    return processed_root(config) / name


# --------------------------------------------------------------------------- #
# Manifest access (D1: reference is a per-SAMPLE attribute)                     #
# --------------------------------------------------------------------------- #


def manifest_samples(config: WGSConfig) -> list[WGSSample]:
    """Every manifest row, in manifest order."""
    return list(load_manifest(config.paths.manifest, reference_root=config.paths.reference_root))


def record_for(config: WGSConfig, sample: str) -> WGSSample:
    """The manifest row for one sample. Raises `KeyError` — loudly — if absent."""
    for record in manifest_samples(config):
        if record.sample == sample:
            return record
    raise KeyError(
        f"sample {sample!r} is not in the manifest ({config.paths.manifest}). "
        "The fan-out unit is a manifest row, so an unknown sample is an error, "
        "not an empty result."
    )


def require_indexed_bam(record: WGSSample) -> Path:
    """Assert this sample's BAM and its index exist, and return the BAM path.

    Stages 1, 3 and 4 all reach the BAM through `fetch`/`pileup`, which require an
    index; TRACE documents neither requirement, so a missing `.bai` surfaces as an
    opaque pysam error deep in a loop instead of here.

    `readers/bam_inventory.py` owns the richer inventory (the `@PG` reference
    cross-check `core.manifest.check_reference_consistency` consumes). This is the
    minimal precondition a step must not run without; see the integration note in
    the module report.
    """
    from ..readers.bam_inventory import BamInventoryError, require_indexed, validate_bam_reference
    try:
        bam = require_indexed([record])[0]
        validate_bam_reference(record)
        return bam
    except BamInventoryError as exc:
        raise StepIncompleteError(str(exc)) from exc


def require_binary(name: str) -> str:
    """Assert an external binary is on PATH, and return its resolved location.

    `sv_plasmid` shells out to `samtools` three times with a bare `"samtools"`
    default. That is not expressible as a Python dependency, so the orchestration
    layer preflights it before doing any work.
    """
    found = shutil.which(name)
    if found is None:
        raise StepIncompleteError(
            f"{name!r} is not on PATH. It is an external binary, not a Python "
            "dependency, so the `magestic[wgs]` extra cannot install it."
        )
    return found


# --------------------------------------------------------------------------- #
# Atomic writes                                                                 #
# --------------------------------------------------------------------------- #


@contextmanager
def atomic(path: Path) -> Iterator[Path]:
    """Yield a temp path in the destination directory; `os.replace` it on success.

    Same filesystem, so the replace is atomic. On any exception the temp file is
    removed and the destination is left untouched — a failed task leaves nothing
    for the next attempt's skip-check to mistake for finished work.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # 🔴 The temp name must be UNIQUE PER WRITER, not `<name>.tmp`.
    # Several samples legitimately write the same shared artifact — the locus sidecar is
    # written by EVERY per-sample task — and under the Snakemake path those tasks run
    # concurrently. With a fixed temp name, writer A's `os.replace` consumes the very
    # file writer B is about to replace, and B dies with
    #   FileNotFoundError: ... 'barcode_locus.sidecar.tmp' -> 'barcode_locus.sidecar'
    # Measured: this took down path (B) at `--cores 8` while the serial path (A) never
    # hit it, which is exactly the §6 array-task concurrency this helper exists to
    # survive. `os.replace` onto the shared destination stays atomic and last-writer-wins,
    # which is correct because the content is deterministic.
    tmp = path.with_name(f"{path.name}.{os.getpid()}.{next(_TMP_COUNTER)}.tmp")
    try:
        yield tmp
        os.replace(tmp, path)
    finally:
        # Guarded: a failing unlink (a concurrent array task, a read-only dir) must
        # not replace the in-flight exception with its own and lose the diagnostic —
        # which is exactly the preempted-job case this helper exists for.
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            logger.warning("could not remove the temp file %s", tmp, exc_info=True)


def write_tsv(path: Path, header: str, rows: Sequence[str]) -> Path:
    """Write `header` + already tab-joined `rows` atomically.

    Byte format matches `core.emit.write_table` (`\\n` line endings, no trailing
    blank line) so a cohort table assembled here is indistinguishable from one
    written row-by-row through emit.
    """
    with atomic(path) as tmp, tmp.open("w", newline="\n") as fh:
        fh.write(header.rstrip("\n") + "\n")
        for row in rows:
            fh.write(row.rstrip("\n") + "\n")
    return path


def publish(src: Path, dest: Path) -> Path:
    """Move a file a TRACE function wrote under a staging name into its final place.

    Keeps "exists implies complete" true for outputs this package does not format
    itself (a cleaned BAM, a `.anno`).
    """
    src, dest = Path(src), Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() == dest.resolve():
        return dest
    try:
        os.replace(src, dest)
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise
        with atomic(dest) as tmp:
            shutil.copy2(src, tmp)
        src.unlink()
    return dest


#: A file whose mtime or ctime is this recent is "racily clean" (git's term) and is hashed
#: afresh rather than from the memo. Oak stores WHOLE-SECOND timestamps (planning/284
#: O-W03 measured a sub-second part of 0), so a same-size rewrite inside one second leaves
#: the stat tuple, and hence the memo key, unchanged.
RACY_WINDOW_NS = 3 * 10**9


def _hash_file(path: str) -> str:
    value = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


@lru_cache(maxsize=256)
def _content_hash(path: str, identity: tuple) -> str:
    return _hash_file(path)


def _file_identity(path: Path) -> dict:
    """Content hashes for small inputs; filesystem identity for large data.

    Size, inode, device, mtime and ctime detect normal replacements and in-place edits
    without rereading every multi-GB BAM on each resume. This is a restart guard, not a
    cryptographic attestation for adversarial filesystems.

    Timestamps are whole seconds on Oak (planning/284 O-W03), so a small file touched
    within `RACY_WINDOW_NS` is re-hashed rather than read from the in-process memo. The
    residual blind spot: a file over 1 MiB (stat only) rewritten in place at the SAME
    size within the same second as the previous stat. A BAM rewrite essentially never
    keeps its size; this is recorded rather than paid for with a hash of every BAM.
    """
    path = Path(path).resolve()
    if not path.is_file():
        return {'path': str(path), 'present': False}
    stat = path.stat()
    identity = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    result = {'path': str(path), 'present': True, 'stat': list(identity)}
    if stat.st_size <= 1024 * 1024:
        racy = max(stat.st_mtime_ns, stat.st_ctime_ns) >= time.time_ns() - RACY_WINDOW_NS
        result['sha256'] = _hash_file(str(path)) if racy else _content_hash(str(path), identity)
    return result


def _configured_paths(value):
    if isinstance(value, Path):
        if not value.is_dir():
            yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from _configured_paths(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _configured_paths(item)


def cache_key(config: WGSConfig, step: str, record: WGSSample | None = None,
              *, inputs: dict[str, Path] | None = None) -> dict:
    """Pin current code, configuration and all declared input dependencies."""
    from trace_crispr import wgs as trace_wgs
    from ..readers.bam_inventory import index_path
    configuration = asdict(config)
    configured_inputs = {k: v for k, v in configuration.items() if k != 'paths'}
    configured_inputs['manifest'] = config.paths.manifest
    paths = set(_configured_paths(configured_inputs))
    paths.update(Path(p) for p in (inputs or {}).values() if p is not None)
    if record is not None:
        paths.update((Path(record.bam), Path(record.reference), Path(str(record.reference)+'.fai')))
        index = index_path(record.bam)
        if index:
            paths.add(index)
    sources = sorted(Path(__file__).resolve().parents[1].rglob('*.py'))
    sources += sorted(Path(trace_wgs.__file__).resolve().parent.rglob('*.py'))
    source_hash = hashlib.sha256()
    for source in sources:
        source_hash.update(str(source).encode())
        source_hash.update(_file_identity(source)['sha256'].encode())
    return {'schema': 1, 'stage': step, 'sample': record.sample if record else None,
            'configuration': json.loads(json.dumps(configuration, default=str, sort_keys=True)),
            'record': json.loads(json.dumps(asdict(record), default=str)) if record else None,
            'source_sha256': source_hash.hexdigest(),
            'inputs': [_file_identity(p) for p in sorted(paths, key=str)]}


def _cache_path(outputs: Sequence[Path]) -> Path:
    if not outputs:
        raise ValueError('Cache requires at least one declared output')
    return Path(str(outputs[0])+'.cache.json')


def cache_valid(outputs: Sequence[Path], key: dict) -> bool:
    """Reuse only a complete output set recorded for this exact input key."""
    try:
        saved = json.loads(_cache_path(outputs).read_text())
        return (saved['key'] == key and all(Path(p).is_file() for p in outputs)
                and saved['outputs'] == [_file_identity(p) for p in outputs])
    except (OSError, ValueError, KeyError, TypeError):
        return False


def cache_commit(outputs: Sequence[Path], key: dict) -> None:
    """Publish the completion record last, after every stage output exists."""
    if not all(Path(p).is_file() for p in outputs):
        raise StepIncompleteError('Cannot record completion with missing stage outputs')
    if key['inputs'] != [_file_identity(Path(row['path'])) for row in key['inputs']]:
        raise StepIncompleteError('Input changed during execution; stage completion is not recorded')
    payload = {'key': key, 'outputs': [_file_identity(p) for p in outputs]}
    with atomic(_cache_path(outputs)) as tmp:
        tmp.write_text(json.dumps(payload, sort_keys=True, indent=2)+'\n')


# --------------------------------------------------------------------------- #
# Reading back                                                                  #
# --------------------------------------------------------------------------- #


def _is_stamp(line: str) -> bool:
    """True for a `# key=value` stamp line (e.g. `# synthetic_error_axis=disabled`).

    A stamp has no tab; a header may legitimately start with `#`
    (`donor_outcome.precall` writes `#sample_name<TAB>...`), so the tab is what
    separates the two.
    """
    return line.startswith("#") and "=" in line and "\t" not in line


def read_tsv(path: Path) -> tuple[str, list[str]]:
    """Read a written table as `(header block, data lines)`.

    Leading stamp lines stay in the header block, so concatenating shards never
    interleaves a stamp with data.
    """
    text = Path(path).read_text()
    if not text.strip():
        raise StepIncompleteError(f"{path} is empty — refusing to treat it as finished")
    lines = text.splitlines()
    block: list[str] = []
    for i, line in enumerate(lines):
        block.append(line)
        if not _is_stamp(line):
            return "\n".join(block), lines[i + 1 :]
    return "\n".join(block), []


def read_table(path: Path) -> list[dict[str, str]]:
    """Read a written table as a list of dicts, tolerating stamps and a `#` header."""
    with Path(path).open(newline="") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    if not rows:
        raise StepIncompleteError(f"{path} is empty")
    start = 0
    for i, row in enumerate(rows):
        if _is_stamp("\t".join(row)):
            continue
        start = i
        break
    columns = [c.lstrip("#") for c in rows[start]]
    if any(not c for c in columns) or len(set(columns)) != len(columns):
        raise StepIncompleteError(f"{path}: blank or duplicate table columns")
    for n, row in enumerate(rows[start + 1:], start=start+2):
        if len(row) != len(columns):
            raise StepIncompleteError(f"{path}:{n}: {len(row)} fields for {len(columns)} columns")
    return [dict(zip(columns, r)) for r in rows[start + 1 :]]


# --------------------------------------------------------------------------- #
# Fan-in                                                                        #
# --------------------------------------------------------------------------- #


def collect_shards(
    records: Sequence[WGSSample],
    shard_path: Callable[[WGSSample], Path],
    *,
    label: str,
) -> tuple[str, list[str]]:
    """Concatenate per-sample shards **in manifest order** into `(header, rows)`.

    Raises `StepIncompleteError` naming the missing samples rather than returning a
    short table, and raises on header drift between shards (which would silently
    mis-align every column after the drift).
    """
    header: str | None = None
    rows: list[str] = []
    missing: list[str] = []
    for record in records:
        path = shard_path(record)
        if not path.exists():
            missing.append(record.sample)
            continue
        shard_header, shard_rows = read_tsv(path)
        columns = shard_header.splitlines()[-1].lstrip('#').split('\t')
        identity = 'sample' if 'sample' in columns else 'sample_name' if 'sample_name' in columns else None
        if len(shard_rows) != 1:
            raise StepIncompleteError(f"{label}: {path} has {len(shard_rows)} rows, expected 1")
        fields = shard_rows[0].split('\t')
        if len(fields) != len(columns) or (identity and fields[columns.index(identity)] != record.sample):
            raise StepIncompleteError(f"{label}: {path} has wrong row shape or sample identity")
        if header is None:
            header = shard_header
        elif shard_header != header:
            raise StepIncompleteError(
                f"{label}: header drift at {path} — shards must share one schema.\n"
                f"  expected: {header!r}\n  found:    {shard_header!r}"
            )
        rows.extend(shard_rows)
    if missing:
        raise StepIncompleteError(
            f"{label}: {len(missing)} of {len(records)} per-sample shards are missing "
            f"(first: {', '.join(missing[:5])}). Refusing to write a short cohort table."
        )
    if header is None:
        raise StepIncompleteError(
            f"{label}: the manifest selected 0 samples. An empty manifest is a "
            "configuration error, not an empty result."
        )
    return header, rows


def require_row_count(path: Path, observed: int, expected: int, *, label: str) -> int:
    """Assert a cohort table has exactly `expected` data rows; return the count.

    Acceptance criterion P1 / fail condition F1: a header-only table with exit code
    0 is the signature failure of BAM-discovery mismatch, so the ROW COUNT — never
    the exit code — is what is checked, against the manifest denominator.
    """
    if observed != expected:
        raise StepIncompleteError(
            f"{label}: {path} has {observed} data rows, expected {expected} "
            f"(one per manifest sample). Not reporting success."
        )
    return observed


def preflight_provenance(step: str) -> None:
    """Refuse a step whose outputs could not carry a truthful receipt, before it computes.

    Called first by every `run_step_N`, so a missing or contradicted code sha fails
    before any table is published (planning/284 X-S1-PARTIAL-PUBLISH, X-PROV-ENV-SHA).
    Imported lazily for the reason `write_provenance` gives.
    """
    from ..qc.provenance import preflight_provenance as _preflight

    _preflight(step)


def write_provenance(
    step: str,
    output: Path,
    *,
    inputs: dict[str, Path] | None = None,
    parameters: dict[str, object] | None = None,
    counts: dict[str, object] | None = None,
    input_md5: dict[str, str] | None = None,
) -> Path | None:
    """Write the §8 `PROVENANCE.json` sidecar beside a cohort table.

    `qc/provenance.py` is imported here rather than at module scope so that the
    import contract of `pipeline/` does not depend on `qc/`'s internals; a name
    mismatch then surfaces at run time as a loud ImportError instead of breaking
    every step's importability. There is deliberately no `except ImportError`
    fallback: a missing sidecar fails acceptance criterion P6.
    """
    from ..qc.provenance import SIDECAR_SUFFIX, write_run_provenance, write_step_provenance

    sidecar = write_step_provenance(
        step,
        output,
        inputs=inputs or {},
        parameters=parameters or {},
        counts=counts or {},
        input_md5=input_md5,
    )

    # Refresh the RUN-level `PROVENANCE.json` indexing every step sidecar in this
    # directory. Done here, on every step-sidecar write, rather than at the end of
    # `run_pipeline.run()`, because only the DIRECT path (A) executes `run()` — the
    # SNAKEMAKE path (B) invokes per-step console scripts and would otherwise never
    # produce a run-level sidecar at all. Acceptance criterion P6 names
    # `PROVENANCE.json` specifically and is evaluated on path B's output, so writing it
    # in the one place both paths reach is what makes the two paths equivalent here.
    # Both sidecars are inside §4's single pre-registered parity carve-out, so this
    # cannot perturb the byte-identity comparison.
    directory = Path(output).parent
    write_run_provenance(
        directory,
        run_name=directory.name,
        step_sidecars=sorted(str(p) for p in directory.glob("*" + SIDECAR_SUFFIX)),
        description="magestic.pipelines.wgs run-level index (refreshed by each step).",
    )
    return sidecar
