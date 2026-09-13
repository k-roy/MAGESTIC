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
import logging
import os
import shutil
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path

from ..config import WGSConfig
from ..core.manifest import WGSSample, load_manifest

logger = logging.getLogger(__name__)


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
    return list(load_manifest(config.paths.manifest))


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
    bam = Path(record.bam)
    if not bam.exists():
        raise StepIncompleteError(f"{record.sample}: BAM not found: {bam}")
    # htslib resolves `<bam>.bai` / `<bam>.csi`, not the picard-style `<stem>.bai`,
    # so accepting the latter here would pass this check and then fail inside fetch.
    if not any(bam.with_suffix(bam.suffix + ext).exists() for ext in (".bai", ".csi")):
        raise StepIncompleteError(
            f"{record.sample}: no {bam.name}.bai or .csi beside {bam}. "
            "Stages 1, 3 and 4 all fetch, and htslib will not find a <stem>.bai."
        )
    if not Path(record.reference).exists():
        raise StepIncompleteError(
            f"{record.sample}: reference not found: {record.reference} "
            "(reference is a per-sample attribute — decision D1)"
        )
    return bam


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
    tmp = path.with_name(path.name + ".tmp")
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
    if src.parent == dest.parent:
        os.replace(src, dest)
    else:
        shutil.move(str(src), str(dest))
    return dest


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
    return [dict(zip(columns, r, strict=False)) for r in rows[start + 1 :]]


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


def write_provenance(
    step: str,
    output: Path,
    *,
    inputs: dict[str, Path] | None = None,
    parameters: dict[str, object] | None = None,
    counts: dict[str, object] | None = None,
) -> Path | None:
    """Write the §8 `PROVENANCE.json` sidecar beside a cohort table.

    `qc/provenance.py` is imported here rather than at module scope so that the
    import contract of `pipeline/` does not depend on `qc/`'s internals; a name
    mismatch then surfaces at run time as a loud ImportError instead of breaking
    every step's importability. There is deliberately no `except ImportError`
    fallback: a missing sidecar fails acceptance criterion P6.
    """
    from ..qc.provenance import write_step_provenance

    return write_step_provenance(
        step,
        output,
        inputs=inputs or {},
        parameters=parameters or {},
        counts=counts or {},
    )
