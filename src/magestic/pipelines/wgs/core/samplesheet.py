"""Build validated donor samplesheets and preserve synthetic-axis provenance.

TRACE requires sample, chrom, donor_start_coord, donor_end_coord, wt_seq,
design_donor, and synthetic_donor. Optional donor_strand ('+'/'-') declares the
supplied sequence orientation relative to the BAM reference and is preserved in
an eighth column. Palindromic WT cannot orient asymmetric design or empirical
sequences without that declaration.

The synthetic-error axis requires an explicit empirical donor sequence. When
synthetic_donor_source='design_donor', the design substitutes for empirical input
and every result is stamped synthetic_error_axis=disabled. A zero in that disabled
axis is not evidence of absent synthesis errors. Choosing 'synthetic_donor' requires
the nonempty empirical column and stamps the axis enabled.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
import os
import tempfile

from . import _trace

__all__ = [
    "AXIS_DISABLED",
    "AXIS_ENABLED",
    "AXIS_STAMP_KEY",
    "SAMPLESHEET_COLUMNS",
    "SYNTHETIC_DONOR_FROM_COLUMN",
    "SYNTHETIC_DONOR_FROM_DESIGN",
    "SamplesheetError",
    "SamplesheetResult",
    "axis_stamp",
    "build_samplesheet",
    "read_axis_stamp",
    "stamp_outcome_table",
    "underscore_contigs",
]

#: The 7 columns `donor_outcome.read_samplesheet` requires, in TRACE's own order.
#: TRACE holds these in a private `_COL` dict, so they are restated here rather than
#: imported — and `build_samplesheet` round-trips every file it writes back through
#: `read_samplesheet`, so a TRACE rename surfaces as a loud failure on the next run
#: rather than as silent column drift.
SAMPLESHEET_COLUMNS: tuple[str, ...] = (
    "sample",
    "chrom",
    "donor_start_coord",
    "donor_end_coord",
    "wt_seq",
    "design_donor",
    "synthetic_donor",
)

#: `synthetic_donor` sources. `FROM_DESIGN` is GATE 1's degenerate input (axis disabled);
#: `FROM_COLUMN` uses each row's own observed step-2 donor (axis enabled).
SYNTHETIC_DONOR_FROM_DESIGN = "design_donor"
SYNTHETIC_DONOR_FROM_COLUMN = "synthetic_donor"

AXIS_STAMP_KEY = "synthetic_error_axis"
AXIS_DISABLED = "disabled"
AXIS_ENABLED = "enabled"

_COORD_COLUMNS = ("donor_start_coord", "donor_end_coord")
_SEQUENCE_COLUMNS = ("wt_seq", "design_donor")


class SamplesheetError(ValueError):
    """A samplesheet row is missing a required field, or the written file failed round-trip."""


@dataclass
class SamplesheetResult:
    """What `build_samplesheet` wrote, and the facts Stage 2 and `qc` need about it.

    Attributes
    ----------
    path : Path
        The samplesheet written.
    n_rows : int
        Data rows written, and the denominator Stage 2 must compare `len(precall(...))`
        against — TRACE's `precall` drops a sample whose VCF is missing with no warning,
        no stderr line and no count, leaving an output that looks complete.
    samples : tuple of str
        Sample names in file order.
    synthetic_error_axis : str
        `AXIS_DISABLED` or `AXIS_ENABLED`. Stamp this into the Stage-2 output and the
        run's `PROVENANCE.json`.
    underscore_contigs : tuple of str
        Contigs containing `_`, retained as compatibility metadata. Structured TRACE
        variants support these names; they are no longer a gate.
    """

    path: Path
    n_rows: int
    samples: tuple[str, ...] = ()
    synthetic_error_axis: str = AXIS_DISABLED
    underscore_contigs: tuple[str, ...] = field(default_factory=tuple)


def underscore_contigs(rows: Iterable[Mapping[str, object]]) -> list[str]:
    """Return underscore contig names for compatibility metadata, without rejecting them."""
    return sorted({str(r.get("chrom", "")) for r in rows if "_" in str(r.get("chrom", ""))})


def build_samplesheet(
    rows: Iterable[Mapping[str, object]],
    out_path: str | Path,
    *,
    synthetic_donor_source: str,
) -> SamplesheetResult:
    """Write TRACE's required donor columns plus supplied strand metadata.

    Parameters
    ----------
    rows : iterable of mapping
        One mapping per sample, carrying at least `sample`, `chrom`, `donor_start_coord`,
        `donor_end_coord`, `wt_seq`, `design_donor` — plus `synthetic_donor` when
        `synthetic_donor_source` is `SYNTHETIC_DONOR_FROM_COLUMN`. Optional
        `donor_strand` ('+'/'-') is retained without inferring it from observed reads.
    out_path : str or Path
        Destination TSV.
    synthetic_donor_source : str
        REQUIRED, no default. `SYNTHETIC_DONOR_FROM_DESIGN` copies `design_donor` into
        `synthetic_donor` and disables the SYNTHETIC_ERROR axis (GATE 1);
        `SYNTHETIC_DONOR_FROM_COLUMN` uses each row's own observed step-2 donor. There is
        deliberately no default: a caller who forgets must get a `TypeError`, never a
        silently-disabled axis with no stamp.

    Returns
    -------
    SamplesheetResult

    Raises
    ------
    SamplesheetError
        Missing or blank required field, non-integer donor coordinate, a duplicate
        sample, or a written file that does not round-trip through
        `donor_outcome.read_samplesheet` at the expected row count.

    Notes
    -----
    Coordinates pass through unchanged: `donor_start_coord`/`donor_end_coord` are
    1-based inclusive on both sides of this boundary (`DonorSpec.donor_start`/`donor_end`),
    so there is no conversion to get wrong.

    Donor orientation is TRACE's problem, not ours: `precall_sample` calls
    `orient_to_genome` itself, so a minus-strand donor needs no pre-orientation here.
    """
    if synthetic_donor_source not in (SYNTHETIC_DONOR_FROM_DESIGN, SYNTHETIC_DONOR_FROM_COLUMN):
        raise SamplesheetError(
            f"synthetic_donor_source must be {SYNTHETIC_DONOR_FROM_DESIGN!r} or "
            f"{SYNTHETIC_DONOR_FROM_COLUMN!r}, got {synthetic_donor_source!r}"
        )
    rows = list(rows)
    if not rows:
        raise SamplesheetError("samplesheet requires at least one donor row")
    columns = SAMPLESHEET_COLUMNS + (("donor_strand",) if any("donor_strand" in row for row in rows) else ())
    required = list(SAMPLESHEET_COLUMNS[:-1])
    if synthetic_donor_source == SYNTHETIC_DONOR_FROM_COLUMN:
        required.append("synthetic_donor")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    samples: list[str] = []
    seen: set[str] = set()
    lines: list[str] = []
    for i, row in enumerate(rows, start=1):
        for column in required:
            if row.get(column) is None or not str(row.get(column, "")).strip():
                raise SamplesheetError(
                    f"{out_path}: row {i} is missing required field {column!r}"
                )
        for column in _COORD_COLUMNS:
            if not str(row[column]).strip().isdigit():
                raise SamplesheetError(
                    f"{out_path}: row {i} has non-integer {column}={row[column]!r} "
                    "(donor coordinates are 1-based inclusive integers)"
                )
        sample = str(row["sample"]).strip()
        if sample in seen:
            raise SamplesheetError(f"{out_path}: duplicate sample {sample!r} on row {i}")
        seen.add(sample)
        samples.append(sample)
        record = {c: str(row[c]).strip() for c in SAMPLESHEET_COLUMNS[:-1]}
        for column in _SEQUENCE_COLUMNS:
            record[column] = record[column].upper()
        record["synthetic_donor"] = (
            record["design_donor"]
            if synthetic_donor_source == SYNTHETIC_DONOR_FROM_DESIGN
            else str(row["synthetic_donor"]).strip().upper()
        )
        if "donor_strand" in columns:
            strand = row.get("donor_strand")
            strand = "" if strand is None else str(strand).strip()
            if strand not in ("", "+", "-"):
                raise SamplesheetError(f"{out_path}: donor_strand must be '+' or '-' on row {i}")
            record["donor_strand"] = strand
        start, end = int(record["donor_start_coord"]), int(record["donor_end_coord"])
        if start < 1 or end < start or len(record["wt_seq"]) != end - start + 1:
            raise SamplesheetError(f"{out_path}: invalid window or wt_seq length on row {i}")
        if any(set(record[column]) - set("ACGTN") for column in ("wt_seq", "design_donor", "synthetic_donor")):
            raise SamplesheetError(f"{out_path}: non-DNA donor sequence on row {i}")
        if any(c in value for value in record.values() for c in "\t\r\n"):
            raise SamplesheetError(f"{out_path}: embedded tab/newline on row {i}")
        if sample in (".", "..") or "/" in sample or "\\" in sample:
            raise SamplesheetError(f"{out_path}: sample is not a safe filename")
        lines.append("\t".join(record[c] for c in columns))

    fd, temporary = tempfile.mkstemp(prefix=out_path.name + ".", suffix=".tmp", dir=out_path.parent)
    try:
        with os.fdopen(fd, "w", newline="\n") as fh:
            fh.write("\t".join(columns) + "\n")
            for line in lines:
                fh.write(line + "\n")
        parsed = _trace.donor_outcome().read_samplesheet(Path(temporary))
        if len(parsed) != len(lines):
            raise SamplesheetError(f"{out_path}: samplesheet round-trip row count mismatch")
        os.replace(temporary, out_path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    if len(parsed) != len(lines):
        raise SamplesheetError(
            f"{out_path}: wrote {len(lines)} rows but donor_outcome.read_samplesheet "
            f"parsed {len(parsed)}"
        )

    axis = (
        AXIS_DISABLED
        if synthetic_donor_source == SYNTHETIC_DONOR_FROM_DESIGN
        else AXIS_ENABLED
    )
    return SamplesheetResult(
        path=out_path,
        n_rows=len(lines),
        samples=tuple(samples),
        synthetic_error_axis=axis,
        underscore_contigs=tuple(underscore_contigs(rows)),
    )


def axis_stamp(synthetic_error_axis: str) -> str:
    """The stamp line, e.g. `# synthetic_error_axis=disabled`."""
    if synthetic_error_axis not in (AXIS_DISABLED, AXIS_ENABLED):
        raise SamplesheetError(
            f"{AXIS_STAMP_KEY} must be {AXIS_DISABLED!r} or {AXIS_ENABLED!r}, "
            f"got {synthetic_error_axis!r}"
        )
    return f"# {AXIS_STAMP_KEY}={synthetic_error_axis}"


def stamp_outcome_table(path: str | Path, *, synthetic_error_axis: str) -> Path:
    """Prepend the axis stamp to a Stage-2 outcome table written by `donor_outcome.precall`.

    TRACE writes its own `#sample_name\\t...` header, so the stamp is inserted as a
    further leading comment line — both lines start with `#`, keeping a comment-skipping
    reader correct. Re-stamping replaces an existing stamp rather than adding a second,
    so array tasks stay idempotent (§6).
    """
    stamp = axis_stamp(synthetic_error_axis)
    path = Path(path)
    existing = path.read_text().splitlines()
    if existing and existing[0].startswith(f"# {AXIS_STAMP_KEY}="):
        existing = existing[1:]
    with path.open("w", newline="\n") as fh:
        fh.write(stamp + "\n")
        for line in existing:
            fh.write(line + "\n")
    return path


def read_axis_stamp(path: str | Path) -> str | None:
    """Read the axis stamp back, or `None` when the file carries none.

    Returns the bare state (`"disabled"` / `"enabled"`). `None` is returned rather than
    raised so `qc` owns the decision to fail the run — which it must, because an
    unstamped Stage-2 output is indistinguishable from one where the axis genuinely ran.
    """
    with Path(path).open() as fh:
        first = fh.readline().strip()
    prefix = f"# {AXIS_STAMP_KEY}="
    return first[len(prefix):] if first.startswith(prefix) else None
