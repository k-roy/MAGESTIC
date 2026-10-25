"""Strict, explicit donor input contracts shared by clonal WGS stages.

All donor coordinates are 1-based inclusive on the sample BAM reference. A donor
sequence table is distinct from a table containing only coordinate windows; empirical
synthetic-donor sequence is never inferred from design data. Optional metadata fields
are retained by the generic reader, including Stage 4 intended_variants.
"""

from __future__ import annotations

import csv
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:  # pragma: no cover - typing only, avoids a config import at runtime
    from ..config import WGSConfig

__all__ = [
    "MissingOligoInputError",
    "donor_samplesheet_rows",
    "load_donor_sequences",
    "read_oligo_table",
]


class MissingOligoInputError(RuntimeError):
    """A designed-oligo input Stage 2 or 3 needs is not configured or does not exist.

    Raised — never swallowed into an empty result. An empty iterable here would produce
    a Stage-2 table with zero rows and exit code 0, the exact silent-empty shape this
    module exists to prevent.
    """


def read_oligo_table(
    path: str | Path,
    *,
    columns: Sequence[str] | None = None,
) -> Iterator[dict[str, str]]:
    """Stream a header'd oligo TSV row by row.

    Streamed, never loaded whole: the harmonized table is ~235 MB / 146,251 rows and
    ARCHITECTURE_POLICY §5 requires bounded memory for any step that scales with row
    count.

    Parameters
    ----------
    path : str or Path
        The oligo table. Pin the `common/annotation_files/` copy (see the module
        docstring) — the `oligo_designs/` copy is a strict subset.
    columns : sequence of str, optional
        Required columns. A missing one raises immediately rather than yielding rows
        with silently absent keys.
    """
    p = Path(path)
    if not p.exists():
        raise MissingOligoInputError(f"oligo table not found: {p}")
    with p.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        if not header or len(header) != len(set(header)) or any(not h for h in header):
            raise MissingOligoInputError(f"{p}: absent, blank, or duplicate TSV header")
        missing = set(columns or ()) - set(header)
        if missing:
            raise MissingOligoInputError(f"{p}: missing required columns {sorted(missing)}")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise MissingOligoInputError(f"{p}:{reader.line_num}: malformed TSV row")
            yield {k: v.strip() for k, v in row.items()}


_DONOR_COLUMNS = ("oligo_name", "chrom", "donor_start_coord", "donor_end_coord", "wt_seq", "design_donor")


def _unique_rows(path: Path, required: Sequence[str], key: str) -> dict[str, dict[str, str]]:
    rows = {}
    for row in read_oligo_table(path, columns=required):
        name = row[key]
        if not name:
            raise MissingOligoInputError(f"{path}: blank {key}")
        if name in rows:
            raise MissingOligoInputError(f"{path}: duplicate {key} {name!r}; explicit deduplication required")
        rows[name] = row
    if not rows:
        raise MissingOligoInputError(f"{path}: no data rows")
    return rows


def donor_samplesheet_rows(config: WGSConfig, samples: Sequence[str]) -> Iterator[dict[str, Any]]:
    """Join exact sample assignments to one explicit donor specification per clone.

    donor_precall.donor_table is preferred. A rich reconcile.oligo_windows table is
    accepted only when it carries every donor column. `reference`, if supplied in
    a donor row, must resolve to that sample's manifest reference. Independently,
    TRACE checks wild-type donor sequence against the actual reference window.
    Optional donor_strand declares '+'/'-' sequence orientation and is retained;
    an asymmetric design/empirical donor with self-reverse-complement WT needs it.
    """
    table = getattr(config.donor_precall, "donor_table", None) or config.reconcile.oligo_windows
    assignments = config.reconcile.oligo_assignments
    if table is None or assignments is None:
        raise MissingOligoInputError("Stage 2 requires donor_precall.donor_table and reconcile.oligo_assignments")
    if not samples or len(samples) != len(set(samples)):
        raise MissingOligoInputError("Stage 2 requested samples must be nonempty and unique")
    required = list(_DONOR_COLUMNS)
    if config.donor_precall.synthetic_donor_source == "synthetic_donor":
        required.append("synthetic_donor")
    by_sample: dict[str, str] = {}
    for row in read_oligo_table(assignments, columns=("sample", "oligo_name")):
        name, oligo = row["sample"], row["oligo_name"]
        if not name or not oligo:
            raise MissingOligoInputError(f"{assignments}: empty sample/oligo assignment")
        if name in by_sample:
            raise MissingOligoInputError(f"{assignments}: Stage 2 requires exactly one assignment per sample: {name}")
        by_sample[name] = oligo
    wanted = {by_sample[sample] for sample in samples if sample in by_sample}
    donors = {}
    for row in read_oligo_table(table, columns=required):
        oligo = row["oligo_name"]
        if oligo not in wanted:
            continue
        if oligo in donors:
            raise MissingOligoInputError(f"{table}: duplicate assigned oligo {oligo!r}")
        donors[oligo] = row
    from ..core.manifest import load_manifest
    records = {record.sample: record for record in load_manifest(
        config.paths.manifest, reference_root=config.paths.reference_root)}
    output = []
    for sample in samples:
        if sample not in records or sample not in by_sample:
            raise MissingOligoInputError(f"Stage 2 missing manifest/oligo assignment for {sample!r}")
        oligo = by_sample[sample]
        if oligo not in donors:
            raise MissingOligoInputError(f"{table}: assigned oligo {oligo!r} for {sample!r} is absent")
        row = dict(donors[oligo], sample=sample)
        for name in required:
            if not row[name]:
                raise MissingOligoInputError(f"{table}: oligo {oligo!r} has blank required {name}")
        if row.get("coordinate_system", "1-based-inclusive") != "1-based-inclusive":
            raise MissingOligoInputError(f"{table}: unsupported donor coordinate_system")
        if row.get("reference"):
            declared = Path(row["reference"])
            if not declared.is_absolute():
                declared = Path(table).parent / declared
            if declared.resolve() != Path(records[sample].reference).resolve():
                raise MissingOligoInputError(f"{table}: donor/reference provenance mismatch for {sample}")
        output.append(row)
    yield from output


def load_donor_sequences(config: WGSConfig) -> dict[str, str]:
    """Read explicit oligo_name/donor_sequence for Stage 3 junction geometry."""
    table = getattr(config.sv_plasmid, "donor_sequences", None)
    if table is None:
        raise MissingOligoInputError("Stage 3 requires sv_plasmid.donor_sequences (oligo_name, donor_sequence)")
    rows = _unique_rows(Path(table), ("oligo_name", "donor_sequence"), "oligo_name")
    result = {}
    for name, row in rows.items():
        sequence = row["donor_sequence"].upper()
        if not sequence or set(sequence) - set("ACGT"):
            raise MissingOligoInputError(f"{table}: invalid or ambiguous donor_sequence for {name!r}")
        result[name] = sequence
    return result
