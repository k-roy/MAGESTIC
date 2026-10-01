"""
magestic.pipelines.wgs.core.emit — TSV writers that take their schema FROM TRACE.

Every writer here reads its column list out of the TRACE module that owns the records —
`bc1_consensus.EMIT_COLUMNS`, `reconcile.EMIT_COLUMNS`, `sv_plasmid.ANNO_HEADER`,
`sv_plasmid.PlasmidIntegrationRow` — and formats each row with TRACE's own
`*_to_row` function. Nothing restates a column name.

That is not fastidiousness. TRACE's reconcile emit schema is **not stable across the
TRACE working trees this program has on disk** (10 columns in one, 13 in another), so a
column list copied into magestic would keep writing the old header against new data and
the drift would be invisible in the file. Taking the schema from the loaded module means
a change surfaces as a different header on the very next run, and
`core._trace.trace_provenance()` records which TRACE produced it.

Two writes live here because TRACE returns rows without writing them, which makes the
write orchestration-owned: `plasmid_integration_summary.tsv` (TRACE's
`summarize_plasmid_integration` returns a sorted list) and the cohort-level `.anno`
aggregate (TRACE writes one `.anno` per sample, never a cohort table).

All writers are deterministic — fixed column order, input row order preserved, explicit
`\\n` newlines, no timestamps — because both orchestration paths call them and §4 requires
byte-identical outputs.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from . import _trace
from .manifest import REFUSED_COLUMNS

if TYPE_CHECKING:  # pragma: no cover - typing only
    from trace_crispr.wgs.bc1_consensus import Bc1Call
    from trace_crispr.wgs.reconcile import CellReconciliation
    from trace_crispr.wgs.sv_plasmid import CleanResult

__all__ = [
    "RECONCILE_COLUMNS",
    "count_data_rows",
    "plasmid_integration_columns",
    "reconcile_row_audited",
    "write_bc1_calls",
    "write_plasmid_integration_summary",
    "write_reconciliations",
    "write_reconciliations_audited",
    "write_refused_table",
    "write_sample_anno_table",
    "write_table",
]

# --------------------------------------------------------------------------- #
# Defect 15 — `confidence` does not mean what a reader assumes                  #
# --------------------------------------------------------------------------- #
#
# `reconcile.reconcile_cell` sets `confidence = "high" if bc_confidence == "high"`,
# downgrading only for `LOW_CONF` or a `NO_EDIT` resting on uncovered sites. It
# therefore reflects the **bc1 barcode call** and never the edit evidence: a 4-read
# site and a 400-read site both ship `confidence=high`. NN 243 measured this — every
# `confidence=high` in its harness was inherited from `bc_confidence`, and in a
# bridge-less cohort `_combine_confidence("high", BRIDGE_NONE)` makes the column
# uniformly "low", carrying no information either way.
#
# Two changes, neither of which touches TRACE:
#   * the bare name `confidence` is shipped as `bc_confidence`, so it cannot be read
#     as edit confidence;
#   * `edit_evidence` is shipped beside it, derived from the scored sites' depth by
#     `core.edit_profile.edit_evidence`.
#
# `write_reconciliations` (TRACE's exact schema) is UNCHANGED, so the dual-path
# parity contract and every existing consumer still hold. The audited schema is
# additive and opt-in.

#: The audited Stage-4 schema: TRACE's columns with `confidence` renamed to
#: `bc_confidence`, plus `edit_evidence`. Built from TRACE's list rather than
#: restated, so a TRACE column addition is picked up automatically.
def _reconcile_columns() -> list[str]:
    tr = _trace.reconcile()
    cols = ["bc_confidence" if c == "confidence" else c for c in tr.EMIT_COLUMNS]
    return cols + ["edit_evidence"]


def __getattr__(name: str):
    """Resolve `RECONCILE_COLUMNS` lazily.

    This module deliberately never imports TRACE at module scope (every other TRACE
    use here is inside a function), so that `core/` stays importable — and testable —
    without it. A module-level constant would break that, so the schema is resolved
    on first attribute access instead.
    """
    if name == "RECONCILE_COLUMNS":
        return _reconcile_columns()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def reconcile_row_audited(cell: CellReconciliation) -> list[str]:
    """One audited Stage-4 row: TRACE's own row, plus the edit-evidence tier.

    The first fields are `reconciliation_to_row`'s verbatim — only the HEADER renames
    `confidence` to `bc_confidence` — so a diff against the unaudited table shows the
    added column and nothing else.
    """
    from .edit_profile import edit_evidence

    tr = _trace.reconcile()
    return list(tr.reconciliation_to_row(cell)) + [edit_evidence(cell.sites)]


def write_reconciliations_audited(
    cells: Iterable[CellReconciliation], path: str | Path
) -> Path:
    """Write the Stage-4 table under `RECONCILE_COLUMNS` (defect 15).

    Use this wherever the table will be read by a person. `write_reconciliations`
    remains the parity writer for TRACE's exact schema and is unchanged.
    """
    return write_table(path, _reconcile_columns(), (reconcile_row_audited(c) for c in cells))


def write_table(
    path: str | Path,
    columns: Sequence[str],
    rows: Iterable[Sequence[Any]],
) -> Path:
    """Write a header'd TSV deterministically, and return the path written.

    Raises `ValueError` on a row whose field count does not match `columns` — a ragged
    row would shift every field after it and still parse.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(columns)
    with path.open("w", newline="\n") as fh:
        fh.write("\t".join(columns) + "\n")
        for i, row in enumerate(rows, start=1):
            values = list(row)
            if len(values) != n:
                raise ValueError(
                    f"{path}: row {i} has {len(values)} fields, header declares {n}"
                )
            fh.write("\t".join("" if v is None else str(v) for v in values) + "\n")
    return path


def count_data_rows(path: str | Path, *, comment_prefix: str = "#") -> int:
    """Count data rows in a written TSV: total lines minus the header and comments.

    Exists because an exit code proves nothing here. TRACE's `bc1_consensus.run` writes a
    header-only file and returns 0 when its BAM glob matches nothing, and `precall` writes
    a header-only file and returns normally when every sample's VCF is missing. Both are
    clean-looking zero-row successes, so every step counts its own output rows against the
    manifest denominator.
    """
    n = 0
    with Path(path).open() as fh:
        for i, line in enumerate(fh):
            if i == 0 or line.startswith(comment_prefix) or not line.strip():
                continue
            n += 1
    return n


def write_bc1_calls(calls: Iterable[Bc1Call], path: str | Path) -> Path:
    """Write the Stage-1 bc1 emit table (`bc1_consensus.EMIT_COLUMNS`, 17 columns today).

    Rows are formatted by `bc1_consensus.call_to_row`, so the barcode written is whatever
    TRACE canonicalized — orientation is TRACE's call, never re-derived here.
    """
    bc1 = _trace.bc1_consensus()
    return write_table(path, bc1.EMIT_COLUMNS, (bc1.call_to_row(c) for c in calls))


def write_reconciliations(cells: Iterable[CellReconciliation], path: str | Path) -> Path:
    """Write the Stage-4 reconciliation table (`reconcile.EMIT_COLUMNS`).

    Column count is version-dependent (see module docstring); it is read from the loaded
    TRACE rather than assumed, and `reconciliation_to_row` formats the `edits` field.
    """
    tr = _trace.reconcile()
    return write_table(path, tr.EMIT_COLUMNS, (tr.reconciliation_to_row(c) for c in cells))


def plasmid_integration_columns() -> list[str]:
    """The `plasmid_integration_summary.tsv` schema, taken from TRACE's TypedDict."""
    return list(_trace.sv_plasmid().PlasmidIntegrationRow.__annotations__)


def write_plasmid_integration_summary(
    rows: Iterable[Mapping[str, Any]],
    path: str | Path,
) -> Path:
    """Write `plasmid_integration_summary.tsv` from `summarize_plasmid_integration` rows.

    TRACE returns the rows already sorted by `total_count` descending but does **not**
    write this file; the write is orchestration-owned. Row order is preserved exactly as
    given, so the sort stays TRACE's.
    """
    columns = plasmid_integration_columns()
    return write_table(path, columns, ([r[c] for c in columns] for r in rows))


def write_sample_anno_table(
    results: Mapping[str, CleanResult],
    path: str | Path,
) -> Path:
    """Write a cohort-level donor-read filtering table: one row per sample.

    TRACE's `clean_donor_reads` writes a two-line `{out_bam}.target_filtering.anno` per
    sample and returns the `CleanResult`; aggregating those into one table is
    orchestration. The header is `sample` plus `sv_plasmid.ANNO_HEADER` verbatim, and each
    row is `CleanResult.anno_row()` verbatim, so the per-sample `.anno` files and this
    table can never disagree on a count.

    Iteration order of `results` is preserved — pass an ordered mapping (a plain `dict`
    built in manifest order) to keep the output byte-identical across paths.
    """
    sv = _trace.sv_plasmid()
    columns = ["sample", *sv.ANNO_HEADER.split("\t")]
    rows = ([sample, *result.anno_row().split("\t")] for sample, result in results.items())
    return write_table(path, columns, rows)


def write_refused_table(rows: Iterable[Mapping[str, Any]], path: str | Path) -> Path:
    """Write the refused cross-reference comparisons table (`manifest.REFUSED_COLUMNS`).

    Same shape as the existing `*_REFUSED_chrVIII.tsv` precedent on Oak: a refusal is
    recorded as data, never silently dropped and never lifted onto another reference.
    """
    return write_table(path, REFUSED_COLUMNS, ([r[c] for c in REFUSED_COLUMNS] for r in rows))
