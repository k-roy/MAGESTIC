"""
magestic.pipelines.wgs.readers.oligo_tables — designed-oligo inputs for Stages 2 and 3.

🔴 **Both public functions currently RAISE on this program's cohort, by design.** They
are the Stage-2/Stage-3 input seam, and the inputs they need do not exist yet. Raising a
`MissingOligoInputError` that names the gate is the honest implementation; returning an
empty iterable would let Stage 2 write a complete-looking table over zero rows, which is
precisely the silent-empty failure the whole module is built to make unreachable.

What is missing, measured rather than assumed (spec §8.1):

* **GATE 1 — no per-oligo step-2 (synthetic) donor sequence exists for the 2021 pool.**
  `donor_outcome.read_samplesheet` requires a `synthetic_donor` column, and
  `bucket_variant` uses it non-optionally (`e3 = lev(tmp, step2)`), so the
  SYNTHETIC_ERROR axis is unreachable without it. It is an EMPIRICAL per-oligo
  observation — the design donor carrying the synthesis error actually present in the
  cloned library — not a designable quantity, so it cannot be derived from design data.
  Every sequence-like candidate column in the harmonized table was checked over all
  40,191 2021-pool rows and ruled out: `donor_changed` is a boolean, `donor_changes` and
  `junction_changes` are empty list literals, `US_syn_codons`/`DS_syn_codons` are codons,
  and `extended_donor` is 1,129 nt against a 117-129 nt donor window (0% filled for the
  2021 pool anyway).
* **`donor_start_coord` / `donor_end_coord` are only 50.2% filled** (20,177/40,191) for
  the 2021 pool, and two of those are required samplesheet columns.

The degenerate unblock TRACE's own tests sanction (`test_precall_partial` sets
`step2 = design`) is available through `config.donor_precall.synthetic_donor_source ==
"design_donor"`, which collapses cleanly to DESIGNED/SPONTANEOUS and stamps the output
`synthetic_error_axis=disabled`. It still needs a real donor table to supply the other
six columns, which is what these functions would read.

**The authoritative table is the `common/annotation_files/` copy**, not the
`common/oligo_designs/` one. They were diffed by CONTENT over all 146,250 comparable
rows: 40,191 differ, in exactly two columns (`subpool_ID`, `subpool_priming_seq`, empty
in the older copy). The annotation_files copy is a strict superset.

Author: Kevin R. Roy
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
        if columns:
            missing = [c for c in columns if c not in (reader.fieldnames or [])]
            if missing:
                raise MissingOligoInputError(
                    f"{p}: missing required column(s) {missing}; "
                    f"found {list(reader.fieldnames or [])[:12]}"
                )
        yield from reader


def donor_samplesheet_rows(
    config: WGSConfig,
    samples: Sequence[str],
) -> Iterator[dict[str, Any]]:
    """Rows for TRACE's 7-column donor samplesheet, one per sample.

    🔴 **Raises today.** Producing these rows requires, per sample, the oligo assigned to
    that clone plus that oligo's `chrom` / `donor_start_coord` / `donor_end_coord` /
    `wt_seq` / `design_donor` / `synthetic_donor`. GATE 1 means `synthetic_donor` does
    not exist for the 2021 pool, and the coordinates are 50.2% filled. No configuration
    key names a table that carries all seven.

    Raises
    ------
    MissingOligoInputError
        Always, until a donor table is configured. The message names what is missing so
        a reader is never left to infer it from an empty output.
    """
    table = getattr(config.reconcile, "oligo_windows", None)
    raise MissingOligoInputError(
        "Stage 2 donor samplesheet rows cannot be built for "
        f"{len(samples)} sample(s): no input supplies all seven columns TRACE's "
        "read_samplesheet requires (sample, chrom, donor_start_coord, donor_end_coord, "
        "wt_seq, design_donor, synthetic_donor).\n"
        "  GATE 1: no per-oligo step-2 (synthetic) donor sequence exists for the 2021 "
        "pool; it is an empirical observation of the cloned library, not a designable "
        "quantity, and every candidate column in the harmonized table was checked and "
        "ruled out over all 40,191 rows.\n"
        "  Coordinates: donor_start_coord/donor_end_coord are 50.2% filled "
        "(20,177/40,191) for that pool.\n"
        f"  reconcile.oligo_windows is currently {table!r}.\n"
        "Returning an empty iterable here would write a complete-looking Stage-2 table "
        "over zero rows, so this raises instead."
    )


def load_donor_sequences(config: WGSConfig) -> dict[str, str]:
    """`{oligo_name: donor sequence}` for Stage 3's per-contig junction derivation.

    Stage 3 Half B derives `junction_distance` per contig via
    `sv_plasmid.integrated_plasmid_junction(donor)` instead of passing TRACE's `1129`
    parity literal, which needs each integrated contig's donor sequence.

    🔴 **Raises until a donor-sequence table is configured.** Under-counting plasmid
    integration because a contig's donor was unknown would look exactly like a real
    biological result, so Stage 3 refuses rather than silently skipping the contig —
    which is what `detect_plasmid_integration` does on its own with a `{contig: distance}`
    mapping that lacks an entry.
    """
    raise MissingOligoInputError(
        "Stage 3 Half B needs {oligo_name: donor_sequence} to derive each integrated "
        "contig's junction distance (sv_plasmid.integrated_plasmid_junction), but no "
        "donor-sequence table is configured. TRACE's detect_plasmid_integration SILENTLY "
        "SKIPS a contig missing from its {contig: distance} mapping, which would "
        "under-count integration and read as a biological result — so this refuses "
        "instead. Configure a donor table, or run Stage 3 Half A only."
    )
