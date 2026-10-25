"""
magestic.pipelines.wgs.qc.pipeline_qc — the loud-failure guards (policy §7).

ARCHITECTURE_POLICY §7: an orchestrator "MUST exit non-zero on a missing input or failed
step, and MUST NOT log 'completed successfully' without having produced outputs."

**TRACE will not fail loudly for us**, and that is not a criticism of TRACE — it is a
parity-faithful port of Perl scripts whose silent-skip behavior is deliberate. Each guard
below maps to one measured way this module's inputs produce a clean-looking wrong answer
(spec §2.7). None of them is hypothetical:

* `bc1_consensus.run` discovers BAMs with `glob(bam_dir + "/*" + bam_suffix)`, default
  suffix `.sort.dedup.recal2.bam`. This program's BAMs are `plate_1B_clone_1.bam`. An
  empty glob writes a HEADER-ONLY TSV and returns 0 — exit code 0 over zero samples.
* `donor_outcome.precall` drops a sample whose VCF is missing with a bare `continue`:
  no warning, no stderr line, no count. The output TSV looks complete.
* `precall_sample` re-parses its own `f"{chrom}_{pos}_{ref}_to_{alt}"` id with
  `split("_")`, so an underscore in a contig name raises `ValueError` mid-loop — and
  every mandated composite reference carries `plasmid_*` / `integrated_plasmid_*` contigs.
* `sv_plasmid` shells out to a bare `"samtools"`, an external binary no Python extra can
  install.

The guards raise `QCFailure`; nothing here logs a warning and continues.

Author: Kevin R. Roy
"""

from __future__ import annotations

import shutil
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

__all__ = [
    "QCFailure",
    "require_axis_stamp",
    "require_manifest_coverage",
    "require_no_underscore_contigs",
    "require_nonempty_output",
    "require_samtools",
    "summarize_bc1_table",
]


class QCFailure(RuntimeError):
    """A §7 guard fired. Always fatal — never downgraded to a warning."""


def require_nonempty_output(path: str | Path, *, label: str, comment_prefix: str = "#") -> int:
    """Assert a written table has at least one DATA row; return the count.

    The §7 anti-pattern is reporting success without having produced output, and the
    shape it takes here is a header-only TSV with exit code 0. Lines beginning with
    `comment_prefix` (TRACE writes `#sample_name\\t...`, and GATE-1 adds a `# key=value`
    stamp) are not data.
    """
    p = Path(path)
    if not p.exists():
        raise QCFailure(f"{label}: {p} was not written")
    # One counting rule for the package: the header is found by content, so a commented
    # TRACE header is not subtracted a second time (audit W14a). Imported here to keep this
    # module free of package imports at load time.
    from ..core.emit import count_data_rows

    n = count_data_rows(p, comment_prefix=comment_prefix)
    if n == 0:
        raise QCFailure(
            f"{label}: {p} has a header but ZERO data rows. Exit code 0 over zero "
            "samples is the signature of a BAM-discovery mismatch (fail condition F1) — "
            "reporting it as success is the §7 anti-pattern."
        )
    return n


def require_manifest_coverage(observed: int, expected: int, *, label: str) -> int:
    """Assert an output covered every manifest row; return `observed`.

    The denominator is always the manifest, never a directory listing — which is what
    makes the silent-empty glob unreachable rather than merely guarded.
    """
    if observed != expected:
        raise QCFailure(
            f"{label}: {observed} rows for {expected} manifest samples. A short table "
            "written with exit code 0 is the most dangerous outcome this module can "
            "produce; refusing to report success."
        )
    return observed


def require_no_underscore_contigs(chroms: Iterable[str]) -> None:
    """Assert no contig name contains `_` (GATE 3, spec §8.3).

    `donor_outcome.precall_sample` rebuilds a variant id as
    `f"{chrom}_{pos}_{ref}_to_{alt}"` and then re-parses it positionally with
    `split("_")`, so `integrated_plasmid_o1` raises
    `ValueError: invalid literal for int() with base 10: 'plasmid'`. It fires only when
    at least one designed or spontaneous variant was recorded, so it is DATA-DEPENDENT,
    not deterministic per sample — which is exactly why it must be checked up front
    rather than discovered as a mid-cohort crash.
    """
    offenders = sorted({c for c in chroms if "_" in c})
    if offenders:
        raise QCFailure(
            "GATE 3: contig name(s) containing '_' reach TRACE's donor_outcome, which "
            "re-parses its own variant id with split('_') and raises mid-loop: "
            f"{', '.join(offenders[:5])}"
            + (f" (+{len(offenders) - 5} more)" if len(offenders) > 5 else "")
            + ". The fix belongs in TRACE (carry structured fields), not a workaround here."
        )


def require_samtools(binary: str = "samtools") -> str:
    """Assert the `samtools` binary is on PATH; return its resolved location.

    `sv_plasmid` shells out to it three times with a bare default string. An external
    binary is not expressible as a Python dependency, so `magestic[wgs]` cannot install
    it and the orchestration layer preflights it instead.
    """
    found = shutil.which(binary)
    if found is None:
        raise QCFailure(
            f"{binary!r} is not on PATH. TRACE's sv_plasmid shells out to it; it is an "
            "external binary, so the magestic[wgs] extra cannot install it."
        )
    return found


def require_axis_stamp(path: str | Path) -> str:
    """Assert a Stage-2 table carries its `synthetic_error_axis` stamp (GATE 1).

    An unstamped Stage-2 output is indistinguishable from one where the axis genuinely
    ran, so a zero in `synthetic_errors` would read as "no synthesis errors found" when
    it means "the axis was never evaluated."
    """
    from ..core.samplesheet import read_axis_stamp

    axis = read_axis_stamp(Path(path))
    if axis is None:
        raise QCFailure(
            f"{path} carries no synthetic_error_axis stamp. GATE 1 disables the "
            "SYNTHETIC_ERROR axis for the 2021 pool; without the stamp a zero in "
            "`synthetic_errors` is a wrong answer waiting to be quoted."
        )
    return axis


def summarize_bc1_table(
    rows: Sequence[Mapping[str, str]],
    *,
    barcode_length: int = 20,
) -> dict[str, object]:
    """Compute the Stage-1 QC numbers, each with its denominator stated.

    Returns keys `n_rows`, `n_called`, `call_rate`, `bad_length` (list of
    `(sample, barcode)`), `contradictions` (list of samples with `support > 0` but an
    empty confidence — fail condition F3), and the P4b partition of the NO-CALL rows by
    `n_reads_seen`: `nocall_reads_seen_gt0` (reads reached the window but no barcode was
    extracted — an extraction/flank problem) vs `nocall_reads_seen_eq0` (no reads reached
    it — mapping across the 20-N run).

    This function does no I/O and takes no config, so it is testable with literal rows.
    Every count is returned with the denominator it was computed over; a rate without a
    stated denominator is the shape of a wrong answer.
    """
    n_rows = len(rows)
    called = [r for r in rows if (r.get("barcode") or "").strip()]
    bad_length = [
        (r.get("sample", ""), r["barcode"])
        for r in called
        if len(r["barcode"].strip()) != barcode_length
    ]
    contradictions = [
        r.get("sample", "")
        for r in rows
        if _as_int(r.get("support")) > 0
        and not ((r.get("caller_confidence") or "").strip() and (r.get("confidence") or "").strip())
    ]
    nocall = [r for r in rows if not (r.get("barcode") or "").strip()]
    seen_gt0 = sum(1 for r in nocall if _as_int(r.get("n_reads_seen")) > 0)
    return {
        "n_rows": n_rows,
        "n_called": len(called),
        "call_rate": (len(called) / n_rows) if n_rows else 0.0,
        "bad_length": bad_length,
        "contradictions": contradictions,
        "n_nocall": len(nocall),
        "nocall_reads_seen_gt0": seen_gt0,
        "nocall_reads_seen_eq0": len(nocall) - seen_gt0,
    }


def _as_int(value: object) -> int:
    """Best-effort int, 0 for anything unparseable (an absent column is not a count)."""
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return 0
