"""
Step 3 — plasmid-integration detection (Half A) and donor-read cleaning (Half B).

Orchestration over TRACE's `trace_crispr.wgs.sv_plasmid`. Half A is a port of
Kevin R. Roy's `screen_bam_for_plasmid_integration.py`; Half B is a bit-for-bit
port of Shengdi Li's `clean_donor_reads.pl`.

Three things this step owns that TRACE deliberately does not:

* **The `samtools` preflight.** `sv_plasmid` shells out to `samtools` three times
  with a bare `"samtools"` default. An external binary is not expressible as a
  Python dependency, so PATH is checked before any work is done.
* **The write of `plasmid_integration_summary.tsv`.** TRACE returns rows; writing
  is orchestration (`core.emit.write_plasmid_integration_summary`). The cohort
  table is assembled by fan-in over shards in manifest order and then STABLY
  sorted by `total_count` descending — byte-for-byte what
  `summarize_plasmid_integration` returns for a manifest-ordered input.
  `summarize_plasmid_integration` itself is a cohort-wide loop; calling it on path
  (A) while path (B) fans in would put two orchestrations over the same numbers,
  and §4 parity is worth more than the one saved call.
* **Per-contig `junction_distance`.** TRACE's `1129` default is, in its own
  docstring, "the parity literal"; production should derive
  `upstream_arm_length + donor_length` per contig. This step derives it with
  `integrated_plasmid_junction` and REFUSES when a contig's donor sequence is
  unavailable — because `detect_plasmid_integration` silently `continue`s past any
  contig missing from a `{contig: distance}` mapping, which under-counts
  integration rather than failing.

🔴 **Correction to the spec's dependency diagram (§3, D2).** Stage 3 is described
as independent of Stage 1. That holds for Half A, which needs only a BAM. It does
NOT hold for Half B: `clean_donor_reads` requires this sample's donor window
`(chrom, donor_start, donor_end)`, and a sample's donor window is only knowable
from its bc1 -> oligo assignment. Half B therefore depends on an oligo assignment
being configured; when none is available (a bridge-less run, e.g. the §7 smoke)
the shard records a REASON per sample and Half A still runs. Never silently
skipped, and never reported as "no chimeric reads found".

Author: Kevin R. Roy
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

import pysam
from trace_crispr.wgs import sv_plasmid

from ..config import SV_PLASMID_SUBDIR, WGSConfig
from ..core.blacklist_adapter import (
    reference_contig_lengths,
    require_blacklist_applicable,
)
from ..core.emit import write_plasmid_integration_summary, write_table
from ..core.manifest import WGSSample
from ..core.sv_status import (
    SV_UNMEASURED,
    cohort_integration_rate,
    integration_measurement,
    require_no_sv_rate,
    require_no_sv_rate_in_columns,
)
from ..core.window_table import load_windows_audited
from ._common import (
    StepIncompleteError,
    atomic,
    cohort_table,
    collect_shards,
    manifest_samples,
    publish,
    read_table,
    record_for,
    require_binary,
    require_indexed_bam,
    require_row_count,
    shard_dir,
    write_provenance,
    write_tsv,
)

logger = logging.getLogger(__name__)

STEP = "step_3_sv_plasmid"
SHARD_SUFFIX = ".sv_plasmid.tsv"
TABLE_NAME = "sv_plasmid_calls.tsv"
INTEGRATION_TABLE_NAME = "plasmid_integration_summary.tsv"

#: Contig-name convention for an integrated guide-donor plasmid. Mirrors TRACE's
#: private `sv_plasmid._integrated_contigs`, named here rather than reaching across
#: the package boundary into a private helper.
INTEGRATED_PREFIX = "integrated_plasmid_"

#: Per-sample shard schema: Half A counts, then Half B status and its `.anno` fields.
#:
#: 🔴 NN 251 defects 7 and 13 add the two status columns. `integration_status`
#: distinguishes a sample that was MEASURED and found negative from one that could
#: not be measured at all (no `integrated_plasmid_` contig) — before NN 251 both
#: wrote `0`, so a cohort combine could print a 0% integration rate that was an
#: absence of measurement. `sv_status` is always `SV_UNMEASURED`, because this
#: pipeline has no inversion or duplication detector and must not let anyone quote a
#: rate for a quantity it never measured.
SHARD_COLUMNS = (
    "sample",
    "upstream_count",
    "downstream_count",
    "total_count",
    "n_integrated_contigs",
    "integration_status",
    "sv_status",
    "half_b_status",
    "read_count_mid_point",
    "frags_genome",
    "frags_plasmid",
    "frags_chimeric",
)

#: The only SV status this pipeline can honestly emit. See `core/sv_status.py` for
#: why the refusal was chosen deliberately over a half-built detector.
SV_STATUS_VALUE = SV_UNMEASURED

HALF_B_NO_WINDOW_TABLE = "not_run_no_oligo_window_table"
HALF_B_NO_ASSIGNMENT = "not_run_no_oligo_assignment"
HALF_B_WINDOW_UNKNOWN = "not_run_oligo_window_unknown"
HALF_B_RUN = "run"


def shard_path(config: WGSConfig, sample: str) -> Path:
    """Per-sample shard — the Snakemake per-sample rule's declared output."""
    return shard_dir(config, SV_PLASMID_SUBDIR) / f"{sample}{SHARD_SUFFIX}"


def cleaned_bam_path(config: WGSConfig, sample: str) -> Path:
    """Half B's cleaned BAM. NOT re-sorted and NOT indexed — TRACE writes it as-is."""
    return shard_dir(config, SV_PLASMID_SUBDIR) / f"{sample}.donor_cleaned.bam"


def table_path(config: WGSConfig) -> Path:
    """Cohort SV/plasmid table (both halves, one row per sample)."""
    return cohort_table(config, TABLE_NAME)


def integration_table_path(config: WGSConfig) -> Path:
    """Cohort `plasmid_integration_summary.tsv` — a Snakemake rule target."""
    return cohort_table(config, INTEGRATION_TABLE_NAME)


def run_step_3(config: WGSConfig, sample: str | None = None) -> list[Path]:
    """Run both halves. Returns every path written (the CLI prints them all)."""
    if sample is not None:
        return [_run_one_sample(config, sample)]
    return _combine(config)


def _run_one_sample(config: WGSConfig, sample: str) -> Path:
    """One sample. Idempotent: the shard is written atomically, so its presence
    means both halves finished (or were recorded as not-run, with a reason)."""
    out = shard_path(config, sample)
    if out.exists():
        # 🔴 NN 251: a shard written BEFORE this change set has 10 columns, not 12
        # (`integration_status` and `sv_status` were added). `_run_one_sample` skips an
        # existing shard, and `collect_shards` only catches drift BETWEEN shards — so a
        # cohort of uniformly-stale shards would be read with the NEW column indices
        # and `half_b_status` would be taken from `frags_genome`. Stage 3 is the LIVE
        # stage (NN 243), so such shards may already exist on Oak. Refuse, and say what
        # to do, rather than resume into a column shift.
        _require_current_shard_schema(out, sample)
        logger.info("step_3: %s already complete — skipping", sample)
        return out

    require_binary(config.sv_plasmid.samtools)
    record = record_for(config, sample)
    bam_path = require_indexed_bam(record)

    contigs = _integrated_contigs(bam_path)
    # 🔴 NN 251 defect 13: `detect_plasmid_integration` returns 0/0 for a BAM with no
    # `integrated_plasmid_` contig, which is numerically identical to a sample that
    # WAS measured and found negative. `integration_measurement` keeps the two apart
    # and writes "" rather than 0 when nothing was measurable.
    measure = integration_measurement(
        str(bam_path),
        sample=sample,
        integrated_contigs=contigs,
        junction_distance=_junction_distances(config, contigs),
    )
    status, anno = _run_half_b(config, record)
    blank = ""

    # 🔴 NN 255 F2: the refusal was applied to the provenance MAPPING only, never to a
    # written HEADER, so an SV-rate COLUMN reached a file unrefused. A column name is
    # published exactly as readily as a count.
    require_no_sv_rate_in_columns(SHARD_COLUMNS, source="step_3 shard header")

    with atomic(out) as tmp:
        write_table(
            tmp,
            SHARD_COLUMNS,
            [
                [
                    sample,
                    blank if not measure.measured else measure.upstream_count,
                    blank if not measure.measured else measure.downstream_count,
                    blank if not measure.measured else measure.total,
                    measure.n_integrated_contigs,
                    measure.status,
                    SV_STATUS_VALUE,
                    status,
                    *anno,
                ]
            ],
        )
    logger.info(
        "step_3: %s integration=%s up=%s down=%s total=%s over %d contig(s); "
        "half_b=%s; sv=%s",
        sample,
        measure.status,
        measure.upstream_count,
        measure.downstream_count,
        measure.total,
        measure.n_integrated_contigs,
        status,
        SV_STATUS_VALUE,
    )
    return out


def _require_current_shard_schema(shard: Path, sample: str) -> None:
    """Refuse a shard written under an older `SHARD_COLUMNS` (NN 251).

    Idempotent resume is the point of skipping an existing shard, but it is only safe
    while the schema is unchanged. A stale shard is a *silent column shift*, not a
    missing column: every field after the insertion point is read under the wrong
    name. Deleting the shard is the right remedy and the message says so.
    """
    header = shard.read_text().splitlines()[:1]
    if not header:
        raise StepIncompleteError(
            f"step_3: {sample}'s shard {shard} is empty. Delete it and re-run."
        )
    found = tuple(header[0].split("\t"))
    if found != SHARD_COLUMNS:
        raise StepIncompleteError(
            f"step_3: {sample}'s existing shard {shard} was written under a DIFFERENT "
            f"schema ({len(found)} columns, expected {len(SHARD_COLUMNS)}).\n"
            f"  found:    {found}\n  expected: {SHARD_COLUMNS}\n"
            "Resuming would read every field after the first difference under the "
            "wrong column name (NN 251 added `integration_status` and `sv_status`). "
            "Delete this shard — and any sibling written before the change — and "
            "re-run the sample."
        )


def _integrated_contigs(bam_path: Path) -> list[str]:
    """Integrated-plasmid contigs present in THIS sample's own reference (D1)."""
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return [ref for ref in bam.references if ref.startswith(INTEGRATED_PREFIX)]


def _junction_distances(config: WGSConfig, contigs: Sequence[str]) -> dict[str, int] | int:
    """Per-contig junction distance, derived — not the `1129` parity literal.

    Returns the scalar `config.sv_plasmid.junction_distance` only when it has been
    set explicitly, which is the documented way to ask for legacy parity behavior.
    Otherwise every integrated contig must resolve to a donor sequence: a PARTIAL
    mapping makes `detect_plasmid_integration` skip the unmapped contigs in
    silence, which reads as "no integration" rather than "not measured".
    """
    if config.sv_plasmid.junction_distance is not None:
        logger.warning(
            "step_3: pinning junction_distance=%s for every contig — TRACE's own "
            "docstring calls 1129 'the parity literal'; production derives it per contig",
            config.sv_plasmid.junction_distance,
        )
        return int(config.sv_plasmid.junction_distance)
    if not contigs:
        return {}

    donors = _donor_sequences(config)
    distances: dict[str, int] = {}
    unresolved: list[str] = []
    for contig in contigs:
        donor = donors.get(contig[len(INTEGRATED_PREFIX) :])
        if donor is None:
            unresolved.append(contig)
            continue
        distances[contig] = sv_plasmid.integrated_plasmid_junction(
            donor, bp_to_extend=config.sv_plasmid.bp_to_extend
        )
    if unresolved:
        raise StepIncompleteError(
            f"step_3: {len(unresolved)} of {len(contigs)} integrated-plasmid contig(s) have "
            f"no donor sequence (first: {', '.join(unresolved[:5])}). "
            "detect_plasmid_integration skips a contig missing from the mapping without a "
            "word, under-counting integration instead of failing. Supply the oligo table, "
            "or set sv_plasmid.junction_distance to opt into the parity literal."
        )
    return distances


def _donor_sequences(config: WGSConfig) -> dict[str, str]:
    """`oligo_name -> donor sequence`, for the per-contig junction derivation.

    Imported inside the function so this module stays importable without
    `readers/`; a name mismatch then fails loudly at run time rather than breaking
    the import contract. See the integration note in the module report.
    """
    from ..readers.oligo_tables import load_donor_sequences

    return dict(load_donor_sequences(config))


def _run_half_b(config: WGSConfig, record: WGSSample) -> tuple[str, list[str]]:
    """Clean donor reads for one sample. Returns `(status, the .anno fields)`."""
    blank = ["", "", "", ""]
    if config.reconcile.oligo_windows is None:
        return HALF_B_NO_WINDOW_TABLE, blank

    oligo = _assigned_oligo(config, record.sample)
    if oligo is None:
        return HALF_B_NO_ASSIGNMENT, blank
    # 🔴 NN 251 defect 12: the audited loader, so Half B's window lookup also has a
    # denominator. `load_oligo_windows` returns a bare dict, so a sample whose window
    # was DROPPED from the table (3001.0, a blank coordinate, a duplicate oligo) was
    # indistinguishable here from a sample whose oligo is genuinely absent.
    audit = load_windows_audited(
        str(config.reconcile.oligo_windows),
        oligo_col=config.columns.oligo_col,
        chrom_col=config.columns.chrom_col,
        start_col=config.columns.start_col,
        end_col=config.columns.end_col,
    )
    window = audit.windows.get(oligo)
    if window is None:
        logger.warning(
            "step_3: %s oligo %s has no window — %s%s",
            record.sample, oligo, audit.summary(),
            " (this oligo's row was DROPPED by the coordinate filter, not absent)"
            if oligo in audit.dropped_oligos else "",
        )
        return HALF_B_WINDOW_UNKNOWN, blank

    blacklist = config.sv_plasmid.blacklist
    if blacklist is None or not Path(blacklist).exists():
        raise StepIncompleteError(
            "step_3: Half B needs the REDI blacklist (REDI_reads_mapped_regions.txt); "
            f"sv_plasmid.blacklist is {blacklist!r}"
        )
    # 🔴 NN 251 defect 2. `exists()` was the ONLY check, and TRACE's loader is
    # deliberately permissive, so a blacklist whose contigs are `chr15`-named against
    # `chrXV` references parsed cleanly and filtered NOTHING — 0 of 9 names present,
    # 100% inert, and the output looked filtered. This refuses that state.
    audit = require_blacklist_applicable(
        blacklist, reference_contig_lengths(record.reference)
    )
    logger.info("step_3: %s blacklist — %s", record.sample, audit.summary())

    chrom, start, end = window
    final_bam = cleaned_bam_path(config, record.sample)
    staged = final_bam.with_name(final_bam.name + ".part")
    staged.parent.mkdir(parents=True, exist_ok=True)
    result = sv_plasmid.clean_donor_reads(
        str(record.bam),
        str(staged),
        chrom,
        start,
        end,
        str(blacklist),
        flexible_bases=config.sv_plasmid.flexible_bases,
        dup_mode=config.sv_plasmid.dup_mode,
        samtools=config.sv_plasmid.samtools,
    )
    publish(staged, final_bam)
    publish(
        staged.with_name(staged.name + ".target_filtering.anno"),
        final_bam.with_name(final_bam.name + ".target_filtering.anno"),
    )

    anno = result.anno_row().strip("\n").split("\t")
    expected = len(sv_plasmid.ANNO_HEADER.split("\t"))
    if len(anno) != expected:
        raise StepIncompleteError(
            f"step_3: clean_donor_reads returned {len(anno)} .anno field(s) for "
            f"{record.sample}, expected {expected} ({sv_plasmid.ANNO_HEADER!r})"
        )
    return HALF_B_RUN, anno


def _assigned_oligo(config: WGSConfig, sample: str) -> str | None:
    """This sample's assigned oligo, from the configured assignment table.

    Deliberately NOT read from Stage 1's bridge output: Half B must not acquire a
    hidden dependency on Stage 1's fan-out completing, which would turn an
    independent stage into a cohort-wide barrier. A cohort with no assignment table
    runs Half A only, and says so per sample.
    """
    table = config.reconcile.oligo_assignments
    if table is None or not Path(table).exists():
        return None
    return _assignment_index(Path(table), config.columns.oligo_col).get(sample)


#: `(resolved path, mtime) -> {sample: oligo}`. A plain per-sample scan is O(n^2)
#: over a fan-in, which is 657 re-reads for the rearray cohort; keying on mtime
#: means an edited table is still picked up within the same process.
_ASSIGNMENTS: dict[tuple[str, float], dict[str, str]] = {}


def _assignment_index(table: Path, oligo_col: str) -> dict[str, str]:
    """`{sample: oligo}` from the assignment table, memoized on path + mtime."""
    key = (str(table.resolve()), table.stat().st_mtime)
    cached = _ASSIGNMENTS.get(key)
    if cached is None:
        cached = {}
        for row in read_table(table):
            sample = row.get("sample", "").strip()
            oligo = row.get(oligo_col, "").strip()
            if sample and oligo:
                cached.setdefault(sample, oligo)
        _ASSIGNMENTS[key] = cached
    return cached


def _combine(config: WGSConfig) -> list[Path]:
    """Fan-in: every sample, then the cohort table and the integration summary."""
    records = manifest_samples(config)
    for record in records:
        _run_one_sample(config, record.sample)

    header, rows = collect_shards(
        records, lambda r: shard_path(config, r.sample), label="step_3 sv/plasmid"
    )
    out = write_tsv(table_path(config), header, rows)
    n = require_row_count(out, len(rows), len(records), label="step_3 sv/plasmid")

    summary = _integration_summary(rows)
    integration_out = integration_table_path(config)
    with atomic(integration_out) as tmp:
        write_plasmid_integration_summary(summary, tmp)

    # 🔴 NN 251 defect 13: the cohort denominator is the MEASURED samples, never the
    # whole cohort. `cohort_integration_rate` refuses outright when any sample was
    # unmeasured, so a 0% that is really an absence of measurement cannot be printed.
    measured = [r for r in summary if r["total_count"] != ""]
    unmeasured = len(summary) - len(measured)
    with_evidence = sum(1 for row in measured if int(row["total_count"]) > 0)
    rate: float | None = None
    if measured and not unmeasured:
        rate = cohort_integration_rate(
            [(str(r["sample_name"]), int(r["total_count"])) for r in measured]
        )
    status_idx = SHARD_COLUMNS.index("half_b_status")
    half_b_run = sum(1 for row in rows if row.split("\t")[status_idx] == HALF_B_RUN)
    logger.info(
        "step_3: %d samples; %d MEASURED for integration, %d NOT MEASURED (no "
        "integrated_plasmid_ contig); %d of the measured with >=1 junction-spanning "
        "read%s; Half B ran on %d of %d; sv_status=%s for every sample",
        n,
        len(measured),
        unmeasured,
        with_evidence,
        f" (rate {rate:.1%} of measured)" if rate is not None
        else " (no cohort rate: the cohort is not fully measured)",
        half_b_run,
        n,
        SV_STATUS_VALUE,
    )

    write_provenance(
        STEP,
        out,
        inputs={"manifest": Path(config.paths.manifest)},
        parameters={
            "bp_to_extend": config.sv_plasmid.bp_to_extend,
            "flexible_bases": config.sv_plasmid.flexible_bases,
            "dup_mode": config.sv_plasmid.dup_mode,
            "junction_distance": (
                "derived per contig via integrated_plasmid_junction"
                if config.sv_plasmid.junction_distance is None
                else int(config.sv_plasmid.junction_distance)
            ),
        },
        counts=_provenance_counts(
            {
                "samples": n,
                "samples_measured_for_integration": len(measured),
                "samples_not_measured_for_integration": unmeasured,
                "samples_with_integration_evidence": with_evidence,
                "half_b_run": half_b_run,
                "sv_status": SV_STATUS_VALUE,
            }
        ),
    )
    return [out, integration_out]


def _provenance_counts(counts: dict[str, object]) -> dict[str, object]:
    """Guard the provenance record against publishing an SV rate (defect 7).

    This is where a rate would actually reach a file, so this is where the refusal
    lives. `require_no_sv_rate` inspects KEYS: a zero is as publishable as a
    non-zero and just as wrong, because there is no detector behind either.

    NN 255 F2: it now normalizes the key (camelCase split, plural-insensitive) and
    recurses into nested values, so `inversions_rate`, `svRate`, `large_deletion_rate`
    and a nested `summary -> sv` + rate key are all refused where they used to pass.
    (Written without the literal quoted key on purpose: a source-text scan in
    `test_251_trace_fixes.py` refuses one anywhere under `pipeline/`, and it is right to.)
    """
    require_no_sv_rate(counts)
    return counts


def _integration_summary(rows: Sequence[str]) -> list[dict[str, object]]:
    """`plasmid_integration_summary.tsv` rows: manifest order, stably sorted desc.

    Stable-sorting a manifest-ordered input by `total_count` descending reproduces
    `summarize_plasmid_integration`'s documented ordering exactly, without adding a
    second cohort-wide code path that path (B) would never take.
    """
    idx = {c: SHARD_COLUMNS.index(c) for c in SHARD_COLUMNS}

    def _count(value: str) -> int | str:
        """`""` for an unmeasured sample, never `0` (NN 251 defect 13).

        The summary table is read by people and by `redi`. Writing `0` where nothing
        was measured is what lets an absence of measurement be averaged into a rate.
        """
        v = value.strip()
        return int(v) if v.lstrip("-").isdigit() else ""

    summary: list[dict[str, object]] = [
        {
            "sample_name": fields[idx["sample"]],
            "upstream_count": _count(fields[idx["upstream_count"]]),
            "downstream_count": _count(fields[idx["downstream_count"]]),
            "total_count": _count(fields[idx["total_count"]]),
        }
        for fields in (row.split("\t") for row in rows)
    ]
    # Unmeasured samples sort last, and are NOT treated as zero.
    summary.sort(
        key=lambda r: (r["total_count"] != "", r["total_count"] or 0), reverse=True
    )
    return summary
