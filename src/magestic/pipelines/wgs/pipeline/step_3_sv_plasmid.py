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
from ..core.emit import write_plasmid_integration_summary, write_table
from ..core.manifest import WGSSample
from ..core.reconcile_driver import load_oligo_windows
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
SHARD_COLUMNS = (
    "sample",
    "upstream_count",
    "downstream_count",
    "total_count",
    "n_integrated_contigs",
    "half_b_status",
    "read_count_mid_point",
    "frags_genome",
    "frags_plasmid",
    "frags_chimeric",
)

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
        logger.info("step_3: %s already complete — skipping", sample)
        return out

    require_binary(config.sv_plasmid.samtools)
    record = record_for(config, sample)
    bam_path = require_indexed_bam(record)

    contigs = _integrated_contigs(bam_path)
    counts = sv_plasmid.detect_plasmid_integration(
        str(bam_path),
        integrated_contigs=contigs,
        junction_distance=_junction_distances(config, contigs),
    )
    total = counts["upstream_count"] + counts["downstream_count"]
    status, anno = _run_half_b(config, record)

    with atomic(out) as tmp:
        write_table(
            tmp,
            SHARD_COLUMNS,
            [
                [
                    sample,
                    counts["upstream_count"],
                    counts["downstream_count"],
                    total,
                    len(contigs),
                    status,
                    *anno,
                ]
            ],
        )
    logger.info(
        "step_3: %s integration up=%d down=%d total=%d over %d contig(s); half_b=%s",
        sample,
        counts["upstream_count"],
        counts["downstream_count"],
        total,
        len(contigs),
        status,
    )
    return out


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
    windows = load_oligo_windows(
        str(config.reconcile.oligo_windows),
        oligo_col=config.columns.oligo_col,
        chrom_col=config.columns.chrom_col,
        start_col=config.columns.start_col,
        end_col=config.columns.end_col,
    )
    window = windows.get(oligo)
    if window is None:
        return HALF_B_WINDOW_UNKNOWN, blank

    blacklist = config.sv_plasmid.blacklist
    if blacklist is None or not Path(blacklist).exists():
        raise StepIncompleteError(
            "step_3: Half B needs the REDI blacklist (REDI_reads_mapped_regions.txt); "
            f"sv_plasmid.blacklist is {blacklist!r}"
        )

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

    with_evidence = sum(1 for row in summary if row["total_count"] > 0)
    status_idx = SHARD_COLUMNS.index("half_b_status")
    half_b_run = sum(1 for row in rows if row.split("\t")[status_idx] == HALF_B_RUN)
    logger.info(
        "step_3: %d samples; %d with >=1 junction-spanning read; Half B ran on %d of %d",
        n,
        with_evidence,
        half_b_run,
        n,
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
        counts={
            "samples": n,
            "samples_with_integration_evidence": with_evidence,
            "half_b_run": half_b_run,
        },
    )
    return [out, integration_out]


def _integration_summary(rows: Sequence[str]) -> list[dict[str, object]]:
    """`plasmid_integration_summary.tsv` rows: manifest order, stably sorted desc.

    Stable-sorting a manifest-ordered input by `total_count` descending reproduces
    `summarize_plasmid_integration`'s documented ordering exactly, without adding a
    second cohort-wide code path that path (B) would never take.
    """
    idx = {c: SHARD_COLUMNS.index(c) for c in SHARD_COLUMNS}
    summary: list[dict[str, object]] = [
        {
            "sample_name": fields[idx["sample"]],
            "upstream_count": int(fields[idx["upstream_count"]]),
            "downstream_count": int(fields[idx["downstream_count"]]),
            "total_count": int(fields[idx["total_count"]]),
        }
        for fields in (row.split("\t") for row in rows)
    ]
    summary.sort(key=lambda r: int(r["total_count"]), reverse=True)
    return summary
