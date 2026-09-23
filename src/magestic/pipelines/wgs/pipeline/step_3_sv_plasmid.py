"""
Step 3 — plasmid-integration detection (Half A) and donor-read cleaning (Half B).

Orchestration over TRACE's `trace_crispr.wgs.sv_plasmid`. Half A is a port of
Kevin R. Roy's `screen_bam_for_plasmid_integration.py`; Half B derives from
Shengdi Li's `clean_donor_reads.pl`, with strict production defaults and explicit
historical parity mode.

Three things this step owns that TRACE deliberately does not:

* **The `samtools` preflight.** `sv_plasmid` shells out to `samtools` three times
  with a bare `"samtools"` default. An external binary is not expressible as a
  Python dependency, so PATH is checked when Half B actually runs.
* **The write of `plasmid_integration_summary.tsv`.** TRACE returns rows; writing
  is orchestration (`core.emit.write_plasmid_integration_summary`). The cohort
  table is assembled by fan-in over shards in manifest order and then STABLY
  sorted by `total_count` descending — byte-for-byte what
  `summarize_plasmid_integration` returns for a manifest-ordered input.
  `summarize_plasmid_integration` itself is a cohort-wide loop; calling it on path
  (A) while path (B) fans in would put two orchestrations over the same numbers,
  and §4 parity is worth more than the one saved call.
* **Per-contig `junction_distance`.** Derive `upstream_arm_length + donor_length`
  and verify both donor copies against the actual sample reference. The historical
  1129 parity literal is no longer an implicit default. Missing distances fail;
  zero eligible integration contigs produce an unmeasured status and blank counts.
* **Reference-specific blacklist and resumability.** Standard Half B requires the
  explicit source reference for blacklist coordinates, exact sequence identity,
  and all output artifacts before a cache completion record is published.

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
import hashlib
import tempfile
from collections.abc import Sequence
from pathlib import Path

import pysam
from trace_crispr.wgs import sv_plasmid

from ..config import SV_PLASMID_SUBDIR, WGSConfig
from ..core.emit import write_plasmid_integration_summary, write_table
from ..core.manifest import WGSSample
from ._common import (
    StepIncompleteError,
    atomic,
    cache_key,
    cache_valid,
    cache_commit,
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
    "half_a_status",
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
HALF_A_NO_CONTIGS = "not_run_no_integrated_contigs"


def shard_path(config: WGSConfig, sample: str) -> Path:
    """Per-sample shard — the Snakemake per-sample rule's declared output."""
    return shard_dir(config, SV_PLASMID_SUBDIR) / f"{sample}{SHARD_SUFFIX}"


def cleaned_bam_path(config: WGSConfig, sample: str) -> Path:
    """Half B's complete, coordinate-order-preserving, indexed cleaned BAM."""
    return shard_dir(config, SV_PLASMID_SUBDIR) / f"{sample}.donor_cleaned.bam"


def cleaned_bam_index_path(config: WGSConfig, sample: str, *, record: WGSSample | None = None) -> Path:
    """Deterministic BAI/CSI companion selected from the source reference lengths."""
    record = record if record is not None else record_for(config, sample)
    with pysam.AlignmentFile(str(record.bam), "rb") as bam:
        suffix = sv_plasmid.bam_index_suffix(bam.lengths)
    return Path(str(cleaned_bam_path(config, sample)) + suffix)


def _remove_cleaned_indexes(path: Path) -> None:
    for suffix in (".bai", ".csi"):
        Path(str(path) + suffix).unlink(missing_ok=True)
        path.with_suffix(suffix).unlink(missing_ok=True)


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


def _run_one_sample(config: WGSConfig, sample: str, *, record: WGSSample | None = None) -> Path:
    """Serialize a complete sample generation, including its cache check/commit."""
    record = record if record is not None else record_for(config, sample)
    if record.sample != sample:
        raise ValueError("Stage 3 record and requested sample disagree")
    with sv_plasmid.output_lock(str(shard_path(config, sample))):
        # Also honor standalone TRACE's destination lock. Holding both locks
        # through cache publication prevents mixed generations across entry points.
        with sv_plasmid.output_lock(str(cleaned_bam_path(config, sample))):
            return _run_one_sample_locked(config, sample, record=record)


def _run_one_sample_locked(config: WGSConfig, sample: str, *, record: WGSSample | None = None) -> Path:
    """One sample. Idempotent: the shard is written atomically, so its presence
    means both halves finished (or were recorded as not-run, with a reason)."""
    out = shard_path(config, sample)
    record = record if record is not None else record_for(config, sample)
    bam_path = require_indexed_bam(record)
    key = cache_key(config, STEP, record)
    outputs = [out]
    if out.exists():
        try:
            existing = list(read_table(out))
        except (OSError, ValueError, StepIncompleteError):
            existing = []
        if len(existing) == 1 and existing[0].get("half_b_status") == HALF_B_RUN:
            cleaned = cleaned_bam_path(config, sample)
            outputs.extend([cleaned, cleaned_bam_index_path(config, sample, record=record), Path(str(cleaned) + ".target_filtering.anno")])
    if cache_valid(outputs, key):
        logger.info("step_3: %s complete for current inputs — skipping", sample)
        return out

    contigs = _integrated_contigs(bam_path)
    half_a_status = "run" if contigs else HALF_A_NO_CONTIGS
    counts = {"upstream_count": None, "downstream_count": None}
    total = None
    if contigs:
        kwargs = {k: getattr(config.sv_plasmid, k) for k in ("min_mapq", "min_anchor") if getattr(config.sv_plasmid, k) is not None}
        counts = sv_plasmid.detect_plasmid_integration(
            str(bam_path), integrated_contigs=contigs,
            junction_distance=_junction_distances(config, contigs, Path(record.reference)), **kwargs,
        )
        total = counts["upstream_count"] + counts["downstream_count"]
    status, anno = _run_half_b(config, record)
    if status != HALF_B_RUN:
        # A changed assignment/config can turn a prior measured result into a
        # not-run result. Remove generated companions so they cannot look current.
        stale = cleaned_bam_path(config, sample)
        stale.unlink(missing_ok=True)
        _remove_cleaned_indexes(stale)
        Path(str(stale) + ".target_filtering.anno").unlink(missing_ok=True)

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
                    half_a_status,
                    status,
                    *anno,
                ]
            ],
        )
    logger.info(
        "step_3: %s integration up=%s down=%s total=%s over %d contig(s); half_b=%s",
        sample,
        counts["upstream_count"],
        counts["downstream_count"],
        total,
        len(contigs),
        status,
    )
    outputs = [out]
    if status == HALF_B_RUN:
        cleaned = cleaned_bam_path(config, sample)
        outputs.extend([cleaned, cleaned_bam_index_path(config, sample, record=record), Path(str(cleaned) + ".target_filtering.anno")])
    cache_commit(outputs, key)
    return out


def _integrated_contigs(bam_path: Path) -> list[str]:
    """Integrated-plasmid contigs present in THIS sample's own reference (D1)."""
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return [ref for ref in bam.references if ref.startswith(INTEGRATED_PREFIX)]


def _junction_distances(config: WGSConfig, contigs: Sequence[str], reference: Path | None = None) -> dict[str, int] | int:
    """Per-contig junction distance, derived — not the `1129` parity literal.

    Returns the scalar `config.sv_plasmid.junction_distance` only when it has been
    set explicitly, which is the documented way to ask for legacy parity behavior.
    Otherwise every integrated contig must resolve to a donor sequence. TRACE
    independently rejects partial mappings; neither layer turns missing inputs
    into an observed integration-negative zero.
    """
    if config.sv_plasmid.junction_distance is not None:
        logger.warning(
            "step_3: explicitly pinning junction_distance=%s for every contig; "
            "the default derives distances from reference-validated donor sequences",
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
            "Supply sv_plasmid.donor_sequences, or an explicit verified uniform junction_distance."
        )
    if reference is not None:
        with pysam.FastaFile(str(reference)) as fasta:
            for contig, distance in distances.items():
                donor = donors[contig[len(INTEGRATED_PREFIX):]].upper()
                length = fasta.get_reference_length(contig)
                left = fasta.fetch(contig, config.sv_plasmid.bp_to_extend, distance).upper()
                right = fasta.fetch(contig, length - distance, length - config.sv_plasmid.bp_to_extend).upper()
                if left != donor or right != donor:
                    raise StepIncompleteError(f"step_3: derived donor boundaries disagree with reference sequence for {contig}")
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
    windows = _donor_windows(config)
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
    if config.sv_plasmid.dup_mode == "standard" and config.sv_plasmid.blacklist_reference is None:
        raise StepIncompleteError("step_3: standard Half B requires sv_plasmid.blacklist_reference to establish coordinate identity")
    require_binary(config.sv_plasmid.samtools)
    final_bam = cleaned_bam_path(config, record.sample)
    final_index = cleaned_bam_index_path(config, record.sample, record=record)
    final_bam.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{record.sample}.cleaned_", dir=final_bam.parent) as temp:
        staged = Path(temp) / "cleaned.bam"
        result = sv_plasmid.clean_donor_reads(
            str(record.bam), str(staged), chrom, start, end, str(blacklist),
            flexible_bases=config.sv_plasmid.flexible_bases,
            dup_mode=config.sv_plasmid.dup_mode,
            samtools=config.sv_plasmid.samtools,
            reference_path=str(record.reference),
            blacklist_reference=str(config.sv_plasmid.blacklist_reference) if config.sv_plasmid.blacklist_reference is not None else None,
            **({"min_mapq": config.sv_plasmid.min_mapq} if config.sv_plasmid.min_mapq is not None else {}),
        )
        anno = result.anno_row().strip("\n").split("\t")
        expected = len(sv_plasmid.ANNO_HEADER.split("\t"))
        if len(anno) != expected:
            raise StepIncompleteError(
                f"step_3: clean_donor_reads returned {len(anno)} .anno field(s) for "
                f"{record.sample}, expected {expected} ({sv_plasmid.ANNO_HEADER!r})"
            )
        with pysam.AlignmentFile(str(staged), "rb") as checked:
            checked.check_index()
        _remove_cleaned_indexes(final_bam)
        publish(staged, final_bam)
        publish(Path(str(staged) + final_index.suffix), final_index)
        publish(Path(str(staged) + ".target_filtering.anno"), Path(str(final_bam) + ".target_filtering.anno"))
    return HALF_B_RUN, anno


_WINDOWS: dict[tuple, dict[str, tuple[str, int, int]]] = {}


def _donor_windows(config: WGSConfig) -> dict[str, tuple[str, int, int]]:
    """Read donor coordinates without requiring Stage 4 intended edit alleles."""
    from ..readers.oligo_tables import read_oligo_table
    path = Path(config.reconcile.oligo_windows)
    columns = (config.columns.oligo_col, config.columns.chrom_col, config.columns.start_col, config.columns.end_col)
    stat = path.stat()
    key = (str(path.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, columns)
    if key not in _WINDOWS:
        windows = {}
        for row in read_oligo_table(path, columns=columns):
            name, chrom = row[columns[0]], row[columns[1]]
            try:
                start, end = int(row[columns[2]]), int(row[columns[3]])
            except ValueError as exc:
                raise StepIncompleteError(f"step_3: invalid donor coordinates for {name}") from exc
            if not name or not chrom or start < 1 or end < start:
                raise StepIncompleteError(f"step_3: invalid donor window for {name!r}")
            window = (chrom, start, end)
            if name in windows and windows[name] != window:
                raise StepIncompleteError(f"step_3: conflicting donor windows for {name}")
            windows[name] = window
        if not windows:
            raise StepIncompleteError(f"step_3: no donor windows in {path}")
        _WINDOWS[key] = windows
    return _WINDOWS[key]


def _assigned_oligo(config: WGSConfig, sample: str) -> str | None:
    """This sample's assigned oligo, from the configured assignment table.

    Deliberately NOT read from Stage 1's bridge output: Half B must not acquire a
    hidden dependency on Stage 1's fan-out completing, which would turn an
    independent stage into a cohort-wide barrier. A cohort with no assignment table
    runs Half A only, and says so per sample.
    """
    table = config.reconcile.oligo_assignments
    if table is None:
        return None
    if not Path(table).is_file():
        raise StepIncompleteError(f"step_3: configured oligo assignments missing: {table}")
    return _assignment_index(Path(table), config.columns.oligo_col).get(sample)


#: Assignment content and configured column are part of the cache identity.
_ASSIGNMENTS: dict[tuple[str, str, str], dict[str, str]] = {}


def _assignment_index(table: Path, oligo_col: str) -> dict[str, str]:
    """`{sample: oligo}` memoized on path, content hash and configured column."""
    key = (str(table.resolve()), hashlib.sha256(table.read_bytes()).hexdigest(), oligo_col)
    cached = _ASSIGNMENTS.get(key)
    if cached is None:
        cached = {}
        from ..readers.oligo_tables import read_oligo_table
        for row in read_oligo_table(table, columns=("sample", oligo_col)):
            sample = row.get("sample", "").strip()
            oligo = row.get(oligo_col, "").strip()
            if not sample or not oligo:
                raise StepIncompleteError(f"step_3: blank sample or oligo in {table}")
            if sample in cached and cached[sample] != oligo:
                raise StepIncompleteError(f"step_3: conflicting oligo assignments for {sample}")
            cached[sample] = oligo
        _ASSIGNMENTS[key] = cached
    return cached


def _combine(config: WGSConfig) -> list[Path]:
    """Fan-in: every sample, then the cohort table and the integration summary."""
    records = manifest_samples(config)
    for record in records:
        _run_one_sample(config, record.sample, record=record)

    header, rows = collect_shards(
        records, lambda r: shard_path(config, r.sample), label="step_3 sv/plasmid"
    )
    out = write_tsv(table_path(config), header, rows)
    n = require_row_count(out, len(rows), len(records), label="step_3 sv/plasmid")

    summary = _integration_summary(rows)
    integration_out = integration_table_path(config)
    with atomic(integration_out) as tmp:
        write_plasmid_integration_summary(summary, tmp)

    with_evidence = sum(1 for row in summary if row["total_count"] is not None and row["total_count"] > 0)
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
            "min_mapq": config.sv_plasmid.min_mapq if config.sv_plasmid.min_mapq is not None else 20,
            "min_anchor": config.sv_plasmid.min_anchor if config.sv_plasmid.min_anchor is not None else 1,
            "integration_count_unit": "qualifying alignment records, including duplicates",
            "blacklist_reference": str(config.sv_plasmid.blacklist_reference) if config.sv_plasmid.blacklist_reference is not None else None,
            "junction_distance": (
                "derived per contig via integrated_plasmid_junction"
                if config.sv_plasmid.junction_distance is None
                else int(config.sv_plasmid.junction_distance)
            ),
        },
        counts={
            "samples": n,
            "samples_with_integration_evidence": with_evidence,
            "half_a_run": sum(row["half_a_status"] == "run" for row in summary),
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
            "upstream_count": int(fields[idx["upstream_count"]]) if fields[idx["upstream_count"]] else None,
            "downstream_count": int(fields[idx["downstream_count"]]) if fields[idx["downstream_count"]] else None,
            "total_count": int(fields[idx["total_count"]]) if fields[idx["total_count"]] else None,
            "half_a_status": fields[idx["half_a_status"]],
        }
        for fields in (row.split("\t") for row in rows)
    ]
    summary.sort(key=lambda r: -1 if r["total_count"] is None else int(r["total_count"]), reverse=True)
    return summary
