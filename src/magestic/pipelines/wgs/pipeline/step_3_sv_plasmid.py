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

from ..config import (
    CONTAINED_AUTO,
    CONTAINED_DISCARD,
    CONTAINED_KEEP,
    SV_PLASMID_SUBDIR,
    WGSConfig,
    parse_interval,
)
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
from ._common import (
    StepIncompleteError,
    atomic,
    cache_key,
    cache_valid,
    cache_commit,
    cohort_table,
    collect_shards,
    manifest_samples,
    preflight_provenance,
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
#: 🔴 NN 251 defects 7 and 13 add two status columns, and the 205 lane a third
#: (`half_a_status`: `run` / `not_run_no_integrated_contigs`, merged in planning/268;
#: it restates `integration_status` in that lane's vocabulary). `integration_status`
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
    "half_a_status",
    "half_b_status",
    "read_count_mid_point",
    "frags_genome",
    "frags_plasmid",
    "frags_chimeric",
    # planning/328 (Kevin, 2026-10-30): the contained-fragment policy applied to this sample
    # and the plasmid-load evidence behind it. Appended, so every earlier index is unchanged.
    "contained_policy",
    "plasmid_backbone_reads",
    "plasmid_expected_reads_per_copy",
    "plasmid_copies_upper95",
    "donor_backbone_pairs",
    "frags_contained_kept",
)

#: The only SV status this pipeline can honestly emit. See `core/sv_status.py` for
#: why the refusal was chosen deliberately over a half-built detector.
SV_STATUS_VALUE = SV_UNMEASURED

HALF_B_NO_WINDOW_TABLE = "not_run_no_oligo_window_table"
HALF_B_NO_ASSIGNMENT = "not_run_no_oligo_assignment"
HALF_B_WINDOW_UNKNOWN = "not_run_oligo_window_unknown"
HALF_B_RUN = "run"
HALF_A_NO_CONTIGS = "not_run_no_integrated_contigs"

#: `contained_policy` values: what Half B did with donor-contained fragments, and why.
POLICY_DISCARD = "discard"                        # the span rule (configured)
POLICY_KEEP_PURGED = "keep_purged"                # auto: this sample passed plasmids_purged
POLICY_DISCARD_NOT_PURGED = "discard_not_purged"  # auto: this sample failed plasmids_purged
POLICY_KEEP_DECLARED = "keep_declared"            # keep: dataset declared purged, no plasmid seen here
POLICY_DISCARD_PLASMID_SEEN = "discard_plasmid_seen"  # keep: but this sample shows plasmid


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
    preflight_provenance(STEP)
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
        # 🔴 NN 251: a shard under an older `SHARD_COLUMNS` would be read with shifted
        # indices. The input/output cache already regenerates one whose inputs or code
        # changed, so this is the invariant check on the path that reuses a shard.
        _require_current_shard_schema(out, sample)
        logger.info("step_3: %s complete for current inputs — skipping", sample)
        return out

    contigs = _integrated_contigs(bam_path)
    # 🔴 NN 251 defect 13: a BAM with no `integrated_plasmid_` contig cannot be measured
    # for integration, and TRACE now refuses that call outright rather than returning a
    # 0/0 that reads like a measured negative. `integration_measurement` does not make
    # the call then, and writes "" rather than 0; `half_a_status` says the same thing in
    # the 205 lane's vocabulary, so both lanes' readers find their column.
    measure = integration_measurement(
        str(bam_path),
        sample=sample,
        integrated_contigs=contigs,
        junction_distance=_junction_distances(config, contigs, Path(record.reference)),
        **{k: getattr(config.sv_plasmid, k) for k in ("min_mapq", "min_anchor")
           if getattr(config.sv_plasmid, k) is not None},
    )
    half_a_status = "run" if measure.measured else HALF_A_NO_CONTIGS
    status, anno, contained = _run_half_b(config, record)
    blank = ""

    # 🔴 NN 255 F2: the refusal was applied to the provenance MAPPING only, never to a
    # written HEADER, so an SV-rate COLUMN reached a file unrefused. A column name is
    # published exactly as readily as a count.
    require_no_sv_rate_in_columns(SHARD_COLUMNS, source="step_3 shard header")
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
                    blank if not measure.measured else measure.upstream_count,
                    blank if not measure.measured else measure.downstream_count,
                    blank if not measure.measured else measure.total,
                    measure.n_integrated_contigs,
                    measure.status,
                    SV_STATUS_VALUE,
                    half_a_status,
                    status,
                    *anno,
                    *contained,
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
    outputs = [out]
    if status == HALF_B_RUN:
        cleaned = cleaned_bam_path(config, sample)
        outputs.extend([cleaned, cleaned_bam_index_path(config, sample, record=record), Path(str(cleaned) + ".target_filtering.anno")])
    cache_commit(outputs, key)
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


def _contained_policy(config: WGSConfig, record: WGSSample, window: tuple[str, int, int]) -> tuple[bool, list[str]]:
    """Decide whether Half B keeps donor-contained fragments for this sample (planning/328).

    Kevin R. Roy, 2026-10-30: "count contained read pairs when the backbone is empty. This is safe when we
    know that we have removed plasmids, so it needs to be a per clone/per dataset option for trace".
    Returns `(keep_contained, [policy, backbone reads, expected reads per copy, copies upper limit,
    donor-backbone pairs])`; `frags_contained_kept` is appended by the caller after cleaning.
    """
    sv = config.sv_plasmid
    if sv.contained_fragments == CONTAINED_DISCARD:
        return False, [POLICY_DISCARD, "", "", "", ""]
    load = sv_plasmid.plasmid_load(
        str(record.bam), sv.donor_plasmid_contigs,
        mask=[parse_interval(m, "sv_plasmid.plasmid_mask") for m in (sv.plasmid_mask or [])],
        donor_windows=[window],
        nuclear_contigs=sv.nuclear_contigs,
        flexible_bases=sv.flexible_bases,
    )
    evidence = [str(load.plasmid_reads), f"{load.expected_reads_per_copy:.3f}",
                f"{load.copies_upper95:.6g}", str(load.donor_backbone_pairs)]
    if sv.contained_fragments == CONTAINED_AUTO:
        keep = sv_plasmid.plasmids_purged(load, max_copies=sv.purged_max_copies)
        return keep, [POLICY_KEEP_PURGED if keep else POLICY_DISCARD_NOT_PURGED, *evidence]
    if sv.contained_fragments == CONTAINED_KEEP:
        seen = load.donor_backbone_pairs > 0 or load.copies > sv.purged_max_copies
        if seen:
            logger.warning(
                "step_3: %s shows plasmid (%d donor-backbone pair(s), %.3g copies per genome) although the "
                "dataset is declared purged; its contained fragments are discarded",
                record.sample, load.donor_backbone_pairs, load.copies,
            )
        return not seen, [POLICY_DISCARD_PLASMID_SEEN if seen else POLICY_KEEP_DECLARED, *evidence]
    raise StepIncompleteError(f"step_3: unknown contained_fragments policy {sv.contained_fragments!r}")


def _run_half_b(config: WGSConfig, record: WGSSample) -> tuple[str, list[str], list[str]]:
    """Clean donor reads for one sample.

    Returns `(status, the .anno fields, the contained-policy fields)`, the last matching the six
    `contained_policy` .. `frags_contained_kept` shard columns.
    """
    blank = ["", "", "", ""]
    blank_policy = ["", "", "", "", "", ""]
    if config.reconcile.oligo_windows is None:
        return HALF_B_NO_WINDOW_TABLE, blank, blank_policy

    oligo = _assigned_oligo(config, record.sample)
    if oligo is None:
        return HALF_B_NO_ASSIGNMENT, blank, blank_policy
    windows = _donor_windows(config)
    window = windows.get(oligo)
    if window is None:
        # `_donor_windows` refuses a malformed or conflicting row, so a missing window
        # here is an oligo genuinely absent from the table, never a dropped row.
        logger.warning("step_3: %s oligo %s has no window in %s",
                       record.sample, oligo, config.reconcile.oligo_windows)
        return HALF_B_WINDOW_UNKNOWN, blank, blank_policy

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
    if config.sv_plasmid.dup_mode == "standard" and config.sv_plasmid.blacklist_reference is None:
        raise StepIncompleteError("step_3: standard Half B requires sv_plasmid.blacklist_reference to establish coordinate identity")
    require_binary(config.sv_plasmid.samtools)
    final_bam = cleaned_bam_path(config, record.sample)
    final_index = cleaned_bam_index_path(config, record.sample, record=record)
    final_bam.parent.mkdir(parents=True, exist_ok=True)
    keep_contained, policy = _contained_policy(config, record, window)
    extra = {}
    if keep_contained:
        extra["keep_contained"] = True
    if config.sv_plasmid.plasmid_mates_are_template:
        extra["plasmid_contigs"] = list(config.sv_plasmid.donor_plasmid_contigs)
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
            **extra,
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
    return HALF_B_RUN, anno, [*policy, str(result.frags_contained_kept) if keep_contained else "0"]


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
    contained = _pooled_contained_evidence(rows)
    if contained:
        logger.info("step_3: contained fragments %s", contained)
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
            "min_mapq": config.sv_plasmid.min_mapq if config.sv_plasmid.min_mapq is not None else 20,
            "min_anchor": config.sv_plasmid.min_anchor if config.sv_plasmid.min_anchor is not None else 1,
            "integration_count_unit": "qualifying alignment records, including duplicates",
            "blacklist_reference": str(config.sv_plasmid.blacklist_reference) if config.sv_plasmid.blacklist_reference is not None else None,
            "junction_distance": (
                "derived per contig via integrated_plasmid_junction"
                if config.sv_plasmid.junction_distance is None
                else int(config.sv_plasmid.junction_distance)
            ),
            "contained_fragments": config.sv_plasmid.contained_fragments,
            "donor_plasmid_contigs": list(config.sv_plasmid.donor_plasmid_contigs),
            "plasmid_mask": list(config.sv_plasmid.plasmid_mask or []),
            "purged_max_copies": config.sv_plasmid.purged_max_copies,
            "plasmid_mates_are_template": config.sv_plasmid.plasmid_mates_are_template,
        },
        counts=_provenance_counts(
            {
                "samples": n,
                "samples_measured_for_integration": len(measured),
                "samples_not_measured_for_integration": unmeasured,
                "samples_with_integration_evidence": with_evidence,
                "half_a_run": sum(row["half_a_status"] == "run" for row in summary),
                "half_b_run": half_b_run,
                "sv_status": SV_STATUS_VALUE,
                **({"contained_fragments": contained} if contained else {}),
            }
        ),
    )
    return [out, integration_out]


def _pooled_contained_evidence(rows: Sequence[str]) -> dict[str, object]:
    """Cohort view of the contained-fragment policy: samples per policy and the POOLED plasmid bound.

    The pooled limit (all backbone reads over all expected reads per copy) is the evidence a dataset-level
    `keep` should rest on: one shallow clone cannot bound its own copy number, a cohort can.
    """
    idx = {name: i for i, name in enumerate(SHARD_COLUMNS)}
    policies: dict[str, int] = {}
    reads = 0
    expected = 0.0
    kept = 0
    for row in rows:
        f = row.split("\t")
        policy = f[idx["contained_policy"]]
        if not policy:
            continue
        policies[policy] = policies.get(policy, 0) + 1
        if f[idx["plasmid_backbone_reads"]]:
            reads += int(f[idx["plasmid_backbone_reads"]])
            expected += float(f[idx["plasmid_expected_reads_per_copy"]])
        if f[idx["frags_contained_kept"]]:
            kept += int(f[idx["frags_contained_kept"]])
    if not policies:
        return {}
    out: dict[str, object] = {"samples_by_policy": policies, "fragments_kept": kept}
    if expected:
        out.update({"pooled_backbone_reads": reads,
                    "pooled_copies_upper95": _poisson_upper(reads) / expected})
    return out


def _poisson_upper(n: int, tail: float = 0.025) -> float:
    """TRACE's `poisson_upper` in closed form, for the POOLED count.

    The lambda at which P(X <= n) = tail is the gamma quantile chi2.ppf(1 - tail, 2(n + 1)) / 2: the same
    exact limit (3.689 for n = 0), constant-time in n. TRACE solves it by bisection over an n-term sum, which
    is right for one sample's backbone reads and quadratic on a plasmid-positive cohort's total.
    """
    from scipy.stats import chi2

    if type(n) is not int or n < 0 or not 0 < tail < 1:
        raise ValueError("n must be a non-negative int and tail in (0, 1)")
    return float(chi2.ppf(1.0 - tail, 2 * (n + 1)) / 2.0)


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
            "half_a_status": fields[idx["half_a_status"]],
        }
        for fields in (row.split("\t") for row in rows)
    ]
    # Unmeasured samples sort last, and are NOT treated as zero.
    summary.sort(
        key=lambda r: (r["total_count"] != "", r["total_count"] or 0), reverse=True
    )
    return summary
