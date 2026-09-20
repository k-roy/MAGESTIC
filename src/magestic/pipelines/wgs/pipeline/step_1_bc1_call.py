"""
Step 1 — de novo bc1 calling from WGS alignments, one task per sample.

Orchestration over TRACE's `trace_crispr.wgs.bc1_consensus` (the algorithm layer,
`ARCHITECTURE_POLICY` §11), reached through `core/` for everything with a schema.
Replaces the ad-hoc `call_barcode_from_bam.pl` / `filter_barcode_region.pl` sweeps.

🔴 **TRACE's glob-based batch entry point is deliberately NOT used** — only the
per-sample library call. That entry point discovers BAMs by
`glob(bam_dir + "/*" + bam_suffix)` with a default suffix (`.sort.dedup.recal2.bam`)
that matches none of this program's BAMs, and an empty glob writes a header-only
TSV and returns 0. It also validates the plasmid track unconditionally, which fails
on a reference that legitimately has no plasmid contig. Driving `call_sample` from
the manifest instead makes both failure modes unreachable rather than merely
guarded — the row-count guard at fan-in stays anyway, because exit code 0 proves
nothing about row count (fail condition F1).

Two artifacts per sample:

* `<sample>.bc1_calls.tsv`  — TRACE's `EMIT_COLUMNS` schema, one row. Declared as
  the Snakemake rule's output and concatenated into the cohort emit table.
* `<sample>.bc1_tracks.tsv` — the per-TRACK diagnostic this module owns. Not a
  declared rule output (like Stage 3's `.anno`), but produced by both paths.

The second file exists because `call_sample` collapses a sample to one winning
track: the plasmid track is consulted only as a *fallback*, and only when the
genome track produced no call. A `Bc1Call` row therefore can never carry both
tracks' barcodes, so the genome-vs-plasmid canonical agreement that acceptance
criterion P5 is defined over is not derivable from the emit table — its
denominator would be empty by construction. Scoring each track independently with
`extract_track` + `top_barcode` costs one extra windowed fetch per track and makes
P5 evaluable, with its denominator stated.

Author: Kevin R. Roy
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

import pysam
from trace_crispr.wgs import bc1_consensus

from ..config import BC1_SUBDIR, WGSConfig
from ..core.emit import write_bc1_calls
from ..core.locus_config import (
    LOCUS_CONFIG_KEYS,
    build_tracks,
    validate_tracks,
    write_locus_sidecar,
)
from ..core.manifest import WGSSample
from ._common import (
    atomic,
    cohort_table,
    collect_shards,
    manifest_samples,
    record_for,
    require_indexed_bam,
    require_row_count,
    shard_dir,
    write_provenance,
    write_tsv,
)

logger = logging.getLogger(__name__)

STEP = "step_1_bc1_call"
SHARD_SUFFIX = ".bc1_calls.tsv"
TRACKS_SHARD_SUFFIX = ".bc1_tracks.tsv"
TABLE_NAME = "bc1_calls.tsv"
TRACKS_TABLE_NAME = "bc1_tracks.tsv"
SIDECAR_NAME = "barcode_locus.sidecar"

#: Per-track diagnostic schema. Deliberately a SEPARATE file rather than extra
#: columns on the emit table: `core/emit.py` takes `EMIT_COLUMNS` from TRACE
#: verbatim so a schema change surfaces as an import error, not column drift.
TRACKS_COLUMNS = (
    "sample",
    "genome_barcode",
    "genome_support",
    "genome_caller_confidence",
    "plasmid_barcode",
    "plasmid_support",
    "plasmid_caller_confidence",
    "both_tracks",
    "canonical_agree",
    "n_genome_bc",
    "n_plasmid_bc",
    "call_track",
    "caller_confidence",
    "confidence",
    "oligos",
)


def shard_path(config: WGSConfig, sample: str) -> Path:
    """Per-sample emit table — the Snakemake per-sample rule's declared output."""
    return shard_dir(config, BC1_SUBDIR) / f"{sample}{SHARD_SUFFIX}"


def tracks_shard_path(config: WGSConfig, sample: str) -> Path:
    """Per-sample per-track diagnostic (the P5 evidence)."""
    return shard_dir(config, BC1_SUBDIR) / f"{sample}{TRACKS_SHARD_SUFFIX}"


def table_path(config: WGSConfig) -> Path:
    """Cohort emit table."""
    return cohort_table(config, TABLE_NAME)


def tracks_table_path(config: WGSConfig) -> Path:
    """Cohort per-track diagnostic table."""
    return cohort_table(config, TRACKS_TABLE_NAME)


def locus_sidecar_path(config: WGSConfig) -> Path:
    """The flat locus sidecar TRACE's `load_locus_config` can actually parse."""
    return config.paths.locus_config / SIDECAR_NAME


def ensure_locus_sidecar(config: WGSConfig) -> Path:
    """Write the locus sidecar (idempotent, deterministic, race-safe).

    Written atomically because every per-sample task calls this, and two array
    tasks writing the same file concurrently must not produce a torn read. The
    bytes are deterministic (fixed key order, no timestamp), so §4 parity holds.

    All six `LOCUS_CONFIG_KEYS` are required even for a genome-only run —
    `tracks_from_config` indexes the plasmid keys unconditionally — so a config
    missing them raises `LocusConfigError` here rather than silently building a
    track from an empty junction sequence.
    """
    if config.locus is None:
        raise ValueError(
            "config.locus is required for Stage 1: the barcode locus and its flanking "
            "junctions cannot be derived from a reference `.fai`, which carries only "
            "contig names and lengths."
        )
    values = {key: getattr(config.locus, key) or "" for key in LOCUS_CONFIG_KEYS}
    out = locus_sidecar_path(config)
    with atomic(out) as tmp:
        write_locus_sidecar(tmp, values)
    return out


def load_bridge(config: WGSConfig) -> bc1_consensus.BridgeIndex | None:
    """Build the bc1 -> oligo bridge index, or `None` for a bridge-less run.

    Bridge-less is a legitimate configuration (the §7 smoke is one). Its
    pre-registered consequence is `bridge_tier == "NONE"` on every row, which
    `_combine_confidence` folds into `confidence == "low"` — expected, not a
    defect, and the reason the per-track diagnostic records `caller_confidence`
    separately from the combined `confidence`.
    """
    table = config.bc1.bridge_table
    if table is None:
        logger.info("step_1: no bridge table configured — running BRIDGE-LESS")
        return None
    path = Path(table)
    if not path.exists():
        raise FileNotFoundError(f"step_1: bridge table not found: {path}")
    bridge = bc1_consensus.BridgeIndex.from_tsv(
        str(path),
        bc1_col=config.bc1.bridge_bc1_col,
        oligo_col=config.bc1.bridge_oligo_col,
        count_col=config.bc1.bridge_count_col,
        table_revcomp=config.bc1.bridge_table_revcomp,
    )
    logger.info("step_1: bridge index loaded (%d barcodes) from %s", len(bridge), path)
    return bridge


def run_step_1(config: WGSConfig, sample: str | None = None) -> Path:
    """Call bc1. With `sample`, one sample; without, every sample + the cohort table.

    Returns the path written — the per-sample shard, or the cohort emit table.
    """
    if sample is not None:
        return _run_one_sample(config, sample)
    return _combine(config)


def _run_one_sample(config: WGSConfig, sample: str) -> Path:
    """One sample. Idempotent: every artifact is written atomically, so a file that
    exists is complete and this skip-check is safe under a requeued array task."""
    out = shard_path(config, sample)
    tracks_out = tracks_shard_path(config, sample)
    if out.exists() and tracks_out.exists():
        logger.info("step_1: %s already complete — skipping", sample)
        return out

    record = record_for(config, sample)
    bam_path = require_indexed_bam(record)

    sidecar = ensure_locus_sidecar(config)
    genome_track, plasmid_track = build_tracks(
        sidecar,
        amplicon_left_window=config.locus.amplicon_left_window,
        amplicon_right_window=config.locus.amplicon_right_window,
    )
    use_plasmid = config.bc1.use_plasmid and config.locus.has_plasmid_track
    tracks = [genome_track, plasmid_track] if use_plasmid else [genome_track]
    validate_tracks(tracks, record.reference)

    bridge = load_bridge(config)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        call = bc1_consensus.call_sample(
            bam, sample, genome_track, plasmid_track if use_plasmid else None, bridge
        )
        genome_call = bc1_consensus.extract_track(bam, genome_track)
        plasmid_call = bc1_consensus.extract_track(bam, plasmid_track) if use_plasmid else None

    genome_top = bc1_consensus.top_barcode(genome_call)
    plasmid_top = bc1_consensus.top_barcode(plasmid_call) if plasmid_call is not None else None
    genome_bcs = bc1_consensus.confident_barcodes(genome_call, bridge)
    plasmid_bcs = (
        bc1_consensus.confident_barcodes(plasmid_call, bridge) if plasmid_call is not None else []
    )

    with atomic(out) as tmp:
        write_bc1_calls([call], tmp)
    write_tsv(
        tracks_out,
        "\t".join(TRACKS_COLUMNS),
        ["\t".join(_tracks_row(sample, call, genome_top, plasmid_top, genome_bcs, plasmid_bcs))],
    )
    logger.info(
        "step_1: %s track=%s barcode=%r support=%d n_genome_bc=%d n_plasmid_bc=%d",
        sample,
        call.track,
        call.barcode,
        call.support,
        len(genome_bcs),
        len(plasmid_bcs),
    )
    return out


def _tracks_row(
    sample: str,
    call: bc1_consensus.Bc1Call,
    genome_top: bc1_consensus.TopBarcode | None,
    plasmid_top: bc1_consensus.TopBarcode | None,
    genome_bcs: Sequence[bc1_consensus.BarcodeCall],
    plasmid_bcs: Sequence[bc1_consensus.BarcodeCall],
) -> list[str]:
    """Assemble the per-track diagnostic row.

    `TopBarcode.barcode` is already canonical on both tracks — `extract_barcode`
    applies `TrackConfig.to_canonical`, which reverse-complements the genome
    integrant and leaves the plasmid cassette alone — so `canonical_agree` is a
    direct string comparison and no orientation is re-derived here. Orientation is
    TRACE's call; this file only records the outcome, which is what P3b asks for.
    """
    genome_bc = genome_top.barcode if genome_top else ""
    plasmid_bc = plasmid_top.barcode if plasmid_top else ""
    both = bool(genome_bc) and bool(plasmid_bc)
    oligos: list[str] = []
    for bc in list(genome_bcs) + list(plasmid_bcs):
        if bc.oligo and bc.oligo not in oligos:
            oligos.append(bc.oligo)
    return [
        sample,
        genome_bc,
        str(genome_top.support if genome_top else 0),
        bc1_consensus.caller_confidence(genome_top) if genome_top else "none",
        plasmid_bc,
        str(plasmid_top.support if plasmid_top else 0),
        bc1_consensus.caller_confidence(plasmid_top) if plasmid_top else "none",
        str(int(both)),
        str(int(both and genome_bc == plasmid_bc)),
        str(len(genome_bcs)),
        str(len(plasmid_bcs)),
        call.track,
        call.caller_confidence,
        call.confidence,
        ";".join(oligos),
    ]


def _combine(config: WGSConfig) -> Path:
    """Fan-in: every sample, then the cohort tables.

    Computes any shard it does not already find, so the DIRECT path (A) is a single
    call while the SNAKEMAKE path (B) — whose per-sample rules already produced them
    — becomes a pure, ordered concatenation. One code path per stage is what makes
    §4 byte-identity structural rather than a convention.
    """
    records = manifest_samples(config)
    for record in records:
        _run_one_sample(config, record.sample)

    header, rows = collect_shards(
        records, lambda r: shard_path(config, r.sample), label="step_1 bc1 emit"
    )
    out = write_tsv(table_path(config), header, rows)
    n = require_row_count(out, len(rows), len(records), label="step_1 bc1 emit")

    tracks_header, tracks_rows = collect_shards(
        records, lambda r: tracks_shard_path(config, r.sample), label="step_1 bc1 tracks"
    )
    tracks_out = write_tsv(tracks_table_path(config), tracks_header, tracks_rows)
    require_row_count(tracks_out, len(tracks_rows), len(records), label="step_1 bc1 tracks")

    barcode_idx = bc1_consensus.EMIT_COLUMNS.index("barcode")
    called = sum(1 for row in rows if row.split("\t")[barcode_idx])
    both_tracks, agree = _track_agreement(tracks_rows)
    logger.info(
        "step_1: %d samples, %d with a bc1 call (%.1f%% of %d); "
        "both-track subset n=%d, canonical agreement %s",
        n,
        called,
        100.0 * called / n if n else 0.0,
        n,
        both_tracks,
        f"{100.0 * agree / both_tracks:.1f}%" if both_tracks else "n/a (empty subset)",
    )

    write_provenance(
        STEP,
        out,
        inputs={
            "manifest": Path(config.paths.manifest),
            "locus_sidecar": locus_sidecar_path(config),
            **_reference_inputs(records),
        },
        parameters={
            "use_plasmid": config.bc1.use_plasmid,
            "bridge_table": str(config.bc1.bridge_table) if config.bc1.bridge_table else None,
            "bam_discovery": "manifest (never a bam_suffix glob)",
        },
        counts={
            "samples": n,
            "bc1_called": called,
            "both_tracks_denominator": both_tracks,
            "both_tracks_canonical_agree": agree,
        },
    )
    return out


def _track_agreement(tracks_rows: Sequence[str]) -> tuple[int, int]:
    """`(both-tracks subset size, canonical agreements)` — the P5 numbers.

    Both are returned, logged and written to PROVENANCE, because an agreement
    percentage over an unstated denominator is the shape of a wrong answer.
    """
    both_idx = TRACKS_COLUMNS.index("both_tracks")
    agree_idx = TRACKS_COLUMNS.index("canonical_agree")
    both = 0
    agree = 0
    for row in tracks_rows:
        fields = row.split("\t")
        if fields[both_idx] == "1":
            both += 1
            agree += int(fields[agree_idx] == "1")
    return both, agree


def _reference_inputs(records: Sequence[WGSSample]) -> dict[str, Path]:
    """One provenance entry per distinct reference in the cohort (D1).

    A cohort does not have "a" reference: the rearray cohort alone spans two, so
    recording a single `genome:` here would be a false provenance claim.
    """
    refs = sorted({Path(r.reference) for r in records})
    return {f"reference_{i}": ref for i, ref in enumerate(refs, start=1)}
