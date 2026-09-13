"""
Step 4 — per-cell edit reconciliation (bc1 multiplicity x designed-site edit VAF).

Orchestration over TRACE's `trace_crispr.wgs.reconcile`, reached through
`core/reconcile_driver.py`. The driver exists because `reconcile.run` is the one
TRACE entry point that is not library-callable with keyword arguments — it takes
an `argparse.Namespace` — and synthesizing a Namespace to call a library would
re-create the inverted dependency `ARCHITECTURE_POLICY` §11 forbids, through the
back door. No threshold is re-implemented here or there: `tau_full` / `tau_wt` /
`min_cov` / `min_base_qual` stay in TRACE, which is why this step passes none of
them.

**Stage 1 -> Stage 4 is the pipeline's one hard dependency, and TRACE's defaults
are not neutral.** `reconcile_cell(..., n_genome_bc=0, n_plasmid_bc=0,
bc_confidence="high")` defaults to a *wrong* answer rather than an absent one:
>=2 genome-track bc1 is the primary contamination signal, so a cell scored with
`n_genome_bc=0` can never report contamination, and `bc_confidence="high"` marks
every call trustworthy. `core.reconcile_driver.reconcile_sample` derives all three
from this sample's own BAM, so this step does NOT read Stage 1's numbers — it
only requires Stage 1's per-track shard to exist, which is the ordering edge the
Snakemake rule declares as `bc1_tracks=BC1_DIR/{sample}.bc1_tracks.tsv`.

🔴 **KNOWN LIMITATION, driver-side, and it moves a number that will be read as a
result.** `reconcile_sample` computes `bc_confidence` unconditionally from
`BarcodeCall.confidence`, which is the COMBINED caller x bridge label. Measured on
the cluster: `_combine_confidence("high", BRIDGE_NONE) == "low"`. So in a
BRIDGE-LESS cohort — the configuration the §7 acceptance smoke uses — every
barcode is `"low"`, `bc_confidence` is therefore `"low"` for every cell, and
`reconcile_cell` stamps `confidence="low"` on the whole cohort because no bridge
was supplied, NOT because the barcode evidence was weak. The fix is a
`bc_confidence=` parameter on `reconcile_sample` (defaulting to today's behavior)
so a caller can pass `caller_confidence` instead; re-deriving the rule here would
put a classification rule in the orchestration layer, which spec D3 forbids. Until
that lands, read a bridge-less Stage-4 `confidence` column as "bridge absent",
never as "call unreliable".

Author: Kevin R. Roy
"""

from __future__ import annotations

import logging
from pathlib import Path

import pysam
from trace_crispr.wgs import reconcile

from ..config import RECONCILE_SUBDIR, WGSConfig
from ..core.emit import write_reconciliations
from ..core.locus_config import build_tracks, validate_tracks
from ..core.reconcile_driver import load_oligo_windows, reconcile_sample
from . import step_1_bc1_call
from ._common import (
    StepIncompleteError,
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

STEP = "step_4_reconcile"
SHARD_SUFFIX = ".reconcile.tsv"
TABLE_NAME = "reconcile_calls.tsv"


def shard_path(config: WGSConfig, sample: str) -> Path:
    """Per-sample shard — the Snakemake per-sample rule's declared output."""
    return shard_dir(config, RECONCILE_SUBDIR) / f"{sample}{SHARD_SUFFIX}"


def table_path(config: WGSConfig) -> Path:
    """Cohort reconciliation table."""
    return cohort_table(config, TABLE_NAME)


def run_step_4(config: WGSConfig, sample: str | None = None) -> Path:
    """Reconcile edits. With `sample`, one cell; without, every cell + the table."""
    if sample is not None:
        return _run_one_sample(config, sample)
    return _combine(config)


def _run_one_sample(config: WGSConfig, sample: str) -> Path:
    """One cell. Idempotent: the shard is written atomically."""
    out = shard_path(config, sample)
    if out.exists():
        logger.info("step_4: %s already complete — skipping", sample)
        return out

    record = record_for(config, sample)
    bam_path = require_indexed_bam(record)
    _require_stage_1(config, sample)

    if config.reconcile.oligo_windows is None:
        raise StepIncompleteError(
            "step_4: reconcile.oligo_windows is required — without a donor-window table "
            "no site is scored, and a cell with zero scored sites classifies LOW_CONF. "
            "That is 'not measured', which must never be written as 'no edit'."
        )
    oligo_windows = load_oligo_windows(
        str(config.reconcile.oligo_windows),
        oligo_col=config.columns.oligo_col,
        chrom_col=config.columns.chrom_col,
        start_col=config.columns.start_col,
        end_col=config.columns.end_col,
    )

    sidecar = step_1_bc1_call.ensure_locus_sidecar(config)
    genome_track, plasmid_track = build_tracks(
        sidecar,
        amplicon_left_window=config.locus.amplicon_left_window,
        amplicon_right_window=config.locus.amplicon_right_window,
    )
    use_plasmid = config.bc1.use_plasmid and config.locus.has_plasmid_track
    validate_tracks([genome_track, plasmid_track] if use_plasmid else [genome_track],
                    record.reference)
    bridge = step_1_bc1_call.load_bridge(config)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
        str(record.reference)
    ) as fasta:
        cell = reconcile_sample(
            bam,
            fasta,
            sample,
            genome_track=genome_track,
            plasmid_track=plasmid_track if use_plasmid else None,
            oligo_windows=oligo_windows,
            bridge=bridge,
        )

    with atomic(out) as tmp:
        write_reconciliations([cell], tmp)
    logger.info(
        "step_4: %s label=%s confidence=%s n_full=%d n_partial=%d n_nocov=%d "
        "n_genome_bc=%d n_plasmid_bc=%d",
        sample,
        cell.label,
        cell.confidence,
        cell.n_full,
        cell.n_partial,
        cell.n_nocov,
        cell.n_genome_bc,
        cell.n_plasmid_bc,
    )
    return out


def _require_stage_1(config: WGSConfig, sample: str) -> Path:
    """Assert Stage 1 ran for this sample, checking the same shard the rule declares.

    `core.reconcile_driver.reconcile_sample` re-derives the barcode calls from the
    BAM, so this is an ORDERING guard rather than a data read: it keeps path (A)
    honest about the edge `wgs_step_4_reconcile` declares as
    `bc1_tracks=BC1_DIR/{sample}.bc1_tracks.tsv`, and it means a cohort whose
    Stage 1 failed cannot quietly produce a full Stage 4 table.
    """
    tracks = step_1_bc1_call.tracks_shard_path(config, sample)
    if not tracks.exists():
        raise StepIncompleteError(
            f"step_4: {sample} has no Stage-1 per-track output ({tracks}). Run step 1 "
            "first — reconcile_cell's n_genome_bc=0 / bc_confidence='high' defaults are "
            "a wrong answer, not a missing one."
        )
    return tracks


def _combine(config: WGSConfig) -> Path:
    """Fan-in: every cell, then the cohort table."""
    records = manifest_samples(config)
    for record in records:
        _run_one_sample(config, record.sample)

    header, rows = collect_shards(
        records, lambda r: shard_path(config, r.sample), label="step_4 reconcile"
    )
    out = write_tsv(table_path(config), header, rows)
    n = require_row_count(out, len(rows), len(records), label="step_4 reconcile")

    label_idx = reconcile.EMIT_COLUMNS.index("label")
    labels: dict[str, int] = {}
    for row in rows:
        label = row.split("\t")[label_idx]
        labels[label] = labels.get(label, 0) + 1
    logger.info(
        "step_4: %d cells — %s", n, ", ".join(f"{k}={v}" for k, v in sorted(labels.items()))
    )

    write_provenance(
        STEP,
        out,
        inputs={
            "manifest": Path(config.paths.manifest),
            "oligo_windows": Path(config.reconcile.oligo_windows),
            "bc1_tracks": step_1_bc1_call.tracks_table_path(config),
        },
        parameters={
            "thresholds": "TRACE defaults (classify_site tau_full/tau_wt/min_cov "
            "and window_edit_vaf min_base_qual are the algorithm's, not ours)",
            "bridge_table": str(config.bc1.bridge_table) if config.bc1.bridge_table else None,
        },
        counts={"cells": n, **labels},
    )
    return out
