"""Stage 4: explicit intended-allele reconciliation with measured bc1 multiplicity.

Windows without intended variants are unassessed compatibility evidence. Explicit
sample assignments take precedence over bridge-derived candidate sites. Shared
windows, assignments and the bridge are loaded once per direct cohort run.
"""
from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from contextlib import ExitStack
from typing import Any
from pathlib import Path

import pysam
from trace_crispr.wgs import reconcile

from ..config import RECONCILE_SUBDIR, WGSConfig
from ..core.emit import write_reconciliations
from ..core.locus_config import build_tracks, validate_tracks
from ..core.reconcile_driver import load_oligo_windows, reconcile_sample
from . import step_1_bc1_call, step_3_sv_plasmid
from ._common import (
    StepIncompleteError,
    cache_key,
    cache_valid,
    cache_commit,
    atomic,
    cohort_table,
    collect_shards,
    manifest_samples,
    record_for,
    require_indexed_bam,
    require_row_count,
    read_table,
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


@dataclass
class _Inputs:
    windows: dict[str, Any]
    assignments: dict[str, list[str]] | None
    bridge: Any


def _load_inputs(config: WGSConfig) -> _Inputs:
    if config.reconcile.oligo_windows is None:
        raise StepIncompleteError("step_4: reconcile.oligo_windows is required; unmeasured sites are not no-edit results")
    windows = load_oligo_windows(str(config.reconcile.oligo_windows),
        oligo_col=config.columns.oligo_col, chrom_col=config.columns.chrom_col,
        start_col=config.columns.start_col, end_col=config.columns.end_col)
    assignments = None
    if config.reconcile.oligo_assignments is not None:
        assignments = {}
        path = config.reconcile.oligo_assignments
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if (not reader.fieldnames or not {"sample", "oligo_name"}.issubset(reader.fieldnames)
                    or len(reader.fieldnames) != len(set(reader.fieldnames))):
                raise ValueError(f"{path}: assignment header requires unique sample and oligo_name columns")
            for lineno, row in enumerate(reader, 2):
                if None in row or not row.get("sample") or not row.get("oligo_name"):
                    raise ValueError(f"{path}:{lineno}: malformed assignment")
                sample, oligo = row["sample"].strip(), row["oligo_name"].strip()
                if not sample or not oligo:
                    raise ValueError(f"{path}:{lineno}: empty assignment")
                if oligo not in windows:
                    raise ValueError(f"{path}:{lineno}: no donor window for assigned oligo {oligo!r}")
                if oligo not in assignments.setdefault(sample, []):
                    assignments[sample].append(oligo)
        expected = {r.sample for r in manifest_samples(config)}
        if set(assignments) != expected:
            raise ValueError(f"{path}: assignment samples do not match manifest; missing={sorted(expected-set(assignments))}, unexpected={sorted(set(assignments)-expected)}")
    return _Inputs(windows, assignments, step_1_bc1_call.load_bridge(config))


def run_step_4(config: WGSConfig, sample: str | None = None) -> Path:
    """One sample shard or complete cohort, using the same scoring code path."""
    if sample is not None:
        return _run_one_sample(config, sample)
    return _combine(config)


def uses_cleaned_reads(config: WGSConfig) -> bool:
    """One policy shared by direct execution and the workflow DAG."""
    source = config.reconcile.read_source
    return source == "donor_cleaned" or (source == "auto" and
        config.sv_plasmid.blacklist is not None and
        config.reconcile.oligo_assignments is not None and
        config.reconcile.oligo_windows is not None)


def _allele_bam_path(config: WGSConfig, sample: str, raw: Path) -> Path:
    if not uses_cleaned_reads(config):
        return raw
    step_3_sv_plasmid.run_step_3(config, sample)
    rows = list(read_table(step_3_sv_plasmid.shard_path(config, sample)))
    if len(rows) != 1 or rows[0].get("half_b_status") != step_3_sv_plasmid.HALF_B_RUN:
        raise StepIncompleteError(f"step_4: {sample}: donor-cleaned evidence requires measured Stage 3 Half B")
    return step_3_sv_plasmid.cleaned_bam_path(config, sample)


def _run_one_sample(config: WGSConfig, sample: str, *, inputs: _Inputs | None = None) -> Path:
    out = shard_path(config, sample)
    record = record_for(config, sample)
    bam_path = require_indexed_bam(record)
    _require_stage_1(config, sample)

    inputs = inputs if inputs is not None else _load_inputs(config)
    allele_path = _allele_bam_path(config, sample, bam_path)
    dependencies = {"bc1_tracks": step_1_bc1_call.tracks_shard_path(config, sample)}
    if uses_cleaned_reads(config):
        dependencies.update(allele_bam=allele_path,
            allele_index=step_3_sv_plasmid.cleaned_bam_index_path(config, sample, record=record),
            filtering=step_3_sv_plasmid.shard_path(config, sample))
    key = cache_key(config, STEP, record, inputs=dependencies)
    if cache_valid([out], key):
        logger.info("step_4: %s verified complete — skipping", sample)
        return out
    sidecar = step_1_bc1_call.ensure_locus_sidecar(config)
    genome_track, plasmid_track = build_tracks(
        sidecar,
        amplicon_left_window=config.locus.amplicon_left_window,
        amplicon_right_window=config.locus.amplicon_right_window,
        min_mapq=config.bc1.min_mapq,
        use_plasmid=config.bc1.use_plasmid,
    )
    use_plasmid = config.bc1.use_plasmid and config.locus.has_plasmid_track
    validate_tracks([genome_track, plasmid_track] if use_plasmid else [genome_track],
                    record.reference, amplicon_window_sequences=(
                        config.locus.amplicon_left_window_sequence,
                        config.locus.amplicon_right_window_sequence))

    with ExitStack() as handles:
        bam = handles.enter_context(pysam.AlignmentFile(str(bam_path), "rb"))
        fasta = handles.enter_context(pysam.FastaFile(str(record.reference)))
        allele_bam = handles.enter_context(pysam.AlignmentFile(str(allele_path), "rb")) if allele_path != bam_path else bam
        allele_bam.check_index()
        cell = reconcile_sample(
            bam,
            fasta,
            sample,
            genome_track=genome_track,
            plasmid_track=plasmid_track if use_plasmid else None,
            oligo_windows=inputs.windows,
            bridge=inputs.bridge,
            min_mapq=getattr(config.reconcile, "min_mapq", None),
            assigned_oligos=inputs.assignments[sample] if inputs.assignments is not None else None,
            allele_bam=allele_bam,
            read_source_flags=("DONOR_CLEANED_READS" if uses_cleaned_reads(config) else "RAW_DONOR_READS",),
        )

    with atomic(out) as tmp:
        write_reconciliations([cell], tmp)
    cache_commit([out], key)
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
    """Verify or refresh Stage 1 through its input/output cache contract."""
    step_1_bc1_call.run_step_1(config, sample)
    tracks = step_1_bc1_call.tracks_shard_path(config, sample)
    if not tracks.exists():
        raise StepIncompleteError(f"step_4: {sample} has no verified Stage-1 per-track output ({tracks})")
    return tracks


def _combine(config: WGSConfig) -> Path:
    """Fan-in: every cell, then the cohort table."""
    records = manifest_samples(config)
    inputs = _load_inputs(config)
    step_1_bc1_call.run_step_1(config)
    for record in records:
        _run_one_sample(config, record.sample, inputs=inputs)

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
            **({"oligo_assignments": Path(config.reconcile.oligo_assignments)}
               if config.reconcile.oligo_assignments else {}),
        },
        parameters={
            "thresholds": "TRACE defaults (classify_site tau_full/tau_wt/min_cov "
            "and score_edit_site min_base_qual are the algorithm's, not ours)",
            "bridge_table": str(config.bc1.bridge_table) if config.bc1.bridge_table else None,
            "allele_read_source_policy": config.reconcile.read_source,
            "allele_reads": "Stage 3 donor-cleaned BAM" if uses_cleaned_reads(config) else "original BAM; donor origin unresolved when plasmid bc1 observed",
            "barcode_reads": "original BAM",
            "variant_min_mapq": (config.reconcile.min_mapq if getattr(config.reconcile, "min_mapq", None) is not None
                                  else "TRACE default 20; MAPQ 255 always rejected"),
        },
        counts={"cells": n, **labels},
    )
    return out
