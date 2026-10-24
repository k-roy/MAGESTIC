"""Stage 4: per-window edit reconciliation with measured bc1 multiplicity.

Windows with intended variants are scored per allele by TRACE; windows without them
are scored from their read evidence. Explicit sample assignments take precedence over
bridge-derived candidate sites. Shared windows, assignments and the bridge are loaded
once per direct cohort run.

Merged in planning/268 with the 222 lane (NN 251/255). Kept from it: the window-table
denominators (defect 12), logged before TRACE's loader refuses a bad row; the reference-
frame check against `wt_seq` when the table carries one (defect 10); the BAM/reference
consistency check (defect 14 / F3b); and the per-sample audit sidecars that make a check
that did not run countable (F10). Scoring is per window (`core.reconcile_driver.score_sites`):
TRACE's allele-specific `score_edit_site` where the window declares `intended_variants`, and
the 222 lane's depth-gated `core.edit_profile` classifier on the read evidence where it does
not, because an intended edit is not required to assign confidence to a WGS call (Kevin,
2026-10-24). Retired: the `bc_confidence`-renamed writer; the shard is TRACE's schema.
"""
from __future__ import annotations

import csv
import json
import logging
from collections.abc import Sequence
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam
from trace_crispr.wgs import reconcile

from ..config import RECONCILE_SUBDIR, WGSConfig
from ..core.blacklist_adapter import reference_contig_lengths
from ..core.emit import write_reconciliations
from ..core.locus_config import build_tracks, validate_tracks
from ..core.manifest import (
    ReferenceUnobservedError,
    WGSSample,
    incomparable_contigs,
    require_reference_consistency,
)
from ..core.reconcile_driver import load_oligo_windows, reconcile_sample
from ..core.window_table import (
    FrameMismatchError,
    WindowAudit,
    designed_positions_by_oligo,
    load_windows_audited,
    verify_window_frame,
)
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


#: Per-sample audit of the two checks whose SKIPPING must not read as PASSING: the
#: reference-frame check (NN 255 F10) and the manifest/BAM reference-NAME check (F3b).
#:
#: 🔴 **Written to a per-sample SIDECAR on disk, not held only in these dicts.** The
#: first NN 255 attempt kept them in module globals and summed them in `_combine` — and
#: that would have shipped **zero on the real path**, twice over: `_run_one_sample`
#: returns early on an existing shard, and the Snakemake fan-out runs every sample in
#: its OWN PROCESS, so `_combine` would read a dict that nothing in that process filled.
#: A detector reading 0 because it was never asked is the exact failure class this change
#: set exists to remove, so the audit is persisted and re-read, and a sample whose
#: sidecar is MISSING is counted as missing rather than as zero.
_FRAME_AUDIT: dict[str, dict] = {}
_REFERENCE_NAME_CHECK: dict[str, str] = {}

AUDIT_SIDECAR_SUFFIX = ".reconcile.audit.json"


def audit_sidecar_path(config: WGSConfig, sample: str) -> Path:
    """Where this sample's frame / reference-name audit is persisted."""
    return shard_dir(config, RECONCILE_SUBDIR) / f"{sample}{AUDIT_SIDECAR_SUFFIX}"


def _write_audit_sidecar(config: WGSConfig, sample: str) -> Path:
    """Persist this sample's audit beside its shard, so `_combine` can read it."""
    path = audit_sidecar_path(config, sample)
    payload = {
        "sample": sample,
        "frame_check": _FRAME_AUDIT.get(sample, {}),
        "reference_name_check": _REFERENCE_NAME_CHECK.get(sample, "not_run"),
    }
    with atomic(path) as tmp:
        tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


def _read_audit_sidecars(
    config: WGSConfig, samples: Sequence[str]
) -> tuple[dict[str, dict], list[str]]:
    """`({sample: payload}, [samples with no sidecar])`.

    The missing list is returned, never silently treated as zero: a shard written before
    NN 255 (or by a task whose sidecar write failed) carries no audit, and summing over
    it would report "0 frame mismatches accepted" for a run in which the check's result
    is simply unknown.
    """
    found: dict[str, dict] = {}
    missing: list[str] = []
    for sample in samples:
        path = audit_sidecar_path(config, sample)
        if not path.exists():
            missing.append(sample)
            continue
        try:
            found[sample] = json.loads(path.read_text())
        except (OSError, ValueError):
            missing.append(sample)
    return found, missing


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
    audit: WindowAudit | None = None


def _load_inputs(config: WGSConfig) -> _Inputs:
    if config.reconcile.oligo_windows is None:
        raise StepIncompleteError("step_4: reconcile.oligo_windows is required; unmeasured sites are not no-edit results")
    # 🔴 NN 251 defect 12: every denominator of the window table, logged BEFORE TRACE's
    # loader below, which refuses a malformed or conflicting row outright. A refusal of
    # the real design table (50.2% filled on coordinates) then arrives with its counts.
    audit = load_windows_audited(
        str(config.reconcile.oligo_windows),
        oligo_col=config.columns.oligo_col,
        chrom_col=config.columns.chrom_col,
        start_col=config.columns.start_col,
        end_col=config.columns.end_col,
        strict=config.reconcile.strict_window_table,
    )
    if audit.n_dropped or audit.n_duplicate_oligos:
        logger.warning("step_4: window table has rows TRACE will not accept — %s", audit.summary())
    else:
        logger.info("step_4: window table — %s", audit.summary())
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
    return _Inputs(windows, assignments, step_1_bc1_call.load_bridge(config), audit)


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
    undeclared = {o for o, w in inputs.windows.items() if not getattr(w, "intended_alleles", ())}
    orientations = _check_reference_frame(config, record, inputs.audit, undeclared=undeclared)
    _check_reference_consistency(config, record)
    _write_audit_sidecar(config, sample)
    # 🔴 NN 251 defects 6 + 11 for windows that declare no intended allele: they are
    # scored from read evidence, gated at the designed base when the frame-verified
    # `wt_seq` + `design_donor` give it, else at the window minimum (planning/268).
    designed = (designed_positions_by_oligo(inputs.audit, orientations=orientations)
                if inputs.audit is not None else {})

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
            designed_positions=designed,
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


def _check_reference_frame(
    config: WGSConfig, record: WGSSample, audit: WindowAudit | None,
    *, undeclared: set[str] | None = None,
) -> dict[str, str]:
    """Verify every scored window against the design's `wt_seq` (NN 251 defect 10).

    🔴 **This is the highest-value guard on the list, because it converts an
    undetectable failure into a detectable one.** NN 243 proved that when the
    reference frame is wrong — the reference already carries the design's ALT at the
    designed base — every per-cell call FLIPS while the apparent aggregate rate lands
    on *exactly truth* (50.0% apparent on a planted 50%, per-cell accuracy **0.0%**),
    and that mixed damage partly cancels so the rate looks BETTER the worse the
    damage is. **The aggregate edit rate is therefore not a detector for frame
    error**, and no rate-level validation can substitute for this check.

    The root cause was that Stage 4 is never handed the design's `wt_seq`: it
    compares reads to whatever `fasta.fetch` returns, and its only reference check,
    `validate_tracks`, validates the bc1 locus alone. `orient_to_genome` already IS
    this comparison and already works (NN 243 confirmed it catches frame error in
    Stage 2); it was simply unreachable from here.

    When the table carries no `wt_seq`, the guard CANNOT run — and says so, loudly,
    rather than passing. A guard that silently does not run is the failure mode this
    whole exercise is about.

    Merge note (planning/268): a window that declares `intended_variants` has every REF
    validated against the FASTA by TRACE, so a wrong frame cannot invert it. A window
    that declares none is scored from read evidence against `fasta.fetch`, where it can.
    So a table with no `wt_seq` column is refused under `require_frame_check` only when
    `undeclared` (the oligos with no intended allele) is non-empty, and is recorded
    otherwise. `undeclared=None` means "assume every window is undeclared".

    Returns
    -------
    dict
        `{oligo: "plus" | "minus"}` for every window whose frame was VERIFIED.
    """
    orientations: dict[str, str] = {}
    if audit is None:
        _FRAME_AUDIT[record.sample] = {"status": "no_window_audit"}
        logger.warning("step_4: NOT VERIFIED: no window audit was supplied for %s, so the "
                       "reference frame was not checked", record.sample)
        return orientations
    if not audit.has_wt_seq:
        exposed = sorted(undeclared) if undeclared is not None else sorted(audit.windows)
        message = (
            f"step_4: the window table {audit.path.name} has no `wt_seq` column, so "
            f"the reference frame CANNOT be verified ({audit.frame_status}), and "
            f"{len(exposed)} window(s) declare no intended allele and are scored from "
            "read evidence against the reference. A wrong frame inverts every such "
            "per-cell call while leaving the aggregate rate on top of truth, so the rate "
            "cannot detect it. Add `wt_seq` (or `intended_variants`) to the window "
            "table, or set reconcile.require_frame_check=False to accept unverified "
            "calls and have the refusal recorded instead."
        )
        if config.reconcile.require_frame_check and exposed:
            raise StepIncompleteError(message)
        logger.warning("%s", message if exposed else (
            f"step_4: the window table {audit.path.name} has no `wt_seq` column; every "
            "window declares REF-validated intended alleles, so the frame is checked by "
            "TRACE instead. Recorded as frame_unverified_no_wt_seq."))
        # 🔴 This branch is the REAL COHORT's case (the design table has no `wt_seq`),
        # and the first NN 255 attempt returned here BEFORE recording anything — so
        # `frame_unverified_no_wt_seq` would have shipped 0 for a run in which not one
        # window was verified. Recorded first, then return.
        _FRAME_AUDIT[record.sample] = {
            "n_windows": audit.n_windows,
            "n_verified": 0,
            "n_frame_mismatch_accepted": 0,
            "n_no_wt_seq": audit.n_windows,
            "frame_mismatch_oligos": [],
            "require_frame_check": config.reconcile.require_frame_check,
            "status": "no_wt_seq_column",
        }
        return orientations

    checked = 0
    unverified: list[str] = []
    mismatched: dict[str, str] = {}
    with pysam.FastaFile(str(record.reference)) as fasta:
        for oligo, (chrom, start, end) in audit.windows.items():
            wt_seq = audit.wt_seqs.get(oligo)
            if not wt_seq:
                if config.reconcile.require_frame_check:
                    raise StepIncompleteError(
                        f"step_4: oligo {oligo} has a window but an empty `wt_seq`, so "
                        "its frame cannot be verified. An unverifiable site must not be "
                        "scored as if it were verified."
                    )
                logger.warning("step_4: oligo %s has no wt_seq; frame unverified", oligo)
                unverified.append(oligo)
                continue
            # 🔴 NN 255 F10. `verify_window_frame` was called UNWRAPPED here, so
            # `require_frame_check=False` suppressed the refusal only for a MISSING or
            # EMPTY `wt_seq`: a `wt_seq` that was present and disagreed with the
            # reference still raised `FrameMismatchError` out of step_4, and the escape
            # hatch the config docstring documents did not exist at any switch setting.
            try:
                orientations[oligo] = verify_window_frame(
                    fasta, chrom, start, end, wt_seq
                ).oriented
            except FrameMismatchError as exc:
                if config.reconcile.require_frame_check:
                    raise
                # Accepted, but NOT silently and NOT as if verified: the oligo is left
                # OUT of `orientations` and named in the audit sidecar.
                mismatched[oligo] = str(exc).split(". ")[0]
                logger.warning(
                    "step_4: NOT VERIFIED: oligo %s FAILED the frame check and "
                    "`reconcile.require_frame_check=False` accepted it — %s. 🔴 The "
                    "aggregate edit rate cannot detect this: NN 243 measured a wrong "
                    "frame flipping every per-cell call while the apparent rate landed "
                    "on exactly truth.",
                    oligo, mismatched[oligo],
                )
                continue
            checked += 1
    logger.info(
        "step_4: reference frame verified for %d of %d window(s) against the design's "
        "wt_seq; %d FAILED and were accepted unverified, %d had no wt_seq",
        checked, audit.n_windows, len(mismatched), len(unverified),
    )
    # Recorded where a reader of the output can find it, not only in the run log: an
    # accepted frame mismatch is the one condition under which the per-cell calls may
    # all be inverted while every rate looks right.
    _FRAME_AUDIT[record.sample] = {
        "n_windows": audit.n_windows,
        "n_verified": checked,
        "n_frame_mismatch_accepted": len(mismatched),
        "n_no_wt_seq": len(unverified),
        "frame_mismatch_oligos": sorted(mismatched),
        "require_frame_check": config.reconcile.require_frame_check,
        "status": "checked",
    }
    n_minus = sum(1 for v in orientations.values() if v == "minus")
    if n_minus:
        logger.info(
            "step_4: %d of %d verified window(s) are MINUS-oriented; their designed "
            "positions are derived after reverse-complementing, not from the raw "
            "string offsets", n_minus, len(orientations),
        )
    return orientations


def _check_reference_consistency(config: WGSConfig, record: WGSSample) -> None:
    """Compare this BAM's header against the manifest's reference (NN 251 defect 14).

    Two guards existed, were documented, and were called by nothing:
    `core.manifest.require_reference_consistency` (BAM `@PG ID:bwa` vs the manifest)
    and `core.manifest.incomparable_contigs` (contig lengths between two references).
    `_common.require_indexed_bam` checks only that the files exist — it compares no
    header field — so a BAM aligned to the wrong reference reached Stage 4 unflagged,
    which is **the frame attack itself**.

    `incomparable_contigs` is wired FIRST because it is the robust one: it compares
    `@SQ` lengths against the reference `.fai`, which every aligner writes, whereas
    `@PG ID:bwa` is absent from a BAM aligned by another tool or merged afterwards.
    An unobservable `@PG` therefore raises `ReferenceUnobservedError`, a different
    class from a mismatch, and is downgraded to a warning unless the screen asks for
    it — because "could not check" must not masquerade as "checked and wrong".
    """
    lengths = reference_contig_lengths(record.reference)
    with pysam.AlignmentFile(str(record.bam), "rb") as bam:
        observed = dict(zip(bam.references, bam.lengths))
        pg_reference = _pg_reference(bam)
    bad = incomparable_contigs(observed, lengths)
    if bad:
        raise StepIncompleteError(
            f"step_4: {record.sample}'s BAM and its declared reference disagree on "
            f"{len(bad)} contig length(s) (first: {', '.join(map(str, bad[:5]))}). "
            "Coordinates are not portable across references (chrVIII is "
            "CUP1-contracted; no single offset exists), so every window coordinate "
            "would be read in the wrong frame."
        )
    try:
        require_reference_consistency(
            [record], {record.sample: pg_reference} if pg_reference else {}
        )
        _REFERENCE_NAME_CHECK[record.sample] = "checked"
    except ReferenceUnobservedError as exc:
        # 🔴 NN 255 F3b: this used to be a passing log line, so a check that did not run
        # read exactly like a check that ran. It is now prefixed `NOT VERIFIED:` and
        # COUNTED into the step's provenance, so a run in which the name check never
        # executed cannot be read as a run in which it passed.
        _REFERENCE_NAME_CHECK[record.sample] = "skipped_pg_unobserved"
        logger.warning(
            "step_4: NOT VERIFIED: %s — contig lengths agree, but no reference could "
            "be read from the BAM's `@PG ID:bwa` command line, so the manifest/BAM "
            "NAME check DID NOT RUN (it is counted as skipped, not as passed). The "
            "residual guard is the `@SQ` length comparison above, which is the stronger "
            "one. Cause: %s",
            record.sample, exc,
        )


#: FASTA suffixes a `bwa` command line can name a reference with, bgzipped forms
#: included. NN 255 F3b: `comp.fa.gz` is a reference `bwa index` accepts and the old
#: tuple did not list, so `_pg_reference` returned `None` for it (2 of 11 real command
#: lines NN 252 tested) and the manifest/BAM name check SILENTLY did not run.
_FASTA_SUFFIXES = (
    ".fa", ".fasta", ".fsa", ".fna",
    ".fa.gz", ".fasta.gz", ".fsa.gz", ".fna.gz",
    ".fa.bgz", ".fasta.bgz",
)


def _pg_reference(bam: pysam.AlignmentFile) -> str | None:
    """The reference path in the BAM's `@PG ID:bwa` command line, or `None`.

    `None` is returned rather than guessed: an absent `@PG` is "unobserved", and
    `require_reference_consistency` raises `ReferenceUnobservedError` for it, which
    is deliberately a different class from a mismatch.

    NN 252 measured this against 11 real `bwa` command lines: **0 returned a WRONG
    token** (including the two traps — a read group whose `ID:` is `comp.fa`, and a
    `.fa`-named FASTQ — because the reference precedes the reads), 9 were correct, and
    **2 returned `None`**: a bgzipped reference and an extensionless `bwa index -p`
    prefix. `None` downgrades to a warning, so those two were a check that quietly did
    not run. The suffix list now covers the gz forms; the index-prefix form is handled
    by `_pg_positional_reference` below, and whatever remains unobserved is reported
    LOUDLY and COUNTED rather than mentioned in passing.
    """
    for pg in bam.header.to_dict().get("PG", []):
        if "bwa" not in str(pg.get("ID", "")).lower():
            continue
        tokens = str(pg.get("CL", "")).split()
        for token in tokens:
            if token.endswith(_FASTA_SUFFIXES):
                return token
        positional = _pg_positional_reference(tokens)
        if positional is not None:
            return positional
    return None


#: `bwa mem` flags that consume the following token, so a positional scan does not read
#: a flag's VALUE as the reference. Deliberately a conservative list: an unknown flag
#: with a value would make the scan return a wrong token, which is worse than `None`,
#: so the scan refuses unless every token it skipped is accounted for.
_BWA_VALUE_FLAGS = frozenset(
    "-t -k -w -d -r -c -A -B -O -E -L -U -R -T -H -I -m -y -h -K -v -G -W -o -f -x"
    .split()
)


def _pg_positional_reference(tokens: Sequence[str]) -> str | None:
    """The `bwa mem` reference argument when it carries no FASTA suffix.

    `bwa index -p <prefix>` makes an extensionless index, and `bwa mem <prefix> r1 r2`
    then names it with no suffix at all. The reference is the FIRST positional argument
    after the subcommand, so it is recoverable — but only if every value-taking flag
    before it is known, because otherwise a flag's value would be returned as the
    reference. A WRONG token is worse than `None` (it would compare the manifest against
    something that is not a reference), so an unrecognized dashed token aborts the scan.
    """
    rest = list(tokens)
    if rest and rest[0].endswith("bwa") or (rest and "bwa" in rest[0]):
        rest = rest[1:]
    if rest and rest[0] in ("mem", "aln", "samse", "sampe", "bwasw"):
        rest = rest[1:]
    else:
        return None
    i = 0
    while i < len(rest):
        token = rest[i]
        if not token.startswith("-"):
            return token
        if token in _BWA_VALUE_FLAGS:
            i += 2
            continue
        if len(token) == 2:  # a known-shape boolean flag, e.g. -p, -a, -C, -M
            i += 1
            continue
        return None  # an unrecognized dashed token: refuse rather than guess
    return None


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

    audits, missing_audits = _read_audit_sidecars(config, [r.sample for r in records])
    if missing_audits:
        logger.warning(
            "step_4: NOT VERIFIED: %d of %d sample(s) have no audit sidecar (%s...), so "
            "the frame-check and reference-name counts below are LOWER BOUNDS, not "
            "totals. A shard written before the sidecar existed is the usual cause.",
            len(missing_audits), len(records), ", ".join(missing_audits[:3]),
        )

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
            "require_frame_check": config.reconcile.require_frame_check,
            "frame_check": {
                s: audits[s].get("frame_check", {}) for s in sorted(audits)
            },
        },
        counts={
            "cells": n,
            **labels,
            # 🔴 NN 255 F10 / F3b, read from the per-sample SIDECARS on disk rather than
            # from a module global, because the fan-out runs each sample in its own
            # process and a cache-valid sample returns before any check. A non-zero
            # `audit_samples_missing` means these sums are lower bounds, not totals.
            "audit_samples_recorded": len(audits),
            "audit_samples_missing": len(missing_audits),
            "frame_mismatch_accepted": sum(
                a.get("frame_check", {}).get("n_frame_mismatch_accepted", 0)
                for a in audits.values()
            ),
            "frame_unverified_no_wt_seq": sum(
                a.get("frame_check", {}).get("n_no_wt_seq", 0) for a in audits.values()
            ),
            "reference_name_check_skipped": sum(
                1 for a in audits.values()
                if a.get("reference_name_check") == "skipped_pg_unobserved"
            ),
        },
    )
    return out
