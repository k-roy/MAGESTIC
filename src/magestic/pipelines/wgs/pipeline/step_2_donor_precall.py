"""Stage 2: explicit donor-window outcomes from sample VCF and gVCF inputs.

Each selected clone needs one donor specification and both variant and reference-depth
inputs. Missing gVCF remains a data gate, not an unedited result. Underscore contigs
are supported by structured TRACE variants. The synthetic-error axis is stamped
DISABLED when design sequence substitutes for an unavailable empirical donor.

Merged in planning/268 with the 222 lane (NN 251/255), whose orthogonal guards are kept:
the outcome-reason columns appended to every shard (defects 4 and 8), the caller
provenance refusal read from each plate VCF's own header (NN 255), and the allele-fraction
annotation from `split_plate_vcf` (defect 1 / F11b). Its workarounds for TRACE `7c4c10f`
are retired, because the pinned TRACE closes those defects itself: the `ContigAlias`
encode/decode (GATE 3; TRACE builds variant ids from structured fields), the gVCF
reference-block expansion (defect 5; TRACE reads `INFO/END` and validates every REF,
which the expander's copied REF failed), and the multi-allelic refusal (defect 14; TRACE
selects the ALT that GT names).
"""

from __future__ import annotations

import hashlib
import logging
import os
from collections.abc import Sequence
from pathlib import Path

from trace_crispr.wgs import donor_outcome

from ..config import DONOR_PRECALL_SUBDIR, WGSConfig
from ..core.donor_precall import audit_vcf_window
from ..core.manifest import WGSSample
from ..core.outcome_reason import explain_outcome
from ..core.samplesheet import (
    build_samplesheet,
    read_axis_stamp,
    stamp_outcome_table,
)
from ..core.vcf_adapter import (
    CALLER_PLOIDY,
    MAG_MINOR_TAG,
    PLOIDY_CAVEAT,
    PRECALL_GVCF_SUBDIR,
    PRECALL_GVCF_SUFFIX,
    PRECALL_VCF_SUBDIR,
    PRECALL_VCF_SUFFIX,
    caller_provenance,
    precall_inputs,
    split_plate_vcf,
    vcf_samples,
)
from ._common import (
    GateError,
    StepIncompleteError,
    cache_key,
    cache_valid,
    cache_commit,
    cohort_table,
    collect_shards,
    manifest_samples,
    preflight_provenance,
    publish,
    record_for,
    require_indexed_bam,
    require_row_count,
    shard_dir,
    write_provenance,
    write_tsv,
)

logger = logging.getLogger(__name__)

STEP = "step_2_donor_precall"
SHARD_SUFFIX = ".donor_outcome.tsv"
TABLE_NAME = "donor_outcome.tsv"

GATE_1 = (
    "The synthetic-error axis requires an explicit empirical donor; substituting "
    "the design disables the axis and is stamped in each output."
)
GATE_2 = (
    "A selected sample is missing its required gVCF reference-depth input; "
    "missing evidence cannot establish an unedited donor."
)


def shard_path(config: WGSConfig, sample: str) -> Path:
    """Per-sample donor-outcome shard."""
    return shard_dir(config, DONOR_PRECALL_SUBDIR) / f"{sample}{SHARD_SUFFIX}"


def table_path(config: WGSConfig) -> Path:
    """Cohort donor-outcome table."""
    return cohort_table(config, TABLE_NAME)


def samplesheet_path(config: WGSConfig, sample: str) -> Path:
    """TRACE's 7-column donor samplesheet for one sample."""
    return config.paths.samplesheets / f"{sample}.samplesheet.tsv"


def split_vcf_path(config: WGSConfig, sample: str) -> Path:
    """The per-sample called-variant VCF that TRACE reads.

    Stage 2 PRODUCES it only when `donor_precall.plate_vcfs` is configured (the plate
    adapter splits it). With `plate_vcfs: []` it is an externally supplied INPUT:
    hashed where consumed, never given a Stage 2 producer receipt.
    """
    return config.paths.donor_precall_workdir / PRECALL_VCF_SUBDIR / f"{sample}{PRECALL_VCF_SUFFIX}"


def gvcf_path(config: WGSConfig, sample: str) -> Path:
    """The per-sample reference-depth gVCF (always an input)."""
    return config.paths.donor_precall_workdir / PRECALL_GVCF_SUBDIR / f"{sample}{PRECALL_GVCF_SUFFIX}"


def _receipt(path: Path) -> Path:
    from ..qc.provenance import SIDECAR_SUFFIX
    return Path(str(path) + SIDECAR_SUFFIX)


def _owned_outputs(config: WGSConfig, sample: str) -> list[Path]:
    """Every data artifact Stage 2 produces for one sample, shard first."""
    owned = [shard_path(config, sample), samplesheet_path(config, sample)]
    if config.donor_precall.plate_vcfs:
        owned.append(split_vcf_path(config, sample))
    return owned


def _completion_set(config: WGSConfig, sample: str) -> list[Path]:
    """What a completed sample consists of: each owned artifact AND its receipt.

    Shared by the per-sample path and the fan-in, so a missing, foreign or damaged
    receipt makes exactly that sample incomplete instead of being skipped.
    """
    owned = _owned_outputs(config, sample)
    return owned + [_receipt(path) for path in owned]


_MD5_MEMO: dict[tuple, str] = {}


def _md5_attested(path: Path) -> str:
    """Raw MD5 of an input, computed once per process per file identity.

    A reference FASTA or plate VCF feeds many sample receipts; rehashing it for each
    would cost samples x bytes. The key is the file's stat identity (device, inode,
    size, nanosecond mtime and ctime), so a replaced or edited file is hashed again.
    """
    path = Path(path)
    stat = path.stat()
    key = (os.path.realpath(path), stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    if key not in _MD5_MEMO:
        digest = hashlib.md5()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 22), b""):
                digest.update(block)
        _MD5_MEMO[key] = digest.hexdigest()
    return _MD5_MEMO[key]


def _write_receipt(output: Path, *, inputs: dict[str, Path], parameters: dict, counts: dict) -> Path:
    """One artifact's own producer receipt, without touching the run-level index.

    Sample workers run concurrently under Snakemake/SLURM; refreshing the shared
    directory index from each of them would let one overwrite another's newer
    membership. The fan-in refreshes it once.
    """
    from ..qc.provenance import write_step_provenance
    return write_step_provenance(STEP, output, inputs=inputs, parameters=parameters, counts=counts,
                                 input_md5={role: _md5_attested(path) for role, path in inputs.items()})


def _donor_inputs(config: WGSConfig) -> dict[str, Path]:
    """The source tables every samplesheet row is derived from."""
    table = getattr(config.donor_precall, "donor_table", None) or config.reconcile.oligo_windows
    return {"manifest": Path(config.paths.manifest), "donor_table": Path(table),
            "oligo_assignments": Path(config.reconcile.oligo_assignments)}


def gate_reason() -> str:
    """The string `run_pipeline` records in place of an artifact path."""
    return f"GATED: {GATE_2} | {GATE_1}"


def run_step_2(config: WGSConfig, sample: str | None = None) -> Path:
    """Pre-call donor-window editing outcomes. Returns the path written.

    Raises `GateError` when a selected sample lacks its required gVCF input.
    """
    preflight_provenance(STEP)
    if sample is not None:
        return _run_one_sample(config, sample)
    return _combine(config)


def _run_one_sample(config: WGSConfig, sample: str, *, record: WGSSample | None = None,
                    prepared_rows: dict | None = None,
                    loaded_reference: donor_outcome.LoadedReference | None = None) -> Path:
    """Run one sample and publish its completed outcome atomically."""
    out = shard_path(config, sample)
    record = record if record is not None else record_for(config, sample)
    workdir = config.paths.donor_precall_workdir
    dependencies = {'vcf': split_vcf_path(config, sample), 'gvcf': gvcf_path(config, sample)}
    completion = _completion_set(config, sample)
    key = cache_key(config, STEP, record, inputs=dependencies)
    if cache_valid(completion, key):
        logger.info("step_2: %s already complete — skipping", sample)
        return out

    require_indexed_bam(record)

    if prepared_rows is None:
        _require_vcf_inputs(config, record, workdir)
    elif sample not in prepared_rows:
        raise StepIncompleteError(f"step_2: prepared inputs omit requested sample {sample!r}")
    sheet, axis = _build_samplesheet(config, sample, row=prepared_rows[sample] if prepared_rows is not None else None)
    _write_receipt(sheet, inputs=_donor_inputs(config),
                   parameters={"sample": sample, "synthetic_donor_source": config.donor_precall.synthetic_donor_source,
                               "synthetic_error_axis": axis},
                   counts={"samplesheet_rows": 1})
    # Adapter publication changes the generated per-sample VCF identity; capture
    # the consumed version only after that step, retaining source plate VCFs too.
    key = cache_key(config, STEP, record, inputs=dependencies)

    raw = out.with_name(out.name + ".raw")
    results = donor_outcome.precall(
        sheet,
        workdir,
        Path(record.reference),
        raw,
        ne_threshold=config.donor_precall.ne_threshold,
        quantification_window_size=config.donor_precall.quantification_window_size,
        quantification_window_cov_perc=config.donor_precall.quantification_window_cov_perc,
        samples={sample},
        loaded_reference=loaded_reference,
    )
    if len(results) != 1:
        raise StepIncompleteError(
            f"step_2: precall returned {len(results)} outcomes for {sample}, expected 1. "
            "The returned sample denominator must match the requested sample."
        )
    _annotate_outcome_reason(
        config, record, workdir, raw, results[0],
        minor_alleles=_minor_alleles_from_vcf(dependencies["vcf"]),
    )
    stamp_outcome_table(raw, synthetic_error_axis=axis)
    publish(raw, out)
    _write_receipt(out, inputs={"samplesheet": sheet, "called_variants_vcf": dependencies["vcf"],
                                "reference_depth_gvcf": dependencies["gvcf"], "reference": Path(record.reference)},
                   parameters={"sample": sample,
                               "called_variants_vcf_origin": ("stage2_plate_adapter" if config.donor_precall.plate_vcfs
                                                              else "external_input"),
                               "bam_preflight": str(record.bam),
                               "ne_threshold": config.donor_precall.ne_threshold,
                               "quantification_window_size": config.donor_precall.quantification_window_size,
                               "quantification_window_cov_perc": config.donor_precall.quantification_window_cov_perc,
                               "synthetic_error_axis": axis},
                   counts={"outcome_rows": len(results)})
    # Completion is recorded last: after every owned artifact AND its receipt exists.
    cache_commit(completion, key)
    logger.info("step_2: %s outcome=%s", sample, results[0].editing_outcome)
    return out


#: Columns NN 251 appends to TRACE's 5-column Stage-2 table (defects 4 and 8).
#: Appended, never inserted, and `step_5` reads `editing_outcome` by NAME
#: (`step_5_join_outcomes.py:139`), so no consumer shifts.
OUTCOME_REASON_COLUMNS = (
    "callable",
    "outcome_reason",
    "qwin_covered",
    "qwin_size",
    "n_variants_discarded_outside_donor_window",
    "minor_allele_variants",
)


def _annotate_outcome_reason(
    config: WGSConfig, record: WGSSample, workdir: Path, table: Path, outcome: object,
    *, minor_alleles: Sequence[str] = (),
) -> None:
    """Append WHY this outcome is what it is to the shipped shard (NN 251 defect 4).

    🔴 **In the TABLE, not only in the log.** `unknown` is the only value in the
    Stage-2 table that means "I could not tell", and TRACE records nothing about
    which cause produced it — no gVCF, unexpanded reference blocks, depth below
    `ne_threshold`, or a genuinely unsequenced locus all collapse to the same five
    letters. A reason that reached only the run log would leave every downstream
    reader unable to tell "uncallable" from "unedited", which is the single most
    load-bearing distinction in this analysis.

    Also carries defect 8's count: how many VCF records fell outside the donor window
    and were therefore discarded in silence, so the VCF reconciles with the outcome.
    """
    sheet = samplesheet_path(config, record.sample)
    try:
        lines = sheet.read_text().splitlines()
        cols = lines[0].split("\t")
        row = lines[1].split("\t")
        chrom = row[cols.index("chrom")].strip()
        dstart = int(row[cols.index("donor_start_coord")])
        dend = int(row[cols.index("donor_end_coord")])
    except (OSError, IndexError, ValueError):
        # Even the failure is explicit in the table: a shard whose reason columns are
        # absent would be indistinguishable from one where the reason was computed.
        _append_columns(table, OUTCOME_REASON_COLUMNS,
                        ["", "samplesheet_unreadable", "", "", "",
                         ";".join(minor_alleles)])
        logger.warning("step_2: %s could not read its samplesheet to explain the "
                       "outcome; the flag stands but carries no reason", record.sample)
        return
    ex = explain_outcome(
        getattr(outcome, "editing_outcome", "unknown"),
        gvcf_path=workdir / PRECALL_GVCF_SUBDIR / f"{record.sample}{PRECALL_GVCF_SUFFIX}",
        chrom=chrom, donor_start=dstart, donor_end=dend, sample=record.sample,
        quantification_window_size=config.donor_precall.quantification_window_size,
        quantification_window_cov_perc=(
            config.donor_precall.quantification_window_cov_perc
        ),
        ne_threshold=config.donor_precall.ne_threshold,
    )
    audit = audit_vcf_window(
        workdir / PRECALL_VCF_SUBDIR / f"{record.sample}{PRECALL_VCF_SUFFIX}",
        chrom=chrom, donor_start=dstart, donor_end=dend,
    )
    _append_columns(
        table,
        OUTCOME_REASON_COLUMNS,
        [
            str(ex.callable_).lower(),
            ex.reason,
            str(ex.n_window_positions_covered),
            str(ex.n_window_positions),
            str(audit.n_discarded),
            ";".join(minor_alleles),
        ],
    )
    logger.info(
        "step_2: %s callable=%s reason=%s qwin=%d-%d covered=%d/%d",
        record.sample, ex.callable_, ex.reason, ex.qwin_start, ex.qwin_end,
        ex.n_window_positions_covered, ex.n_window_positions,
    )
    if audit.n_discarded:
        # NN 251 defect 8: distal unintended edits are invisible to Stage 2. Counted
        # here rather than discarded, so the VCF's records reconcile with the outcome.
        logger.info("step_2: %s %s", record.sample, audit.summary())


def _append_columns(table: Path, columns: Sequence[str], values: Sequence[str]) -> None:
    """Append fixed columns to TRACE's Stage-2 table, header and data rows alike.

    TRACE writes `#sample_name\\t...` as its header, so the header is the `#`-prefixed
    line and every other non-blank line is data. Appending (never inserting) keeps
    `step_5`'s `precall[sample]["editing_outcome"]` name-based lookup correct, and
    re-running is safe because a shard that already carries the columns is skipped.
    """
    lines = table.read_text().splitlines()
    if any(line.startswith("#") and columns[0] in line.split("\t") for line in lines):
        return  # idempotent: already annotated
    out: list[str] = []
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("#"):
            out.append(line + "\t" + "\t".join(columns))
        else:
            out.append(line + "\t" + "\t".join(values))
    with table.open("w", newline="\n") as fh:
        for line in out:
            fh.write(line + "\n")


def _minor_alleles_from_vcf(path: Path) -> tuple[str, ...]:
    """Variant ids `split_plate_vcf` tagged `MAG_MINOR` in this sample's own VCF.

    Read back from the file rather than carried in memory, so the per-sample SLURM path,
    the cohort path and a cache-reused split all report the same set. An externally
    supplied VCF (`plate_vcfs: []`) was never annotated and so yields `()`, exactly as
    the 222 lane's in-memory route did.
    """
    if not Path(path).exists():
        return ()
    found: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) > 7 and MAG_MINOR_TAG in cols[7].split(";"):
                found.append(f"{cols[0]}_{cols[1]}_{cols[3]}_to_{cols[4]}")
    return tuple(dict.fromkeys(found))


def _require_vcf_inputs(config: WGSConfig, record: WGSSample, workdir: Path) -> None:
    """Use the same preparation path for a single scheduler sample."""
    _prepare_vcf_inputs(config, [record], workdir)


def _prepare_vcf_inputs(config: WGSConfig, records: list[WGSSample], workdir: Path) -> None:
    """Validate reference-depth inputs, then scan each selected plate at most once."""
    samples = [record.sample for record in records]
    inputs = precall_inputs(workdir, samples)
    if inputs.gvcf_missing:
        raise GateError(f"step_2: {GATE_2} (no gVCF for {', '.join(inputs.gvcf_missing)})")
    sources = {sample: [] for sample in samples}
    groups = []
    for plate in config.donor_precall.plate_vcfs:
        names = set(vcf_samples(plate))
        selected = [sample for sample in samples if sample in names]
        if selected:
            groups.append((plate, selected))
        for sample in selected:
            sources[sample].append(plate)
    if config.donor_precall.plate_vcfs:
        for sample, matching in sources.items():
            if not matching:
                raise StepIncompleteError(f"step_2: sample {sample!r} absent from every configured plate VCF")
            if len(matching) != 1:
                raise StepIncompleteError(f"step_2: sample {sample!r} occurs in multiple plate VCFs: {matching}")
    for plate, selected in groups:
        _require_recoverable_caller(config, plate)
        result = split_plate_vcf(
            plate, workdir, samples=selected, ad_slot=config.donor_precall.ad_format_slot_out,
            min_allele_fraction=config.donor_precall.min_allele_fraction,
            min_alt_reads=config.donor_precall.min_alt_reads,
        )
        if (
            result.n_rescued_by_af or result.n_no_ad or result.n_dropped_below_af
            or result.n_flagged_minor or result.n_dropped_below_min_alt_reads
            or result.n_no_ad_carrying_alt
        ):
            logger.info(
                "step_2: %s allele fraction — %d of %d record(s) had a usable AD; "
                "flagged MAG_MINOR %d (of which %d ALREADY CARRIED THE ALT IN GT, the "
                "NN 252 F11b over-call class); rescued %d, below fraction %d, below "
                "%s alt read(s) %d, no AD %d, no AD WITH an alt GT %d",
                Path(plate).name, result.n_af_computed, result.n_records_in,
                result.n_flagged_minor, result.n_flagged_minor_carrying_alt,
                result.n_rescued_by_af, result.n_dropped_below_af,
                config.donor_precall.min_alt_reads if config.donor_precall.min_alt_reads is not None
                else "the depth-conditional",
                result.n_dropped_below_min_alt_reads, result.n_no_ad,
                result.n_no_ad_carrying_alt,
            )
        if result.n_multiallelic:
            # Not refused since planning/268: the pinned TRACE selects the ALT that GT
            # names. Counted, and recorded in each sample's split receipt below.
            logger.info("step_2: %s carries %d multi-allelic record(s); TRACE selects "
                        "the allele GT names", Path(plate).name, result.n_multiallelic)
        for sample in selected:
            produced = Path(result.vcf_paths[sample])
            if produced.resolve() != split_vcf_path(config, sample).resolve():
                raise StepIncompleteError(f"step_2: adapter wrote {produced}, expected {split_vcf_path(config, sample)}")
            _write_receipt(produced, inputs={"plate_vcf": Path(plate)},
                           parameters={"sample": sample, "ad_format_slot_out": result.ad_slot,
                                       "adapter": "magestic.pipelines.wgs.core.vcf_adapter.split_plate_vcf"},
                           counts={"plate_records_in": result.n_records_in,
                                   "sample_records_out": result.n_records_out[sample],
                                   "plate_multiallelic_records": result.n_multiallelic})
    inputs = precall_inputs(workdir, samples)
    if inputs.vcf_missing:
        raise StepIncompleteError(
            f"step_2: no per-sample VCF for {', '.join(inputs.vcf_missing)} under {workdir}. "
            "A called-variant input is required independently of reference coverage."
        )
    for sample in samples:
        for path in (workdir / "vcf" / f"{sample}.fil.vcf", workdir / "gvcf" / f"{sample}.g.vcf"):
            names = vcf_samples(path)
            if names != (sample,):
                raise StepIncompleteError(f"step_2: {path} sample header {names!r} does not match {sample!r}")


def _require_recoverable_caller(config: WGSConfig, plate_vcf: Path) -> None:
    """Read the plate VCF's OWN caller command line and refuse an unrecoverable combination.

    🔴 **Detected, not declared — NN 255.** Which ploidy these VCFs were called at is
    the premise the whole of defect 1 rests on; NN 251 asserted it in a docstring and
    NN 252 measured the opposite. `bcftools` writes `##bcftools_callCommand=` into every
    VCF it calls and `bcftools merge` preserves it, so the question is answerable from
    the file and is answered here, every run.

    Two measured facts drive the two refusals (NN 255 floor sweep, bcftools 1.21, flat
    40x, alt counts 0..40):

    * At `--ploidy 1` every minority allele is `GT 0`, so `carries_alt` drops it and the
      ALLELE-FRACTION RESCUE is the only thing that can see it. `--ploidy 1` with
      `min_allele_fraction=None` therefore DISCARDS every partial edit while looking
      like a clean haploid run — refused.
    * With `-v`, a site at which no sample was called non-reference is not emitted at
      all, so a mixture private to one well is unrecoverable at any threshold. Measured
      private-mixture emission floor at 40x: `-mv` 6/40 = 15%, `-mv --ploidy 1` 34/40 =
      85%, `-m --ploidy 1 -A` 1/40 = 2.5% (magnitudes are sequence-context dependent;
      the ORDERING is what is stable). **86.6% of wells are alone at their designed
      position within their Tn5 plate** — measured on the 20231213 plate-2348/2349
      cohort (n=768), the only cohort whose per-well table carries an oligo assignment;
      `spg_20240702_POST`'s keyfile names only a yL sublibrary, which pools many oligos,
      so its private fraction can only be HIGHER. Either way the PRIVATE floor is the
      one that applies — warned loudly, with the number, because it bounds the
      correct-edit rate and no code change can lift it.
    """
    prov = caller_provenance(plate_vcf)
    if prov.n_call_commands > 1:
        raise StepIncompleteError(
            f"step_2: {Path(plate_vcf).name} carries {prov.n_call_commands} different "
            "`##bcftools_callCommand=` lines, so its records were not all called the "
            "same way and no single detection floor describes the file. Re-call the "
            "plate in one pass."
        )
    if prov.ploidy == CALLER_PLOIDY and config.donor_precall.min_allele_fraction is None:
        raise StepIncompleteError(
            f"step_2: {Path(plate_vcf).name} was called at `--ploidy {prov.ploidy}` "
            "(read from its own header) but `donor_precall.min_allele_fraction` is "
            "None. At haploid ploidy EVERY minority allele is `GT 0`, so the "
            "allele-fraction rescue is the only thing that can see one, and with it "
            "disabled every partial edit is discarded while the run looks clean. Set "
            "`min_allele_fraction` (0.05 with the default "
            f"min_alt_reads={config.donor_precall.min_alt_reads or 'depth-conditional'} floor), or re-call "
            "the plate at the diploid default."
        )
    if prov.ploidy is None:
        logger.warning(
            "step_2: NOT VERIFIED: %s records no `##bcftools_callCommand=` (or an "
            "unparseable `--ploidy`), so the caller's ploidy is UNOBSERVED, not "
            "diploid. %s",
            Path(plate_vcf).name, PLOIDY_CAVEAT,
        )
    elif prov.variants_only:
        logger.warning(
            "step_2: %s was called with `-v` (`%s`), so a site at which NO sample was "
            "called non-reference carries no record and a minority allele PRIVATE to "
            "one well is unrecoverable at any threshold. Measured private-mixture "
            "emission floor at 40x: 6/40 (15%%) at the diploid default, 34/40 (85%%) at "
            "`--ploidy 1`; `-m --ploidy 1 -A` (no `-v`) reaches 1/40 (2.5%%). 86.6%% of "
            "wells are alone at their designed position within their Tn5 plate "
            "(measured on the 20231213 plate-2348/2349 cohort, n=768), so this is the "
            "floor that applies. Regenerate with `MPILEUP_CALL_WINDOWS_COMMAND` "
            "(all-sites RESTRICTED to the donor windows) to lift it; the unrestricted "
            "all-sites form is ~1.2e7 records per plate and is re-read once per sample.",
            Path(plate_vcf).name, prov.command,
        )
    else:
        logger.info(
            "step_2: %s caller provenance ploidy=%s variants_only=%s keep_alts=%s — "
            "private minor alleles ARE recoverable from this file",
            Path(plate_vcf).name, prov.ploidy, prov.variants_only, prov.keep_alts,
        )


def _build_samplesheet(config: WGSConfig, sample: str, *, row: dict | None = None) -> tuple[Path, str]:
    """Build this sample's explicit donor specification and axis metadata."""
    from ..readers.oligo_tables import donor_samplesheet_rows

    rows = list(donor_samplesheet_rows(config, [sample])) if row is None else [row]
    if any(r.get("sample") != sample for r in rows):
        raise StepIncompleteError(f"step_2: donor row identity does not match {sample!r}")
    result = build_samplesheet(
        rows,
        samplesheet_path(config, sample),
        synthetic_donor_source=config.donor_precall.synthetic_donor_source,
    )
    if result.n_rows != 1:
        raise StepIncompleteError(
            f"step_2: samplesheet for {sample} has {result.n_rows} rows, expected 1"
        )
    return result.path, result.synthetic_error_axis


def prepare_step_2(config: WGSConfig, records: list[WGSSample]) -> dict[str, dict]:
    """Prepare selected inputs once for direct-cohort execution.

    Returns in-memory donor rows keyed by sample after validating/splitting source
    VCFs. This is a preparation seam, not a cross-process cache or completion flag;
    scheduler preparation requires separately persisted, fingerprinted artifacts.
    """
    from ..readers.oligo_tables import donor_samplesheet_rows
    if not records:
        return {}
    samples = [record.sample for record in records]
    if len(samples) != len(set(samples)):
        raise StepIncompleteError("step_2: duplicate requested preparation sample")
    _prepare_vcf_inputs(config, records, config.paths.donor_precall_workdir)
    rows = list(donor_samplesheet_rows(config, samples))
    prepared = {row["sample"]: row for row in rows}
    if set(prepared) != set(samples) or len(rows) != len(samples):
        raise StepIncompleteError("step_2: donor preparation sample denominator mismatch")
    return prepared


def _combine(config: WGSConfig) -> Path:
    """Fan-in: every sample, then the cohort table, then re-read the GATE-1 stamp."""
    records = manifest_samples(config)
    groups: dict[Path, list[WGSSample]] = {}
    workdir = config.paths.donor_precall_workdir
    pending = []
    for record in records:
        sample = record.sample
        key = cache_key(config, STEP, record, inputs={"vcf": split_vcf_path(config, sample),
                                                    "gvcf": gvcf_path(config, sample)})
        if not cache_valid(_completion_set(config, sample), key):
            pending.append(record)
            groups.setdefault(Path(record.reference).resolve(), []).append(record)
    prepared = prepare_step_2(config, pending)
    for reference, group in groups.items():
        # Release each loaded reference before the next group: memory scales with
        # one reference, not the number of distinct backgrounds in a cohort.
        loaded = donor_outcome.load_reference(reference)
        for record in group:
            _run_one_sample(config, record.sample, record=record, prepared_rows=prepared,
                            loaded_reference=loaded)
        del loaded

    header, rows = collect_shards(
        records, lambda r: shard_path(config, r.sample), label="step_2 donor outcome"
    )
    out = write_tsv(table_path(config), header, rows)
    n = require_row_count(out, len(rows), len(records), label="step_2 donor outcome")
    _require_axis_stamp(out)

    unknown = sum(1 for row in rows if row.split("\t")[1].startswith("unknown"))
    logger.info(
        "step_2: %d samples; %d (%.1f%% of %d) outcome=unknown. unknown is NOT unedited "
        "(reference coverage or genotype evidence is insufficient).",
        n,
        unknown,
        100.0 * unknown / n if n else 0.0,
        n,
    )

    # The cohort receipt binds its full source closure and the shards it concatenated.
    inputs = dict(_donor_inputs(config))
    for index, plate in enumerate(config.donor_precall.plate_vcfs):
        inputs[f"plate_vcf:{index}"] = Path(plate)
    for index, reference in enumerate(dict.fromkeys(Path(r.reference) for r in records)):
        inputs[f"reference:{index}"] = reference
    for record in records:
        sample = record.sample
        inputs[f"called_variants_vcf:{sample}"] = split_vcf_path(config, sample)
        inputs[f"reference_depth_gvcf:{sample}"] = gvcf_path(config, sample)
        inputs[f"samplesheet:{sample}"] = samplesheet_path(config, sample)
        inputs[f"shard:{sample}"] = shard_path(config, sample)
    write_provenance(
        STEP,
        out,
        inputs=inputs,
        input_md5={role: _md5_attested(path) for role, path in inputs.items()},
        parameters={
            "synthetic_error_axis": read_axis_stamp(out),
            "synthetic_donor_source": config.donor_precall.synthetic_donor_source,
            "ne_threshold": config.donor_precall.ne_threshold,
            "quantification_window_size": config.donor_precall.quantification_window_size,
            "quantification_window_cov_perc": config.donor_precall.quantification_window_cov_perc,
            "reference_depth_inputs": "validated_per_sample_gvcf",
            "unedited_coverage_policy": (
                "central_window_fraction_and_all_diagnostic_positions; "
                "length_changing_design_requires_whole_donor_reference_coverage"
            ),
        },
        counts={"samples": n, "outcome_unknown": unknown},
    )
    return out


def _require_axis_stamp(path: Path) -> None:
    """Re-read the written table and confirm the GATE-1 stamp survived the fan-in."""
    axis = read_axis_stamp(path)
    if axis is None:
        raise StepIncompleteError(
            f"step_2: {path} carries no synthetic_error_axis stamp. An unstamped Stage-2 "
            "output is indistinguishable from one where the axis genuinely ran, so a zero "
            "in `synthetic_errors` would read as 'no synthesis errors found'."
        )
