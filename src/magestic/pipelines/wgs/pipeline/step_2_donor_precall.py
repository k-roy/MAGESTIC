"""Stage 2: explicit donor-window outcomes from sample VCF and gVCF inputs.

Each selected clone needs one donor specification and both variant and reference-depth
inputs. Missing gVCF remains a data gate, not an unedited result. Underscore contigs
are supported by structured TRACE variants. The synthetic-error axis is stamped
DISABLED when design sequence substitutes for an unavailable empirical donor.
"""

from __future__ import annotations

import hashlib
import logging
import os
from pathlib import Path

from trace_crispr.wgs import donor_outcome

from ..config import DONOR_PRECALL_SUBDIR, WGSConfig
from ..core.manifest import WGSSample
from ..core.samplesheet import (
    build_samplesheet,
    read_axis_stamp,
    stamp_outcome_table,
)
from ..core.vcf_adapter import (
    PRECALL_GVCF_SUBDIR,
    PRECALL_GVCF_SUFFIX,
    PRECALL_VCF_SUBDIR,
    PRECALL_VCF_SUFFIX,
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
        result = split_plate_vcf(plate, workdir, samples=selected, ad_slot=config.donor_precall.ad_format_slot_out)
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
