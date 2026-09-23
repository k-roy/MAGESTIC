"""Stream a plate VCF into atomic per-sample variant and unresolved-genotype VCFs.

GT is selected by FORMAT name. Multiallelic records retain their allele ordering and
GT/AD fields so TRACE can select the actual genotype allele; the adapter never derives
an alternate genotype. Missing/incomplete genotypes remain explicit evidence rather
than disappearing as presumed reference calls. AD reordering is a compatibility representation only, because
the repaired TRACE parser reads named fields. No reference-depth gVCF is fabricated.
"""

from __future__ import annotations

import gzip
import os
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "PRECALL_GVCF_SUBDIR",
    "PRECALL_GVCF_SUFFIX",
    "PRECALL_VCF_SUBDIR",
    "PRECALL_VCF_SUFFIX",
    "PrecallInputs",
    "VcfAdapterError",
    "VcfSplitResult",
    "carries_alt",
    "precall_inputs",
    "reorder_format_field",
    "split_plate_vcf",
    "vcf_samples",
]

#: The layout `donor_outcome.precall` hardcodes. Owned here, not left implicit.
PRECALL_VCF_SUBDIR = "vcf"
PRECALL_VCF_SUFFIX = ".fil.vcf"
PRECALL_GVCF_SUBDIR = "gvcf"
PRECALL_GVCF_SUFFIX = ".g.vcf"

#: Legacy-compatible 1-based output slot for AD. TRACE now reads named FORMAT fields.
TRACE_AD_SLOT = 2

_MISSING = "."


class VcfAdapterError(ValueError):
    """The source VCF is malformed, or a requested sample is not in it."""


@dataclass
class VcfSplitResult:
    """What `split_plate_vcf` produced, with every count's denominator stated.

    Attributes
    ----------
    workdir : Path
        The precall workdir; per-sample VCFs are under `workdir/vcf/`.
    vcf_paths : dict
        `{sample: path}` for every file written.
    n_records_in : int
        Data records read from the source VCF (the denominator for everything below).
    n_records_out : dict
        `{sample: records written}`. A sample with 0 is written as a valid header-only
        VCF, which is a real "no variants called" result — distinct from an absent file,
        which strict `precall` rejects.
    n_multiallelic : int
        Source records whose ALT carried more than one allele. See Notes: these are
        retained with allele-indexed GT/AD for TRACE genotype selection.
    ad_slot : int
        1-based FORMAT slot `AD` was placed in.
    """

    workdir: Path
    vcf_paths: dict[str, Path] = field(default_factory=dict)
    n_records_in: int = 0
    n_records_out: dict[str, int] = field(default_factory=dict)
    n_multiallelic: int = 0
    ad_slot: int = TRACE_AD_SLOT


@dataclass
class PrecallInputs:
    """Filesystem readiness of the required VCF/gVCF input pairs."""

    workdir: Path
    vcf_present: tuple[str, ...] = ()
    vcf_missing: tuple[str, ...] = ()
    gvcf_present: tuple[str, ...] = ()
    gvcf_missing: tuple[str, ...] = ()

    @property
    def ready(self) -> tuple[str, ...]:
        """Samples with BOTH a VCF and a gVCF — the only ones `precall` will process."""
        gvcf = set(self.gvcf_present)
        return tuple(s for s in self.vcf_present if s in gvcf)

    @property
    def silently_skipped(self) -> tuple[str, ...]:
        """Compatibility name for samples that fail strict precall input validation."""
        ready = set(self.ready)
        ordered = list(self.vcf_present) + [s for s in self.vcf_missing if s not in ready]
        return tuple(s for s in dict.fromkeys(ordered) if s not in ready)


def carries_alt(genotype: str) -> bool:
    """True when a sample's GT names at least one non-reference allele.

    This is a representation-level per-sample filter based on the called genotype.
    Repaired TRACE independently validates GT before selecting a clonal allele.

    `./.`, `.`, and an empty field are not calls and return False.
    """
    gt = genotype.split(":", 1)[0].strip()
    if not gt or gt in (_MISSING, "./.", ".|."):
        return False
    alleles = gt.replace("|", "/").split("/")
    for allele in alleles:
        allele = allele.strip()
        if allele != _MISSING and (not allele.isdigit() or int(allele) < 0):
            raise VcfAdapterError(f"invalid genotype {genotype!r}")
    return any(allele != _MISSING and int(allele) > 0 for allele in alleles)


def reorder_format_field(
    format_keys: Sequence[str],
    sample_fields: Sequence[str],
    *,
    key: str = "AD",
    slot: int = TRACE_AD_SLOT,
) -> tuple[list[str], list[str]]:
    """Move one FORMAT key (and its value) to a given position, keeping the rest in order.

    This pure list permutation retains all values and allele indexing. The default
    preserves historical output layout; the repaired algorithm reads FORMAT by name.

    Parameters
    ----------
    format_keys : sequence of str
        The record's FORMAT keys, e.g. `["GT", "PL", "DP", "AD"]`.
    sample_fields : sequence of str
        That sample's `:`-split values. VCF allows trailing fields to be dropped, so a
        short list is padded with `.` to the FORMAT length.
    key : str, default "AD"
        The key to move.
    slot : int, default `TRACE_AD_SLOT`
        Destination **1-based** FORMAT slot, preserving the historical default of 2.

    Returns
    -------
    tuple
        `(reordered_keys, reordered_values)`.

    Raises
    ------
    VcfAdapterError
        `key` is absent, the slot is below 1, keys repeat, or sample fields exceed
        the declared FORMAT width.
    """
    keys = list(format_keys)
    if len(keys) != len(set(keys)) or len(sample_fields) > len(keys):
        raise VcfAdapterError("duplicate FORMAT keys or excess sample fields")
    if key not in keys:
        raise VcfAdapterError(
            f"FORMAT {keys} does not contain {key!r}; cannot place it at slot {slot}"
        )
    if slot < 1:
        raise VcfAdapterError(f"FORMAT slot is 1-based and must be >= 1, got {slot}")
    values = list(sample_fields) + [_MISSING] * (len(keys) - len(sample_fields))
    values = values[: len(keys)]
    at = keys.index(key)
    moved_key, moved_value = keys.pop(at), values.pop(at)
    index = min(slot - 1, len(keys))
    keys.insert(index, moved_key)
    values.insert(index, moved_value)
    return keys, values


def _open_text(path: Path):
    """Open a plain or bgzipped/gzipped VCF as text. bgzip is gzip-compatible."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def _sample_header(line: str, path: str | Path) -> list[str]:
    header = line.rstrip("\r\n").split("\t")
    fixed = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    if len(header) < 10 or header[:9] != fixed:
        raise VcfAdapterError(f"{path}: invalid VCF sample header")
    samples = header[9:]
    if any(not sample for sample in samples) or len(samples) != len(set(samples)):
        raise VcfAdapterError(f"{path}: absent or duplicate VCF sample IDs in header")
    return header


def vcf_samples(path: str | Path) -> tuple[str, ...]:
    """Read exact source sample IDs without opening or overwriting any output."""
    with _open_text(Path(path)) as fh:
        for line in fh:
            if line.startswith("#CHROM\t"):
                return tuple(_sample_header(line, path)[9:])
            if not line.startswith("#") and line.strip():
                break
    raise VcfAdapterError(f"{path}: no #CHROM sample header")


def split_plate_vcf(
    vcf_path: str | Path,
    workdir: str | Path,
    *,
    samples: Iterable[str] | None = None,
    sample_map: Mapping[str, str] | None = None,
    ad_slot: int = TRACE_AD_SLOT,
) -> VcfSplitResult:
    """Split a multi-sample plate VCF into the per-sample `.fil.vcf` layout TRACE reads.

    Takes primitive arguments, not a `WGSConfig`: `core/` never imports `..config`, so
    every function here stays callable from a test with no config at all. That follows
    the house exemplar — `redi.core.wgs_check.compare_to_wgs` and `redi.core.purity` take
    DataFrames and primitives, and only `pipeline/step_*.py` unpacks the config object.
    The call is therefore `split_plate_vcf(config.donor.plate_vcfs[i], workdir, ...)`, with the
    unpacking living in `pipeline/_common.py` alongside `step_dir` and `record_for`.

    Parameters
    ----------
    vcf_path : str or Path
        The plate VCF. Plain, gzipped or bgzipped (bgzip is gzip-compatible, so both read
        through `gzip`); the **outputs are always plain text**, because TRACE opens them
        with a bare `open()`.
    workdir : str or Path
        Precall workdir. Files land in `workdir/vcf/{sample}.fil.vcf`.
    samples : iterable of str, optional
        Restrict output to these sample names (post-mapping). A requested name absent
        from the VCF raises — that is the silent-empty failure this guards.
    sample_map : mapping, optional
        `{vcf_column_name: output_sample_name}`. **Sample identity is the trap here**: a
        bcftools merge names its columns from `@RG SM` or the BAM path, which need not
        equal the manifest's `sample`. Columns absent from the map keep their own name.
    ad_slot : int, default `TRACE_AD_SLOT`
        Destination 1-based FORMAT slot for `AD`.

    Returns
    -------
    VcfSplitResult

    Raises
    ------
    VcfAdapterError
        No `#CHROM` header, a requested sample missing from the VCF, a record with fewer
        columns than the header, or missing/invalid genotype fields.

    Notes
    -----
    The source is streamed once with all output handles open, so peak memory is bounded by
    the number of samples (~96 handles for a plate), not by the record count (§5).

    Multiallelic ALT/GT/AD are retained consistently. TRACE chooses the allele named by
    GT and reports ambiguous genotypes as unknown for a clonal donor. A source
    header with duplicate identities or a many-to-one sample map is rejected.
    """
    vcf_path = Path(vcf_path)
    workdir = Path(workdir)
    out_dir = workdir / PRECALL_VCF_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = set(samples) if samples is not None else None
    result = VcfSplitResult(workdir=workdir, ad_slot=ad_slot)

    temporary_paths: dict[str, Path] = {}
    with ExitStack() as stack:
        def cleanup():
            for path in temporary_paths.values():
                path.unlink(missing_ok=True)
        stack.callback(cleanup)
        source = stack.enter_context(_open_text(vcf_path))
        meta: list[str] = []
        columns: list[str] | None = None
        handles: dict[str, object] = {}
        emit_names: list[str] = []

        for raw in source:
            line = raw.rstrip("\r\n")
            if line.startswith("##"):
                meta.append(line)
                continue
            if line.startswith("#CHROM"):
                if columns is not None:
                    raise VcfAdapterError(f"{vcf_path}: repeated #CHROM header")
                header = _sample_header(line, vcf_path)
                columns = header[9:]
                names = [(sample_map or {}).get(c, c) for c in columns]
                if len(names) != len(set(names)):
                    raise VcfAdapterError(f"{vcf_path}: duplicate or many-to-one sample names")
                if any(not n or n in (".", "..") or any(c in n for c in "/\\\t\r\n") for n in names):
                    raise VcfAdapterError(f"{vcf_path}: sample names must be safe filenames")
                if wanted is not None:
                    unknown = sorted(wanted - set(names))
                    if unknown:
                        raise VcfAdapterError(
                            f"{vcf_path}: requested sample(s) {unknown} are not in the VCF "
                            f"(it carries {len(names)} samples)"
                        )
                emit_names = [n for n in names if wanted is None or n in wanted]
                for name in emit_names:
                    path = out_dir / f"{name}{PRECALL_VCF_SUFFIX}"
                    fd, temp = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=out_dir)
                    temporary_paths[name] = Path(temp)
                    fh = stack.enter_context(os.fdopen(fd, "w", newline="\n"))
                    for m in meta:
                        fh.write(m + "\n")
                    fh.write("\t".join(header[:9] + [name]) + "\n")
                    handles[name] = fh
                    result.vcf_paths[name] = path
                    result.n_records_out[name] = 0
                continue
            if not line.strip():
                continue
            if columns is None:
                raise VcfAdapterError(f"{vcf_path}: data record before the #CHROM header line")

            cols = line.split("\t")
            if len(cols) != 9 + len(columns):
                raise VcfAdapterError(
                    f"{vcf_path}: record {cols[:2]!r} has {len(cols)} columns, "
                    f"header declares {9 + len(columns)}"
                )
            result.n_records_in += 1
            if "," in cols[4]:
                result.n_multiallelic += 1
            format_keys = cols[8].split(":")
            if "GT" not in format_keys or len(format_keys) != len(set(format_keys)):
                raise VcfAdapterError(f"{vcf_path}: missing GT or duplicate FORMAT key")
            for i, column in enumerate(columns):
                name = (sample_map or {}).get(column, column)
                fh = handles.get(name)
                if fh is None:
                    continue
                field_value = cols[9 + i]
                values = field_value.split(":")
                if len(values) > len(format_keys):
                    raise VcfAdapterError(f"{vcf_path}: too many sample FORMAT values")
                values += ["."] * (len(format_keys) - len(values))
                gt = values[format_keys.index("GT")]
                # Missing alleles are unresolved evidence, even when no called
                # ALT remains (0/., ./., .). Retain them for TRACE's explicit
                # unknown outcome; omitting them can manufacture an unedited call.
                has_alt = carries_alt(gt)
                incomplete = "." in gt.replace("|", "/").split("/")
                if not has_alt and not incomplete:
                    continue
                alleles = [int(a) for a in gt.replace("|", "/").split("/") if a != "."]
                if any(a > len(cols[4].split(",")) for a in alleles):
                    raise VcfAdapterError(f"{vcf_path}: GT exceeds ALT allele count")
                if "AD" in format_keys:
                    keys, values = reorder_format_field(format_keys, values, key="AD", slot=ad_slot)
                else:
                    keys = format_keys
                fh.write(
                    "\t".join(cols[:8] + [":".join(keys), ":".join(values)]) + "\n"
                )
                result.n_records_out[name] += 1

        if columns is None:
            raise VcfAdapterError(f"{vcf_path}: no #CHROM header line found")
        for name, fh in handles.items():
            fh.flush()
            fh.close()
        for name, path in temporary_paths.items():
            os.replace(path, result.vcf_paths[name])
    return result


def precall_inputs(workdir: str | Path, samples: Iterable[str]) -> PrecallInputs:
    """Enumerate required VCF/gVCF files before calling TRACE.

    This never substitutes called-variant depth for absent reference-depth evidence.
    The compatibility property `silently_skipped` identifies unready inputs; repaired
    TRACE raises for them rather than silently skipping.
    """
    workdir = Path(workdir)
    vcf_dir = workdir / PRECALL_VCF_SUBDIR
    gvcf_dir = workdir / PRECALL_GVCF_SUBDIR
    names = list(dict.fromkeys(samples))
    vcf_present = tuple(s for s in names if (vcf_dir / f"{s}{PRECALL_VCF_SUFFIX}").exists())
    gvcf_present = tuple(s for s in names if (gvcf_dir / f"{s}{PRECALL_GVCF_SUFFIX}").exists())
    return PrecallInputs(
        workdir=workdir,
        vcf_present=vcf_present,
        vcf_missing=tuple(s for s in names if s not in set(vcf_present)),
        gvcf_present=gvcf_present,
        gvcf_missing=tuple(s for s in names if s not in set(gvcf_present)),
    )
