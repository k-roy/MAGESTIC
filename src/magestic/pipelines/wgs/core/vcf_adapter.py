"""
magestic.pipelines.wgs.core.vcf_adapter — plate VCF -> the per-sample layout TRACE hardcodes.

`donor_outcome.precall` reads `{workdir}/vcf/{sample}.fil.vcf` and
`{workdir}/gvcf/{sample}.g.vcf` — relative subpaths and suffixes hardcoded in its body,
opened with plain `open()` (so **plain text only**, no bgzip). What this program actually
has is six per-plate **96-sample bcftools merges** (`mpileup -q 20 -Q 20 -d 500
-a FORMAT/AD,FORMAT/DP | call -mv`, `FORMAT = GT:PL:DP:AD`).

This module does the **re-representation** half of the bridge, and only that half:

* 96-sample merge -> one plain-text `.fil.vcf` per sample;
* `AD` moved from FORMAT slot 4 (bcftools) to slot 2, which is where TRACE's *positional*
  gVCF parser reads reference depth (`cols[9].split(":")[1].split(",")[0]`);
* records restricted to the samples that actually carry a non-reference allele.

🔴 **It writes nothing into `{workdir}/gvcf/`, deliberately.** No gVCF exists anywhere in
the cohort (GATE 2), and a re-keyed bcftools merge dropped into that directory would be a
*fabricated* gVCF from which TRACE would happily read reference depth. Synthesizing
per-position reference depth from a BAM pileup is **deriving a new quantity**, which by
this module's scope test is an algorithm and belongs in TRACE, not in a magestic
workaround. Until TRACE grows that primitive, every in-window-uncovered call collapses to
`unknown` — and `unknown` must never be read as `unedited`.

The consequence of a missing gVCF is not an error: `precall` hits
`if not vcf_path.exists() or not gvcf_path.exists(): continue` and drops the sample with
no warning, no stderr line and no count. Miss every sample and it writes a header-only
file and returns normally. `precall_inputs()` enumerates exactly that set so `qc` can
fail the run loudly instead.
"""

from __future__ import annotations

import gzip
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
]

#: The layout `donor_outcome.precall` hardcodes. Owned here, not left implicit.
PRECALL_VCF_SUBDIR = "vcf"
PRECALL_VCF_SUFFIX = ".fil.vcf"
PRECALL_GVCF_SUBDIR = "gvcf"
PRECALL_GVCF_SUFFIX = ".g.vcf"

#: 1-based FORMAT slot TRACE's positional gVCF parser reads `AD` from.
#: `sample_field.split(":")[1]` is the SECOND field, so slot 2. Stated 1-based to match
#: the spec's language and `config.donor.ad_format_slot_out`, which removes an off-by-one
#: at the boundary; conversion to a list index happens inside `reorder_format_field`.
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
        which `precall` would silently skip.
    n_multiallelic : int
        Source records whose ALT carried more than one allele. See Notes: these are
        passed through unchanged and are a known hazard, not a handled case.
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
    """Which samples `precall` will actually process, and which it will silently drop."""

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
        """Requested samples `precall` would drop with no warning and no count."""
        ready = set(self.ready)
        ordered = list(self.vcf_present) + [s for s in self.vcf_missing if s not in ready]
        return tuple(s for s in dict.fromkeys(ordered) if s not in ready)


def carries_alt(genotype: str) -> bool:
    """True when a sample's GT names at least one non-reference allele.

    This is the per-sample record filter, and it is the equivalent of
    `bcftools view -s <sample> -c1`: the information is already in the GT field, so
    selecting on it re-represents the merge rather than deriving anything new. It matters
    because TRACE's `_iter_vcf_records` reads only columns 0/1/3/4 — it never looks at
    GT, and it deliberately ignores FILTER — so **without this filter every sample would
    inherit every variant called anywhere on the plate.**

    `./.`, `.`, and an empty field are not calls and return False.
    """
    gt = genotype.split(":", 1)[0].strip()
    if not gt or gt in (_MISSING, "./.", ".|."):
        return False
    for allele in gt.replace("|", "/").split("/"):
        allele = allele.strip()
        if allele and allele != _MISSING and allele != "0":
            return True
    return False


def reorder_format_field(
    format_keys: Sequence[str],
    sample_fields: Sequence[str],
    *,
    key: str = "AD",
    slot: int = TRACE_AD_SLOT,
) -> tuple[list[str], list[str]]:
    """Move one FORMAT key (and its value) to a given position, keeping the rest in order.

    bcftools emits `GT:PL:DP:AD`, putting `AD` in slot 4; TRACE's gVCF parser reads
    reference depth *positionally* from slot 2 (`split(":")[1].split(",")[0]`). This is a
    pure list permutation — no value is recomputed — so it is orchestration, and it is
    unit-testable without a VCF.

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
        Destination **1-based** FORMAT slot: `slot=2` puts `AD` in the field TRACE's
        positional parser reads. 1-based because that is how the spec, TRACE's own
        docstring and `config.donor.ad_format_slot_out` all state it.

    Returns
    -------
    tuple
        `(reordered_keys, reordered_values)`.

    Raises
    ------
    VcfAdapterError
        `key` is absent from `format_keys`, or `slot` is below 1. Raising is deliberate:
        a record silently written without `AD` in the expected slot would make TRACE read
        some other field's first comma element as the reference depth.

    Notes
    -----
    Only the gVCF-shaped file is parsed positionally — TRACE's *VCF* parser reads columns
    0/1/3/4 and never touches FORMAT. The permutation is applied to both outputs anyway so
    the layout is uniform, and so the helper is exercised on the path that exists today.
    """
    keys = list(format_keys)
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
        columns than the header, or a record whose FORMAT has no `AD`.

    Notes
    -----
    The source is streamed once with all output handles open, so peak memory is bounded by
    the number of samples (~96 handles for a plate), not by the record count (§5).

    🔴 **Multi-allelic records are passed through unchanged and counted, not resolved.**
    TRACE reads column 4 whole, so an ALT of `"A,G"` reaches `bucket_variant` as a literal
    allele. Choosing one allele per sample would be deriving a call, so this module
    refuses to; normalize upstream (`bcftools norm -m-`) and check
    `VcfSplitResult.n_multiallelic` is 0 before trusting Stage 2.
    """
    vcf_path = Path(vcf_path)
    workdir = Path(workdir)
    out_dir = workdir / PRECALL_VCF_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = set(samples) if samples is not None else None
    result = VcfSplitResult(workdir=workdir, ad_slot=ad_slot)

    with ExitStack() as stack:
        source = stack.enter_context(_open_text(vcf_path))
        meta: list[str] = []
        columns: list[str] | None = None
        handles: dict[str, object] = {}
        emit_names: list[str] = []

        for raw in source:
            line = raw.rstrip("\n")
            if line.startswith("##"):
                meta.append(line)
                continue
            if line.startswith("#CHROM"):
                header = line.split("\t")
                if len(header) < 10:
                    raise VcfAdapterError(f"{vcf_path}: #CHROM line carries no sample columns")
                columns = header[9:]
                names = [(sample_map or {}).get(c, c) for c in columns]
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
                    fh = stack.enter_context(path.open("w", newline="\n"))
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
            if len(cols) < 9 + len(columns):
                raise VcfAdapterError(
                    f"{vcf_path}: record {cols[0]}:{cols[1]} has {len(cols)} columns, "
                    f"header declares {9 + len(columns)}"
                )
            result.n_records_in += 1
            if "," in cols[4]:
                result.n_multiallelic += 1
            format_keys = cols[8].split(":")
            for i, column in enumerate(columns):
                name = (sample_map or {}).get(column, column)
                fh = handles.get(name)
                if fh is None:
                    continue
                field_value = cols[9 + i]
                if not carries_alt(field_value):
                    continue
                keys, values = reorder_format_field(
                    format_keys, field_value.split(":"), key="AD", slot=ad_slot
                )
                fh.write(
                    "\t".join(cols[:8] + [":".join(keys), ":".join(values)]) + "\n"
                )
                result.n_records_out[name] += 1

        if columns is None:
            raise VcfAdapterError(f"{vcf_path}: no #CHROM header line found")
    return result


def precall_inputs(workdir: str | Path, samples: Iterable[str]) -> PrecallInputs:
    """Enumerate which samples `precall` will process and which it will silently skip.

    TRACE drops a sample whose `.fil.vcf` OR `.g.vcf` is missing with a bare `continue` —
    no warning, no count — so an output TSV missing half a cohort looks exactly like a
    complete one. `qc` compares `len(precall(...))` against the samplesheet row count and
    uses `silently_skipped` to say *which* samples went missing and why.

    No gVCF exists in this program today (GATE 2), so `gvcf_missing` is expected to be the
    full sample list until TRACE grows a reference-depth primitive. That is a gate, not a
    bug, and it must be reported rather than worked around.
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
