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
    "CALLER_PLOIDY",
    "DEFAULT_FULL_ALLELE_FRACTION",
    "DEFAULT_MIN_ALT_READS",
    "MAG_AF_TAG",
    "MAG_MINOR_TAG",
    "MPILEUP_CALL_ALL_SITES_COMMAND",
    "MPILEUP_CALL_COMMAND",
    "MPILEUP_CALL_WINDOWS_COMMAND",
    "PLOIDY_CAVEAT",
    "PRECALL_GVCF_SUBDIR",
    "PRECALL_GVCF_SUFFIX",
    "PRECALL_VCF_SUBDIR",
    "PRECALL_VCF_SUFFIX",
    "CallerProvenance",
    "PrecallInputs",
    "VcfAdapterError",
    "VcfSplitResult",
    "allele_fraction",
    "alt_depth",
    "caller_provenance",
    "carries_alt",
    "precall_inputs",
    "reorder_format_field",
    "split_plate_vcf",
    "vcf_samples",
]

# --------------------------------------------------------------------------- #
# Defect 1 — the caller is haploid, and the GT is not the evidence             #
# --------------------------------------------------------------------------- #

#: Saccharomyces cerevisiae MAGESTIC clones are HAPLOID. The plate VCFs this adapter
#: reads were produced at bcftools' DIPLOID default, so a 5-33% minor allele is
#: emitted as `GT 0/0`; `carries_alt()` then drops the record, Stage 2 sees no
#: variant, and — Stage 2 having NO allele-fraction threshold of its own — the cell
#: is written as a bare `designed_variant_unedited`. **A real partial edit is
#: silently converted into a no-edit.** Magnitude measured at ~1.1 pp, so it does not
#: explain the 27%, but it is a correctness bug.
CALLER_PLOIDY = 1

#: ⚠️ **SUPERSEDED by `MPILEUP_CALL_WINDOWS_COMMAND` for the designed-site call.**
#: Kept, and kept accurate, because it documents what the EXISTING plate VCFs were made
#: with (minus the ploidy flag they do not actually carry — measured) and because the
#: genome-wide variants-only call is still the right input for the DISTAL spontaneous
#: and unintended-variant counts, which an `-R`-restricted call cannot see at all.
#:
#: 🔴 NN 255 measured that `--ploidy 1` **with `-v`** is the worst of the three options
#: for a minority allele private to one well: `-v` drops a site no sample was called
#: variant at, and the haploid model calls `GT 0` far further up the allele-fraction
#: range than the diploid one does. Use this for the distal axis; use
#: `MPILEUP_CALL_WINDOWS_COMMAND` for the designed site. See `PLOIDY_CAVEAT`.
MPILEUP_CALL_COMMAND = (
    "bcftools mpileup -q 20 -Q 20 -d 500 -a FORMAT/AD,FORMAT/DP "
    "-f {reference} {bams} | bcftools call -mv --ploidy 1 -Ov -o {out}"
)

#: 🔴 Stated plainly, because it is a limitation of this fix and not a footnote.
#:
#: NN 255: the old text opened "UNTESTED AGAINST THE REAL BINARY — `bcftools` is absent
#: from this toolchain". **That premise was false** (NN 252 F13): bcftools 1.21 / htslib
#: 1.21 is installed and `MPILEUP_CALL_COMMAND` runs clean (exit 0). The premise is
#: retired; its two real limits survive, now as MEASURED numbers rather than
#: expectations, plus a third that the measurement found and nobody predicted.
PLOIDY_CAVEAT = (
    "MEASURED against bcftools 1.21 (NN 252 F11/F13, NN 255 floor sweep over a flat 40x "
    "fixture, alt counts 0..40). Three limits, all measured, none to be assumed away. "
    "(1) `--ploidy 1` makes the caller choose between GT 0 and GT 1, so a 33% and a 10% "
    "minor allele are both called GT 0 — the ploidy flag makes the genotype CORRECT for "
    "a haploid, it does not recover minor alleles, which is why `min_allele_fraction` "
    "ignores GT entirely and reads AD. (2) `call -v` emits no record at a site where NO "
    "sample was called non-reference, so the EMISSION floor depends on whether the "
    "mixture is private to one well or co-occurs with a called well at the same "
    "position: measured at 40x in one sequence context, `-mv` (diploid default) emits a "
    "private mixture from 6/40 = 15% alt and `-mv --ploidy 1` only from 34/40 = 85%, "
    "while either emits a co-occurring one from 1/40 = 2.5%. 🔴 Those thresholds are "
    "CONTEXT-DEPENDENT, not constants: re-measured in a second 40x context (a different "
    "reference sequence and read length) the diploid private floor fell to <=4/40 = 10% "
    "and the `-mv --ploidy 1` floor to <=13/40 = 32.5%. What is stable across both "
    "contexts is the ORDERING and its sign, which is the part to rely on. (3) "
    "Consequence of (2) that inverts the naive "
    "advice: 🔴 **adopting `--ploidy 1` WITH `-v` makes the private-mixture floor ~5.7x "
    "WORSE than the diploid default it replaces.** `-A`/`--keep-alts` does not help, "
    "because `-v` drops the site before alts are kept (measured: identical 34/40 floor). "
    "The command that fixes it is `call -m --ploidy 1 -A` — WITHOUT `-v` — which "
    "measured a 1/40 = 2.5% private floor with the correct haploid GT and a real ALT "
    "base at GT 0. Use it RESTRICTED TO THE DONOR WINDOWS "
    "(`MPILEUP_CALL_WINDOWS_COMMAND`): all-sites emits one record per covered base, "
    "~1.2e7 per plate genome-wide, and step_2 re-reads each plate once per sample, so "
    "the unrestricted form is infeasible rather than merely large. Keep the genome-wide "
    "variants-only call for the DISTAL unintended-variant axis, which an -R VCF cannot "
    "see at all."
)

#: The command that actually recovers a private minor allele, measured. `-v` is the
#: problem, not the ploidy: see `PLOIDY_CAVEAT` limit (3).
#:
#: ⚠️ It emits EVERY COVERED POSITION. Measured: one record per covered base (100 records
#: over a 100 bp covered stretch), so genome-wide on S. cerevisiae that is ~1.2e7 records
#: per plate. `step_2._require_vcf_inputs` re-reads each plate VCF in full ONCE PER
#: SAMPLE, so at 96 samples that is ~1.2e9 record parses per plate in Python — not an
#: engineering inconvenience, an infeasibility. **Use the `-R` form below instead**, and
#: keep `MPILEUP_CALL_COMMAND`'s genome-wide variants-only call for the distal axis.
MPILEUP_CALL_ALL_SITES_COMMAND = (
    "bcftools mpileup -q 20 -Q 20 -d 500 -a FORMAT/AD,FORMAT/DP "
    "-f {reference} {bams} | bcftools call -m --ploidy 1 -A -Ov -o {out}"
)

#: 🔴 **The command to regenerate the designed-site plate VCFs with.** All-sites inside
#: the donor windows only, which is where the designed-variant call is made:
#:
#: * `-R {regions}` — a BED of the donor windows. ~480 windows x ~60 bp is ~29 kb, so the
#:   all-sites output is ~3e4 records rather than ~1.2e7, and the per-sample re-read is
#:   trivial again.
#: * no `-v` — so a site at which NO sample on the plate was called variant is still
#:   emitted, which is the only way a minority allele PRIVATE to one well survives
#:   (measured private floor 1/40 = 2.5% vs 6/40 and 34/40 for the `-v` forms).
#: * `--ploidy 1` — correct for a haploid clone; the minority allele then arrives as
#:   `GT 0` + a real ALT base + AD, which `min_allele_fraction` reads.
#: * `| bcftools norm -m-` — still recommended. The pinned TRACE now selects the ALT that
#:   GT names, so `step_2` no longer refuses a multi-allelic record (it did while TRACE
#:   `7c4c10f` read `"A,G"` as one literal allele; planning/268). But the `MAG_MINOR`
#:   variant id below is built from the whole ALT column, so a minority allele on a
#:   multi-allelic record cannot match a designed variant id in `step_5`. NN 255 measured
#:   `-A` NOT producing a multi-allelic record even at a genuine 12-G/12-T site, on one
#:   fixture — reassuring and not a guarantee, and `norm -m-` costs nothing.
#:
#: ⚠️ **What this command CANNOT see, stated so it is not discovered later:** anything
#: outside the donor windows. NN 251 defect 8's
#: `n_variants_discarded_outside_donor_window` becomes structurally 0 on an `-R` VCF, so
#: the DISTAL unintended-variant axis must keep coming from the genome-wide
#: `MPILEUP_CALL_COMMAND` output. Two VCFs, two questions.
MPILEUP_CALL_WINDOWS_COMMAND = (
    "bcftools mpileup -q 20 -Q 20 -d 500 -a FORMAT/AD,FORMAT/DP -R {regions} "
    "-f {reference} {bams} | bcftools call -m --ploidy 1 -A -Ou "
    "| bcftools norm -m- -f {reference} -Ov -o {out}"
)

#: INFO tags this adapter writes. Owned as constants so a consumer can grep for them.
MAG_AF_TAG = "MAG_AF"
MAG_MINOR_TAG = "MAG_MINOR"

#: 🔴 A fraction threshold with no COUNT floor selects for low coverage, and the depth
#: this program works at makes that decisive. Measured per-well depth, **named by
#: dataset because they differ**: `spg_20240702_POST` (n=480, `depth_chrIV_baseline`)
#: median 20.1x, q25 9.1x, 26.7% of wells below 10x; the 20231213 plate-2348/2349 cohort
#: (n=768, `ave_cov`) median 14.4x, q25 6.2x. At either, a 5% allele-fraction floor
#: rounds to a ONE-read minimum, and one read is a ~1-2% per-position false positive
#: from sequencing error alone (binomial, q=1e-3 to 5e-3). Two reads puts the
#: per-position false-positive rate below 1e-3 at 40x while sensitivity to a true 10%
#: allele is still 0.92. So a RESCUE also requires this many alt reads.
#: (Project precedent: an |log2FC| threshold with no baseMean floor inverted a +0.765
#: correlation to -0.076.)
DEFAULT_MIN_ALT_READS = 2

#: Below this allele fraction a RESCUED record is flagged `MAG_MINOR`. It matters
#: because TRACE's `_iter_vcf_records` reads only columns 0/1/3/4 — it never looks at
#: GT, FILTER, INFO or FORMAT — so any record this adapter writes is treated by Stage
#: 2 as a 100% call. Rescuing a 10% allele without flagging it would turn defect 1
#: into the opposite error: the live over-call class of 242 section 7a (`AD=29,1` read
#: as `designed_variant_installed`).
DEFAULT_FULL_ALLELE_FRACTION = 0.85

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
    n_rescued_by_af : int
        (sample, record) pairs emitted because their AD-derived allele fraction
        cleared `min_allele_fraction` despite a reference GT. Zero unless
        `min_allele_fraction` was passed, so the default path is unchanged.
    n_dropped_below_af : int
        Pairs with a non-reference AD but a fraction below the threshold. Counted,
        because "below threshold" and "no evidence" are different results.
    n_no_ad : int
        Pairs where the allele fraction was needed and `AD` was absent or
        unparseable. These are DROPPED, and counted so the drop is never silent.
    n_dropped_below_min_alt_reads : int
        Pairs whose fraction cleared `min_allele_fraction` but whose ALT read count was
        below `min_alt_reads`. Counted separately because "too few reads" and "too low a
        fraction" are different results at this cohort's 14.4x median depth, where a 5%
        fraction floor rounds to one read.
    n_af_computed : int
        (sample, record) pairs where an allele fraction COULD be computed from `AD` —
        **the denominator for every `MAG_AF` number below.** Without it, a zero
        `n_flagged_minor` cannot be distinguished from "AD was never usable".
    n_flagged_minor : int
        Pairs tagged `MAG_MINOR` (fraction below `full_allele_fraction`), on either
        branch.
    n_flagged_minor_carrying_alt : int
        🔴 The subset of the above whose **GT already named the alt** — i.e. records the
        old code passed through completely untagged. This is NN 252's F11b: at the
        diploid default a 10% allele is `GT 0/1`, `carries_alt` is True, and the whole
        annotation block was skipped. A non-zero value here is the over-call that
        defect 1 exists to catch, finally visible.
    n_no_ad_carrying_alt : int
        Pairs with a non-reference GT whose `AD` was unusable, so no fraction could be
        computed and the record is emitted UNFLAGGED. Counted, because this is the one
        residual route by which a minority allele still reaches Stage 2 silently.
    minor_allele_variants : dict
        `{sample: (variant_id, ...)}` for every record below `full_allele_fraction`,
        **on either branch**. 🔴 Stage 2 cannot see a fraction, so this is the only
        record that a call it will treat as 100% was in fact a minority allele.
    """

    workdir: Path
    vcf_paths: dict[str, Path] = field(default_factory=dict)
    n_records_in: int = 0
    n_records_out: dict[str, int] = field(default_factory=dict)
    n_multiallelic: int = 0
    ad_slot: int = TRACE_AD_SLOT
    n_rescued_by_af: int = 0
    n_dropped_below_af: int = 0
    n_dropped_below_min_alt_reads: int = 0
    n_no_ad: int = 0
    n_af_computed: int = 0
    n_flagged_minor: int = 0
    n_flagged_minor_carrying_alt: int = 0
    n_no_ad_carrying_alt: int = 0
    minor_allele_variants: dict[str, tuple[str, ...]] = field(default_factory=dict)


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


def allele_fraction(
    format_keys: Sequence[str], sample_fields: Sequence[str]
) -> float | None:
    """Non-reference allele fraction from `AD`, or `None` when `AD` is unusable.

    🔴 **This deliberately ignores GT.** GT is the thing defect 1 is about: at the
    diploid default a haploid 33% alt is `0/0`, and even at `--ploidy 1` it is `0`.
    The evidence is the allele depths, so the evidence is what this reads.

    `AD` is `ref,alt1[,alt2...]`; the fraction is `sum(alts) / sum(all)`. Returns
    `None` for a missing key, a missing value (`.`), a non-integer, or total depth 0
    — all of which are "cannot compute", never "zero".
    """
    keys = list(format_keys)
    if "AD" not in keys:
        return None
    index = keys.index("AD")
    values = list(sample_fields)
    if len(values) <= index:
        return None
    raw = values[index].strip()
    if not raw or raw == _MISSING:
        return None
    parts = raw.split(",")
    try:
        depths = [int(p) for p in parts]
    except ValueError:
        return None
    total = sum(depths)
    if total <= 0 or len(depths) < 2:
        return None
    return sum(depths[1:]) / total


def alt_depth(
    format_keys: Sequence[str], sample_fields: Sequence[str]
) -> int | None:
    """Total non-reference read count from `AD`, or `None` when `AD` is unusable.

    The COUNT, not the fraction. A fraction threshold alone selects for low coverage —
    see `DEFAULT_MIN_ALT_READS` — so the count is read separately rather than
    back-derived from `af * DP`, which would round.
    """
    keys = list(format_keys)
    if "AD" not in keys:
        return None
    index = keys.index("AD")
    values = list(sample_fields)
    if len(values) <= index:
        return None
    raw = values[index].strip()
    if not raw or raw == _MISSING:
        return None
    parts = raw.split(",")
    try:
        depths = [int(p) for p in parts]
    except ValueError:
        return None
    if len(depths) < 2:
        return None
    return sum(depths[1:])


@dataclass(frozen=True)
class CallerProvenance:
    """What the VCF's own header says about how it was called.

    🔴 **Detected, not declared.** `vcf_adapter`'s whole defect-1 story turns on which
    ploidy the plate VCFs were produced at, and NN 251 asserted it in a docstring while
    NN 252 measured the opposite. `bcftools` writes `##bcftools_callCommand=` into every
    VCF it calls, and `bcftools merge` preserves it (measured, NN 255), so the question
    has a factual answer in the file.

    Attributes
    ----------
    command : str or None
        The `##bcftools_callCommand=` value, or `None` when the header carries none
        (a hand-written VCF, or one from another caller) — which is `ploidy=None`, i.e.
        "unobserved", never "diploid".
    ploidy : int or None
        `1` for `--ploidy 1`; `2` for a `bcftools call` with no `--ploidy` (its
        documented default); `None` when unobservable or when `--ploidy` named a file or
        an assembly rather than a number.
    variants_only : bool or None
        True when `-v`/`--variants-only` was used. 🔴 This, not the ploidy, is what sets
        the private-mixture emission floor — see `PLOIDY_CAVEAT` limit (3).
    keep_alts : bool or None
        True when `-A`/`--keep-alts` was used.
    n_call_commands : int
        How many `##bcftools_callCommand=` lines the header carries. `bcftools merge`
        keeps only the FIRST input's, so a plate assembled from inputs called different
        ways reports 1 and hides the disagreement; a value above 1 means the header
        itself disagrees and must be refused rather than reconciled.
    """

    command: str | None = None
    ploidy: int | None = None
    variants_only: bool | None = None
    keep_alts: bool | None = None
    n_call_commands: int = 0

    @property
    def recovers_private_minor_alleles(self) -> bool:
        """True only for a call that emits a site no sample was called variant at.

        This is the property that decides whether a minority allele private to one well
        can be recovered AT ALL, and it is `not variants_only` — the ploidy is
        irrelevant to it. Measured in the NN 255 floor sweep.
        """
        return self.variants_only is False


def caller_provenance(vcf_path: str | Path) -> CallerProvenance:
    """Read the caller's own command line out of a VCF header.

    Scans the `##` block only (stops at `#CHROM`), so it is O(header) on a plate VCF.
    """
    commands: list[str] = []
    with _open_text(Path(vcf_path)) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.startswith("##"):
                break
            if line.startswith("##bcftools_callCommand="):
                commands.append(line.split("=", 1)[1])
    if not commands:
        return CallerProvenance(n_call_commands=0)
    command = commands[0]
    tokens = command.split()
    # 🔴 `--ploidy=1` and `--ploidy 1` are both valid getopt_long syntax. Reading only
    # the space-separated form would parse the `=` form as the DIPLOID DEFAULT, which
    # bypasses `step_2`'s ploidy-1 refusal — the same shape of hole as F2's `svRate`.
    ploidy: int | None = 2
    value: str | None = None
    for i, token in enumerate(tokens):
        if token == "--ploidy":
            value = tokens[i + 1] if i + 1 < len(tokens) else ""
            break
        if token.startswith("--ploidy="):
            value = token.split("=", 1)[1]
            break
    if value is not None:
        # A non-numeric value names a ploidy FILE or a species preset, which this
        # function cannot resolve: that is "unobserved", never "diploid".
        ploidy = int(value) if value.isdigit() else None
    # 🔴 `-O` is the OUTPUT-TYPE flag and its value is usually attached (`-Ov`, `-Oz`,
    # `-Ob`, `-Ou`), so a naive scan for the letter `v` in every short cluster reads
    # `-Ov` as `--variants-only`. Measured: that mis-read `call -m --ploidy 1 -A -Ov`
    # as variants-only, i.e. as UNABLE to recover a private minor allele when it is in
    # fact the one command that can. `-O` clusters are therefore dropped before the
    # letters are collected.
    letters: set[str] = set()
    for token in tokens:
        if not token.startswith("-") or token.startswith("--"):
            continue
        if token.startswith("-O"):
            continue
        letters.update(token[1:])
    variants_only = ("--variants-only" in tokens) or ("v" in letters)
    keep_alts = ("--keep-alts" in tokens) or ("A" in letters)
    return CallerProvenance(
        command=command, ploidy=ploidy, variants_only=variants_only,
        keep_alts=keep_alts, n_call_commands=len(commands),
    )


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
    """Open a plain or bgzipped/gzipped VCF as text. bgzip is gzip-compatible.

    🔴 A **BCF** is refused by name rather than attempted. This cohort's plate files are
    produced with `call ... -Ob` (measured in their own `##bcftools_callCommand`) and
    converted to `.vcf.gz` afterwards, so a `.bcf` reaching here is a wiring mistake;
    opening one as text raises a `UnicodeDecodeError` from deep inside `gzip`, which
    names neither the file nor the cause.
    """
    if path.suffix == ".bcf":
        raise VcfAdapterError(
            f"{path} is a BCF (binary). This adapter is text-only, because TRACE's "
            "`donor_outcome.precall` opens its inputs with a bare `open()`. Convert "
            f"first: `bcftools view -Ov -o {path.with_suffix('.vcf')} {path}`."
        )
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
    min_allele_fraction: float | None = None,
    min_alt_reads: int = DEFAULT_MIN_ALT_READS,
    full_allele_fraction: float = DEFAULT_FULL_ALLELE_FRACTION,
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
    min_allele_fraction : float, optional
        When set, a record whose GT is reference but whose **AD-derived allele
        fraction** is `>= min_allele_fraction` is emitted anyway, with
        `MAG_AF=<fraction>` appended to its INFO and `MAG_MINOR` added when the
        fraction is below `full_allele_fraction`. `None` (the default) leaves the
        GT-only behaviour byte-identical, so no existing output can move.

        🔴 This is the half of defect 1 that actually recovers a partial edit.
        `--ploidy 1` (see `MPILEUP_CALL_COMMAND`) makes the GENOTYPE correct for a
        haploid; it does not make a 33% allele a `1`. See `PLOIDY_CAVEAT`.

        It governs the RESCUE only. **Annotation of a minority allele is unconditional**
        — see `min_alt_reads` and the F11b comment in the body.
    min_alt_reads : int, default `DEFAULT_MIN_ALT_READS`
        Minimum ALT read count for a RESCUE. A fraction floor with no count floor
        selects for low coverage: at this cohort's 14.4x median depth a 5% fraction
        rounds to one read, whose per-position false-positive rate from sequencing error
        alone is ~1-2%.
    full_allele_fraction : float, default `DEFAULT_FULL_ALLELE_FRACTION`
        Below this, a record is flagged `MAG_MINOR` and listed in
        `VcfSplitResult.minor_allele_variants` — **whether or not its GT carried the
        alt**. TRACE reads only columns 0/1/3/4, so without this flag a 10% allele
        reaches Stage 2 as a full call. A record at or above this fraction is written
        byte-identically to before, so the only outputs that move are the wrong ones.

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
                info = cols[7]
                # Missing alleles are unresolved evidence, even when no called
                # ALT remains (0/., ./., .). Retain them for TRACE's explicit
                # unknown outcome; omitting them can manufacture an unedited call.
                has_alt = carries_alt(gt)
                incomplete = "." in gt.replace("|", "/").split("/")
                af = allele_fraction(format_keys, values)
                n_alt = alt_depth(format_keys, values)
                if af is not None:
                    result.n_af_computed += 1
                rescued = False
                if not has_alt and not incomplete:
                    # GT says reference. Defect 1: at `--ploidy 1` that is exactly what
                    # a 5-33% haploid minor allele looks like, so fall back to the
                    # ALLELE DEPTHS, which are the actual evidence. RESCUE adds a record
                    # TRACE will treat as a real variant, so it stays opt-in.
                    if min_allele_fraction is None:
                        continue
                    if af is None:
                        result.n_no_ad += 1
                        continue
                    if af < min_allele_fraction:
                        result.n_dropped_below_af += 1
                        continue
                    if n_alt is not None and n_alt < min_alt_reads:
                        # A fraction floor with no COUNT floor selects for low coverage;
                        # see DEFAULT_MIN_ALT_READS for the measured false-positive rate.
                        result.n_dropped_below_min_alt_reads += 1
                        continue
                    result.n_rescued_by_af += 1
                    rescued = True
                alleles = [int(a) for a in gt.replace("|", "/").split("/") if a != "."]
                if any(a > len(cols[4].split(",")) for a in alleles):
                    raise VcfAdapterError(f"{vcf_path}: GT exceeds ALT allele count")
                # ------------------------------------------------------------------ #
                # 🔴 NN 255 F11b — ANNOTATION IS UNCONDITIONAL, on BOTH branches.     #
                # ------------------------------------------------------------------ #
                # Before NN 255 this block lived inside `if not carries_alt`, so on a
                # DIPLOID-DEFAULT plate VCF — which `vcf_adapter`'s own defect-1
                # docstring says is what this cohort has — a 10% allele is `GT 0/1`,
                # `carries_alt` is True, and the record was passed through COMPLETELY
                # UNTAGGED. Measured by NN 252 on a real `call -mv` plate VCF carrying a
                # 10% and a 33% minor allele: `n_rescued_by_af=0`,
                # `minor_allele_variants={}`, no INFO tag. The over-call defect 1 exists
                # to catch was therefore unflagged on precisely the input defect 1
                # describes — a guard that only fired on a non-default branch.
                #
                # Annotation is safe to run unconditionally because it ADDS INFO and
                # drops nothing: TRACE's genotype selection reads GT/AD from FORMAT, not
                # INFO, so no INFO tag can change what Stage 2 does. Only the RESCUE
                # above changes which records exist, and that is still opt-in.
                #
                # A clean, high-fraction record is left byte-identical (no tag at all),
                # so the only outputs that move are the ones that were wrong.
                if af is not None and (af < full_allele_fraction or rescued):
                    tags = [f"{MAG_AF_TAG}={af:.3f}"]
                    if af < full_allele_fraction:
                        tags.append(MAG_MINOR_TAG)
                        result.n_flagged_minor += 1
                        if has_alt:
                            result.n_flagged_minor_carrying_alt += 1
                        v_id = f"{cols[0]}_{cols[1]}_{cols[3]}_to_{cols[4]}"
                        result.minor_allele_variants[name] = (
                            result.minor_allele_variants.get(name, ()) + (v_id,)
                        )
                    info = ";".join(
                        ([] if info in ("", _MISSING) else [info]) + tags
                    )
                elif af is None and has_alt:
                    # The one residual silent route: a non-reference GT whose AD is
                    # unusable cannot be evaluated, so it ships unflagged. Counted.
                    result.n_no_ad_carrying_alt += 1
                if "AD" in format_keys:
                    keys, values = reorder_format_field(format_keys, values, key="AD", slot=ad_slot)
                else:
                    keys = format_keys
                fh.write(
                    "\t".join(
                        cols[:7] + [info] + [":".join(keys), ":".join(values)]
                    )
                    + "\n"
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
