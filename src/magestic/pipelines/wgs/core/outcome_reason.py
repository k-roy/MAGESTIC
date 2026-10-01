"""magestic.pipelines.wgs.core.outcome_reason — WHY a Stage-2 outcome is what it is.

THE DISTINCTION THIS MODULE EXISTS FOR (defect 4)
--------------------------------------------------
`donor_outcome.precall_sample` ends in a four-way `if`, and its last branch is a bare
string::

    elif flag_covered >= quantification_window_size * quantification_window_cov_perc:
        flag = "designed_variant_unedited"
    else:
        flag = "unknown"

`unknown` therefore means **"I could not tell"**, and it is the only value in the
table that does. It records nothing about *why*: no gVCF at all, a gVCF whose
reference blocks were not expanded (defect 5), depth below `ne_threshold`, a
quantification window off the end of the contig, and a genuinely unsequenced locus
all collapse to the same five letters. Downstream, a reader who tabulates
`designed_variant_unedited` against everything else silently folds "uncallable" into
"no edit" — **the single most load-bearing distinction in this analysis.**

Worse, a sample whose `.fil.vcf` or `.g.vcf` is missing gets **no row at all**:
`precall` hits `if not vcf_path.exists() or not gvcf_path.exists(): continue`. An
output table missing half a cohort looks exactly like a complete one.

WHAT THIS MODULE GUARANTEES
---------------------------
* Every outcome string maps to exactly one machine-readable `reason`.
* Every explanation carries `callable_`, a boolean that is **False for every
  `unknown`** and True for every determinate call. "Uncallable" is never
  indistinguishable from "unedited".
* Every sample the driver did not process gets a row (`rows_for_missing_inputs`),
  with the reason naming which input was absent.
* The coverage numbers quoted in the reason come from
  `gvcf_adapter.window_reference_coverage`, which is the same arithmetic TRACE's
  `flag_covered` attempts. The test suite asserts the two never disagree: a reason
  that contradicted the flag would be a new silent defect, not a fix.

This module derives nothing about the biology. It re-states, in a stable vocabulary,
a decision TRACE already made plus the inputs TRACE already read — so it is an
adapter, and it is magestic's while TRACE is pinned.

THE TRACE-SIDE FIX
------------------
`SampleOutcome` should carry `flag_covered`, `n_window_positions` and a `reason`
field, and `precall` should emit a row for a skipped sample instead of `continue`.
Then this module collapses to a vocabulary mapping.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from .gvcf_adapter import window_reference_coverage

__all__ = [
    "REASON_BLOCKS_NOT_EXPANDED",
    "REASON_GVCF_ABSENT",
    "REASON_INSTALLED",
    "REASON_PARTIAL",
    "REASON_QWIN_UNCOVERED",
    "REASON_UNEDITED_COVERED",
    "REASON_UNRECOGNISED_FLAG",
    "REASON_VCF_ABSENT",
    "REASON_VCF_AND_GVCF_ABSENT",
    "OutcomeExplanation",
    "explain_outcome",
    "quantification_window",
    "rows_for_missing_inputs",
]

#: Determinate outcomes (`callable_ is True`).
REASON_INSTALLED = "installed_designed_variant_reconstructs_design"
REASON_PARTIAL = "partial_designed_variant_closer_to_design_than_wt"
REASON_UNEDITED_COVERED = "unedited_quantification_window_covered"

#: Indeterminate outcomes (`callable_ is False`). Each names a DIFFERENT cause.
REASON_QWIN_UNCOVERED = "uncallable_quantification_window_uncovered"

#: 🔴 The reason that keeps this module from contradicting its own numbers.
#: TRACE's `flag_covered` compares only POS, so an UNEXPANDED compressed gVCF yields
#: `unknown` even when its reference blocks cover the whole window (defect 5). Without
#: this value, `explain_outcome` would print "window uncovered" beside "covered 10/10"
#: — a reason that argues with its own evidence, which is a new silent defect rather
#: than a fix. The remedy is `gvcf_adapter.expand_reference_blocks`, named in the
#: string so the reader knows what to do.
REASON_BLOCKS_NOT_EXPANDED = "uncallable_gvcf_reference_blocks_not_expanded"
REASON_GVCF_ABSENT = "uncallable_gvcf_absent"
REASON_VCF_ABSENT = "uncallable_vcf_absent"
REASON_VCF_AND_GVCF_ABSENT = "uncallable_vcf_and_gvcf_absent"
REASON_UNRECOGNISED_FLAG = "uncallable_unrecognised_outcome_flag"

_DETERMINATE = {
    "designed_variant_installed": REASON_INSTALLED,
    "designed_variant_partially_installed": REASON_PARTIAL,
    "designed_variant_unedited": REASON_UNEDITED_COVERED,
}


@dataclass(frozen=True)
class OutcomeExplanation:
    """One Stage-2 outcome, with its reason and the coverage that produced it.

    Attributes
    ----------
    sample : str
        May be `""` when explaining a bare flag rather than a row.
    flag : str
        The outcome flag WITHOUT its `;synthetic_error_installed` /
        `;spontaneous_mutations` overlays, which do not affect callability.
    reason : str
        One of this module's `REASON_*` constants. Machine-readable and stable.
    callable_ : bool
        `True` only for a determinate call. **Trailing underscore** because
        `callable` is a builtin; the name is deliberate, not a typo.
    n_window_positions : int
        Size of the quantification window (its denominator).
    n_window_positions_covered : int
        Distinct positions in it at depth >= `ne_threshold`, honouring gVCF `END`.
    qwin_start, qwin_end : int
        The window, so a reader can check it against the contig.
    overlays : tuple of str
        The `;`-separated overlays stripped from `flag`, preserved rather than lost.
    """

    flag: str
    reason: str
    callable_: bool
    sample: str = ""
    n_window_positions: int = 0
    n_window_positions_covered: int = 0
    qwin_start: int = 0
    qwin_end: int = 0
    overlays: tuple[str, ...] = ()

    @property
    def window_fraction_covered(self) -> float:
        """Covered fraction, or 0.0 for an empty window (never a ZeroDivisionError)."""
        if self.n_window_positions <= 0:
            return 0.0
        return self.n_window_positions_covered / self.n_window_positions


def quantification_window(
    donor_start: int, donor_end: int, size: int = 10
) -> tuple[int, int]:
    """TRACE's own quantification window, restated so a caller can report it.

    Copied arithmetic, deliberately::

        qwin_start = math.ceil(dend - (dend - dstart) / 2 - quantification_window_size / 2)
        qwin_end   = qwin_start + quantification_window_size - 1

    Restating it here rather than re-deriving a "better" centre is the point: the
    reason must describe the window TRACE actually used, not one this module prefers.
    """
    start = math.ceil(donor_end - (donor_end - donor_start) / 2 - size / 2)
    return start, start + size - 1


def explain_outcome(
    flag: str,
    *,
    gvcf_path: str | Path,
    chrom: str,
    donor_start: int,
    donor_end: int,
    sample: str = "",
    quantification_window_size: int = 10,
    quantification_window_cov_perc: float = 1.0,
    ne_threshold: int = 1,
    ad_slot: int = 2,
) -> OutcomeExplanation:
    """Explain one Stage-2 outcome flag.

    The flag is taken as given — this function never second-guesses TRACE's call. It
    only (a) names the reason in a stable vocabulary and (b) attaches the coverage
    numbers, so that `unknown` stops being opaque.

    `gvcf_path` not existing is reported as `REASON_GVCF_ABSENT` regardless of the
    flag, because in that state no flag TRACE produced can be trusted.
    """
    base, _, rest = str(flag).partition(";")
    overlays = tuple(o for o in rest.split(";") if o)
    qs, qe = quantification_window(donor_start, donor_end, quantification_window_size)
    path = Path(gvcf_path)

    if not path.exists():
        return OutcomeExplanation(
            flag=base, reason=REASON_GVCF_ABSENT, callable_=False, sample=sample,
            n_window_positions=quantification_window_size,
            n_window_positions_covered=0, qwin_start=qs, qwin_end=qe, overlays=overlays,
        )

    covered = window_reference_coverage(
        path, chrom=chrom, window_start=qs, window_end=qe,
        ne_threshold=ne_threshold, ad_slot=ad_slot,
    )
    common = {
        "sample": sample,
        "n_window_positions": quantification_window_size,
        "n_window_positions_covered": covered,
        "qwin_start": qs,
        "qwin_end": qe,
        "overlays": overlays,
    }
    if base in _DETERMINATE:
        return OutcomeExplanation(
            flag=base, reason=_DETERMINATE[base], callable_=True, **common
        )
    if base == "unknown":
        # 🔴 The reason is derived from the COVERAGE, not restated from the flag.
        # Deriving it from the flag would make this function unable to disagree with
        # TRACE, and a test that it "agrees" would then be a tautology proving nothing.
        # The two genuinely differ in exactly one case: an unexpanded compressed gVCF,
        # where TRACE says `unknown` (it compares only POS) while the END-aware count
        # says the window IS covered. That case gets its own reason and names the fix.
        required = quantification_window_size * quantification_window_cov_perc
        if covered >= required and required > 0:
            return OutcomeExplanation(
                flag=base, reason=REASON_BLOCKS_NOT_EXPANDED, callable_=False, **common
            )
        return OutcomeExplanation(
            flag=base, reason=REASON_QWIN_UNCOVERED, callable_=False, **common
        )
    return OutcomeExplanation(
        flag=base, reason=REASON_UNRECOGNISED_FLAG, callable_=False, **common
    )


def rows_for_missing_inputs(
    samples: list[str] | tuple[str, ...],
    *,
    processed: list[str] | tuple[str, ...],
    workdir: str | Path,
    vcf_subdir: str = "vcf",
    vcf_suffix: str = ".fil.vcf",
    gvcf_subdir: str = "gvcf",
    gvcf_suffix: str = ".g.vcf",
) -> tuple[OutcomeExplanation, ...]:
    """One explicit row per requested sample `precall` did NOT return an outcome for.

    `precall` drops such a sample with a bare `continue`. This turns that silence into
    a row whose `reason` names which input was absent, so the cohort table's row count
    always equals the manifest's and a missing half-cohort is visible.
    """
    workdir = Path(workdir)
    done = set(processed)
    out: list[OutcomeExplanation] = []
    for s in samples:
        if s in done:
            continue
        has_vcf = (workdir / vcf_subdir / f"{s}{vcf_suffix}").exists()
        has_gvcf = (workdir / gvcf_subdir / f"{s}{gvcf_suffix}").exists()
        if not has_vcf and not has_gvcf:
            reason = REASON_VCF_AND_GVCF_ABSENT
        elif not has_gvcf:
            reason = REASON_GVCF_ABSENT
        elif not has_vcf:
            reason = REASON_VCF_ABSENT
        else:
            # Both inputs present yet no outcome: `orient_to_genome` raised, or the
            # driver filtered the sample. Not a missing input, and not silent either.
            reason = REASON_UNRECOGNISED_FLAG
        out.append(
            OutcomeExplanation(
                flag="unknown", reason=reason, callable_=False, sample=s
            )
        )
    return tuple(out)
