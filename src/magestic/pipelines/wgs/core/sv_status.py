"""magestic.pipelines.wgs.core.sv_status — refuse to emit a rate with no detector.

DEFECT 7 — THE DELIBERATE CHOICE, STATED ONCE
----------------------------------------------
The brief offered two options: implement inversion/duplication detection, or make
the pipeline refuse to emit an SV rate with an explicit "SV unmeasured" status.

🔴 **This module implements the REFUSAL.** The reasoning, so a reviewer does not
"complete" it later by adding a half-built detector:

* There is no inversion or tandem-duplication detector anywhere in TRACE
  (asserted by `test_repro_no_detector_exists`, which searches `dir(sv_plasmid)`).
* The translocation and large-deletion helpers that DO exist
  (`mate_pair_large_deletion`, `within_read_large_deletions`,
  `one_read_inside_large_deletion`) are called by nothing in the pipeline, so there
  is no path from a BAM to an SV number of any kind, and those paths have never been
  exercised on real data.
* An inversion caller built tonight would be unvalidated, and its output would be a
  *number*, which is exactly the artefact that cannot be withdrawn once quoted. A
  half-built detector is the worse failure: silence here would let someone quote a
  rate that does not exist, and a bad detector lets them quote one that is wrong
  while looking measured.

So: `sv_status()` is always `SV_UNMEASURED`, step 3 ships it as a column, and
`require_no_sv_rate` refuses any output mapping that contains an SV rate key.

DEFECT 13 — "0" IS NOT "NOT MEASURED"
--------------------------------------
`sv_plasmid.detect_plasmid_integration` returns `{"upstream_count": 0,
"downstream_count": 0}` for a BAM with no `integrated_plasmid_` contig — which is
**numerically identical** to a sample that was measured and found negative (asserted
by `test_repro_zero_is_indistinguishable_from_unmeasured`). A cohort combine over
those rows can print a 0% integration rate that is an absence of measurement.

`integration_measurement` makes the distinction explicit: `measured=False`,
`total=None`, `status="not_measured_no_integrated_contig"`. `cohort_integration_rate`
then REFUSES a cohort containing an unmeasured sample rather than treating its `None`
as a zero. Same class as GATE 2, same remedy.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

__all__ = [
    "INTEGRATION_MEASURED",
    "INTEGRATION_NOT_MEASURED",
    "IntegrationMeasurement",
    "IntegrationNotMeasured",
    "SV_UNMEASURED",
    "SvRateRefused",
    "cohort_integration_rate",
    "integration_measurement",
    "is_sv_rate_name",
    "require_no_sv_rate",
    "require_no_sv_rate_in_columns",
    "sv_status",
]

#: The only SV status this pipeline can honestly emit.
SV_UNMEASURED = "SV_UNMEASURED"

INTEGRATION_MEASURED = "measured"
INTEGRATION_NOT_MEASURED = "not_measured_no_integrated_contig"

#: Contig-name convention for an integrated guide-donor plasmid (mirrors
#: `step_3_sv_plasmid.INTEGRATED_PREFIX` and TRACE's private `_integrated_contigs`).
INTEGRATED_PREFIX = "integrated_plasmid_"

# --------------------------------------------------------------------------- #
# NN 255 F2 — the refusal READ as comprehensive and was not                     #
# --------------------------------------------------------------------------- #
#
# NN 252 measured **10 of 12** plausible SV-rate keys passing: `large_deletion_rate`,
# `translocation_rate`, `translocation_frequency`, `inversions_rate` (**plural defeats a
# `[^a-z]|$` anchor**), `n_inversions_per_cell`, `svRate` (**camelCase lowers to
# `svrate`, and `(sv)([^a-z]|$)` then fails on the following `r`**),
# `duplications_fraction`, `rearrangement_rate`, `cnv_rate`, `indel_gt50bp_rate`. It
# also never looked inside a nested dict, and was never applied to a TSV header at all.
#
# `large_deletion_rate` is the name a well-intentioned author reaches for, and it is
# exactly the class TRACE ships three UNCALLED helpers for (`mate_pair_large_deletion`,
# `within_read_large_deletions`, `one_read_inside_large_deletion`), which is what makes
# the hole realistic rather than adversarial.
#
# The fix is a NORMALIZED TOKEN match, not a longer regex: split camelCase, split on any
# non-alphanumeric, lowercase, and singularize. Matching tokens rather than substrings
# is also what keeps it from over-firing — `sv_plasmid_calls` has no rate token and
# `upstream_count` has neither.
#
# 🔴 **The direction of the asymmetry is deliberate.** A false refusal is a loud error an
# author fixes by renaming a key; a false pass is a published rate with no detector
# behind it, which cannot be withdrawn once quoted. So where a word is ambiguous the
# refusal wins.
_CAMEL = re.compile(r"(?<=[a-z0-9])(?=[A-Z])")
_NON_ALNUM = re.compile(r"[^a-z0-9]+")

#: Always an SV class: no detector exists for any of these anywhere in the pipeline.
_SV_TOKENS = frozenset({
    "sv", "svs", "structural", "inversion", "inverted", "duplication", "dup",
    "translocation", "rearrangement", "cnv", "cna", "aneuploidy", "aneuploid",
    "breakpoint", "chimera", "chimeric", "fusion", "transposition",
})

#: An SV class ONLY when a size/evidence qualifier is also present. A deletion or an
#: insertion inside the donor window IS measurable (`edit_profile.variant_class` records
#: DELETION / INSERTION), so refusing a bare `deletion_rate` would refuse a real
#: quantity; `large_deletion_rate` or `indel_gt50bp_rate` is an SV by size.
_SV_IF_QUALIFIED = frozenset({"deletion", "insertion", "indel", "del", "ins"})
_SIZE_TOKENS = frozenset({
    "large", "big", "huge", "long", "kb", "mb", "megabase", "kilobase",
    "matepair", "mate", "gt50bp", "gt50", "gt100", "gt1kb", "genomic", "span",
})

_RATE_TOKENS = frozenset({
    "rate", "rates", "freq", "freqs", "frequency", "frequencies", "fraction",
    "fractions", "pct", "percent", "percentage", "proportion", "prevalence",
    "burden", "density", "incidence",
})
#: Multi-token rate idioms, matched as adjacent token pairs.
_RATE_PAIRS = frozenset({("per", "cell"), ("per", "sample"), ("per", "clone"),
                         ("per", "genome"), ("per", "well")})


class SvRateRefused(ValueError):
    """An output tried to carry an SV rate for a quantity with no detector."""


class IntegrationNotMeasured(ValueError):
    """A cohort rate was asked for over at least one unmeasured sample."""


def sv_status() -> str:
    """Always `SV_UNMEASURED`. A function, not a bare constant, so that a future
    detector has one place to change and every caller picks it up."""
    return SV_UNMEASURED


def _tokens(key: object) -> list[str]:
    """`svRate` -> `['sv', 'rate']`; `inversions_rate` -> `['inversion', 'rate']`.

    camelCase is split BEFORE lowercasing (that is what `svRate` needed), then every
    run of non-alphanumerics, then a crude singularization. `-es` is stripped before
    `-s` so `frequencies`/`fractions` and `inversions` both land on a stem the token
    sets carry.
    """
    text = _CAMEL.sub("_", str(key))
    out: list[str] = []
    for raw in _NON_ALNUM.split(text.lower()):
        if not raw:
            continue
        if raw.endswith("ies") and len(raw) > 4:
            raw = raw[:-3] + "y"
        elif raw.endswith("es") and len(raw) > 3 and raw[-3] in "sxzh":
            raw = raw[:-2]
        elif raw.endswith("s") and not raw.endswith("ss") and len(raw) > 2:
            raw = raw[:-1]
        out.append(raw)
    return out


def is_sv_rate_name(key: object) -> bool:
    """True when a name would publish a rate for a quantity with no detector.

    Normalized token match on BOTH halves — an SV class and a rate-like quantity —
    because matching one half alone would refuse `sv_status` (legitimate: it is the
    refusal itself) or `read_count_mid_point`.
    """
    tokens = _tokens(key)
    if not tokens:
        return False
    tset = set(tokens)
    has_sv = bool(tset & _SV_TOKENS) or bool(
        (tset & _SV_IF_QUALIFIED) and (tset & _SIZE_TOKENS)
    )
    if not has_sv:
        return False
    pairs = set(zip(tokens, tokens[1:]))
    return bool(tset & _RATE_TOKENS) or bool(pairs & _RATE_PAIRS)


def _walk(value: object, prefix: str = "") -> Iterable[str]:
    """Every key path inside a nested mapping/sequence, as `a.b.c` strings.

    NN 252: `{"summary": {"sv_rate": 0.07}}` passed, because only top-level keys were
    inspected. A provenance dict is nested by construction (`write_provenance` takes
    `parameters=` and `counts=` sub-mappings), so this is the realistic shape.
    """
    if isinstance(value, Mapping):
        for k, v in value.items():
            path = f"{prefix}.{k}" if prefix else str(k)
            yield path
            yield from _walk(v, path)
    elif isinstance(value, (list, tuple, set, frozenset)):
        for i, v in enumerate(value):
            yield from _walk(v, f"{prefix}[{i}]" if prefix else f"[{i}]")


def require_no_sv_rate(outputs: Mapping[str, object]) -> None:
    """Refuse any mapping whose keys would publish an SV rate, at ANY nesting depth.

    Call sites, stated accurately because the old docstring claimed one that did not
    exist (NN 252 B4, defect-14 class — a guard documented as wired where it was not):

    * `pipeline/step_3_sv_plasmid.py` `_provenance_counts` (the original site), and
    * `pipeline/step_5_join_outcomes.py` `_counts`, plus
      `require_no_sv_rate_in_columns` over step 5's written header and matrix header.

    It inspects NAMES, not values: a zero is as publishable as a non-zero and just as
    wrong. The last token path is reported so a nested offender can be found.
    """
    offenders = sorted(p for p in _walk(outputs) if is_sv_rate_name(p.split(".")[-1]))
    if offenders:
        raise SvRateRefused(
            f"refusing to emit {offenders}: there is no inversion, duplication, "
            "translocation, rearrangement or large-deletion detector in this pipeline, "
            f"so no SV rate exists to report. Emit {SV_UNMEASURED} instead (see "
            "core/sv_status.py for why the refusal was chosen over a half-built "
            "detector). A name refused here is not a name to work around: either the "
            "quantity has a validated detector or it does not."
        )


def require_no_sv_rate_in_columns(columns: Iterable[str], *, source: str = "") -> None:
    """Refuse a TSV HEADER that would publish an SV rate as a column.

    NN 252: `require_no_sv_rate` was applied to provenance mappings only, never to
    `step_3.SHARD_COLUMNS` nor to any written header, so an SV-rate COLUMN reached a
    file unrefused. A column name is published exactly as readily as a count.
    """
    offenders = sorted(c for c in columns if is_sv_rate_name(c))
    if offenders:
        raise SvRateRefused(
            f"refusing to write the column(s) {offenders}"
            + (f" in {source}" if source else "")
            + f": no detector exists for them, so the column would carry an "
            f"unmeasured quantity. Emit {SV_UNMEASURED} instead."
        )


@dataclass(frozen=True)
class IntegrationMeasurement:
    """Plasmid-integration evidence for one sample, with measurement made explicit.

    `total is None` **iff** `measured is False`. That invariant is what lets a
    cohort rate refuse rather than silently average a zero that was never measured.
    """

    sample: str
    measured: bool
    status: str
    n_integrated_contigs: int = 0
    upstream_count: int | None = None
    downstream_count: int | None = None

    @property
    def total(self) -> int | None:
        if not self.measured:
            return None
        return int(self.upstream_count or 0) + int(self.downstream_count or 0)


def integration_measurement(
    bam_path: str,
    *,
    sample: str = "",
    integrated_contigs: Iterable[str] | None = None,
    junction_distance: int | Mapping[str, int] = 1129,
    window: int = 200,
) -> IntegrationMeasurement:
    """Measure plasmid integration, or say plainly that it was not measurable.

    A BAM with no `integrated_plasmid_` contig cannot be measured for integration —
    the junctions this assay counts do not exist in its coordinate space. Returning
    `0` for that case, as TRACE does, is the defect.
    """
    import pysam

    from . import _trace

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = list(bam.references)
    contigs = (
        list(integrated_contigs)
        if integrated_contigs is not None
        else [r for r in refs if r.startswith(INTEGRATED_PREFIX)]
    )
    present = [c for c in contigs if c in refs]
    if not present:
        return IntegrationMeasurement(
            sample=sample, measured=False, status=INTEGRATION_NOT_MEASURED,
            n_integrated_contigs=0, upstream_count=None, downstream_count=None,
        )
    counts = _trace.sv_plasmid().detect_plasmid_integration(
        bam_path, present, junction_distance, window
    )
    return IntegrationMeasurement(
        sample=sample, measured=True, status=INTEGRATION_MEASURED,
        n_integrated_contigs=len(present),
        upstream_count=int(counts["upstream_count"]),
        downstream_count=int(counts["downstream_count"]),
    )


def cohort_integration_rate(samples: Iterable[tuple[str, int | None]]) -> float:
    """Fraction of samples with integration evidence, refusing unmeasured ones.

    Parameters
    ----------
    samples : iterable of (sample, total)
        `total` is `IntegrationMeasurement.total`: `None` means not measured.

    Raises
    ------
    IntegrationNotMeasured
        Any `total is None`. Averaging over it would publish an absence of
        measurement as a 0%, which is defect 13.
    """
    rows = list(samples)
    unmeasured = sorted(s for s, t in rows if t is None)
    if unmeasured:
        raise IntegrationNotMeasured(
            f"{len(unmeasured)} of {len(rows)} sample(s) were never measured for "
            f"plasmid integration (first: {', '.join(unmeasured[:5])}). A cohort rate "
            "over them would report an absence of measurement as 0%. Report the "
            f"measured subset with its own denominator, or {INTEGRATION_NOT_MEASURED}."
        )
    if not rows:
        raise IntegrationNotMeasured("no samples: a rate over an empty cohort is not 0%")
    return sum(1 for _, t in rows if (t or 0) > 0) / len(rows)
