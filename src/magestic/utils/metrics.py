"""
magestic.utils.metrics — representation and evenness statistics that carry their own depth.

Gini, Shannon entropy, Pielou evenness, Simpson and inverse Simpson, Hill numbers, top-N share
and richness, for any count vector a MAGESTIC pipeline produces: oligos per library, barcodes per
oligo, reads per guide, clones per position.

🔴 THE TRAP THIS MODULE EXISTS TO ENCODE
========================================
**Every statistic here is depth-dependent. Computed on a sample, it describes the sample, not the
library.** Sequence the same plasmid pool twice at different depths and every number below moves.
Sequence two libraries at different depths and the difference you measure is partly the depth.

Two consequences are built into the API rather than left to the reader:

1. **No function returns a bare float.** They return :class:`Metric`, which carries the depth the
   statistic was computed at, the count floor applied, how many observations survived that floor,
   the denominator universe, and the grain. ``Metric`` deliberately has no ``__float__``, so
   ``float(m)`` raises and ``f"{m}"`` prints the depth alongside the value. Quoting the number
   without its depth takes a deliberate ``.value``.

2. **Rarefaction is offered, and the un-rarefied form is labelled sample-specific.**
   :func:`rarefied_metric` subsamples reads *without* replacement to a common depth so two
   libraries can be compared honestly, and records the depth, the replicate count and the seed.
   Anything not rarefied says so in its label.

**The count floor is not cosmetic.** In the 2026-11 AsLOV2/POL3 guide-donor-bc0 data, 65–84 % of
distinct barcodes are singletons, so the distinct-key denominator is dominated by the synthesis and
sequencing error tail. An unfloored Gini over distinct keys is close to meaningless: it mostly
measures how many error variants the run happened to see. Hence ``min_count`` is an explicit
argument everywhere, :func:`representation_metrics` reports at several floors by default rather
than one, and ``min_count=0`` (keep the designed zeros) is a distinct, labelled choice.

**The same floor choice already differs inside this package**, which is why the universe is
recorded too: ``guide_donor_bc0.core.library_representation.gini_index`` documents that zero-count
entries are *included* so missing oligos drag the index toward 1, while its own caller
``plot_gini_indices`` passes ``x[x > 0]`` and so silently applies a floor of 1. Those are two
different statistics of the same library. ``Metric.context.universe`` and ``min_count`` name which
one you have.

WORKED EXAMPLE — the distinct-key vs read-weighted denominator
---------------------------------------------------------------
From the 2026-11 rerun of the POL3 guide-donor-bc0 libraries (unit 2168 §6.0.2 and §5.5a;
job 53132655), scoring one arm frame against Alex's designed donor table — **the same comparison,
two denominators**:

====================================  =====================  ======================
frame / source                        distinct-key hit rate  read/link-weighted rate
====================================  =====================  ======================
``(donor + SPS[:5])[-40:]`` (§6.0.2)  **0.90 %**             **89.31 %**
``revcomp(donor[-35:] + SPS[:5])``    **0.59 %**             **90.91 %**
(§5.5a, link-weighted)
====================================  =====================  ======================

A factor of ~100–150 between two correct numbers for one comparison. The distinct-key rate is low
only because most distinct keys are low-abundance error variants; it is almost never the number
anyone wants, and reading one as the other was the single most repeated mistake in that unit.
:func:`coverage_fractions` returns both from one call so neither can be quoted alone.

Author: Kevin R. Roy
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np

__all__ = [
    "WEIGHTING_READ",
    "WEIGHTING_DISTINCT",
    "UNIVERSE_OBSERVED",
    "UNIVERSE_DESIGNED",
    "DEFAULT_MIN_COUNTS",
    "DEFAULT_HILL_ORDERS",
    "DEFAULT_TOP_N",
    "DEFAULT_SEED",
    "MetricContext",
    "Metric",
    "MetricReport",
    "gini",
    "shannon_entropy",
    "pielou_evenness",
    "simpson_index",
    "inverse_simpson",
    "hill_number",
    "effective_number_of_species",
    "top_n_share",
    "richness",
    "coverage_fractions",
    "rarefy_counts",
    "rarefied_metric",
    "representation_metrics",
    "STATISTICS",
]

# Vocabulary, so a sidecar column is comparable across pipelines.
WEIGHTING_READ = "read_weighted"
WEIGHTING_DISTINCT = "distinct_key"
UNIVERSE_OBSERVED = "observed"      # only keys the run actually saw
UNIVERSE_DESIGNED = "designed"      # designed members absent from the run are present as zeros

DEFAULT_MIN_COUNTS: Tuple[int, ...] = (1, 2, 5, 10)
DEFAULT_HILL_ORDERS: Tuple[float, ...] = (0.0, 1.0, 2.0, math.inf)
DEFAULT_TOP_N: Tuple[int, ...] = (1, 10, 100)
DEFAULT_SEED = 20262174
DEFAULT_REPLICATES = 20


# ============================================================================
# Arithmetic kernels — private on purpose: they return bare floats, which is
# exactly what the public API refuses to do.
# ============================================================================

def _as_counts(counts: Any) -> np.ndarray:
    """Coerce to a 1-D float64 array of non-negative finite counts."""
    arr = np.asarray(counts, dtype=float).ravel()
    if arr.size and not np.all(np.isfinite(arr)):
        raise ValueError("counts contain NaN or inf; drop or impute them before measuring")
    if arr.size and np.any(arr < 0):
        raise ValueError("counts contain negative values; these statistics are undefined there")
    return arr


def _proportions(x: np.ndarray) -> np.ndarray:
    total = x.sum()
    if total <= 0:
        return np.empty(0, dtype=float)
    p = x[x > 0] / total
    return p


def _gini(x: np.ndarray) -> float:
    """Gini coefficient (0 = every member equal, 1 = all mass on one member).

    Zeros present in *x* are included, so a vector padded with designed-but-absent members
    scores closer to 1. Matches ``library_representation.gini_index`` exactly.
    """
    n = x.size
    total = x.sum()
    if n == 0 or total == 0:
        return float("nan")
    s = np.sort(x)
    idx = np.arange(1, n + 1)
    return float((2.0 * np.sum(idx * s) / total - (n + 1)) / n)


def _shannon(x: np.ndarray) -> float:
    """Shannon entropy H in nats. 0 for a single occupied member."""
    p = _proportions(x)
    if p.size == 0:
        return float("nan")
    return float(-np.sum(p * np.log(p)))


def _pielou(x: np.ndarray) -> float:
    """Pielou evenness J = H / ln(S), normalized to [0, 1].

    Undefined (nan) for S < 2: one occupied member has no evenness to speak of.
    """
    p = _proportions(x)
    s = p.size
    if s == 0:
        return float("nan")
    if s == 1:
        return float("nan")
    return float(_shannon(x) / math.log(s))


def _simpson(x: np.ndarray) -> float:
    """Simpson index lambda = sum(p^2), the plug-in (biased) estimator.

    The plug-in form is used deliberately, not the unbiased ``sum n(n-1)/N(N-1)``: it is the one
    for which ``1/lambda`` equals the Hill number of order 2 exactly, which is a tested invariant
    of this module. lambda is the probability two reads drawn with replacement are the same member.
    """
    p = _proportions(x)
    if p.size == 0:
        return float("nan")
    return float(np.sum(p ** 2))


def _inv_simpson(x: np.ndarray) -> float:
    lam = _simpson(x)
    if not np.isfinite(lam) or lam <= 0:
        return float("nan")
    return float(1.0 / lam)


def _hill(x: np.ndarray, q: float) -> float:
    """Hill number of order q — the effective number of equally abundant members.

    q=0 richness, q=1 exp(Shannon) (the limit, special-cased), q=2 inverse Simpson,
    q=inf 1/max(p) (the Berger-Parker reciprocal). Units are members, which is what makes it the
    interpretable evenness measure for a library: "this pool behaves like N equally abundant
    designs".
    """
    p = _proportions(x)
    if p.size == 0:
        return float("nan")
    if q == 1.0:
        return float(math.exp(_shannon(x)))
    if math.isinf(q):
        if q > 0:
            return float(1.0 / p.max())
        return float(1.0 / p.min())
    return float(np.sum(p ** q) ** (1.0 / (1.0 - q)))


def _top_n_share(x: np.ndarray, n: int) -> float:
    """Fraction of total counts held by the n most abundant members.

    n greater than the number of occupied members returns 1.0 (they hold everything).
    """
    if n <= 0:
        raise ValueError("n must be >= 1")
    total = x.sum()
    if x.size == 0 or total == 0:
        return float("nan")
    s = np.sort(x)[::-1]
    return float(s[:n].sum() / total)


def _richness(x: np.ndarray) -> float:
    """Number of members with a non-zero count. The floor is applied before this is reached."""
    if x.size == 0:
        return float("nan")
    return float(np.count_nonzero(x > 0))


# ============================================================================
# The context every value travels with
# ============================================================================

@dataclass(frozen=True)
class MetricContext:
    """Everything needed to say what a statistic's value is a statistic *of*.

    Attributes
    ----------
    depth : int
        Total counts the statistic was computed at, after the floor. **The number that makes the
        value comparable or not.**
    n_observations : int
        Members retained after the floor.
    min_count : int
        Count floor applied. ``0`` keeps zero-count members (the designed universe);
        ``1`` is the conventional "observed at all".
    n_excluded_by_floor : int
        Members dropped by the floor.
    counts_excluded_by_floor : int
        Counts dropped with them, so the discarded mass is visible.
    weighting : str
        ``read_weighted`` or ``distinct_key`` — which denominator the value is on.
    universe : str
        ``designed`` if absent designed members were present as zeros before the floor,
        ``observed`` otherwise.
    grain : str
        What one member is: ``oligo``, ``guide_bc0``, ``bc1``, ``position`` …
    n_designed : int, optional
        Size of the designed universe when known, for a coverage denominator.
    rarefied_to : int, optional
        Common depth subsampled to; ``None`` means **not rarefied, hence sample-specific**.
    rarefaction_replicates, rarefaction_seed : int, optional
        Set when rarefied, so the value is reproducible.
    note : str, optional
        Anything that qualifies the value (for example a library shallower than the target depth).
    """

    depth: int
    n_observations: int
    min_count: int
    n_excluded_by_floor: int = 0
    counts_excluded_by_floor: int = 0
    weighting: str = WEIGHTING_READ
    universe: str = UNIVERSE_OBSERVED
    grain: str = "unspecified"
    n_designed: Optional[int] = None
    rarefied_to: Optional[int] = None
    rarefaction_replicates: Optional[int] = None
    rarefaction_seed: Optional[int] = None
    note: Optional[str] = None

    @property
    def is_rarefied(self) -> bool:
        return self.rarefied_to is not None

    @property
    def zeros_included(self) -> bool:
        """True when the floor kept zero-count members, so absences drag the statistic."""
        return self.min_count <= 0

    def label(self) -> str:
        """One-line provenance for the value. Always mentions the depth."""
        bits = [
            f"depth={self.depth:,}",
            f"floor={self.min_count}",
            f"n={self.n_observations:,} {self.grain}",
            f"universe={self.universe}",
        ]
        if self.is_rarefied:
            bits.append(
                f"RAREFIED to {self.rarefied_to:,} x{self.rarefaction_replicates} "
                f"seed={self.rarefaction_seed}"
            )
        else:
            bits.append("NOT rarefied: sample-specific")
        if self.note:
            bits.append(self.note)
        return "; ".join(bits)

    def to_dict(self) -> Dict[str, Any]:
        d = {
            "depth": self.depth,
            "n_observations": self.n_observations,
            "min_count": self.min_count,
            "n_excluded_by_floor": self.n_excluded_by_floor,
            "counts_excluded_by_floor": self.counts_excluded_by_floor,
            "weighting": self.weighting,
            "universe": self.universe,
            "grain": self.grain,
            "n_designed": self.n_designed,
            "zeros_included": self.zeros_included,
            "rarefied": self.is_rarefied,
            "rarefied_to": self.rarefied_to,
            "rarefaction_replicates": self.rarefaction_replicates,
            "rarefaction_seed": self.rarefaction_seed,
            "note": self.note,
        }
        return d


@dataclass(frozen=True)
class Metric:
    """A statistic bound to the depth it was computed at.

    🔴 There is **no** ``__float__`` and no ``__index__`` on purpose: ``float(metric)`` raises, and
    ``str(metric)`` / ``f"{metric}"`` always print the depth and floor. Reaching the bare number
    takes an explicit ``.value``, which is the point — a value quoted without its depth is the
    error this module exists to prevent.
    """

    name: str
    value: float
    context: MetricContext
    sd: Optional[float] = None          # across rarefaction replicates
    n_replicates: Optional[int] = None

    def __str__(self) -> str:
        v = "nan" if self.value is None or not np.isfinite(self.value) else f"{self.value:.6g}"
        sd = "" if self.sd is None or not np.isfinite(self.sd) else f" ± {self.sd:.3g}"
        return f"{self.name}={v}{sd} [{self.context.label()}]"

    __repr__ = __str__

    def to_dict(self) -> Dict[str, Any]:
        row: Dict[str, Any] = {"metric": self.name, "value": self.value, "sd": self.sd,
                               "n_replicates": self.n_replicates}
        row.update(self.context.to_dict())
        return row


@dataclass
class MetricReport:
    """A batch of :class:`Metric` values for one count vector, plus the labels identifying it."""

    metrics: List[Metric] = field(default_factory=list)
    labels: Dict[str, Any] = field(default_factory=dict)

    def __iter__(self):
        return iter(self.metrics)

    def __len__(self) -> int:
        return len(self.metrics)

    def extend(self, other: Union["MetricReport", Iterable[Metric]]) -> "MetricReport":
        self.metrics.extend(other.metrics if isinstance(other, MetricReport) else list(other))
        return self

    def rows(self) -> List[Dict[str, Any]]:
        out = []
        for m in self.metrics:
            row = dict(self.labels)
            row.update(m.to_dict())
            out.append(row)
        return out

    def to_dataframe(self):
        """Long-format table, one row per (label set, metric, floor). Needs pandas."""
        import pandas as pd
        return pd.DataFrame(self.rows())

    def to_dict(self) -> Dict[str, Any]:
        return {"labels": dict(self.labels), "metrics": [m.to_dict() for m in self.metrics]}

    def get(self, name: str, *, min_count: Optional[int] = None) -> Optional[Metric]:
        for m in self.metrics:
            if m.name != name:
                continue
            if min_count is not None and m.context.min_count != min_count:
                continue
            return m
        return None


# ============================================================================
# Floor + context assembly
# ============================================================================

def _prepare(
    counts: Any,
    *,
    min_count: int,
    universe: str,
    grain: str,
    n_designed: Optional[int],
    weighting: str,
    note: Optional[str] = None,
) -> Tuple[np.ndarray, MetricContext]:
    """Apply the floor and build the context. One place, so every statistic agrees."""
    if universe not in (UNIVERSE_OBSERVED, UNIVERSE_DESIGNED):
        raise ValueError(f"universe must be {UNIVERSE_OBSERVED!r} or {UNIVERSE_DESIGNED!r}")
    raw = _as_counts(counts)
    if min_count <= 0:
        kept = raw
    else:
        kept = raw[raw >= min_count]
    notes = [note] if note else []
    if min_count > 0 and universe == UNIVERSE_DESIGNED:
        n_zeros = int(np.count_nonzero(raw == 0))
        if n_zeros:
            notes.append(
                f"floor {min_count} dropped {n_zeros:,} designed zero-count members, so the "
                f"effective universe is 'observed at floor {min_count}'"
            )
    ctx = MetricContext(
        depth=int(kept.sum()),
        n_observations=int(kept.size),
        min_count=int(min_count),
        n_excluded_by_floor=int(raw.size - kept.size),
        counts_excluded_by_floor=int(raw.sum() - kept.sum()),
        weighting=weighting,
        universe=universe,
        grain=grain,
        n_designed=n_designed,
        note="; ".join(notes) if notes else None,
    )
    return kept, ctx


def _single(
    kernel,
    name: str,
    counts: Any,
    *,
    min_count: int,
    universe: str,
    grain: str,
    n_designed: Optional[int],
    weighting: str,
    **kernel_kw: Any,
) -> Metric:
    kept, ctx = _prepare(
        counts, min_count=min_count, universe=universe, grain=grain,
        n_designed=n_designed, weighting=weighting,
    )
    return Metric(name=name, value=kernel(kept, **kernel_kw), context=ctx)


# ============================================================================
# Public statistics — each returns a Metric, never a bare number
# ============================================================================

def gini(counts, *, min_count: int = 1, universe: str = UNIVERSE_OBSERVED,
         grain: str = "unspecified", n_designed: Optional[int] = None) -> Metric:
    """Gini coefficient of a count vector, with its depth.

    0 = perfectly even. ``min_count=0`` keeps zero-count members, so designed-but-absent members
    push the value up; ``min_count=1`` measures only what was seen. **Both are legitimate and they
    are different statistics** — the returned context says which.

    🔴 The upper bound is ``(n-1)/n``, **not 1**: all the mass on one of n members scores exactly
    (n-1)/n. So 0.90 over ten members is maximal inequality, while 0.90 over ten thousand is not,
    and two libraries with different member counts are scored against different ceilings. Read the
    value against ``context.n_observations``, which is why it travels with it.
    """
    return _single(_gini, "gini", counts, min_count=min_count, universe=universe,
                   grain=grain, n_designed=n_designed, weighting=WEIGHTING_READ)


def shannon_entropy(counts, *, min_count: int = 1, universe: str = UNIVERSE_OBSERVED,
                    grain: str = "unspecified", n_designed: Optional[int] = None) -> Metric:
    """Shannon entropy in nats, with its depth. Rises with both richness and evenness."""
    return _single(_shannon, "shannon_entropy", counts, min_count=min_count, universe=universe,
                   grain=grain, n_designed=n_designed, weighting=WEIGHTING_READ)


def pielou_evenness(counts, *, min_count: int = 1, universe: str = UNIVERSE_OBSERVED,
                    grain: str = "unspecified", n_designed: Optional[int] = None) -> Metric:
    """Pielou evenness J = H / ln(S) in [0, 1], with its depth.

    🔴 J divides out richness, so it is **more** depth-sensitive than it looks: deeper sequencing
    finds more rare members, which raises ln(S) faster than H. Compare J across libraries only
    after :func:`rarefied_metric`.
    """
    return _single(_pielou, "pielou_evenness", counts, min_count=min_count, universe=universe,
                   grain=grain, n_designed=n_designed, weighting=WEIGHTING_READ)


def simpson_index(counts, *, min_count: int = 1, universe: str = UNIVERSE_OBSERVED,
                  grain: str = "unspecified", n_designed: Optional[int] = None) -> Metric:
    """Simpson index lambda = sum(p^2): the chance two reads hit the same member."""
    return _single(_simpson, "simpson_index", counts, min_count=min_count, universe=universe,
                   grain=grain, n_designed=n_designed, weighting=WEIGHTING_READ)


def inverse_simpson(counts, *, min_count: int = 1, universe: str = UNIVERSE_OBSERVED,
                    grain: str = "unspecified", n_designed: Optional[int] = None) -> Metric:
    """1/lambda, in members. Identical to :func:`hill_number` at q=2."""
    return _single(_inv_simpson, "inverse_simpson", counts, min_count=min_count, universe=universe,
                   grain=grain, n_designed=n_designed, weighting=WEIGHTING_READ)


def hill_number(counts, q: float = 1.0, *, min_count: int = 1,
                universe: str = UNIVERSE_OBSERVED, grain: str = "unspecified",
                n_designed: Optional[int] = None) -> Metric:
    """Hill number of order *q* — the effective number of equally abundant members.

    The most interpretable evenness measure for a library, because its unit is members:
    "this pool behaves like 1,840 equally abundant designs out of 2,744 designed".
    q=0 is richness, q=1 is exp(Shannon), q=2 is inverse Simpson, q=inf is 1/max(p).
    Rising q weights abundant members more, so the sequence over q *is* the skew profile.
    """
    return _single(_hill, f"hill_q{q:g}", counts, min_count=min_count, universe=universe,
                   grain=grain, n_designed=n_designed, weighting=WEIGHTING_READ, q=q)


def effective_number_of_species(counts, q: float = 1.0, **kw) -> Metric:
    """Alias of :func:`hill_number` under its ecology name."""
    return hill_number(counts, q, **kw)


def top_n_share(counts, n: int, *, min_count: int = 1, universe: str = UNIVERSE_OBSERVED,
                grain: str = "unspecified", n_designed: Optional[int] = None) -> Metric:
    """Fraction of counts held by the *n* most abundant members, with its depth.

    🔴 Strongly depth- and floor-dependent: the floor removes members from the denominator, which
    *raises* the share. Compare only at the same floor and depth.
    """
    return _single(_top_n_share, f"top_{n}_share", counts, min_count=min_count, universe=universe,
                   grain=grain, n_designed=n_designed, weighting=WEIGHTING_READ, n=n)


def richness(counts, *, min_count: int = 1, universe: str = UNIVERSE_OBSERVED,
             grain: str = "unspecified", n_designed: Optional[int] = None) -> Metric:
    """Number of members at or above *min_count*, with the depth that produced it.

    🔴 Richness is the most depth-dependent statistic in this module — it only ever rises with
    depth — so "richness at a stated count floor" is the only form worth quoting, and the floor is
    required reading. ``min_count=1`` at one depth is not comparable to ``min_count=1`` at another.
    """
    return _single(_richness, f"richness_min{min_count}", counts, min_count=min_count,
                   universe=universe, grain=grain, n_designed=n_designed,
                   weighting=WEIGHTING_DISTINCT)


def coverage_fractions(
    counts,
    selected,
    *,
    grain: str = "unspecified",
    label: str = "selected",
    n_designed: Optional[int] = None,
) -> Dict[str, Metric]:
    """Return the distinct-key AND read-weighted fraction of a selection, from one call.

    This is the worked example in the module docstring expressed as code. *selected* is a boolean
    mask over the same members as *counts* — "this key matched a design", "this oligo is present",
    "this barcode is pure". The two fractions answer different questions and in this data they
    differ by a factor of ~100:

    - ``distinct_key_fraction`` — what share of distinct members are selected. Dominated by the
      singleton error tail, so usually not the number anyone wants.
    - ``read_weighted_fraction`` — what share of the reads are selected. Usually the right one.

    Returns
    -------
    dict[str, Metric]
        Keys ``distinct_key_fraction`` and ``read_weighted_fraction``. Both are returned so
        neither can be reported without the other being in hand.
    """
    x = _as_counts(counts)
    mask = np.asarray(selected).ravel().astype(bool)
    if mask.size != x.size:
        raise ValueError(f"selected has {mask.size} entries but counts has {x.size}")
    total = x.sum()
    n = x.size
    ctx_common = dict(min_count=0, universe=UNIVERSE_OBSERVED, grain=grain,
                      n_designed=n_designed, n_excluded_by_floor=0,
                      counts_excluded_by_floor=0)
    distinct_ctx = MetricContext(
        depth=int(total), n_observations=int(n), weighting=WEIGHTING_DISTINCT,
        note=f"share of distinct {grain} keys that are {label}", **ctx_common)
    read_ctx = MetricContext(
        depth=int(total), n_observations=int(n), weighting=WEIGHTING_READ,
        note=f"share of counts on {label} {grain} keys", **ctx_common)
    distinct_val = float(mask.sum() / n) if n else float("nan")
    read_val = float(x[mask].sum() / total) if total > 0 else float("nan")
    return {
        "distinct_key_fraction": Metric("distinct_key_fraction", distinct_val, distinct_ctx),
        "read_weighted_fraction": Metric("read_weighted_fraction", read_val, read_ctx),
    }


# ============================================================================
# Rarefaction — the honest way to compare two libraries
# ============================================================================

def rarefy_counts(counts, depth: int, *, rng: Optional[np.random.Generator] = None,
                  seed: Optional[int] = None) -> np.ndarray:
    """Subsample *depth* reads **without replacement** from a count vector.

    Without replacement is the right model: a sequencing run draws distinct molecules from a pool,
    so the subsample is multivariate hypergeometric, not multinomial. Returns a new count vector of
    the same length summing to *depth*.

    Raises
    ------
    ValueError
        If *depth* exceeds the available total — you cannot rarefy upward. Callers that want a
        flagged nan instead should use :func:`rarefied_metric`.
    """
    x = _as_counts(counts)
    total = int(x.sum())
    if depth < 0:
        raise ValueError("depth must be >= 0")
    if depth > total:
        raise ValueError(f"cannot rarefy to {depth:,} from a total of {total:,}")
    if rng is None:
        rng = np.random.default_rng(seed if seed is not None else DEFAULT_SEED)
    colors = np.rint(x).astype(np.int64)
    # multivariate_hypergeometric needs a positive-length color vector; preserve the zeros.
    nz = colors > 0
    out = np.zeros_like(colors)
    if nz.any() and depth > 0:
        out[nz] = rng.multivariate_hypergeometric(colors[nz], depth, method="marginals")
    return out.astype(float)


def rarefied_metric(
    statistic,
    counts,
    target_depth: int,
    *,
    replicates: int = DEFAULT_REPLICATES,
    seed: int = DEFAULT_SEED,
    min_count: int = 1,
    universe: str = UNIVERSE_OBSERVED,
    grain: str = "unspecified",
    n_designed: Optional[int] = None,
    **statistic_kw: Any,
) -> Metric:
    """Compute *statistic* at a common subsampled depth, averaged over replicates.

    *statistic* is any public function in this module (``gini``, ``pielou_evenness``,
    ``hill_number`` …). The floor is applied **after** subsampling, which is the point: a member
    seen 3 times in the full run may be seen 0 times in the subsample, and that is precisely the
    depth effect being controlled for.

    A library whose total is below *target_depth* returns ``nan`` with the shortfall in
    ``context.note``. It never silently returns the un-rarefied value.

    Returns
    -------
    Metric
        ``value`` is the mean over replicates, ``sd`` their standard deviation (the sampling
        spread, which is itself worth reporting), and the context records the target depth,
        replicate count and seed.
    """
    x = _as_counts(counts)
    total = int(x.sum())
    probe = statistic(x, min_count=min_count, universe=universe, grain=grain,
                      n_designed=n_designed, **statistic_kw)
    name = probe.name

    if total < target_depth:
        ctx = MetricContext(
            depth=total, n_observations=int(x.size), min_count=int(min_count),
            weighting=probe.context.weighting, universe=universe, grain=grain,
            n_designed=n_designed, rarefied_to=int(target_depth),
            rarefaction_replicates=int(replicates), rarefaction_seed=int(seed),
            note=(f"NOT COMPUTED: total depth {total:,} is below the rarefaction target "
                  f"{target_depth:,}; no value is comparable here"),
        )
        return Metric(name=name, value=float("nan"), context=ctx, sd=None,
                      n_replicates=0)

    rng = np.random.default_rng(seed)
    values: List[float] = []
    depths: List[int] = []
    n_obs: List[int] = []
    for _ in range(max(1, int(replicates))):
        sub = rarefy_counts(x, target_depth, rng=rng)
        m = statistic(sub, min_count=min_count, universe=universe, grain=grain,
                      n_designed=n_designed, **statistic_kw)
        values.append(m.value)
        depths.append(m.context.depth)
        n_obs.append(m.context.n_observations)

    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]
    ctx = MetricContext(
        depth=int(round(float(np.mean(depths)))),
        n_observations=int(round(float(np.mean(n_obs)))),
        min_count=int(min_count),
        weighting=probe.context.weighting,
        universe=universe,
        grain=grain,
        n_designed=n_designed,
        rarefied_to=int(target_depth),
        rarefaction_replicates=int(max(1, int(replicates))),
        rarefaction_seed=int(seed),
        note=(f"mean over {max(1, int(replicates))} subsamples of {target_depth:,} from "
              f"{total:,}; floor applied after subsampling"),
    )
    return Metric(
        name=name,
        value=float(np.mean(finite)) if finite.size else float("nan"),
        context=ctx,
        sd=float(np.std(finite, ddof=1)) if finite.size > 1 else None,
        n_replicates=int(finite.size),
    )


# Registry so a caller (or a config file) can name the statistics it wants.
STATISTICS = {
    "gini": gini,
    "shannon_entropy": shannon_entropy,
    "pielou_evenness": pielou_evenness,
    "simpson_index": simpson_index,
    "inverse_simpson": inverse_simpson,
    "hill_number": hill_number,
    "top_n_share": top_n_share,
    "richness": richness,
}


# ============================================================================
# The batch entry point
# ============================================================================

def representation_metrics(
    counts,
    *,
    min_counts: Sequence[int] = DEFAULT_MIN_COUNTS,
    hill_orders: Sequence[float] = DEFAULT_HILL_ORDERS,
    top_n: Sequence[int] = DEFAULT_TOP_N,
    universe: str = UNIVERSE_OBSERVED,
    grain: str = "unspecified",
    n_designed: Optional[int] = None,
    rarefy_to: Optional[int] = None,
    replicates: int = DEFAULT_REPLICATES,
    seed: int = DEFAULT_SEED,
    labels: Optional[Dict[str, Any]] = None,
) -> MetricReport:
    """Every statistic in this module, at several count floors, for one count vector.

    Reporting at several floors is the default rather than a option because a single floor hides
    the one thing a reader needs: whether the value is a property of the library or of the error
    tail. If gini moves from 0.93 at floor 1 to 0.41 at floor 10, the first number was mostly
    singletons.

    Parameters
    ----------
    counts : array-like
        Counts per member at one grain.
    min_counts : sequence of int
        Floors to report at. ``0`` keeps zero-count members (the designed universe).
    hill_orders : sequence of float
        Orders q for the Hill profile. ``inf`` is accepted.
    top_n : sequence of int
        N values for top-N share.
    universe, grain, n_designed
        Recorded in every context; see :class:`MetricContext`.
    rarefy_to : int, optional
        If given, every statistic is **also** computed at this common depth, labelled rarefied.
    replicates, seed : int
        Rarefaction settings, recorded in the context.
    labels : dict, optional
        Identifying columns (library, sample, pool) copied onto every output row.

    Returns
    -------
    MetricReport
    """
    report = MetricReport(labels=dict(labels or {}))
    x = _as_counts(counts)

    for floor in min_counts:
        common = dict(min_count=int(floor), universe=universe, grain=grain,
                      n_designed=n_designed)
        report.metrics.append(gini(x, **common))
        report.metrics.append(shannon_entropy(x, **common))
        report.metrics.append(pielou_evenness(x, **common))
        report.metrics.append(simpson_index(x, **common))
        report.metrics.append(inverse_simpson(x, **common))
        report.metrics.append(richness(x, **common))
        for q in hill_orders:
            report.metrics.append(hill_number(x, float(q), **common))
        for n in top_n:
            report.metrics.append(top_n_share(x, int(n), **common))

        if rarefy_to is not None:
            for fn, kw in (
                (gini, {}), (shannon_entropy, {}), (pielou_evenness, {}),
                (simpson_index, {}), (inverse_simpson, {}), (richness, {}),
            ):
                report.metrics.append(rarefied_metric(
                    fn, x, int(rarefy_to), replicates=replicates, seed=seed, **common, **kw))
            for q in hill_orders:
                report.metrics.append(rarefied_metric(
                    hill_number, x, int(rarefy_to), replicates=replicates, seed=seed,
                    q=float(q), **common))
            for n in top_n:
                report.metrics.append(rarefied_metric(
                    top_n_share, x, int(rarefy_to), replicates=replicates, seed=seed,
                    n=int(n), **common))

    return report
