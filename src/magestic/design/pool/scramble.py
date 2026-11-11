"""GC **distribution**-matched non-targeting ("scrambled") controls.

Lifted from the screen packer ``467aa_pack_pool.py``, which is where this logic was written
and where both of its correctness bugs were paid for. Two things here are not obvious and
were each a failed build:

1. 🔴 **The match is to the GC *distribution*, not the GC *mean*.** The real pool's
   ``gc_frac`` is binned (2 % by default), scrambles are drawn per bin in proportion, and a
   two-sample KS statistic is reported every build. A mean-matched set passes a test of means
   and fails this one, which is the point of Kevin's wording.

2. 🔴 **An exact N needs a top-up *and* a trim.** The per-bin targets round independently, so
   they sum to one or two **less** than N (every early build landed short and was chased by
   hand) *and* sometimes one **more** — the v5 pool's GC histogram differed from v4's by six
   designs, rounded to 1001, and the job died with the whole pack discarded (NN 496/497).
   Both corrections draw from the modal GC bin, where the real designs are densest, so the
   distribution match is perturbed least.

Making exactly N of a *generated* arm is not "padding to fit": the count is the arm's
declared target, and the generator raises if it cannot reach it. G-COMPOSITION's refusal to
pad is about the **pool total**, which is never reached by inventing oligos.
"""

from __future__ import annotations

import random
from collections import Counter
from dataclasses import dataclass
from typing import Any, Sequence

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

#: SapI and its reverse complement — the internal cloning site that must appear exactly once.
DEFAULT_RESTRICTION_SITES = ("GCTCTTC", "GAAGAGC")


def rev_comp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def ks_statistic(x: Sequence[float], y: Sequence[float]) -> float:
    """Two-sample Kolmogorov-Smirnov statistic, reported for every generated arm."""
    xs, ys = sorted(x), sorted(y)
    if not xs or not ys:
        return 0.0
    i = j = 0
    d = 0.0
    while i < len(xs) and j < len(ys):
        if xs[i] <= ys[j]:
            i += 1
        else:
            j += 1
        d = max(d, abs(i / len(xs) - j / len(ys)))
    return d


@dataclass(frozen=True)
class OligoArchitecture:
    """The fixed parts of an oligo, around which the variable region is drawn.

    ``forward_priming + variable[:guide_length] + internal_cloning_site +
    variable[guide_length:] + reverse_tail`` makes a ``total_length``-mer, exactly as the
    designed arms are built, so a scrambled control is indistinguishable in architecture
    from a real design and amplifies with the same primers.
    """

    forward_priming: str
    internal_cloning_site: str
    reverse_tail: str
    total_length: int = 200
    guide_length: int = 20

    @property
    def variable_length(self) -> int:
        n = (
            self.total_length
            - len(self.forward_priming)
            - len(self.internal_cloning_site)
            - len(self.reverse_tail)
        )
        if n <= 0:
            raise ValueError(
                f"architecture leaves no variable region: total_length "
                f"{self.total_length} minus fixed parts is {n}"
            )
        return n

    def assemble(self, variable: str) -> str:
        g = self.guide_length
        return (
            self.forward_priming
            + variable[:g]
            + self.internal_cloning_site
            + variable[g:]
            + self.reverse_tail
        )


def _draw(rng: random.Random, target_gc: float, length: int) -> str:
    return "".join(
        rng.choice("GC") if rng.random() < target_gc else rng.choice("AT")
        for _ in range(length)
    )


def make_nontargeting_controls(
    n: int,
    reference_gc: Sequence[float],
    architecture: OligoArchitecture,
    *,
    seed: int = 20261029,
    bin_width: float = 0.02,
    restriction_sites: Sequence[str] = DEFAULT_RESTRICTION_SITES,
    id_prefix: str = "scrambled",
    id_width: int = 5,
    max_tries_per_oligo: int = 4000,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Make exactly *n* non-targeting controls matching the *reference_gc* distribution.

    Parameters
    ----------
    n
        The arm's declared target. The function makes exactly this many or raises.
    reference_gc
        ``gc_frac`` of every real designed oligo — the distribution being matched.
    architecture
        See :class:`OligoArchitecture`.
    seed
        Fixed so a rebuild from the same recipe reproduces the same controls.
    restriction_sites
        A candidate is kept only if these appear exactly once in total, matching the
        designed arms' single internal cloning site.

    Returns
    -------
    (rows, report)
        *rows* carry ``oligo_id``, ``guide_seq``, ``oligo_seq`` and ``gc_frac``; the caller
        maps them onto its own schema. *report* carries the KS statistic against
        *reference_gc*, both medians, and the top-up / trim counts.

    Raises
    ------
    ValueError
        If *n* cannot be reached within the try budget — a short arm is never silently
        shipped.
    """
    if n < 0:
        raise ValueError(f"n must be >= 0, got {n}")
    if n == 0:
        return [], {"requested": 0, "made": 0, "topped_up": 0, "trimmed": 0,
                    "ks_statistic_vs_reference": None}
    if not reference_gc:
        raise ValueError(
            "reference_gc is empty: there is no distribution to match, so a 'GC "
            "distribution-matched' control set cannot be made"
        )

    var_len = architecture.variable_length
    rng = random.Random(seed)

    hist = Counter(int(g / bin_width) for g in reference_gc)
    total = sum(hist.values())
    want_per_bin = {b: int(round(n * c / total)) for b, c in hist.items()}

    def _acceptable(oligo: str) -> bool:
        return sum(oligo.count(s) for s in restriction_sites) == 1

    rows: list[dict[str, Any]] = []

    def _make_in_bin(b: int, k: int, budget: int) -> int:
        """Append up to *k* oligos whose GC lands in bin *b*; return how many were made."""
        lo, hi = b * bin_width, (b + 1) * bin_width
        target = (lo + hi) / 2
        made = tries = 0
        while made < k and tries < budget:
            tries += 1
            variable = _draw(rng, target, var_len)
            oligo = architecture.assemble(variable)
            g = gc_fraction(oligo)
            if not (lo <= g < hi) or not _acceptable(oligo):
                continue
            made += 1
            rows.append(
                {
                    "guide_seq": variable[: architecture.guide_length],
                    "oligo_seq": oligo,
                    "gc_frac": round(g, 5),
                }
            )
        return made

    for b, k in sorted(want_per_bin.items()):
        _make_in_bin(b, k, k * max_tries_per_oligo + 20_000)

    modal_bin = max(hist, key=lambda k: hist[k])

    topped_up = 0
    if len(rows) < n:
        topped_up = n - len(rows)
        _make_in_bin(modal_bin, topped_up, topped_up * 10 * max_tries_per_oligo + 20_000)

    trimmed = 0
    if len(rows) > n:
        trimmed = len(rows) - n
        lo, hi = modal_bin * bin_width, (modal_bin + 1) * bin_width
        drop = [i for i, r in enumerate(rows) if lo <= r["gc_frac"] < hi][-trimmed:]
        if len(drop) < trimmed:  # modal bin too small; trim from the end
            extra = [i for i in range(len(rows) - 1, -1, -1) if i not in set(drop)]
            drop = drop + extra[: trimmed - len(drop)]
        keep = set(range(len(rows))) - set(drop)
        rows = [rows[i] for i in sorted(keep)]

    if len(rows) != n:
        raise ValueError(
            f"asked for {n:,} non-targeting controls, made {len(rows):,}. The GC "
            f"distribution could not be matched within the try budget; widen bin_width or "
            f"raise max_tries_per_oligo rather than accepting a short arm."
        )

    for i, r in enumerate(rows, start=1):
        r["oligo_id"] = f"{id_prefix}_{i:0{id_width}d}"

    made_gc = [r["gc_frac"] for r in rows]
    ref_sorted = sorted(reference_gc)
    report = {
        "requested": n,
        "made": len(rows),
        "topped_up": topped_up,
        "trimmed": trimmed,
        "bin_width": bin_width,
        "seed": seed,
        "ks_statistic_vs_reference": round(ks_statistic(reference_gc, made_gc), 5),
        "reference_gc_median": round(ref_sorted[len(ref_sorted) // 2], 5),
        "generated_gc_median": round(sorted(made_gc)[len(made_gc) // 2], 5),
    }
    return rows, report
