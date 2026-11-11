"""A per-build diff against the previous pool, by category and by oligo segment.

NN 506 §5, NN 507 §6. The model is NN 505 §4: of 6,964 rebuilt designs, **6,603 differed only
in the reverse priming site and 394 in the donor**. That split turns "6,964 oligos changed"
into two numbers with different meanings and different costs, and it takes one pass to
compute. A pool rebuild that reports only a single changed-row count hides which of those two
it was.

Segments are declared in the recipe::

    [segments]
    forward_priming = [0, 20]
    guide           = [20, 40]
    cloning_site    = [40, 51]
    donor           = [51, -20]
    reverse_priming = [-20, 0]     # 0 as the end means "to the 3' end"

A negative index counts from the 3' end, as in Python slicing.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from .recipe import PoolRecipe


@dataclass
class DiffReport:
    """What changed between two pools."""

    n_old: int = 0
    n_new: int = 0
    added: int = 0
    removed: int = 0
    changed: int = 0
    unchanged: int = 0
    by_arm: dict[str, dict[str, int]] = field(default_factory=dict)
    #: ``{segment: n}`` — rows whose ONLY difference is that segment. The NN 505 split.
    changed_only_in: dict[str, int] = field(default_factory=dict)
    #: ``{segment: n}`` — rows in which that segment differs, however many others do too.
    changed_in: dict[str, int] = field(default_factory=dict)
    #: Rows that changed in more than one segment at once.
    changed_in_multiple: int = 0
    #: Rows that changed but no declared segment explains it (sequence lengths differ, or
    #: the declared segments do not cover the difference).
    changed_unattributed: int = 0
    added_ids: list[Any] = field(default_factory=list)
    removed_ids: list[Any] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_old": self.n_old,
            "n_new": self.n_new,
            "added": self.added,
            "removed": self.removed,
            "changed": self.changed,
            "unchanged": self.unchanged,
            "by_arm": self.by_arm,
            "changed_only_in": self.changed_only_in,
            "changed_in": self.changed_in,
            "changed_in_multiple": self.changed_in_multiple,
            "changed_unattributed": self.changed_unattributed,
            "added_ids": self.added_ids[:50],
            "removed_ids": self.removed_ids[:50],
        }

    def summary(self) -> str:
        lines = [
            f"old {self.n_old:,} -> new {self.n_new:,}: "
            f"{self.added:,} added, {self.removed:,} removed, "
            f"{self.changed:,} changed, {self.unchanged:,} unchanged",
        ]
        for arm, counts in sorted(self.by_arm.items()):
            lines.append(
                f"  {arm}: old {counts.get('old', 0):,} new {counts.get('new', 0):,} "
                f"(+{counts.get('added', 0):,} -{counts.get('removed', 0):,} "
                f"~{counts.get('changed', 0):,})"
            )
        if self.changed_only_in:
            lines.append("  changed in exactly one segment:")
            for seg, n in sorted(
                self.changed_only_in.items(), key=lambda kv: -kv[1]
            ):
                lines.append(f"    {seg}: {n:,}")
        if self.changed_in_multiple:
            lines.append(f"  changed in more than one segment: {self.changed_in_multiple:,}")
        if self.changed_unattributed:
            lines.append(
                f"  changed but not explained by a declared segment: "
                f"{self.changed_unattributed:,}"
            )
        return "\n".join(lines)


def _slice(seq: str, span: tuple[int, int | None]) -> str:
    start, end = span
    if end is None or end == 0:
        return seq[start:]
    return seq[start:end]


def diff_pools(
    old: Sequence[Mapping[str, Any]],
    new: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe | None = None,
    *,
    id_column: str | None = None,
    sequence_column: str | None = None,
    arm_column: str | None = None,
    segments: Mapping[str, tuple[int, int | None]] | None = None,
) -> DiffReport:
    """Diff two pools by id, reporting by category and by oligo segment.

    *recipe* supplies the column names and segment definitions; any of them can be
    overridden with the keyword arguments, which is what lets this run against a legacy pool
    written before the recipe existed.
    """
    id_column = id_column or (recipe.id_column if recipe else "oligo_id")
    sequence_column = sequence_column or (recipe.sequence_column if recipe else "oligo_seq")
    arm_column = arm_column or (recipe.arm_column if recipe else "arm")
    if segments is None:
        segments = dict(recipe.segments) if recipe else {}

    old_by_id = {r.get(id_column): r for r in old}
    new_by_id = {r.get(id_column): r for r in new}

    rep = DiffReport(n_old=len(old), n_new=len(new))
    by_arm: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for r in old:
        by_arm[str(r.get(arm_column))]["old"] += 1
    for r in new:
        by_arm[str(r.get(arm_column))]["new"] += 1

    only_in: Counter = Counter()
    in_seg: Counter = Counter()

    for oid, nrow in new_by_id.items():
        orow = old_by_id.get(oid)
        arm = str(nrow.get(arm_column))
        if orow is None:
            rep.added += 1
            by_arm[arm]["added"] += 1
            rep.added_ids.append(oid)
            continue
        oseq = str(orow.get(sequence_column, "") or "")
        nseq = str(nrow.get(sequence_column, "") or "")
        if oseq == nseq:
            rep.unchanged += 1
            continue
        rep.changed += 1
        by_arm[arm]["changed"] += 1

        differing = [
            name
            for name, span in segments.items()
            if _slice(oseq, span) != _slice(nseq, span)
        ]
        if not differing:
            rep.changed_unattributed += 1
            continue
        for name in differing:
            in_seg[name] += 1
        if len(differing) == 1:
            only_in[differing[0]] += 1
        else:
            rep.changed_in_multiple += 1

    for oid, orow in old_by_id.items():
        if oid not in new_by_id:
            rep.removed += 1
            by_arm[str(orow.get(arm_column))]["removed"] += 1
            rep.removed_ids.append(oid)

    rep.by_arm = {k: dict(v) for k, v in sorted(by_arm.items())}
    rep.changed_only_in = dict(only_in)
    rep.changed_in = dict(in_seg)
    return rep
