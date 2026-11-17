"""The gates the composition module owns — G-COMPOSITION, G-SCHEMA, G-PARENT, G-COAMPLIFY.

NN 506 §4, NN 507 §3. The package's existing suite is strong on *is this pool well-formed*
(200 nt, one SapI site, no duplicate id, every donor reconstructs from the reference) and had
no notion of *is this the pool we meant*. These four gates ask the second question, by
re-reading the same ``pool_recipe.toml`` the build read.

🔴 **G-COMPOSITION never pads to fit.** If the arms do not sum to ``total`` it fails. The
whole incident is a pool that hit its target by letting one arm absorb the others, so a total
that matches must never again be read as evidence that the composition is right.

Every gate returns a :class:`GateResult` rather than raising, so a report can carry all four
verdicts; :func:`assert_gates` turns failures into one loud exception. A gate that cannot
apply says so explicitly (``status="not_applicable"`` with a reason) — a silent skip is the
failure mode this package exists to end.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Sequence

from .recipe import PoolRecipe

#: Gate verdicts.
PASS = "pass"
FAIL = "fail"
NOT_APPLICABLE = "not_applicable"

#: Separator for a list-valued parent key on one row (see :func:`gate_parent`).
PARENT_KEY_SEP = ";"


class PoolGateError(AssertionError):
    """One or more gates failed. Carries the failing :class:`GateResult` objects."""

    def __init__(self, message: str, results: Sequence["GateResult"]):
        super().__init__(message)
        self.results = tuple(results)

    @property
    def failures(self) -> tuple["GateResult", ...]:
        return tuple(r for r in self.results if r.failed)


@dataclass(frozen=True)
class GateResult:
    """One gate's verdict."""

    gate: str
    status: str
    detail: str = ""
    rows: tuple[tuple[Any, ...], ...] = ()
    header: tuple[str, ...] = ()

    @property
    def failed(self) -> bool:
        return self.status == FAIL

    def __str__(self) -> str:  # pragma: no cover - formatting only
        mark = {PASS: "PASS", FAIL: "FAIL", NOT_APPLICABLE: "n/a "}[self.status]
        return f"[{mark}] {self.gate}: {self.detail}"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _columns_of(pool: Sequence[Mapping[str, Any]]) -> list[str]:
    """Union of column names across the pool, order of first appearance."""
    seen: dict[str, None] = {}
    for row in pool:
        for k in row:
            seen.setdefault(k, None)
    return list(seen)


def _counts_by_arm(pool: Sequence[Mapping[str, Any]], arm_column: str) -> Counter:
    return Counter(row.get(arm_column) for row in pool)


# ---------------------------------------------------------------------------
# G-COMPOSITION
# ---------------------------------------------------------------------------

def gate_composition(
    pool: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe,
) -> GateResult:
    """🔴 Every arm's built count equals its target, and the sum equals ``total``.

    Catches all three NN 505 regressions:

    * a whole arm absent — its built count is 0 against a target of 252;
    * an arm partially lost — built < target;
    * composition drifting while the total holds — the remainder arm's target is derived
      from ``total`` minus the **declared** targets, so a shrunken arm still fails its own
      target and an over-supplied arm pushes the sum past ``total``.

    It **never pads to fit**: a mismatch is a failure, not something to correct.
    """
    arm_column = recipe.arm_column
    targets = recipe.resolved_targets()
    built = _counts_by_arm(pool, arm_column)

    header = ("arm", "kind", "target", "built", "delta")
    rows: list[tuple[Any, ...]] = []
    problems: list[str] = []

    for spec in recipe.arms:
        want = targets[spec.name]
        have = built.get(spec.name, 0)
        delta = have - want
        rows.append((spec.name, spec.kind, want, have, delta))
        if delta != 0:
            if have == 0:
                problems.append(
                    f"arm '{spec.name}' ({spec.kind}) is ABSENT from the pool: "
                    f"target {want:,}, built 0"
                )
            else:
                problems.append(
                    f"arm '{spec.name}' ({spec.kind}): target {want:,}, built {have:,} "
                    f"({delta:+,})"
                )

    undeclared = {
        str(name): n
        for name, n in built.items()
        if name not in targets
    }
    for name, n in sorted(undeclared.items()):
        rows.append((name, "UNDECLARED", 0, n, n))
        problems.append(
            f"arm '{name}' is in the pool but is not declared in the recipe ({n:,} oligos)"
        )

    total_built = len(pool)
    if total_built != recipe.total:
        problems.append(
            f"pool total: declared {recipe.total:,}, built {total_built:,} "
            f"({total_built - recipe.total:+,}). G-COMPOSITION never pads to fit."
        )

    if problems:
        return GateResult(
            gate="G-COMPOSITION",
            status=FAIL,
            detail="; ".join(problems),
            rows=tuple(rows),
            header=header,
        )
    return GateResult(
        gate="G-COMPOSITION",
        status=PASS,
        detail=(
            f"{len(recipe.arms)} arms, {total_built:,} oligos, every arm on target"
        ),
        rows=tuple(rows),
        header=header,
    )


# ---------------------------------------------------------------------------
# G-SCHEMA
# ---------------------------------------------------------------------------

def gate_schema(
    pool: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe,
    columns: Sequence[str] | None = None,
) -> GateResult:
    """Every column and column-group the recipe names is present in the output.

    Catches the NN 505 loss of 24 of 51 annotation columns, including all 15 re-cut scores —
    a 47 % loss of the schema that was found only by diffing two pools months later.

    *columns* lets a caller pass the writer's intended column order when the pool rows
    themselves are sparse; by default the union of keys across the pool is used.
    """
    if not recipe.schema.columns and not recipe.schema.require_groups:
        return GateResult(
            gate="G-SCHEMA",
            status=NOT_APPLICABLE,
            detail=(
                "the recipe declares no [schema] columns or require_groups, so there is "
                "nothing to check — add one to pin the output schema"
            ),
        )

    present = list(columns) if columns is not None else _columns_of(pool)
    missing_cols, missing_groups = recipe.schema.missing(present)

    if not missing_cols and not missing_groups:
        return GateResult(
            gate="G-SCHEMA",
            status=PASS,
            detail=(
                f"{len(recipe.schema.columns)} declared columns and "
                f"{len(recipe.schema.require_groups)} column-group(s) all present "
                f"({len(present)} columns in the output)"
            ),
        )

    problems: list[str] = []
    rows: list[tuple[Any, ...]] = []
    if missing_cols:
        problems.append(
            f"{len(missing_cols)} declared column(s) missing: {', '.join(missing_cols)}"
        )
        rows.extend(("column", c, "missing") for c in missing_cols)
    for gname, gaps in missing_groups.items():
        problems.append(f"column-group '{gname}' incomplete: {', '.join(gaps)}")
        rows.extend(("group", f"{gname}:{g}", "missing") for g in gaps)

    return GateResult(
        gate="G-SCHEMA",
        status=FAIL,
        detail="; ".join(problems),
        rows=tuple(rows),
        header=("what", "name", "status"),
    )


# ---------------------------------------------------------------------------
# G-PARENT
# ---------------------------------------------------------------------------

def gate_populated(
    pool: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe,
) -> GateResult:
    """🔴 Does every column the recipe says must CARRY A VALUE actually carry one?

    G-SCHEMA asks whether a column is present. This asks whether it is populated, on the rows
    where it means something. The two are not the same question, and the TRIUMPH v6 pool is the
    proof: it shipped `aa_change` present and empty for all 1,821 of its natural missense designs
    (NN 512 finding 1). The column was in the header, so the schema diff and G-SCHEMA both passed
    it, and a collaborator holding the file could not say which substitution any natural missense
    variant makes.
    """
    rules = list(recipe.schema.populated.values())
    if not rules:
        return GateResult(
            gate="G-POPULATED",
            status=NOT_APPLICABLE,
            detail="the recipe declares no [schema.populated.*] rule",
        )

    header = ("rule", "column", "in_scope", "populated", "fraction", "required")
    rows: list[tuple[Any, ...]] = []
    problems: list[str] = []
    for rule in rules:
        n, got = rule.check(pool)
        if n == 0:
            rows.append((rule.name, ",".join(rule.columns), 0, 0, "-", rule.min_fraction))
            problems.append(
                f"rule '{rule.name}' matched NO rows (where={dict(rule.where)}) — a population "
                f"rule that selects nothing checks nothing, so it is a defect in the rule"
            )
            continue
        for col, k in got.items():
            frac = k / n
            rows.append((rule.name, col, n, k, f"{frac:.4f}", rule.min_fraction))
            if frac < rule.min_fraction:
                problems.append(
                    f"'{col}' is populated on {k:,} of {n:,} rows in scope "
                    f"({100 * frac:.2f}%), below the declared {100 * rule.min_fraction:.2f}%"
                    + (f" (where={dict(rule.where)})" if rule.where else "")
                )

    if problems:
        return GateResult(
            gate="G-POPULATED",
            status=FAIL,
            detail="; ".join(problems),
            rows=tuple(rows),
            header=header,
        )
    return GateResult(
        gate="G-POPULATED",
        status=PASS,
        detail=f"{len(rules)} population rule(s) satisfied over "
               f"{sum(1 for _ in rules)} scope(s)",
        rows=tuple(rows),
        header=header,
    )


def gate_parent(
    pool: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe,
) -> GateResult:
    """Every control with a ``parent`` points at a design present in the pool.

    A recoding helper whose recoded-missense designs are absent controls for nothing. When
    the child arm declares a ``parent_key`` the check is row by row; otherwise it is arm by
    arm, which is the weaker check and says so in its detail.
    """
    children = [a for a in recipe.arms if a.parent is not None]
    if not children:
        return GateResult(
            gate="G-PARENT",
            status=NOT_APPLICABLE,
            detail="no arm declares a parent",
        )

    arm_column = recipe.arm_column
    by_arm: dict[Any, list[Mapping[str, Any]]] = defaultdict(list)
    for row in pool:
        by_arm[row.get(arm_column)].append(row)

    def parent_arms(child: Any) -> list[str]:
        """🔴 ONE control arm can control SEVERAL parent arms.

        The helper arm holds the silent change for every recoded-alternative category, so a
        single parent name here reports every helper of the other categories as an orphan. A
        string stays a string so every existing recipe is unchanged.
        """
        p = child.parent
        return [p] if isinstance(p, str) else list(p)

    problems: list[str] = []
    rows: list[tuple[Any, ...]] = []
    row_level = 0
    for child in children:
        kids = by_arm.get(child.name, [])
        pnames = parent_arms(child)
        pname = ", ".join(pnames)
        parents = [r for n in pnames for r in by_arm.get(n, [])]
        if not kids:
            # An absent child is G-COMPOSITION's finding, not this gate's.
            rows.append((child.name, pname, 0, len(parents), "child absent"))
            continue
        if not parents:
            problems.append(
                f"arm '{child.name}' ({len(kids):,} oligos) declares parent "
                f"'{pname}', which is ABSENT from the pool — it controls for nothing"
            )
            rows.append((child.name, pname, len(kids), 0, "parent absent"))
            continue
        if child.parent_key:
            key = child.parent_key
            parent_keys = {p.get(key) for p in parents}
            # 🔴 A child may legitimately carry SEVERAL parent keys: one molecule can serve two
            # parents once duplicate designs are collapsed onto a single row, which is why the
            # recipe already allows `parent` to be a list. Resolving the value ATOMICALLY read such
            # a row as an orphan and failed a correct pool (measured 2026-11-17: "1 of 1,344 rows
            # have a recode_pair_key absent from parent"). Every element must resolve, so this is
            # strictly stronger than the old check, not a relaxation of it.
            def _elements(value: Any) -> list[Any]:
                if isinstance(value, str) and PARENT_KEY_SEP in value:
                    return [v.strip() for v in value.split(PARENT_KEY_SEP) if v.strip()]
                return [value]

            orphans = [
                k.get(recipe.id_column)
                for k in kids
                if any(e not in parent_keys for e in _elements(k.get(key)))
            ]
            row_level += 1
            if orphans:
                problems.append(
                    f"arm '{child.name}': {len(orphans):,} of {len(kids):,} rows have a "
                    f"{key} absent from parent '{pname}', e.g. {orphans[:3]}"
                )
                rows.append(
                    (child.name, pname, len(kids), len(parents),
                     f"{len(orphans)} orphan rows")
                )
            else:
                rows.append(
                    (child.name, pname, len(kids), len(parents), "every row matched")
                )
        else:
            rows.append(
                (child.name, pname, len(kids), len(parents), "arm-level only")
            )

    header = ("arm", "parent", "child_rows", "parent_rows", "note")
    if problems:
        return GateResult(
            gate="G-PARENT",
            status=FAIL,
            detail="; ".join(problems),
            rows=tuple(rows),
            header=header,
        )
    arm_only = len(children) - row_level
    note = (
        f"{len(children)} parented arm(s) checked"
        + (f", {row_level} row-by-row" if row_level else "")
        + (
            f"; {arm_only} checked arm-level only (declare parent_key for a row-by-row check)"
            if arm_only
            else ""
        )
    )
    return GateResult(gate="G-PARENT", status=PASS, detail=note, rows=tuple(rows), header=header)


# ---------------------------------------------------------------------------
# G-COAMPLIFY
# ---------------------------------------------------------------------------

def gate_coamplify(
    pool: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe,
) -> GateResult:
    """Every control shares a reverse priming site with its parent, or the pool has one.

    A control amplified by a primer its parent does not share is measured in a different
    reaction from the thing it controls, which is the defect an earlier library shipped.

    A control with a declared ``parent`` is checked against that parent. A control with no
    parent — a non-targeting scramble is the clear case — controls for the pool as a whole
    rather than for one arm, so it is checked against the **designed** (non-control) arms.
    An arm may be exempted only by declaring ``coamplify = false`` in the recipe; an
    exemption is never inferred.
    """
    primer_col = recipe.primer_column
    present_cols = set(_columns_of(pool))
    if primer_col not in present_cols:
        return GateResult(
            gate="G-COAMPLIFY",
            status=NOT_APPLICABLE,
            detail=(
                f"the pool carries no '{primer_col}' column; set [pool] primer_column to "
                f"the column naming each oligo's reverse priming site to enable this gate"
            ),
        )

    primers = {row.get(primer_col) for row in pool}
    if len(primers) <= 1:
        return GateResult(
            gate="G-COAMPLIFY",
            status=PASS,
            detail=(
                f"the whole pool shares one reverse priming site "
                f"({next(iter(primers), None)!r}), so every control co-amplifies"
            ),
        )

    arm_column = recipe.arm_column
    by_arm: dict[Any, set] = defaultdict(set)
    for row in pool:
        by_arm[row.get(arm_column)].add(row.get(primer_col))

    problems: list[str] = []
    rows: list[tuple[Any, ...]] = []
    controls = [
        a for a in recipe.arms
        if a.kind.startswith("control_") and a.coamplify and by_arm.get(a.name)
    ]
    if not controls:
        return GateResult(
            gate="G-COAMPLIFY",
            status=NOT_APPLICABLE,
            detail=(
                f"the pool has {len(primers)} reverse priming sites but no non-exempt "
                f"control arm is present to check"
            ),
        )

    # A control with no parent controls for the pool as a whole rather than for one arm
    # (a non-targeting scramble is the clear case), so it must co-amplify with the DESIGNED
    # arms — the things it will be compared against.
    designed_primers: set = set()
    for spec in recipe.arms:
        if not spec.kind.startswith("control_"):
            designed_primers |= by_arm.get(spec.name, set())

    for arm in controls:
        mine = by_arm.get(arm.name, set())
        if arm.parent is None:
            theirs = designed_primers or (
                set().union(*(v for k, v in by_arm.items() if k != arm.name)) or set()
            )
            shared = mine & theirs
            rows.append(
                (arm.name, "<designed arms>", sorted(mine), sorted(theirs),
                 "shared" if shared else "DISJOINT")
            )
            if not shared:
                problems.append(
                    f"control arm '{arm.name}' ({sorted(mine)!r}) shares no reverse priming "
                    f"site with the designed arms ({sorted(theirs)!r}), so it is amplified "
                    f"in a different reaction from everything it controls for — give it a "
                    f"shared priming site, declare a parent, or set coamplify = false"
                )
            continue
        theirs = by_arm.get(arm.parent, set())
        shared = mine & theirs
        rows.append(
            (arm.name, arm.parent, sorted(mine), sorted(theirs),
             "shared" if shared else "DISJOINT")
        )
        if not shared:
            problems.append(
                f"control arm '{arm.name}' ({sorted(mine)!r}) shares no reverse priming "
                f"site with its parent '{arm.parent}' ({sorted(theirs)!r})"
            )

    header = ("arm", "parent", "arm_primers", "parent_primers", "note")
    if problems:
        return GateResult(
            gate="G-COAMPLIFY",
            status=FAIL,
            detail="; ".join(problems),
            rows=tuple(rows),
            header=header,
        )
    return GateResult(
        gate="G-COAMPLIFY",
        status=PASS,
        detail=(
            f"{len(controls)} control arm(s) each share a reverse priming site with their "
            f"parent across {len(primers)} sites in the pool"
        ),
        rows=tuple(rows),
        header=header,
    )


# ---------------------------------------------------------------------------
# runner
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# G-RESITE and G-RECUT
#
# 🔴 WHY THESE EXIST (2026-11-17). The design code already refuses a bad oligo as it builds one:
# ``remove_bspqi_sites_from_donor`` drops a design whose extra restriction site cannot be removed,
# and ``guide_disruption`` drops a guide the repaired donor would still match. Both are correct and
# both are per-oligo, at design time. **Nothing re-asked either question of the assembled pool.**
# So when a LATER stage rebuilt a donor — which is exactly what the codon unit rule does when it
# splits a linked variant — the rebuilt molecule was never re-checked, and a 24,000-oligo order was
# assembled carrying 63 donors with an internal SapI site (the cloning enzyme would cut the insert)
# and 127 whose repair product is still a perfect substrate for its own guide (the edit is installed
# and then cut back out). The pool that shipped before the rebuild had ZERO of either.
#
# A guard that runs only at design time protects only designs that nothing touches afterwards.
# These gates ask the same two questions of whatever actually ends up in the pool.
# ---------------------------------------------------------------------------

def _donor_bounds(recipe: PoolRecipe) -> tuple[int, int] | None:
    seg = recipe.segments.get("donor") if recipe.segments else None
    if not seg:
        return None
    lo, hi = seg
    if hi is None:
        hi = recipe.oligo_nt
    if lo is None or hi is None:
        return None
    return int(lo), int(hi)


def gate_restriction_sites(
    pool: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe,
) -> GateResult:
    """G-RESITE — the cloning enzyme must cut the linker and nothing else.

    Fails if any oligo carries a declared recognition site INSIDE the donor window, or if the
    per-oligo site count differs from ``restriction_sites_per_oligo``.
    """
    sites = recipe.restriction_sites
    if not sites:
        return GateResult(
            gate="G-RESITE", status=NOT_APPLICABLE,
            detail="no [cloning] restriction_sites declared, so nothing to check",
        )
    seqcol = recipe.sequence_column
    bounds = _donor_bounds(recipe)
    want = recipe.restriction_sites_per_oligo
    internal: list[tuple[Any, ...]] = []
    miscount: list[tuple[Any, ...]] = []
    for row in pool:
        seq = str(row.get(seqcol) or "").upper()
        if not seq:
            continue
        offsets: list[int] = []
        for site in sites:
            start = seq.find(site)
            while start != -1:
                offsets.append(start)
                start = seq.find(site, start + 1)
        if want is not None and len(offsets) != want:
            miscount.append((row.get(recipe.id_column), row.get(recipe.arm_column),
                             f"{len(offsets)} site(s), expected {want}"))
        if bounds is not None:
            lo, hi = bounds
            inside = sorted(o for o in offsets if lo <= o < hi)
            if inside:
                internal.append((row.get(recipe.id_column), row.get(recipe.arm_column),
                                 ",".join(str(i) for i in inside)))
    if internal or miscount:
        bits = []
        if internal:
            bits.append(
                f"{len(internal):,} oligo(s) carry a restriction site INSIDE the donor window "
                f"{bounds} — the cloning enzyme would cut the insert"
            )
        if miscount:
            bits.append(
                f"{len(miscount):,} oligo(s) do not carry exactly {want} site(s)"
            )
        return GateResult(
            gate="G-RESITE", status=FAIL, detail="; ".join(bits),
            header=("oligo_id", "arm", "detail"),
            # both lists carry (oligo_id, arm, detail); keeping one shape is what makes this
            # safe to concatenate -- an earlier version mixed a 4-tuple in and raised here.
            rows=tuple((a, b, c) for a, b, c in (internal + miscount)[:20]),
        )
    detail = f"every oligo carries exactly {want} site(s)" if want is not None else "sites counted"
    if bounds is not None:
        detail += f" and none inside the donor window {bounds}"
    return GateResult(gate="G-RESITE", status=PASS, detail=detail)


def gate_recut(
    pool: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe,
    guide_column: str = "guide_seq",
) -> GateResult:
    """G-RECUT — the repair product must not still be a perfect target for its own guide.

    Scans each oligo's donor window for ``guide + PAM`` on either strand. A perfect match means
    the installed edit is cut out again and the lineage is lost, so it fails.
    """
    pam = recipe.recut_pam
    if not pam:
        return GateResult(
            gate="G-RECUT", status=NOT_APPLICABLE,
            detail="no [cloning] recut_pam declared, so nothing to check",
        )
    bounds = _donor_bounds(recipe)
    if bounds is None:
        return GateResult(
            gate="G-RECUT", status=NOT_APPLICABLE,
            detail="no donor segment declared in [segments], so the donor cannot be located",
        )
    lo, hi = bounds
    seqcol = recipe.sequence_column
    comp = str.maketrans("ACGTN", "TGCAN")

    def matches(window: str, guide: str) -> bool:
        if len(window) != len(guide) + len(pam):
            return False
        if window[: len(guide)] != guide:
            return False
        tail = window[len(guide):]
        return all(p == "N" or p == t for p, t in zip(pam, tail))

    hits: list[tuple[Any, ...]] = []
    checked = 0
    for row in pool:
        guide = str(row.get(guide_column) or "").upper()
        seq = str(row.get(seqcol) or "").upper()
        if not guide or not seq:
            continue
        donor = seq[lo:hi]
        checked += 1
        span = len(guide) + len(pam)
        rc_donor = donor.translate(comp)[::-1]
        retained = any(
            matches(d[i: i + span], guide)
            for d in (donor, rc_donor)
            for i in range(max(0, len(d) - span + 1))
        )
        if retained:
            hits.append((row.get(recipe.id_column), row.get(recipe.arm_column), guide))
    if hits:
        return GateResult(
            gate="G-RECUT", status=FAIL,
            detail=(
                f"{len(hits):,} of {checked:,} design(s) whose donor still contains "
                f"guide+{pam} on one strand, so the repaired site is re-cut and the edit lost"
            ),
            header=("oligo_id", "arm", "guide_seq"),
            rows=tuple(hits[:20]),
        )
    return GateResult(
        gate="G-RECUT", status=PASS,
        detail=f"no donor of {checked:,} design(s) retains guide+{pam} on either strand",
    )


def run_gates(
    pool: Sequence[Mapping[str, Any]],
    recipe: PoolRecipe,
    columns: Sequence[str] | None = None,
) -> list[GateResult]:
    """Run every composition gate and return their verdicts in a fixed order."""
    return [
        gate_composition(pool, recipe),
        gate_schema(pool, recipe, columns=columns),
        # G-POPULATED runs right after G-SCHEMA because it asks the question G-SCHEMA cannot:
        # presence is not population (NN 512 finding 1).
        gate_populated(pool, recipe),
        gate_parent(pool, recipe),
        gate_coamplify(pool, recipe),
        # G-RESITE and G-RECUT re-ask at POOL level what the design code already enforces per
        # oligo, because a later stage can rebuild a donor and nothing else would notice.
        gate_restriction_sites(pool, recipe),
        gate_recut(pool, recipe),
    ]


def format_gate_results(results: Iterable[GateResult]) -> str:
    """A compact multi-line rendering, used in reports and in the exception message."""
    lines = []
    for r in results:
        lines.append(str(r))
        if r.failed and r.rows:
            if r.header:
                lines.append("        " + " | ".join(str(h) for h in r.header))
            for row in r.rows:
                lines.append("        " + " | ".join(str(c) for c in row))
    return "\n".join(lines)


def assert_gates(results: Sequence[GateResult]) -> None:
    """Raise :class:`PoolGateError` if any gate failed, with every verdict in the message."""
    failures = [r for r in results if r.failed]
    if not failures:
        return
    names = ", ".join(r.gate for r in failures)
    raise PoolGateError(
        f"{len(failures)} composition gate(s) failed: {names}\n"
        + format_gate_results(results),
        results,
    )
