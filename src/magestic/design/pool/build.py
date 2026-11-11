"""``build_pool`` — assemble a pool from its declared recipe.

NN 506 §3, NN 507 §1.

🔴 **The one rule that matters.** ``build_pool`` **enumerates the categories itself from the
recipe** rather than taking whatever arms a caller happens to pass. A ``kind`` the recipe
names but the caller does not supply is an **error, not an omission** — that single rule is
what would have caught the loss of the 252-oligo ``syn_recoded`` control arm the moment the
pipeline changed, instead of a collaborator finding it by eye after the pool shipped.

The symmetric rule holds too: an arm the caller supplies that the recipe does **not** declare
is an error, so a category cannot be added to a pool without being written down either.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from .gates import GateResult, assert_gates, run_gates
from .recipe import PoolRecipe
from .scramble import OligoArchitecture, gc_fraction, make_nontargeting_controls


class PoolBuildError(ValueError):
    """The pool could not be built as declared."""


class MissingCategoryError(PoolBuildError):
    """🔴 The recipe names a category the caller did not supply.

    This is the NN 505 ``syn_recoded`` regression, caught at build time.
    """


class UndeclaredArmError(PoolBuildError):
    """The caller supplied an arm the recipe does not declare."""


class GeneratedArmSuppliedError(PoolBuildError):
    """The caller supplied rows for an arm the module generates.

    Two things would then decide that arm's composition, which is the hodge-podge this
    module exists to end.
    """


@dataclass
class ArmReport:
    """Per-arm accounting: target, built, delta, and the reason for every drop."""

    name: str
    kind: str
    target: int
    supplied: int
    built: int
    generated: bool = False
    trimmed: int = 0
    select: str = "error"
    drops: Mapping[str, int] = field(default_factory=dict)
    extra: Mapping[str, Any] = field(default_factory=dict)

    @property
    def delta(self) -> int:
        return self.built - self.target

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "kind": self.kind,
            "target": self.target,
            "supplied": self.supplied,
            "built": self.built,
            "delta": self.delta,
            "generated": self.generated,
            "trimmed": self.trimmed,
            "select": self.select,
            "drops": dict(self.drops),
            **({"extra": dict(self.extra)} if self.extra else {}),
        }


@dataclass
class BuildReport:
    """Everything a reader needs to see what was intended and what was built."""

    recipe_name: str
    total_target: int
    total_built: int
    arms: list[ArmReport] = field(default_factory=list)
    gates: list[GateResult] = field(default_factory=list)
    reference: str | None = None
    reference_md5: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not any(g.failed for g in self.gates)

    def to_dict(self) -> dict[str, Any]:
        return {
            "recipe": self.recipe_name,
            "total_target": self.total_target,
            "total_built": self.total_built,
            "total_delta": self.total_built - self.total_target,
            "reference": self.reference,
            "reference_md5": self.reference_md5,
            "arms": [a.to_dict() for a in self.arms],
            "gates": [
                {"gate": g.gate, "status": g.status, "detail": g.detail} for g in self.gates
            ],
            "passed": self.passed,
            "notes": list(self.notes),
        }


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_reference(
    recipe: PoolRecipe, reference: str | Path | None, report: BuildReport
) -> None:
    """Refuse a build against a reference the recipe did not declare."""
    if reference is None:
        if recipe.reference:
            report.notes.append(
                f"recipe declares reference {recipe.reference!r} but the build passed none, "
                f"so it was not verified"
            )
        return

    path = Path(reference)
    report.reference = path.name
    if recipe.reference and path.name != recipe.reference:
        raise PoolBuildError(
            f"reference mismatch: the recipe declares {recipe.reference!r} and the build "
            f"was given {path.name!r}. Every design is made against one reference; "
            f"building against another silently changes every coordinate."
        )
    if path.exists() and path.is_file():
        digest = _md5(path)
        report.reference_md5 = digest
        if recipe.reference_md5 and digest != recipe.reference_md5:
            raise PoolBuildError(
                f"reference md5 mismatch for {path.name}: recipe declares "
                f"{recipe.reference_md5}, file is {digest}. The file under the declared "
                f"name is not the declared file."
            )
    elif recipe.reference_md5:
        report.notes.append(
            f"reference {path} is not a readable file, so its md5 was not checked against "
            f"the declared {recipe.reference_md5}"
        )


def _gc_of(row: Mapping[str, Any], sequence_column: str) -> float | None:
    g = row.get("gc_frac")
    if isinstance(g, (int, float)):
        return float(g)
    seq = row.get(sequence_column)
    if isinstance(seq, str) and seq:
        return gc_fraction(seq)
    return None


def build_pool(
    recipe: PoolRecipe,
    arms: Mapping[str, Sequence[Mapping[str, Any]]],
    reference: str | Path | None = None,
    *,
    drops: Mapping[str, Mapping[str, int]] | None = None,
    architecture: OligoArchitecture | Mapping[str, OligoArchitecture] | None = None,
    row_template: Mapping[str, Any] | None = None,
    columns: Sequence[str] | None = None,
    on_gate_failure: str = "raise",
) -> tuple[list[dict[str, Any]], BuildReport]:
    """Assemble the pool the *recipe* declares from the caller's *arms*.

    Parameters
    ----------
    recipe
        The declared composition, from :func:`~magestic.design.pool.recipe.load_recipe`.
    arms
        ``{arm_name: rows}``. 🔴 Every arm the recipe declares and does not generate must
        appear here; a missing one raises :class:`MissingCategoryError`. An arm here that the
        recipe does not declare raises :class:`UndeclaredArmError`.
    reference
        The reference the designs were made against. Checked by name, and by md5 when the
        file is readable, against the recipe's declaration.
    drops
        ``{arm_name: {reason: count}}`` from the caller's own filtering, recorded verbatim in
        the report so every drop has a reason beside it.
    architecture
        An :class:`~magestic.design.pool.scramble.OligoArchitecture` for generated arms —
        one for all, or ``{arm_name: architecture}``.
    row_template
        Defaults merged **under** each generated row, so a generated control carries the
        caller's full schema (empty coordinate columns and so on).
    columns
        The writer's intended column order, passed to G-SCHEMA when the rows themselves are
        sparse.
    on_gate_failure
        ``"raise"`` (default) raises :class:`~magestic.design.pool.gates.PoolGateError`;
        ``"report"`` returns with the failures recorded in ``report.gates``.

    Returns
    -------
    (pool, report)
    """
    if on_gate_failure not in {"raise", "report"}:
        raise ValueError(f"on_gate_failure must be 'raise' or 'report', got {on_gate_failure!r}")

    targets = recipe.resolved_targets()
    report = BuildReport(
        recipe_name=recipe.name,
        total_target=recipe.total,
        total_built=0,
        reference=recipe.reference,
        reference_md5=recipe.reference_md5,
    )
    _check_reference(recipe, reference, report)

    supplied_names = set(arms)
    declared_names = set(recipe.arm_names)

    # ---- 🔴 the rule: enumerate the recipe's categories, never the caller's --------
    undeclared = sorted(supplied_names - declared_names)
    if undeclared:
        raise UndeclaredArmError(
            f"arm(s) {undeclared} were supplied but are not declared in recipe "
            f"'{recipe.name}'"
            + (f" ({recipe.source})" if recipe.source else "")
            + f". Declared arms: {sorted(declared_names)}. Add an [arm.<name>] block with "
            f"its target, or stop passing it — a category that is not written down cannot "
            f"be checked."
        )

    missing: list[str] = []
    for spec in recipe.arms:
        if spec.is_generated:
            continue
        if spec.name not in supplied_names:
            missing.append(spec.name)
    if missing:
        detail = ", ".join(
            f"'{n}' (kind {recipe.arm(n).kind}, target {targets[n]:,})" for n in missing
        )
        raise MissingCategoryError(
            f"recipe '{recipe.name}' declares {len(missing)} categor"
            f"{'y' if len(missing) == 1 else 'ies'} the build did not supply: {detail}. "
            f"🔴 A declared category that is absent is an ERROR, not an omission — this is "
            f"exactly how a 252-oligo control arm left a shipped pool unnoticed (NN 505)."
        )

    for spec in recipe.arms:
        if spec.is_generated and spec.name in supplied_names:
            raise GeneratedArmSuppliedError(
                f"arm '{spec.name}' is generated by the module (generate="
                f"{spec.generate!r}) but the build also supplied {len(arms[spec.name]):,} "
                f"rows for it. Two things would then decide its composition; remove one."
            )

    # ---- assemble, in recipe order ------------------------------------------------
    pool: list[dict[str, Any]] = []
    designed_gc: list[float] = []
    arm_reports: dict[str, ArmReport] = {}
    drops = drops or {}

    # Pass 1: the supplied arms, so the generated arms can match their GC distribution.
    for spec in recipe.arms:
        if spec.is_generated:
            continue
        rows = [dict(r) for r in arms[spec.name]]
        supplied_n = len(rows)
        target = targets[spec.name]
        trimmed = 0
        if spec.select == "head" and supplied_n > target:
            trimmed = supplied_n - target
            rows = rows[:target]
        elif spec.select == "tail" and supplied_n > target:
            trimmed = supplied_n - target
            rows = rows[-target:]

        for r in rows:
            r[recipe.arm_column] = spec.name
            g = _gc_of(r, recipe.sequence_column)
            if g is not None:
                designed_gc.append(g)
        pool.extend(rows)
        arm_reports[spec.name] = ArmReport(
            name=spec.name,
            kind=spec.kind,
            target=target,
            supplied=supplied_n,
            built=len(rows),
            trimmed=trimmed,
            select=spec.select,
            drops=dict(drops.get(spec.name, {})),
        )

    # Pass 2: the generated arms.
    for spec in recipe.arms:
        if not spec.is_generated:
            continue
        target = targets[spec.name]
        arch = (
            architecture.get(spec.name)
            if isinstance(architecture, Mapping)
            else architecture
        )
        if arch is None:
            raise PoolBuildError(
                f"arm '{spec.name}' is generated (generate={spec.generate!r}) but no "
                f"OligoArchitecture was passed for it; build_pool cannot invent the fixed "
                f"priming and cloning sequences."
            )
        rows, gen_report = make_nontargeting_controls(
            target,
            designed_gc,
            arch,
            id_prefix=spec.name,
        )
        full: list[dict[str, Any]] = []
        for r in rows:
            row = dict(row_template or {})
            row.update(r)
            row[recipe.arm_column] = spec.name
            full.append(row)
        pool.extend(full)
        arm_reports[spec.name] = ArmReport(
            name=spec.name,
            kind=spec.kind,
            target=target,
            supplied=0,
            built=len(full),
            generated=True,
            select=spec.select,
            drops=dict(drops.get(spec.name, {})),
            extra=gen_report,
        )

    report.arms = [arm_reports[s.name] for s in recipe.arms]
    report.total_built = len(pool)

    unused_drops = sorted(set(drops) - declared_names)
    if unused_drops:
        report.notes.append(
            f"drop counts were passed for undeclared arm(s) {unused_drops} and were ignored"
        )

    # ---- the gates re-read the same recipe ----------------------------------------
    report.gates = run_gates(pool, recipe, columns=columns)
    if on_gate_failure == "raise":
        assert_gates(report.gates)
    return pool, report
