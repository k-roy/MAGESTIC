"""Parse and validate ``pool_recipe.toml`` — the declared composition of a pool.

NN 506 §2, NN 507 §2. The recipe is read by the **build** and by the **gate**, so the two
cannot drift: the build makes what the recipe says, and the gate re-reads the same file and
refuses a pool that disagrees.

A worked recipe::

    [pool]
    name            = "TRIUMPH_yKR888"
    total           = 24000       # the gate FAILS if the built pool differs; it never pads
    reference       = "467_scaffold_composite.fasta"
    reference_md5   = "4dd4ec2024c1c754f1d386c49582214e"
    oligo_nt        = 200
    reverse_primers = 1
    arm_column      = "subpool"   # this screen's schema spells the category column "subpool"

    [arm.natural]
    kind   = "natural_variant"
    target = "remainder"          # exactly one arm may be the remainder
    donors = ["BY", "RM11", "M22", "Y10", "PW5", "YJM981"]

    [arm.ptc]
    kind   = "control_ptc"
    target = 2447

    [arm.syn_recoded]
    kind   = "control_recoding_helper"
    target = 252
    parent = "missense_recoded"   # a helper whose parent is absent controls for nothing

    [arm.scrambled]
    kind     = "control_nontargeting"
    target   = 1000
    gc_match = "designed_pool"    # KS against the real GC distribution, reported every build

    [schema]
    columns        = ["oligo_id", "subpool", "control_class", "gene"]
    require_groups = ["recut_cfd"]

    [schema.groups.recut_cfd]
    prefix    = "recut_cfd_"
    min_count = 15

🔴 **Vocabulary.** An *arm* is a **category** (natural, ptc, syn, …). It is **not** the
package's PAM-preference synthesis ``subpool_idx``. See the package docstring.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

try:  # Python >= 3.11
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on 3.10 only
    import tomli as tomllib  # type: ignore[no-redef]


#: ``target = "remainder"`` — the one arm whose count is derived rather than declared.
REMAINDER = "remainder"

#: ``target = "fill"`` — the one arm that absorbs the slots left over, and that MAY be trimmed
#: by its declared ranking to reach them exactly.
#:
#: 🔴 The difference from :data:`REMAINDER` is the whole point. A *remainder* arm's count is
#: derived and then checked: if it does not match, G-COMPOSITION fails. A *fill* arm is declared
#: to be the flexible one, so trimming it is its stated job rather than a silent correction, and
#: the rule it is trimmed by is written in the recipe. Every OTHER arm still has to match its
#: own target exactly, so a loss elsewhere still fails the gate — a fill arm cannot absorb one
#: the way ``natural`` absorbed 454 oligos in NN 505. An UNDER-supplied fill arm fails too.
FILL = "fill"

#: The two derived target forms. At most one arm in a recipe may use either.
DERIVED = frozenset({REMAINDER, FILL})

#: Every category the module understands. A recipe naming anything else is an error, so a
#: typo in a category name cannot create a silently-empty arm.
KNOWN_KINDS = frozenset(
    {
        "natural_variant",
        "control_ptc",
        "control_synonymous",
        "recoded_missense",
        "control_recoding_helper",
        "control_nontargeting",
        # Extra guides for a variant whose first-choice design is most at risk of being re-cut —
        # the "second and third shots on goal". A distinct KIND because it is a distinct purpose:
        # not another variant, but another attempt at one already in the pool.
        "natural_extra_shots",
    }
)

#: Categories the module can make for itself. Every OTHER category must be supplied by the
#: caller, and omitting one is an error rather than an absence.
GENERATED_KINDS = frozenset({"control_nontargeting"})

#: Selection policy when a supplied arm does not match its target. ``"error"`` keeps every
#: row untouched and lets G-COMPOSITION fail, which is the default and the safe answer.
#: ``"head"``/``"tail"`` trim to the target and are recorded in the report — trimming is
#: possible only because a human declared it in the recipe.
#:
#: ``"rank"`` trims to the target by a DECLARED ordering (``rank_by``), which is the only
#: honest way to cut an over-supplied arm: ``head``/``tail`` cut by whatever order the caller
#: happened to pass, so the same recipe over the same designs could drop different oligos.
#: ``rank`` always appends the id column as the final tie-break, so a rebuild reproduces the
#: same pool exactly.
SELECT_POLICIES = frozenset({"error", "head", "tail", "rank"})


class RecipeError(ValueError):
    """The recipe is malformed, self-contradictory, or names an unknown category."""


@dataclass(frozen=True)
class ColumnGroup:
    """A named set of columns the output must carry.

    Either an explicit ``columns`` list, or a ``prefix`` with a ``min_count``. The prefix
    form is what catches the NN 505 loss of all 15 re-cut score columns: the group says
    "at least 15 columns beginning ``recut_cfd_``" without pinning their exact names.
    """

    name: str
    columns: tuple[str, ...] = ()
    prefix: str | None = None
    min_count: int = 1

    def missing_from(self, present: Sequence[str]) -> tuple[str, ...]:
        """Return what this group is missing from *present*, empty if satisfied."""
        present_set = set(present)
        if self.columns:
            return tuple(c for c in self.columns if c not in present_set)
        if self.prefix is not None:
            n = sum(1 for c in present_set if c.startswith(self.prefix))
            if n < self.min_count:
                return (f"{self.prefix}* ({n} present, {self.min_count} required)",)
        return ()


@dataclass(frozen=True)
class SchemaSpec:
    """The intended output schema: exact columns plus named column-groups."""

    columns: tuple[str, ...] = ()
    groups: Mapping[str, ColumnGroup] = field(default_factory=dict)
    require_groups: tuple[str, ...] = ()

    def missing(self, present: Sequence[str]) -> tuple[tuple[str, ...], dict[str, tuple[str, ...]]]:
        """``(missing_columns, {group_name: missing_members})`` against *present*."""
        present_set = set(present)
        missing_cols = tuple(c for c in self.columns if c not in present_set)
        missing_groups: dict[str, tuple[str, ...]] = {}
        for gname in self.require_groups:
            group = self.groups.get(gname)
            if group is None:
                missing_groups[gname] = (f"group '{gname}' is required but not defined",)
                continue
            gaps = group.missing_from(present)
            if gaps:
                missing_groups[gname] = gaps
        return missing_cols, missing_groups


@dataclass(frozen=True)
class ArmSpec:
    """One declared category of the pool.

    Attributes
    ----------
    name
        The arm's key in the recipe and the value written into the pool's category column.
    kind
        One of :data:`KNOWN_KINDS`.
    target
        The declared count, or ``None`` for the single remainder arm.
    parent
        For a control that only means something beside another arm (a recoding helper beside
        its recoded missense designs). Checked by G-PARENT and G-COAMPLIFY.
    parent_key
        Optional column joining each child row to a parent row. When set, G-PARENT checks
        row by row instead of only arm by arm.
    generate
        ``"scrambled"`` for a module-generated arm; ``None`` means the caller must supply it.
    select
        See :data:`SELECT_POLICIES`.
    coamplify
        Set ``false`` to declare that this arm is deliberately not co-amplifiable with a
        parent, which exempts it from G-COAMPLIFY. An exemption must be declared; it is
        never inferred.
    """

    name: str
    kind: str
    target: int | None = None
    parent: str | None = None
    parent_key: str | None = None
    donors: tuple[str, ...] = ()
    gc_match: str | None = None
    generate: str | None = None
    select: str = "error"
    coamplify: bool = True
    #: ``"remainder"``, ``"fill"``, or ``None`` when the target is declared outright.
    derive: str | None = None
    #: Columns to order by when ``select == "rank"``. A leading ``-`` means descending.
    rank_by: tuple[str, ...] = ()
    #: ``pct`` as written in the recipe, kept for the report; ``target`` already holds the
    #: resolved integer.
    pct: float | None = None
    options: Mapping[str, Any] = field(default_factory=dict)

    @property
    def is_remainder(self) -> bool:
        """True only for a ``remainder`` arm, whose count must match once derived."""
        return self.derive == REMAINDER

    @property
    def is_fill(self) -> bool:
        """True only for a ``fill`` arm, which may be trimmed by its declared ranking."""
        return self.derive == FILL

    @property
    def is_derived(self) -> bool:
        """True for either derived form — the arm whose target comes from ``total``."""
        return self.derive is not None

    @property
    def is_generated(self) -> bool:
        return self.generate is not None


@dataclass(frozen=True)
class PoolRecipe:
    """The whole declared pool."""

    name: str
    total: int
    arms: tuple[ArmSpec, ...]
    reference: str | None = None
    reference_md5: str | None = None
    oligo_nt: int | None = None
    reverse_primers: int = 1
    schema: SchemaSpec = field(default_factory=SchemaSpec)
    segments: Mapping[str, tuple[int, int | None]] = field(default_factory=dict)
    #: Arms the recipe declares but switches off with ``enabled = false``. Reported, never silent.
    disabled_arms: tuple[str, ...] = ()
    arm_column: str = "arm"
    id_column: str = "oligo_id"
    sequence_column: str = "oligo_seq"
    primer_column: str = "primer_pair"
    source: Path | None = None

    # -- lookups ---------------------------------------------------------------
    def arm(self, name: str) -> ArmSpec:
        for a in self.arms:
            if a.name == name:
                return a
        raise KeyError(name)

    @property
    def arm_names(self) -> tuple[str, ...]:
        return tuple(a.name for a in self.arms)

    @property
    def remainder_arm(self) -> ArmSpec | None:
        for a in self.arms:
            if a.is_remainder:
                return a
        return None

    @property
    def fill_arm(self) -> ArmSpec | None:
        for a in self.arms:
            if a.is_fill:
                return a
        return None

    @property
    def derived_arm(self) -> ArmSpec | None:
        """The one arm whose target comes from ``total``, in either derived form."""
        for a in self.arms:
            if a.is_derived:
                return a
        return None

    def resolved_targets(self) -> dict[str, int]:
        """Every arm's intended count, with the derived arm's taken from ``total``.

        🔴 This is why a derived arm cannot absorb the others: its target is
        ``total`` minus the **declared** targets, so if another arm shrinks the derived arm's
        target grows to match and the shrunken arm still fails its own target. The 454-oligo
        drift of NN 505 fails here regardless of which arm it moved into. That holds for a
        ``fill`` arm exactly as it does for a ``remainder`` arm — what ``fill`` changes is only
        whether an OVER-supply of that one arm may be trimmed to its target, not whether the
        other arms are checked.
        """
        declared = {a.name: a.target for a in self.arms if not a.is_derived}
        out: dict[str, int] = {k: int(v) for k, v in declared.items()}  # type: ignore[arg-type]
        der = self.derived_arm
        if der is not None:
            out[der.name] = self.total - sum(out.values())
        return out


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_target(raw: Any, arm_name: str) -> int | None:
    if isinstance(raw, str):
        if raw not in DERIVED:
            raise RecipeError(
                f"arm '{arm_name}': target must be an integer or one of "
                f"{sorted(DERIVED)}, got {raw!r}"
            )
        return None
    if isinstance(raw, bool) or not isinstance(raw, int):
        raise RecipeError(f"arm '{arm_name}': target must be an integer, got {raw!r}")
    if raw < 0:
        raise RecipeError(f"arm '{arm_name}': target must be >= 0, got {raw}")
    return raw


def _parse_arm(name: str, body: Mapping[str, Any]) -> ArmSpec:
    if not isinstance(body, Mapping):
        raise RecipeError(f"arm '{name}' must be a table")
    kind = body.get("kind")
    if kind is None:
        raise RecipeError(f"arm '{name}': 'kind' is required")
    if kind not in KNOWN_KINDS:
        raise RecipeError(
            f"arm '{name}': unknown kind {kind!r}. Known kinds: "
            f"{', '.join(sorted(KNOWN_KINDS))}"
        )
    # `pct` is an alternative to an absolute count, for an arm sized as a share of the pool —
    # the non-targeting floor is specified that way ("1 to 5% of the pool"), so the recipe can
    # now say so directly instead of carrying a number whose relation to the total is implicit.
    pct = body.get("pct")
    if pct is not None:
        if "target" in body:
            raise RecipeError(
                f"arm '{name}': give either 'target' or 'pct', not both — two numbers that "
                f"can disagree is how an intended composition stops being readable."
            )
        if isinstance(pct, bool) or not isinstance(pct, (int, float)):
            raise RecipeError(f"arm '{name}': pct must be a number, got {pct!r}")
        if not 0 < float(pct) <= 100:
            raise RecipeError(f"arm '{name}': pct must be in (0, 100], got {pct}")
        # resolved against [pool] total by load_recipe, which is where total is known
        return _finish_arm(name, body, kind, None, float(pct))

    if "target" not in body:
        raise RecipeError(
            f"arm '{name}': 'target' is required. Declare a count, a 'pct', or one of "
            f"{sorted(DERIVED)} for the single derived arm — an undeclared count is how "
            f"composition drifted."
        )
    target = _parse_target(body["target"], name)

    select = body.get("select", "error")
    if select not in SELECT_POLICIES:
        raise RecipeError(
            f"arm '{name}': select must be one of {sorted(SELECT_POLICIES)}, got {select!r}"
        )

    generate = body.get("generate")
    if generate is None and kind in GENERATED_KINDS:
        generate = "scrambled"
    if generate is not None and kind not in GENERATED_KINDS:
        raise RecipeError(
            f"arm '{name}': kind {kind!r} cannot be generated by the module; only "
            f"{sorted(GENERATED_KINDS)} can."
        )

    return _finish_arm(name, body, kind, target, None, select=select, generate=generate)


def _finish_arm(name: str, body: Mapping[str, Any], kind: str, target: int | None,
                pct: float | None, *, select: str | None = None,
                generate: str | None = None) -> ArmSpec:
    if select is None:
        select = body.get("select", "error")
        if select not in SELECT_POLICIES:
            raise RecipeError(
                f"arm '{name}': select must be one of {sorted(SELECT_POLICIES)}, got {select!r}"
            )
    if generate is None:
        generate = body.get("generate")
        if generate is None and kind in GENERATED_KINDS:
            generate = "scrambled"

    raw_target = body.get("target")
    derive = raw_target if isinstance(raw_target, str) and raw_target in DERIVED else None

    rank_by = tuple(body.get("rank_by", ()) or ())
    if select == "rank" and not rank_by:
        raise RecipeError(
            f"arm '{name}': select='rank' needs 'rank_by'. Trimming an arm without saying "
            f"what to keep is the same as trimming it arbitrarily."
        )
    if rank_by and select != "rank":
        raise RecipeError(
            f"arm '{name}': 'rank_by' is set but select is {select!r}; it would be ignored."
        )

    known_keys = {
        "kind", "target", "pct", "parent", "parent_key", "donors", "gc_match",
        "generate", "select", "coamplify", "rank_by", "enabled",
    }
    options = {k: v for k, v in body.items() if k not in known_keys}

    return ArmSpec(
        name=name,
        kind=kind,
        target=target,
        parent=body.get("parent"),
        parent_key=body.get("parent_key"),
        donors=tuple(body.get("donors", ()) or ()),
        gc_match=body.get("gc_match"),
        generate=generate,
        select=select,
        coamplify=bool(body.get("coamplify", True)),
        derive=derive,
        rank_by=rank_by,
        pct=pct,
        options=options,
    )


def _parse_schema(body: Mapping[str, Any]) -> SchemaSpec:
    groups_raw = body.get("groups", {}) or {}
    groups: dict[str, ColumnGroup] = {}
    for gname, gbody in groups_raw.items():
        if isinstance(gbody, (list, tuple)):
            groups[gname] = ColumnGroup(name=gname, columns=tuple(gbody))
            continue
        if not isinstance(gbody, Mapping):
            raise RecipeError(f"schema group '{gname}' must be a list or a table")
        cols = tuple(gbody.get("columns", ()) or ())
        prefix = gbody.get("prefix")
        if not cols and prefix is None:
            raise RecipeError(
                f"schema group '{gname}': give either 'columns' or 'prefix'"
            )
        groups[gname] = ColumnGroup(
            name=gname,
            columns=cols,
            prefix=prefix,
            min_count=int(gbody.get("min_count", 1)),
        )

    require = tuple(body.get("require_groups", ()) or ())
    for gname in require:
        if gname not in groups:
            raise RecipeError(
                f"schema: require_groups names '{gname}' but no [schema.groups.{gname}] "
                f"defines it"
            )
    return SchemaSpec(
        columns=tuple(body.get("columns", ()) or ()),
        groups=groups,
        require_groups=require,
    )


def _parse_segments(body: Mapping[str, Any]) -> dict[str, tuple[int, int | None]]:
    out: dict[str, tuple[int, int | None]] = {}
    for name, span in body.items():
        if not isinstance(span, (list, tuple)) or len(span) != 2:
            raise RecipeError(
                f"segment '{name}' must be a two-element [start, end] list; use a negative "
                f"index to count from the 3' end"
            )
        start, end = span
        if not isinstance(start, int):
            raise RecipeError(f"segment '{name}': start must be an integer")
        if end is not None and not isinstance(end, int):
            raise RecipeError(f"segment '{name}': end must be an integer or absent")
        out[name] = (int(start), None if end is None else int(end))
    return out


def _validate(recipe: PoolRecipe) -> PoolRecipe:
    if recipe.total <= 0:
        raise RecipeError(f"[pool] total must be > 0, got {recipe.total}")
    if not recipe.arms:
        raise RecipeError("a recipe with no [arm.*] declares nothing and can check nothing")

    names = [a.name for a in recipe.arms]
    dupes = {n for n in names if names.count(n) > 1}
    if dupes:
        raise RecipeError(f"duplicate arm names: {sorted(dupes)}")

    derived = [a.name for a in recipe.arms if a.is_derived]
    if len(derived) > 1:
        raise RecipeError(
            f"exactly one arm may be derived ({sorted(DERIVED)}), got {len(derived)}: "
            f"{sorted(derived)}. With two derived arms neither one is pinned and a loss "
            f"in one is absorbed by the other, which is the NN 505 failure."
        )
    fills = [a for a in recipe.arms if a.is_fill]
    for a in fills:
        if a.select != "rank":
            raise RecipeError(
                f"arm '{a.name}': target='{FILL}' requires select='rank'. A fill arm is the one "
                f"arm allowed to be trimmed to fit, so the rule it is trimmed by has to be "
                f"declared — otherwise the pool depends on the order rows happened to arrive in."
            )

    declared_sum = sum(a.target for a in recipe.arms if not a.is_derived)  # type: ignore[misc]
    if derived:
        if declared_sum > recipe.total:
            raise RecipeError(
                f"declared arm targets sum to {declared_sum:,}, which already exceeds "
                f"[pool] total {recipe.total:,}; the derived arm "
                f"'{derived[0]}' would need a negative count"
            )
    elif declared_sum != recipe.total:
        raise RecipeError(
            f"declared arm targets sum to {declared_sum:,} but [pool] total is "
            f"{recipe.total:,}. With no derived arm the two must agree exactly — "
            f"G-COMPOSITION never pads to fit."
        )

    for a in recipe.arms:
        if a.parent is not None and a.parent not in names:
            raise RecipeError(
                f"arm '{a.name}': parent '{a.parent}' is not a declared arm "
                f"(declared: {sorted(names)})"
            )
        if a.parent == a.name:
            raise RecipeError(f"arm '{a.name}': cannot be its own parent")
        if a.parent_key is not None and a.parent is None:
            raise RecipeError(
                f"arm '{a.name}': parent_key is set but no parent is declared"
            )
    return recipe


def _from_mapping(data: Mapping[str, Any], source: Path | None = None) -> PoolRecipe:
    pool = data.get("pool")
    if not isinstance(pool, Mapping):
        raise RecipeError("a recipe needs a [pool] table")
    if "total" not in pool:
        raise RecipeError("[pool] total is required — it is what G-COMPOSITION checks")

    arms_raw = data.get("arm", {}) or {}
    if not isinstance(arms_raw, Mapping):
        raise RecipeError("[arm.<name>] tables expected")
    total = int(pool["total"])

    # 🔴 `enabled = false` is how a category is left OUT of a pool. It is deliberately not the
    # same as deleting the block: the block stays, diffable, saying what was excluded and why,
    # and the build report names every disabled arm out loud. That is the difference between
    # this and NN 505 — there the missing arm existed nowhere a reader or a checker could see it.
    all_arms = []
    disabled = []
    for n, b in arms_raw.items():
        if isinstance(b, Mapping) and b.get("enabled") is False:
            disabled.append(str(n))
            continue
        all_arms.append(_parse_arm(n, b))

    # resolve `pct` now that [pool] total is known
    arms = tuple(
        (a if a.pct is None else
         ArmSpec(**{**a.__dict__, "target": int(round(total * a.pct / 100.0))}))
        for a in all_arms
    )

    schema = _parse_schema(data.get("schema", {}) or {})
    segments = _parse_segments(data.get("segments", {}) or {})

    recipe = PoolRecipe(
        name=str(pool.get("name", "unnamed_pool")),
        total=total,
        arms=arms,
        disabled_arms=tuple(disabled),
        reference=pool.get("reference"),
        reference_md5=pool.get("reference_md5"),
        oligo_nt=pool.get("oligo_nt"),
        reverse_primers=int(pool.get("reverse_primers", 1)),
        schema=schema,
        segments=segments,
        arm_column=str(pool.get("arm_column", "arm")),
        id_column=str(pool.get("id_column", "oligo_id")),
        sequence_column=str(pool.get("sequence_column", "oligo_seq")),
        primer_column=str(pool.get("primer_column", "primer_pair")),
        source=source,
    )
    return _validate(recipe)


def load_recipe(path: str | Path) -> PoolRecipe:
    """Read and validate ``pool_recipe.toml``.

    Raises :class:`RecipeError` on anything self-contradictory, so a recipe that could not
    have described a real pool fails before a build consumes hours.
    """
    path = Path(path)
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    return _from_mapping(data, source=path)


def load_recipe_string(text: str, source: Path | None = None) -> PoolRecipe:
    """Same as :func:`load_recipe` for a recipe already in memory (used by the selftests)."""
    return _from_mapping(tomllib.loads(text), source=source)
