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
    }
)

#: Categories the module can make for itself. Every OTHER category must be supplied by the
#: caller, and omitting one is an error rather than an absence.
GENERATED_KINDS = frozenset({"control_nontargeting"})

#: Selection policy when a supplied arm does not match its target. ``"error"`` keeps every
#: row untouched and lets G-COMPOSITION fail, which is the default and the safe answer.
#: ``"head"``/``"tail"`` trim to the target and are recorded in the report — trimming is
#: possible only because a human declared it in the recipe.
SELECT_POLICIES = frozenset({"error", "head", "tail"})


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
    options: Mapping[str, Any] = field(default_factory=dict)

    @property
    def is_remainder(self) -> bool:
        return self.target is None

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

    def resolved_targets(self) -> dict[str, int]:
        """Every arm's intended count, with the remainder arm's derived from ``total``.

        🔴 This is why a remainder arm cannot absorb the others: its target is
        ``total`` minus the **declared** targets, so if another arm shrinks the remainder's
        target grows to match and the shrunken arm still fails its own target. The 454-oligo
        drift of NN 505 fails here regardless of which arm it moved into.
        """
        declared = {a.name: a.target for a in self.arms if not a.is_remainder}
        out: dict[str, int] = {k: int(v) for k, v in declared.items()}  # type: ignore[arg-type]
        rem = self.remainder_arm
        if rem is not None:
            out[rem.name] = self.total - sum(out.values())
        return out


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_target(raw: Any, arm_name: str) -> int | None:
    if isinstance(raw, str):
        if raw != REMAINDER:
            raise RecipeError(
                f"arm '{arm_name}': target must be an integer or the string "
                f"'{REMAINDER}', got {raw!r}"
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
    if "target" not in body:
        raise RecipeError(
            f"arm '{name}': 'target' is required. Declare a count, or '{REMAINDER}' for "
            f"the single derived arm — an undeclared count is how composition drifted."
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

    known_keys = {
        "kind", "target", "parent", "parent_key", "donors", "gc_match",
        "generate", "select", "coamplify",
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

    remainders = [a.name for a in recipe.arms if a.is_remainder]
    if len(remainders) > 1:
        raise RecipeError(
            f"exactly one arm may be '{REMAINDER}', got {len(remainders)}: "
            f"{sorted(remainders)}. With two derived arms neither one is pinned and a loss "
            f"in one is absorbed by the other, which is the NN 505 failure."
        )

    declared_sum = sum(a.target for a in recipe.arms if not a.is_remainder)  # type: ignore[misc]
    if remainders:
        if declared_sum > recipe.total:
            raise RecipeError(
                f"declared arm targets sum to {declared_sum:,}, which already exceeds "
                f"[pool] total {recipe.total:,}; the remainder arm "
                f"'{remainders[0]}' would need a negative count"
            )
    elif declared_sum != recipe.total:
        raise RecipeError(
            f"declared arm targets sum to {declared_sum:,} but [pool] total is "
            f"{recipe.total:,}. With no '{REMAINDER}' arm the two must agree exactly — "
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
    arms = tuple(_parse_arm(n, b) for n, b in arms_raw.items())

    schema = _parse_schema(data.get("schema", {}) or {})
    segments = _parse_segments(data.get("segments", {}) or {})

    recipe = PoolRecipe(
        name=str(pool.get("name", "unnamed_pool")),
        total=int(pool["total"]),
        arms=arms,
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
