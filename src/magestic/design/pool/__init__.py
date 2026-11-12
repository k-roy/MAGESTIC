"""
magestic.design.pool — declared pool composition, read by the build AND by the gate.

**Why this package exists (NN 505, NN 506, NN 507).** The TRIUMPH guide-donor pool was
rebuilt onto a new reference and three regressions reached a pool that was uploaded,
announced to collaborators and had its checksum quoted in email:

===========================================  ========  ====================================
regression                                   size      how it escaped
===========================================  ========  ====================================
the ``syn_recoded`` control arm disappeared   252       the arm was never named anywhere a
                                                        checker could read
24 of 51 annotation columns disappeared,     47 % of   nothing compared the shipped schema
including all 15 re-cut scores                the        against an intended one
                                              schema
composition drifted, ``natural`` absorbing   454       🔴 **the total stayed at exactly
454 slots from every other category                     24,000**
===========================================  ========  ====================================

The common root is one sentence: *the intended composition existed only in prose and in a
sequence of ``--arm`` flags typed into an sbatch script, so nothing in the build could read
it and nothing could check it.*

This package makes composition a **declared input**. One file, ``pool_recipe.toml``, is the
single source of truth for what the pool should contain, and it is an input to the build and
to the gate alike. The build makes what it says; the gate re-reads it and refuses a pool that
disagrees.

🔴 **The rule that would have caught the ``syn_recoded`` loss the moment the pipeline
changed:** :func:`~magestic.design.pool.build.build_pool` enumerates the categories **itself**
from the recipe rather than taking whatever arms a caller happens to pass. A ``kind`` the
recipe names but the caller does not supply is an **error, not an omission**.

🔴 **G-COMPOSITION never pads to fit.** If the arms do not sum to ``total``, it fails. The
whole incident is a pool that hit its target by letting one arm absorb the others, so a total
that matches its target must never again be read as evidence that the composition is right.

Vocabulary warning
------------------
🔴 **The word "subpool" is overloaded and this package does not use it.** In
:mod:`magestic.design.oligo_generation` a *subpool* is a **PAM-preference synthesis subpool**
(``subpool_idx``, an index into ``rev_primer_sequences``, NGNG before NGNH). In a screen's
pool table "subpool" has meant the **category** (natural, ptc, syn, missense_recoded,
scrambled). Those are different things and wiring them together would be a bug.

This package says **arm** (equivalently *category*) and never ``subpool_idx``. A caller whose
output schema spells the category column ``subpool`` says so once, in the recipe, via
``[pool] arm_column = "subpool"``.

Modules
-------
recipe
    :class:`PoolRecipe`, :class:`ArmSpec`, :class:`SchemaSpec` and :func:`load_recipe` —
    parsing and validation of ``pool_recipe.toml``.
build
    :func:`build_pool` — the core function; enumerates the recipe's categories itself.
gates
    G-COMPOSITION, G-SCHEMA, G-PARENT, G-COAMPLIFY and :func:`run_gates`.
scramble
    GC **distribution**-matched non-targeting controls (lifted from the screen packer, whose
    exact-N top-up and trim were each paid for by a failed build).
diff
    :func:`diff_pools` — a per-build diff against the previous pool, by category and by oligo
    segment.
"""

from .recipe import (
    PoolRecipe,
    ArmSpec,
    SchemaSpec,
    ColumnGroup,
    load_recipe,
    load_recipe_string,
    KNOWN_KINDS,
    GENERATED_KINDS,
    REMAINDER,
    RecipeError,
)
from .build import (
    build_pool,
    BuildReport,
    ArmReport,
    MissingCategoryError,
    UndeclaredArmError,
    GeneratedArmSuppliedError,
    PoolBuildError,
)
from .gates import (
    GateResult,
    PoolGateError,
    gate_composition,
    gate_schema,
    gate_parent,
    gate_populated,
    gate_coamplify,
    run_gates,
    assert_gates,
    format_gate_results,
)
from .scramble import make_nontargeting_controls, OligoArchitecture, ks_statistic
from .diff import diff_pools, DiffReport

__all__ = [
    # recipe
    "PoolRecipe",
    "ArmSpec",
    "SchemaSpec",
    "ColumnGroup",
    "load_recipe",
    "load_recipe_string",
    "KNOWN_KINDS",
    "GENERATED_KINDS",
    "REMAINDER",
    "RecipeError",
    # build
    "build_pool",
    "BuildReport",
    "ArmReport",
    "MissingCategoryError",
    "UndeclaredArmError",
    "GeneratedArmSuppliedError",
    "PoolBuildError",
    # gates
    "GateResult",
    "PoolGateError",
    "gate_composition",
    "gate_schema",
    "gate_parent",
    "gate_populated",
    "gate_coamplify",
    "run_gates",
    "assert_gates",
    "format_gate_results",
    # scramble
    "make_nontargeting_controls",
    "OligoArchitecture",
    "ks_statistic",
    # diff
    "diff_pools",
    "DiffReport",
]
