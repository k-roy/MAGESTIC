"""Shared fixtures for the pool-composition selftests (NN 505 / 506 / 507).

Not a test module (the leading underscore keeps pytest from collecting it). It builds a
small pool with the same *shape* as the TRIUMPH pool — a remainder natural arm, four
designed control arms, one module-generated non-targeting arm, one of the controls parented
to another arm — so each gate's selftest can plant one defect into a pool that is otherwise
correct.
"""

from __future__ import annotations

import random

from magestic.design.pool import OligoArchitecture, load_recipe_string

#: Total 100 = natural 45 (remainder) + ptc 10 + syn 10 + missense_recoded 20
#: + syn_recoded 5 + scrambled 10.
RECIPE_TOML = """
[pool]
name            = "SELFTEST_POOL"
total           = 100
oligo_nt        = 200
reverse_primers = 1
arm_column      = "subpool"

[arm.natural]
kind   = "natural_variant"
target = "remainder"
donors = ["BY", "RM11"]

[arm.ptc]
kind   = "control_ptc"
target = 10

[arm.syn]
kind   = "control_synonymous"
target = 10

[arm.missense_recoded]
kind   = "recoded_missense"
target = 20

[arm.syn_recoded]
kind   = "control_recoding_helper"
target = 5
parent = "missense_recoded"

[arm.scrambled]
kind     = "control_nontargeting"
target   = 10
gc_match = "designed_pool"

[schema]
columns        = ["oligo_id", "subpool", "gene", "oligo_seq", "primer_pair"]
require_groups = ["recut_cfd"]

[schema.groups.recut_cfd]
prefix    = "recut_cfd_"
min_count = 3

[segments]
forward_priming = [0, 20]
guide           = [20, 40]
cloning_site    = [40, 51]
donor           = [51, -20]
reverse_priming = [-20, 0]
"""

FORWARD_PRIMING = "CTGCGATTGGCAGGCGCGCC"
CLONING_SITE = "GTTTGAAGAGC"
REVERSE_TAIL = "CGATAGACGACTGGACAGCA"

ARCHITECTURE = OligoArchitecture(
    forward_priming=FORWARD_PRIMING,
    internal_cloning_site=CLONING_SITE,
    reverse_tail=REVERSE_TAIL,
    total_length=200,
    guide_length=20,
)

#: The declared targets of the designed (caller-supplied) arms.
DESIGNED_TARGETS = {
    "natural": 45,
    "ptc": 10,
    "syn": 10,
    "missense_recoded": 20,
    "syn_recoded": 5,
}


def recipe():
    """A fresh parsed recipe (never share a mutated one between tests)."""
    return load_recipe_string(RECIPE_TOML)


def _oligo(rng: random.Random) -> str:
    """A 200-mer with the same architecture as a real design."""
    var_len = 200 - len(FORWARD_PRIMING) - len(CLONING_SITE) - len(REVERSE_TAIL)
    variable = "".join(rng.choice("ACGT") for _ in range(var_len))
    return (
        FORWARD_PRIMING
        + variable[:20]
        + CLONING_SITE
        + variable[20:]
        + REVERSE_TAIL
    )


def make_rows(arm: str, n: int, *, seed: int = 7, primer: str = "revP1") -> list[dict]:
    """*n* rows for *arm*, carrying the full declared schema."""
    rng = random.Random(f"{arm}:{seed}")
    rows = []
    for i in range(n):
        rows.append(
            {
                "oligo_id": f"{arm}_{i:05d}",
                "gene": f"Y{arm[:2].upper()}{i % 7:03d}W",
                "oligo_seq": _oligo(rng),
                "primer_pair": primer,
                "recut_cfd_max": 0.1,
                "recut_cfd_sum": 0.2,
                "recut_cfd_n": 3,
            }
        )
    return rows


def designed_arms(**overrides: int) -> dict[str, list[dict]]:
    """The caller-supplied arms at their declared targets.

    ``designed_arms(ptc=0)`` plants a loss; ``designed_arms(natural=50)`` plants a drift.
    """
    counts = dict(DESIGNED_TARGETS)
    counts.update(overrides)
    return {arm: make_rows(arm, n) for arm, n in counts.items()}


def built_pool(**overrides: int) -> list[dict]:
    """A finished pool WITHOUT running build_pool, so a gate can be tested on its own.

    The scrambled arm is included as plain rows here; its real generation is covered by
    ``test_pool_scramble.py``.
    """
    rows: list[dict] = []
    counts = dict(DESIGNED_TARGETS)
    counts["scrambled"] = 10
    counts.update(overrides)
    for arm, n in counts.items():
        for r in make_rows(arm, n):
            r["subpool"] = arm
            rows.append(r)
    return rows
