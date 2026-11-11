"""The gates fire on the REAL shipped TRIUMPH pool, with the real NN 505 numbers.

**Why this test exists.** A gate that has never been shown to fire is not evidence, and a
gate shown to fire only on a fixture its author wrote is weak evidence. This pins the two
new gates against the actual incident: the pool that was uploaded, announced to collaborators
at EMBL and had its checksum quoted in email.

**Provenance of every number below**, measured 2026-11-11 and not invented here:

* shipped pool — Oak ``…/20260930_TRIUMPH_EMBL_demo_library/n467_polished_frame/pool_v6/
  496v6_pool.tsv``, 24,000 rows x 27 columns, md5 ``62fe24af21758e97682dce08d5541164``;
* intended composition and the per-category deltas — ``505_category_reconciliation.tsv``;
* the 24 lost columns — ``505_columns_lost.tsv``.

The gates' verdicts depend only on the per-category counts and the column names, so this
test reconstructs a pool at those exact counts rather than vendoring a 9.7 MB TSV. When the
real file happens to be reachable it is used instead, and the two must agree.

🔴 **The trap this pins:** ``len(pool) == 24000`` exactly. The shipped pool hit its target
while four categories drifted and a fifth vanished. Any gate that reads a matching total as
evidence of a correct composition passes this pool. G-COMPOSITION must not.
"""

from __future__ import annotations

import csv
import os
import pathlib

import pytest

from magestic.design.pool import gate_composition, gate_schema, load_recipe_string
from magestic.design.pool.gates import FAIL

#: The shipped v6 pool, by category. Sums to exactly 24,000.
SHIPPED = {
    "natural": 16320,
    "ptc": 2447,
    "syn": 2436,
    "missense_recoded": 1797,
    "scrambled": 1000,
    # "syn_recoded": 0 — the whole category is gone
}

#: What the pool was meant to contain (Library A). Also sums to exactly 24,000.
INTENDED = {
    "natural": 15866,
    "ptc": 2491,
    "syn": 2494,
    "missense_recoded": 1897,
    "syn_recoded": 252,
    "scrambled": 1000,
}

#: NN 505's reconciliation, shipped minus intended.
EXPECTED_DELTAS = {
    "natural": 454,
    "ptc": -44,
    "syn": -58,
    "missense_recoded": -100,
    "syn_recoded": -252,
    "scrambled": 0,
}

#: The 27 columns the shipped pool carries.
SHIPPED_COLUMNS = [
    "oligo_id", "subpool", "control_class", "gene", "common_name", "chrom", "pam_coord",
    "pam_strand", "guide_seq", "guide_start", "guide_end", "donor_start", "donor_end",
    "donor_strand", "edit_pos", "ref_allele", "alt_allele", "donor_strains",
    "variant_class", "consequence", "aa_change", "pam_disruption", "recoded",
    "recoded_changes", "primer_pair", "oligo_seq", "gc_frac",
]

#: The 24 that vanished, 15 of them re-cut scores (8 ``recut_*`` + 7 ``moff_*``).
LOST_COLUMNS = [
    "cds_gene", "cds_consequence", "cds_aa_change", "consequence_source",
    "disruption_position", "disruption_position_recomputed", "protection_class",
    "recut_hazard", "recut_mm_prox", "recut_mm_total", "recut_cfd", "recut_cfd_pam",
    "recut_cfd_n_mm", "recut_cfd_band", "recut_cfd_indel_window",
    "moff_mde", "moff_ce", "moff_gop", "moff_gmt", "moff_raw", "moff_recut", "moff_n_mm",
    "recoded_tract_nt", "recoded_gap_nt",
]

RECIPE_TOML = """
[pool]
name            = "TRIUMPH_yKR888"
total           = 24000
oligo_nt        = 200
reverse_primers = 1
arm_column      = "subpool"

[arm.natural]
kind   = "natural_variant"
target = "remainder"

[arm.ptc]
kind   = "control_ptc"
target = 2491

[arm.syn]
kind   = "control_synonymous"
target = 2494

[arm.missense_recoded]
kind   = "recoded_missense"
target = 1897

[arm.syn_recoded]
kind   = "control_recoding_helper"
target = 252
parent = "missense_recoded"

[arm.scrambled]
kind     = "control_nontargeting"
target   = 1000
gc_match = "designed_pool"

[schema]
columns = %(columns)r
require_groups = ["recut", "moff", "consequence_provenance"]

[schema.groups.recut]
prefix    = "recut_"
min_count = 8

[schema.groups.moff]
prefix    = "moff_"
min_count = 7

[schema.groups.consequence_provenance]
columns = ["cds_gene", "cds_consequence", "cds_aa_change", "consequence_source"]
""" % {"columns": SHIPPED_COLUMNS + LOST_COLUMNS}


def recipe():
    return load_recipe_string(RECIPE_TOML)


def shipped_pool() -> list[dict]:
    """A pool at the shipped per-category counts, carrying the shipped 27 columns."""
    rows = []
    for arm, n in SHIPPED.items():
        for i in range(n):
            row = {c: "" for c in SHIPPED_COLUMNS}
            row["oligo_id"] = f"{arm}_{i:06d}"
            row["subpool"] = arm
            row["primer_pair"] = "revP2"
            rows.append(row)
    return rows


def test_the_shipped_pool_hits_its_target_total_exactly():
    """🔴 The trap. Every existing gate was satisfied by this number."""
    assert sum(SHIPPED.values()) == 24000
    assert sum(INTENDED.values()) == 24000
    assert len(shipped_pool()) == 24000


def test_the_recipe_derives_library_a_s_natural_count():
    """The remainder arm's target is arithmetic, so it reproduces Library A exactly."""
    assert recipe().resolved_targets() == INTENDED


def test_g_composition_fires_on_the_shipped_pool_with_the_nn505_deltas():
    """🔴 The decisive test: real pool, real recipe, the published reconciliation."""
    result = gate_composition(shipped_pool(), recipe())

    assert result.status == FAIL, (
        "the shipped TRIUMPH pool must fail G-COMPOSITION; if this passes, the gate has "
        "not been shown to catch the incident it was written for"
    )
    deltas = {row[0]: row[4] for row in result.rows}
    assert deltas == EXPECTED_DELTAS, (
        "the gate's per-arm deltas must reproduce 505_category_reconciliation.tsv"
    )


def test_the_vanished_control_arm_is_named_as_absent():
    result = gate_composition(shipped_pool(), recipe())
    assert "syn_recoded" in result.detail
    assert "ABSENT" in result.detail
    assert "target 252" in result.detail


def test_g_schema_fires_and_names_exactly_the_24_lost_columns():
    pool = shipped_pool()
    result = gate_schema(pool, recipe(), columns=SHIPPED_COLUMNS)

    assert result.status == FAIL
    named = {row[1] for row in result.rows if row[0] == "column"}
    assert named == set(LOST_COLUMNS), (
        "G-SCHEMA must name exactly the 24 columns 505_columns_lost.tsv measured"
    )
    assert "24 declared column(s) missing" in result.detail


def test_the_fifteen_recut_scores_are_reported_as_two_empty_groups():
    """'All 15 re-cut scores' = 8 ``recut_*`` + 7 ``moff_*``."""
    result = gate_schema(shipped_pool(), recipe(), columns=SHIPPED_COLUMNS)

    assert "recut_* (0 present, 8 required)" in result.detail
    assert "moff_* (0 present, 7 required)" in result.detail
    assert len([c for c in LOST_COLUMNS if c.startswith(("recut_", "moff_"))]) == 15


def test_a_pool_built_to_the_recipe_passes_both_gates():
    """The positive control: without it, "it fails" is not evidence the gate discriminates."""
    rows = []
    for arm, n in INTENDED.items():
        for i in range(n):
            row = {c: "" for c in SHIPPED_COLUMNS + LOST_COLUMNS}
            row["oligo_id"] = f"{arm}_{i:06d}"
            row["subpool"] = arm
            row["primer_pair"] = "revP2"
            rows.append(row)

    assert gate_composition(rows, recipe()).status == "pass"
    assert gate_schema(rows, recipe(), columns=SHIPPED_COLUMNS + LOST_COLUMNS).status == "pass"


# ---------------------------------------------------------------------------
# The same checks against the real file, when it is reachable.
# ---------------------------------------------------------------------------

REAL_POOL_ENV = "MAGESTIC_TRIUMPH_V6_POOL"
REAL_POOL_MD5 = "62fe24af21758e97682dce08d5541164"


def test_the_real_shipped_tsv_agrees_with_the_counts_pinned_here():
    """Runs only where the 9.7 MB pool is reachable; set ``MAGESTIC_TRIUMPH_V6_POOL``."""
    path = os.environ.get(REAL_POOL_ENV)
    if not path or not pathlib.Path(path).exists():
        pytest.skip(f"set {REAL_POOL_ENV} to the shipped 496v6_pool.tsv to run this")

    import hashlib

    digest = hashlib.md5(pathlib.Path(path).read_bytes()).hexdigest()
    assert digest == REAL_POOL_MD5, f"{path} is not the shipped pool ({digest})"

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = list(reader)

    assert header == SHIPPED_COLUMNS
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["subpool"]] = counts.get(r["subpool"], 0) + 1
    assert counts == SHIPPED
    assert len(rows) == 24000

    assert gate_composition(rows, recipe()).status == FAIL
    assert gate_schema(rows, recipe(), columns=header).status == FAIL
