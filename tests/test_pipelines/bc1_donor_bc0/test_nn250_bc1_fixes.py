"""NN 250 -- regression tests for the defects the bc1 audits (NN 239/240/241/249) unearthed.

Every test in this file was written to FAIL against the code as it stood on 2026-10-01 and to pass
only after the matching fix. The defect ids are the ones in
``planning/222_triumph_embl_library_design/250_bc1_pipeline_fixes.md``: ``D*`` from the original
audit inventory, ``C*`` from the NN 240 known-answer roundtrip.

🔴 TWO HOUSE RULES THIS FILE OBEYS, because the audits were produced by breaking them:

1. **Fixtures use REAL column names.** An earlier ``collapse_gdb_bridge`` test suite was green and
   vacuous because its fixture column was named ``count`` -- a name that does not exist in any gdb
   table (the real ones are ``grouped_guide_bc0_counts`` / ``total_counts`` / ``total_fwd_counts``).
   A fixture built from plausible names tests nothing.
2. **Tie-breaks are tested by PERMUTATION INVARIANCE, not by reading the sort call.** Asserting
   that ``kind="mergesort"`` is present proves nothing: a stable sort *inherits* input order, so the
   defect it is meant to fix survives. Every tie-break test here shuffles the input rows under
   several fixed seeds and asserts the output is identical.
"""

import importlib.util
import itertools
from pathlib import Path

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core import matching as M

_SRC = Path(__file__).resolve().parents[3] / "src/magestic/pipelines/bc1_donor_bc0"
STEP05 = _SRC / "snakemake/05_assign_sublibrary.py"
STEP01 = _SRC / "snakemake/01_parse_bc1_donor_bc0.py"


def _load(path, name):
    """Import one of the digit-prefixed snakemake step scripts by path."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Shared fixture: a two-design library whose donors are perfect-lookup sized
# and whose guides are distinct, built through the REAL builder so the lookup
# structure under test is the production one.
# ---------------------------------------------------------------------------

ANCHOR = M.DONOR_VARIABLE_REGION_ANCHORS[0]      # GTTTGAAGAGC, SpCas9/SpG/SaCas9
# A 98-nt constant 5' flank puts the anchor at DEFAULT_ANCHOR_POS, as a real read does.
FLANK = "A" * M.DEFAULT_ANCHOR_POS

GUIDE_A = "ACGTACGTACGTACGTACGT"
GUIDE_B = "TTGCATTGCATTGCATTGCA"

# 129-nt designed donors differing across the whole length, so neither the perfect lookup nor any
# 60-nt sliding window can confuse them.
DONOR_A = ("ACGT" * 33)[:129]
DONOR_B = ("TTGCA" * 26)[:129]


def _oligo_df():
    """Two designs, in the harmonized table's real column names."""
    return pd.DataFrame(
        {
            "oligo_name": ["design_A", "design_B"],
            "guide": [GUIDE_A, GUIDE_B],
            "donor": [DONOR_A, DONOR_B],
            M.VARIANT_COLUMN: ["chrI:1000:A/T", "chrII:2000:C/G"],
        }
    )


@pytest.fixture()
def tables():
    return M.build_lookup_tables(_oligo_df())


def _read(donor_variable_region):
    """An observed donor as stage 01 emits it: constant flank + anchor + variable region."""
    return FLANK + ANCHOR + donor_variable_region


# ===========================================================================
# D12 -- guide_status == CONFLICT must force is_unambiguous = False
# ===========================================================================
#
# THE ATTACK (NN 241): one sequencing error turns a read into another design's exact donor while
# the read still carries its own design's guide. The donor answer wins, the disagreement is
# recorded in guide_status, and the row nevertheless ships `match_method='perfect'` with
# `is_unambiguous=True` -- so every consumer that filters on is_unambiguous accepts a call the
# matcher itself knows is internally contradicted.


def test_d12_perfect_donor_with_conflicting_guide_is_not_unambiguous(tables):
    """A perfect donor hit on design A carrying design B's guide is a CONFLICT, so not unambiguous."""
    res = M.match_sequence(_read(DONOR_A), GUIDE_B, tables)

    # The donor answer is kept and the conflict is recorded (both pre-existing, both correct).
    assert res.perfect_donor is True
    assert res.oligo_name == "design_A"
    assert res.guide_status == M.GUIDE_CONFLICT

    # THE DEFECT: an internally contradicted call must not advertise itself as unambiguous.
    assert res.is_unambiguous is False, (
        "guide_status=conflict but is_unambiguous=True: the matcher knows the guide names a "
        "different design and still tells every downstream filter the call is clean"
    )


def test_d12_sliding_window_with_conflicting_guide_is_not_unambiguous(tables):
    """Same guard on the sliding-window branch, which has its own copy of the conflict logic."""
    # A 1-nt 3' truncation of DONOR_A: no perfect hit, but every 60-nt window still matches A.
    res = M.match_sequence(_read(DONOR_A[:-1]), GUIDE_B, tables)

    assert res.perfect_donor is False
    assert res.match_method.startswith("sliding_window")
    assert res.oligo_name == "design_A"
    assert res.guide_status == M.GUIDE_CONFLICT
    assert res.is_unambiguous is False


def test_d12_conflict_guard_holds_through_match_dataframe(tables):
    """Through the PRODUCTION entry point, not the inner helper.

    match_dataframe has a vectorized perfect-donor fast path with its own copy of the result
    tuple. It excludes rows whose guide IS a designed guide, so a conflict row should fall through
    to match_sequence -- but that is a property to TEST, not to assume: a fix made only in
    match_sequence would be silently bypassed if the mask ever changed.
    """
    df = pd.DataFrame(
        {
            "donor": [_read(DONOR_A), _read(DONOR_A), _read(DONOR_A)],
            "guide": [GUIDE_A, GUIDE_B, None],
        }
    )
    out = M.match_dataframe(df, tables)

    # row 0: guide agrees -> clean unambiguous match
    assert out.loc[0, "guide_status"] == M.GUIDE_MATCH
    assert bool(out.loc[0, "is_unambiguous"]) is True
    # row 1: guide names the other design -> conflict, and NOT unambiguous
    assert out.loc[1, "guide_status"] == M.GUIDE_CONFLICT
    assert bool(out.loc[1, "is_unambiguous"]) is False
    # row 2: no guide at all -> nothing to conflict with, still unambiguous
    assert out.loc[2, "guide_status"] == M.GUIDE_NO_GUIDE
    assert bool(out.loc[2, "is_unambiguous"]) is True


def test_d12_guard_does_not_touch_the_clean_cases(tables):
    """The guard must fire ONLY on conflict: no-guide, undesigned-guide and agreeing-guide rows
    keep the unambiguity they had, or the fix buys correctness with a loss of sensitivity."""
    # no guide
    assert M.match_sequence(_read(DONOR_A), None, tables).is_unambiguous is True
    # a guide that is not a designed guide at all -- not evidence against the donor
    assert M.match_sequence(_read(DONOR_A), "GGGGGGGGGGGGGGGGGGGG", tables).is_unambiguous is True
    # the right guide
    assert M.match_sequence(_read(DONOR_A), GUIDE_A, tables).is_unambiguous is True
    # guide-only fallback: nothing to conflict WITH
    assert M.match_sequence(None, GUIDE_A, tables).is_unambiguous is True


# ===========================================================================
# D6 -- top_oligo_read_share was computed against a DIFFERENT oligo than the row reports
# ===========================================================================
#
# Step 05 reduces a bc1's many links to one reference row under a four-key ordering
# (counts, num_PCR_replicates, donor, bc0). But `top_oligo_read_share` named "the bc1's top link"
# using a counts-ONLY ordering. Whenever two links tie on counts and differ on
# num_PCR_replicates, the two orderings pick different links, so the share shipped on the row
# describes an oligo the row does not report. Reproduced below: the row reports
# `design_HIGH_rep` and carries share 0.5, which is `design_LOW_rep`'s share.


def _two_links_tied_on_counts():
    """One bc1, two links TIED on counts, differing on num_PCR_replicates and oligo_name.

    Column names are step 05's real ones. The rows are ordered so that input order and the
    canonical tie order DISAGREE -- a stable sort left alone would keep the wrong one.
    """
    return pd.DataFrame(
        {
            "bc1": ["A" * 20] * 2,
            "oligo_name": ["design_LOW_rep", "design_HIGH_rep"],
            "counts": [10, 10],
            "num_PCR_replicates": [1, 3],
            "donor": ["ACGT", "TTGC"],
            "bc0": ["BC0X", "BC0Y"],
        }
    )


def test_d6_the_old_counts_only_ordering_disagreed_with_the_reference_ordering():
    """Documents the defect numerically. The two expressions are FROZEN COPIES of step 05 as it
    stood on 2026-10-01; this test pins the disagreement that motivated the fix, so a future
    reader can see the bug was real rather than inferred from reading the sort call."""
    df = _two_links_tied_on_counts()
    cnt = "counts"

    old_top = (
        df.sort_values([cnt], ascending=False, kind="mergesort")
        .drop_duplicates("bc1")
        .set_index("bc1")["oligo_name"]
    )
    tie_order = {cnt: False, "num_PCR_replicates": False, "donor": True, "bc0": True}
    tie_cols = [c for c in tie_order if c in df.columns]
    ref = df.sort_values(
        tie_cols, ascending=[tie_order[c] for c in tie_cols], kind="mergesort"
    ).drop_duplicates(subset=["bc1"], keep="first")

    assert old_top.iloc[0] == "design_LOW_rep"
    assert ref["oligo_name"].iloc[0] == "design_HIGH_rep"
    assert old_top.iloc[0] != ref["oligo_name"].iloc[0]


def test_d6_one_canonical_link_ordering_exists_and_both_consumers_use_it():
    """THE FIX: a single exported ordering, so the share and the reported oligo cannot diverge."""
    step05 = _load(STEP05, "bdb_step05_d6")
    df = _two_links_tied_on_counts()

    cols, asc = step05.canonical_bc1_link_order(df, "counts")
    # The ordering must consider replicate support, not counts alone -- that is the whole defect.
    assert "num_PCR_replicates" in cols, cols
    assert cols[0] == "counts" and asc[0] is False

    top = step05.top_link_oligo(df, "counts")
    ref = df.sort_values(cols, ascending=asc, kind="mergesort").drop_duplicates(
        subset=["bc1"], keep="first"
    )
    assert top.iloc[0] == ref["oligo_name"].iloc[0] == "design_HIGH_rep", (
        "top_oligo_read_share is computed against a different oligo than the reference row reports"
    )


@pytest.mark.parametrize("seed", [0, 1, 7, 42, 1234])
def test_d6_top_link_is_permutation_invariant(seed):
    """PERMUTATION INVARIANCE, the property that actually matters.

    `kind="mergesort"` alone does not make a tie-break deterministic -- a stable sort inherits the
    input order, so the same bc1 arriving in a different row order would get a different top link.
    The ordering must be total over the tied rows.
    """
    step05 = _load(STEP05, "bdb_step05_d6_perm")
    df = _two_links_tied_on_counts()
    shuffled = df.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    assert (
        step05.top_link_oligo(shuffled, "counts").iloc[0]
        == step05.top_link_oligo(df, "counts").iloc[0]
        == "design_HIGH_rep"
    )


def test_d6_ordering_is_total_over_every_permutation_of_a_full_tie():
    """A tie on EVERY support key still has to resolve to one answer, over all permutations.

    donor/bc0 ascending are the terminal keys, so the answer is the lexicographically smallest
    (donor, bc0) -- deterministic rather than input-order-dependent.
    """
    step05 = _load(STEP05, "bdb_step05_d6_total")
    rows = [
        {"bc1": "A" * 20, "oligo_name": f"d{i}", "counts": 5, "num_PCR_replicates": 2,
         "donor": d, "bc0": b}
        for i, (d, b) in enumerate([("AAAA", "B2"), ("AAAA", "B1"), ("TTTT", "B1")])
    ]
    answers = {
        step05.top_link_oligo(pd.DataFrame(list(p)), "counts").iloc[0]
        for p in itertools.permutations(rows)
    }
    assert len(answers) == 1, f"top link depends on input order: {answers}"


# ===========================================================================
# C17 -- the gdb fallback fragment key disagreed with stage 01's, so the bridge
#        silently recovered ZERO and the run read as a biological null
# ===========================================================================
#
# Stage 01 writes donor_bc0_fragment = donor[-30:] + SPS(20) + bc0(10) = 60 nt
# (01_parse_bc1_donor_bc0.py:568, DONOR_BC0_FRAGMENT_LENGTH = 60). Two fallbacks in
# data_loader.py rebuilt the key when the column was absent and left the SPS out:
#   :585  else-branch (no SPS column)        -> donor[-30:] + bc0            = 40 nt
#   :1019 the gdb path, with NO SPS branch    -> donor[-30:] + bc0            = 40 nt
# A 40-nt key cannot join a 60-nt key, so a gdb table lacking the column bridged zero bc1s --
# and `recover_full_donors` logged "0/N (0.0%)" at INFO and returned normally. That is the
# defect shape that has produced wrong conclusions all night: a pipeline failure that reads as
# a measured biological zero.

from magestic.pipelines.bc1_donor_bc0 import config as CFG            # noqa: E402
from magestic.pipelines.bc1_donor_bc0.readers import data_loader as DL  # noqa: E402

GDB_DONOR = "".join("ACGT"[i % 4] for i in range(129))
GDB_SPS = "A" * CFG.SPS_LENGTH
GDB_BC0 = "C" * CFG.BC0_LENGTH


def _stage01_fragment(donor=GDB_DONOR, sps=GDB_SPS, bc0=GDB_BC0):
    """A frozen copy of stage 01's builder (01_parse_bc1_donor_bc0.py:568)."""
    return donor[-30:] + sps + bc0


def test_c17_stage01_key_is_60nt_and_the_old_sps_less_fallback_was_40nt():
    """Documents the defect numerically, so the length claim is measured rather than asserted."""
    assert len(_stage01_fragment()) == CFG.DONOR_BC0_FRAGMENT_LENGTH == 60
    old_fallback = GDB_DONOR[-30:] + GDB_BC0           # frozen copy of data_loader.py:1019
    assert len(old_fallback) == 40
    assert old_fallback != _stage01_fragment(), "a 40-nt key cannot join a 60-nt key"


def test_c17_canonical_builder_is_shared_and_reproduces_stage01():
    """THE FIX, part 1: ONE builder, so a producer and a fallback cannot drift apart."""
    assert CFG.build_donor_bc0_fragment(GDB_DONOR, GDB_SPS, GDB_BC0) == _stage01_fragment()
    assert len(CFG.build_donor_bc0_fragment(GDB_DONOR, GDB_SPS, GDB_BC0)) == 60


def test_c17_builder_refuses_to_pad_when_the_sps_is_absent():
    """A joinable key cannot be built without the SPS, so the builder must REFUSE rather than
    emit a short key that silently joins nothing."""
    with pytest.raises(ValueError, match="SPS"):
        CFG.build_donor_bc0_fragment(GDB_DONOR, None, GDB_BC0)
    with pytest.raises(ValueError, match="SPS"):
        CFG.build_donor_bc0_fragment(GDB_DONOR, "", GDB_BC0)


def test_c17_gdb_fallback_builds_a_key_that_joins_stage01s():
    """THE FIX, part 2: the gdb fallback (data_loader.py:1019) must use the shared builder.

    Real gdb column names: donor / bc0 / SPS / oligo_name / total_counts.
    """
    gdb = pd.DataFrame(
        {
            "donor": [GDB_DONOR],
            "bc0": [GDB_BC0],
            "SPS": [GDB_SPS],
            "oligo_name": ["design_A"],
            "total_counts": [42],
        }
    )
    out = DL.ensure_donor_bc0_fragment(gdb, source="gdb")
    assert out["donor_bc0_fragment"].iloc[0] == _stage01_fragment()
    assert len(out["donor_bc0_fragment"].iloc[0]) == 60


def test_c17_fallback_refuses_a_gdb_table_with_no_sps_column():
    """Without an SPS column no joinable key exists; building a 40-nt one is worse than failing."""
    gdb = pd.DataFrame({"donor": [GDB_DONOR], "bc0": [GDB_BC0], "total_counts": [42]})
    with pytest.raises(ValueError, match="SPS"):
        DL.ensure_donor_bc0_fragment(gdb, source="gdb")


def test_c17_a_zero_bridge_on_nonempty_input_is_an_error_not_a_quiet_zero():
    """THE FIX, part 3: zero recovery from two non-empty tables is a FAILURE, never a result."""
    bdb = pd.DataFrame({"bc1": ["A" * 20], "donor_bc0_fragment": [_stage01_fragment()]})
    # a lookup keyed on the OLD 40-nt key -- exactly the production failure
    bad_lookup = pd.DataFrame(
        {"donor_bc0_fragment": [GDB_DONOR[-30:] + GDB_BC0], "donor": [GDB_DONOR]}
    )
    with pytest.raises(DL.BridgeRecoveredNothing):
        DL.recover_full_donors(bdb, bad_lookup)


def test_c17_a_matching_lookup_still_bridges_and_is_not_broken_by_the_guard():
    """The guard must not fire on a healthy bridge, or it would break NN 240's 48/48 roundtrip."""
    bdb = pd.DataFrame({"bc1": ["A" * 20], "donor_bc0_fragment": [_stage01_fragment()]})
    good = pd.DataFrame({"donor_bc0_fragment": [_stage01_fragment()], "donor": [GDB_DONOR]})
    out = DL.recover_full_donors(bdb, good)
    assert out["donor"].notna().sum() == 1


def test_c17_guard_stays_quiet_when_there_is_genuinely_nothing_to_bridge():
    """An EMPTY lookup is a different condition from a non-joining one: NN 240's non-vacuity
    control deliberately runs without a bridge, and must keep working."""
    bdb = pd.DataFrame({"bc1": ["A" * 20], "donor_bc0_fragment": [_stage01_fragment()]})
    out = DL.recover_full_donors(bdb, pd.DataFrame())
    assert len(out) == 1


# ===========================================================================
# D5 -- the donor-bc0 lookup collapse was support-blind, uncounted, and could emit
#       a (donor, bc0) pair that never co-occurred in any read
# ===========================================================================
#
# data_loader.py:849-854 collapsed the donor-bc0 table with
#   .groupby('donor_bc0_fragment').agg({'donor':'first','bc0':'first','counts':'sum', ...})
# Three measured failure modes:
#   1. SUPPORT-BLIND  -- 'first' takes whichever row came first in the file. Reproduced: a
#      fragment whose donor_WRONG had 1 read and whose donor_RIGHT had 180 returned donor_WRONG.
#   2. UNCOUNTED      -- nothing in the output says the fragment was ambiguous at all, so the
#      rate cannot be audited. merge.py already picked its best gdb match by max(counts), so a
#      support-blind collapse DEMOTED a correct tie-break to an arbitrary one.
#   3. CHIMERIC PAIR  -- groupby.agg('first') skips NaN INDEPENDENTLY PER COLUMN, so donor can
#      come from one row and bc0 from another. Reproduced exactly: rows (NaN, BC0_A) and
#      (donor_B, NaN) emitted ('donor_B', 'BC0_A'), a pair present in no read.

import numpy as np  # noqa: E402

LOOKUP_FRAG = "F" * 60


def _support_blind_fixture():
    """Worst-supported row FIRST, so 'first' and 'max support' disagree.

    Real column names from the loader's own select list (`counts`, `donor`, `bc0`,
    `donor_bc0_fragment`, `sample_number`).
    """
    return pd.DataFrame(
        {
            "donor_bc0_fragment": [LOOKUP_FRAG] * 3,
            "donor": ["donor_WRONG", "donor_RIGHT", "donor_RIGHT"],
            "bc0": ["BC0_WRONG", "BC0_RIGHT", "BC0_RIGHT"],
            "counts": [1, 100, 80],
            "sample_number": ["s1", "s2", "s3"],
        }
    )


def test_d5_the_old_first_collapse_was_support_blind():
    """Frozen copy of data_loader.py:849-854, pinning the defect as measured."""
    old = (
        _support_blind_fixture()
        .groupby("donor_bc0_fragment")
        .agg({"donor": "first", "bc0": "first", "counts": "sum"})
        .reset_index()
    )
    assert old["donor"].iloc[0] == "donor_WRONG"      # 1 read beat 180
    assert "n_distinct_donors" not in old.columns     # and nothing recorded the ambiguity


def test_d5_the_old_first_collapse_could_emit_a_pair_that_never_co_occurred():
    """Frozen copy. 'first' skips NaN per column, so donor and bc0 come from different rows."""
    c = pd.DataFrame(
        {
            "donor_bc0_fragment": [LOOKUP_FRAG] * 2,
            "donor": [np.nan, "donor_B"],
            "bc0": ["BC0_A", np.nan],
            "counts": [5, 5],
        }
    )
    old = (
        c.groupby("donor_bc0_fragment")
        .agg({"donor": "first", "bc0": "first", "counts": "sum"})
        .reset_index()
    )
    assert (old["donor"].iloc[0], old["bc0"].iloc[0]) == ("donor_B", "BC0_A")
    assert ("donor_B", "BC0_A") not in set(zip(c["donor"], c["bc0"]))


def test_d5_collapse_picks_the_best_supported_donor():
    """THE FIX: support decides, by summed reads across samples -- not file order."""
    out, stats = DL.collapse_donor_bc0_lookup(_support_blind_fixture())
    assert len(out) == 1
    assert out["donor"].iloc[0] == "donor_RIGHT"
    assert out["bc0"].iloc[0] == "BC0_RIGHT"
    assert int(out["counts"].iloc[0]) == 181          # total support preserved
    assert int(out["num_samples_donor_bc0"].iloc[0]) == 3


def test_d5_collapse_counts_what_it_collapsed():
    """THE FIX: the ambiguity is COUNTED and returned, never silently resolved."""
    out, stats = DL.collapse_donor_bc0_lookup(_support_blind_fixture())
    assert int(out["n_distinct_donors"].iloc[0]) == 2
    assert int(out["donor_support"].iloc[0]) == 180    # support behind the WINNER
    assert stats["ambiguous_fragments"] == 1
    assert stats["unique_fragments"] == 1
    assert stats["rows_in"] == 3
    # no silent drop anywhere
    assert stats["rows_in"] - stats["dropped_null_key"] >= stats["rows_out"]


def test_d5_donor_and_bc0_always_come_from_the_same_observation():
    """THE FIX: the whole winning row is taken, so the pair is one that was really seen."""
    c = pd.DataFrame(
        {
            "donor_bc0_fragment": [LOOKUP_FRAG] * 2,
            "donor": [np.nan, "donor_B"],
            "bc0": ["BC0_A", np.nan],
            "counts": [5, 5],
        }
    )
    out, _ = DL.collapse_donor_bc0_lookup(c)
    pair = (out["donor"].iloc[0], out["bc0"].iloc[0])
    assert pair != ("donor_B", "BC0_A"), "emitted a donor/bc0 pair that never co-occurred"
    assert pair in set(zip(c["donor"], c["bc0"])) or pd.isna(pair[0]) or pd.isna(pair[1])


@pytest.mark.parametrize("seed", [0, 1, 7, 42, 1234])
def test_d5_collapse_is_permutation_invariant(seed):
    """Support-based, therefore independent of the order rows arrive in."""
    df = _support_blind_fixture()
    a, _ = DL.collapse_donor_bc0_lookup(df)
    b, _ = DL.collapse_donor_bc0_lookup(df.sample(frac=1.0, random_state=seed).reset_index(drop=True))
    assert a["donor"].iloc[0] == b["donor"].iloc[0] == "donor_RIGHT"
    assert int(a["counts"].iloc[0]) == int(b["counts"].iloc[0])


def test_d5_collapse_resolves_a_total_support_tie_deterministically():
    """Equal support must still give ONE answer over every permutation."""
    rows = [
        {"donor_bc0_fragment": LOOKUP_FRAG, "donor": d, "bc0": b, "counts": 7}
        for d, b in [("donor_A", "BC0_1"), ("donor_B", "BC0_2"), ("donor_C", "BC0_3")]
    ]
    answers = {
        DL.collapse_donor_bc0_lookup(pd.DataFrame(list(p)))[0]["donor"].iloc[0]
        for p in itertools.permutations(rows)
    }
    assert len(answers) == 1, f"collapse depends on input order under a full tie: {answers}"


def test_d5_collapse_handles_the_real_gdb_support_column_names():
    """The fixture-naming trap that made an earlier suite vacuous.

    Verified real gdb header (Oak, V455_guide_donor_bc0_counts.tsv):
      library_ID total_counts total_fwd_counts total_rev_counts leader guide donor SPS bc0
      msd donor_bc0_fragment seq
    A table carrying `total_counts` and no `counts` must still be collapsed BY SUPPORT.
    """
    df = pd.DataFrame(
        {
            "donor_bc0_fragment": [LOOKUP_FRAG] * 2,
            "donor": ["donor_WRONG", "donor_RIGHT"],
            "bc0": ["BC0_WRONG", "BC0_RIGHT"],
            "total_counts": [1, 100],
        }
    )
    out, stats = DL.collapse_donor_bc0_lookup(df)
    assert out["donor"].iloc[0] == "donor_RIGHT", (
        "collapse fell back to first-in-order because it only looked for a column named 'counts'"
    )
    assert stats["support_column"] == "total_counts"


# ===========================================================================
# C18 -- a non-N20 barcode is FABRICATED into a 20-mer and emitted at tier `confirmed`
# ===========================================================================
#
# REPRODUCED 2026-10-01, exactly as reported. `extract_bc1` takes the barcode as a fixed-width
# slice s[b0 : b0+20]. When the true barcode is 17 or 18 nt, that slice runs 2-3 nt INTO the
# post-barcode spacer, so the emitted 20-mer is part real barcode and part spacer. Measured on a
# layout-A read (18-nt spacer, LAYOUT_BY_OFFSET[18]):
#
#   true barcode   tier         anchor_offset   layout         emitted bc1
#   20 nt (right)  confirmed    18              A_1.0/2.5      the real barcode
#   18 nt          confirmed    16              None           barcode + 2 nt of spacer
#   17 nt          confirmed    15              None           barcode + 3 nt of spacer
#
# The evidence is right there on the call -- the confirming offset is NOT the sample's
# established offset -- and the tier says `confirmed` anyway.
#
# A read that truly ENDS short is already handled correctly (tier `length_mismatch`, bc1=None),
# so the module docstring's "nothing is ever padded" is true of truncation and false of this.
#
# THE FIX routes these reads into the EXISTING `confirmed_indel` tier -- the vocabulary already
# has a name for "the anchor is intact but at a shifted position". That means no new column, no
# schema change, and no silent drop: `n_confirmed_indel` is already in TIER_COUNT_COLUMNS, is
# already written by stage 01, and is already read by its consumers. The read keeps its bc1 and
# `counts == sum(n_<tier>)` still holds.

from magestic.utils import bc1_extraction as B  # noqa: E402

_U, _D = B.UPSTREAM_ANCHOR, B.UNIVERSAL_DOWNSTREAM
_SPACER_A = "TTGGCCAATTGGCCAATT"        # 18 nt -> the layout-A offset
_BC20 = "ACGTACGTACGTACGTACGT"
_BC18, _BC17 = _BC20[:18], _BC20[:17]


def _layout_a_read(barcode):
    """A layout-A read: upstream anchor, barcode, 18-nt spacer, universal downstream anchor."""
    return "AAAA" + _U + barcode + _SPACER_A + _D + "GGGGGGGGGG"


def test_c18_a_short_barcode_is_fabricated_from_the_spacer():
    """Documents the defect: the emitted 20-mer is NOT the true barcode."""
    call = B.extract_bc1(_layout_a_read(_BC18), established_offset=18,
                         search_both_strands=False)
    assert call.bc1 is not None and len(call.bc1) == 20
    assert call.bc1 != _BC18
    assert call.bc1 == _BC18 + _SPACER_A[:2], "2 nt of spacer was read as barcode"
    # the evidence that it is fabricated: the confirming offset is not the established one
    assert call.anchor_offset == 16 != 18
    assert call.layout is None


@pytest.mark.parametrize("barcode,expect_offset", [(_BC18, 16), (_BC17, 15)])
def test_c18_fabricated_barcode_must_not_carry_a_bare_confirmed_tier(barcode, expect_offset):
    """THE FIX: when the confirming offset disagrees with the sample's ESTABLISHED offset, the
    barcode's 3' boundary is shifted, so the call is not plain `confirmed`."""
    call = B.extract_bc1(_layout_a_read(barcode), established_offset=18,
                         search_both_strands=False)
    assert call.anchor_offset == expect_offset
    assert call.anchor_tier != "confirmed", (
        f"a {len(barcode)}-nt barcode padded out of the spacer still reports tier 'confirmed'"
    )
    assert call.anchor_tier == "confirmed_indel"
    # no silent drop: the read still yields a barcode, in an EMITTING tier that is already counted
    assert call.bc1 is not None
    assert call.anchor_tier in B.EMITTING_TIERS
    assert f"n_{call.anchor_tier}" in B.TIER_COUNT_COLUMNS


def test_c18_the_correct_20mer_is_untouched():
    """The guard must not demote a correct read, or it costs sensitivity for nothing."""
    call = B.extract_bc1(_layout_a_read(_BC20), established_offset=18,
                         search_both_strands=False)
    assert call.anchor_tier == "confirmed"
    assert call.bc1 == _BC20
    assert call.anchor_offset == 18
    assert call.layout == "A_1.0/2.5"


def test_c18_no_established_offset_means_no_demotion():
    """With no established offset there is nothing to disagree with, so behaviour is unchanged.

    This is the guard's scope limit and it is deliberate: inventing a demotion rule for samples
    whose layout was never established would reclassify reads on no evidence.
    """
    call = B.extract_bc1(_layout_a_read(_BC20), search_both_strands=False)
    assert call.anchor_tier == "confirmed"


def test_c18_offset_agreement_is_exposed_per_read():
    """THE EVIDENCE must be readable, not merely computable: `as_row()` is what the per-read
    debug path writes, and it now carries the offset-agreement verdict."""
    good = B.extract_bc1(_layout_a_read(_BC20), established_offset=18, search_both_strands=False)
    bad = B.extract_bc1(_layout_a_read(_BC18), established_offset=18, search_both_strands=False)
    assert good.barcode_boundary_trusted is True
    assert bad.barcode_boundary_trusted is False
    assert "barcode_boundary_trusted" in good.as_row()
    # and it VARIES -- a constant column names nothing
    assert len({good.barcode_boundary_trusted, bad.barcode_boundary_trusted}) == 2


# ===========================================================================
# C19 -- bc1_donor_bc0 stage 01 weights every collapsed record as ONE read
# ===========================================================================
#
# `bc1_from_screen/02_parse_bc1_reads.py` was fixed to use the canonical
# `magestic.utils.collapsed_fastq` reader. `bc1_donor_bc0/01_parse_bc1_donor_bc0.py` was NOT: it
# has its own `read_fastq_entries`, never looks at the read name, and `_tally` does
# `vec[tier] += 1`. So a collapsed FASTQ is undercounted by exactly its collapse factor
# (measured elsewhere at 1.4-18.9x), and the distortion SATURATES on the abundant lineages
# rather than scaling, so the top barcodes are hit hardest.
#
# Auto-detection is safe here: both recognized patterns are whole-name anchored
# (`^@?(\d+)\s*$` and `count=(\d+)\s*$`), so a plain Illumina read name matches neither and
# keeps weight 1 -- i.e. non-collapsed input behaves exactly as before.

from magestic.utils import collapsed_fastq as CF  # noqa: E402


def test_c19_stage01_does_not_use_the_canonical_reader_by_accident():
    """Baseline: the shared parser resolves both producers' formats."""
    assert CF.parse_collapsed_header_count("@seq_1_count=208") == (208, "count_kv")
    assert CF.parse_collapsed_header_count("@12") == (12, "bare_int")
    with pytest.raises(CF.CollapsedHeaderError):
        CF.parse_collapsed_header_count("@A00123:45:HXXX:1:1101:1000:1000 1:N:0:ATCG")


@pytest.mark.parametrize(
    "header,expect_count,expect_fmt",
    [
        ("@seq_1_count=208", 208, "count_kv"),
        ("@12", 12, "bare_int"),
        # a plain Illumina read name is NOT collapsed and must stay weight 1
        ("@A00123:45:HXXX:1:1101:1000:1000 1:N:0:ATCG", 1, CF.UNKNOWN_FORMAT),
    ],
)
def test_c19_stage01_resolves_a_read_weight_from_the_header(header, expect_count, expect_fmt):
    """THE FIX, part 1: stage 01 can say how many reads a record is worth."""
    step01 = _load(STEP01, "bdb_step01_c19")
    assert step01.read_weight_from_header(header) == (expect_count, expect_fmt)


def test_c19_tally_weights_by_the_collapsed_count():
    """THE FIX, part 2: the weight reaches the counts table.

    Without it a 208-read record contributes 1, so `counts` undercounts by 207 on that record.
    """
    step01 = _load(STEP01, "bdb_step01_c19_tally")
    from collections import Counter

    parsed = {
        "anchor_tier": "confirmed",
        "read_orientation": "fwd",
        "bc1": "A" * 20,
        "donor": "ACGT",
        "bc0": "C" * 10,
        "SPS": "G" * 20,
        "donor_bc0_fragment": "F" * 60,
        "bc1_length": 20,
    }
    results, tier_acct, length_stats = {}, Counter(), Counter()
    assert step01._tally(parsed, results, tier_acct, length_stats, weight=208) is True
    vec = next(iter(results.values()))
    assert sum(vec) == 208, "a 208-read collapsed record was counted as 1 read"
    # the per-tier accounting and the length histogram must be weighted the SAME way,
    # or `counts == sum(n_<tier>)` and the length stats disagree with each other
    assert tier_acct[("confirmed", "fwd")] == 208
    assert length_stats[20] == 208


def test_c19_tally_defaults_to_weight_one():
    """Non-collapsed input is unchanged -- the fix must not silently rescale ordinary runs."""
    step01 = _load(STEP01, "bdb_step01_c19_default")
    from collections import Counter

    parsed = {"anchor_tier": "confirmed", "read_orientation": "fwd", "bc1": "A" * 20,
              "donor": "", "bc0": "", "SPS": "", "donor_bc0_fragment": "", "bc1_length": 20}
    results, tier_acct, length_stats = {}, Counter(), Counter()
    step01._tally(parsed, results, tier_acct, length_stats)
    assert sum(next(iter(results.values()))) == 1


def test_c19_undercount_is_exactly_the_collapse_factor():
    """Quantifies the defect: records vs reads over a mixed-format collapsed file."""
    step01 = _load(STEP01, "bdb_step01_c19_factor")
    headers = ["@seq_1_count=208", "@12", "@seq_2_count=5"]
    reads = sum(step01.read_weight_from_header(h)[0] for h in headers)
    assert len(headers) == 3 and reads == 225
    assert reads / len(headers) == 75.0      # the collapse factor the old code discarded


def test_c19_tier_accounting_balances_against_weighted_reads_not_records(tmp_path):
    """🔴 A REGRESSION THE WEIGHTING WOULD HAVE CAUSED, now guarded.

    `_write_sidecars` asserts that the tier accounting balances. The accounting is in READS, so
    once records are collapse-weighted it must balance against the WEIGHTED total; comparing it
    against the RECORD count would raise on every collapsed run. Both are checked here.
    """
    from collections import Counter

    step01 = _load(STEP01, "bdb_step01_c19_sidecar")
    out = tmp_path / "s1_parsed.tsv"
    out.write_text("")
    tier_acct = Counter({("confirmed", "fwd"): 225})     # 3 records worth 208 + 12 + 5 reads

    # balances against the weighted total
    step01._write_sidecars(str(out), "r1.fq", "r2.fq", None, tier_acct,
                           3, 3, 0, reads_weighted=225,
                           collapse_fmts=Counter({"count_kv": 2, "bare_int": 1}))

    # and would NOT balance against the record count -- the bug this guards
    with pytest.raises(RuntimeError, match="tier accounting"):
        step01._write_sidecars(str(out), "r1.fq", "r2.fq", None, tier_acct, 3, 3, 0)

    import json
    meta = json.loads((tmp_path / "s1_parsed.tsv.parse_meta.json").read_text())
    assert meta["reads_weighted"] == 225 and meta["reads_processed"] == 3
    assert meta["collapse_factor"] == 75.0
    assert meta["collapsed_header_formats"] == {"count_kv": 2, "bare_int": 1}


def test_c19_non_collapsed_input_balances_exactly_as_before(tmp_path):
    """Weight 1 everywhere: records == reads, so the pre-existing invariant is unchanged."""
    from collections import Counter

    step01 = _load(STEP01, "bdb_step01_c19_plain")
    out = tmp_path / "s2_parsed.tsv"
    out.write_text("")
    step01._write_sidecars(str(out), "r1.fq", "r2.fq", None,
                           Counter({("confirmed", "fwd"): 10}), 10, 10, 0,
                           reads_weighted=10,
                           collapse_fmts=Counter({CF.UNKNOWN_FORMAT: 10}))
    import json
    meta = json.loads((tmp_path / "s2_parsed.tsv.parse_meta.json").read_text())
    assert meta["collapse_factor"] == 1.0


# ===========================================================================
# D11 / D7 -- SPS scoping: the shipped mapping is rejected, and absent scoping is silent
# ===========================================================================
#
# The shipped mapping is /oak/.../common/keyfiles/sps_library_mapping.tsv, whose real header is
#   sps  V_library  oligo_pool  project  trafo_context  description
# Three defects, all reproduced 2026-10-01 against that exact structure:
#
#  D11a  `load_sps_mapping` demands `library_id`, so the shipped file raises outright.
#  D11b  one SPS legitimately maps to SEVERAL V_libraries (V629/V545 share one SPS, V631/V546 the
#        other), and the dict comprehension collapsed them LAST-WINS, silently losing V629/V631.
#  D11c  worse, the v_libraries filter ran AFTER that collapse, so asking to scope to V629
#        raised "empty after filtering" -- a legitimate request failing with a misleading
#        message, because V629 had already been overwritten by V545.
#
# Raising rather than silently unscoping is CORRECT and is kept (D11). What changes is which
# inputs are legitimate.

SHIPPED_SPS_ROWS = [
    ("CTTGACTCGACAGTGAGAGC", "V629"),
    ("CTTGACTCGACAGTGAGAGC", "V545"),
    ("CGTATGGAGATGACGTCAGG", "V631"),
    ("CGTATGGAGATGACGTCAGG", "V546"),
]


def _sps_file(tmp_path, col="V_library", rows=None):
    p = tmp_path / f"sps_mapping_{col}.tsv"
    pd.DataFrame(rows or SHIPPED_SPS_ROWS, columns=["sps", col]).to_csv(p, sep="\t", index=False)
    return str(p)


def test_d11_accepts_the_shipped_V_library_column(tmp_path):
    """THE FIX: the column the shipped file actually uses is accepted, by name, not by position."""
    m = M.load_sps_mapping(_sps_file(tmp_path, "V_library"), ["V545"])
    assert m == {"CTTGACTCGACAGTGAGAGC": "V545"}


def test_d11_still_accepts_library_id(tmp_path):
    """The canonical name keeps working -- the alias is additive."""
    m = M.load_sps_mapping(_sps_file(tmp_path, "library_id"), ["V546"])
    assert m == {"CGTATGGAGATGACGTCAGG": "V546"}


def test_d11_scoping_to_a_library_that_is_not_last_in_the_file_works(tmp_path):
    """D11c: V629 appears BEFORE V545 on the same SPS. Filtering must precede the collapse."""
    m = M.load_sps_mapping(_sps_file(tmp_path, "V_library"), ["V629"])
    assert m == {"CTTGACTCGACAGTGAGAGC": "V629"}, (
        "the v_libraries filter ran after a last-wins collapse, so V629 was already gone"
    )


def test_d11_an_sps_ambiguous_WITHIN_scope_is_dropped_loudly_not_guessed(tmp_path, caplog):
    """Silently picking one of two libraries for one SPS is how 7% cross-subpool naming happens.

    With both V629 and V545 in scope that SPS names neither uniquely, so it is DROPPED from the
    scoping map (its rows fall back to the full scope, which is reported per run) rather than
    resolved by file order. The other SPS, which V546 names uniquely in this scope, survives.
    """
    import logging

    with caplog.at_level(logging.WARNING, logger=M.logger.name):
        m = M.load_sps_mapping(_sps_file(tmp_path, "V_library"), ["V629", "V545", "V546"])
    assert m == {"CGTATGGAGATGACGTCAGG": "V546"}
    assert "CTTGACTCGACAGTGAGAGC" not in m
    assert any("MORE THAN ONE" in r.message or "more than one" in r.message
               for r in caplog.records), [r.message for r in caplog.records]


def test_d11_unscoped_load_does_not_silently_lose_libraries(tmp_path):
    """Called with no scope, the SHIPPED file is ambiguous for EVERY SPS, so nothing can be
    scoped and the loader must refuse rather than return a map that quietly answers with
    whichever row came last. The message must name ambiguity as the cause, not just emptiness."""
    with pytest.raises(ValueError, match="shared by several libraries|uniquely names"):
        M.load_sps_mapping(_sps_file(tmp_path, "V_library"))


def test_d11_a_genuinely_unambiguous_unscoped_file_still_loads(tmp_path):
    """One SPS, one library: nothing to refuse."""
    rows = [("CTTGACTCGACAGTGAGAGC", "V545"), ("CGTATGGAGATGACGTCAGG", "V546")]
    m = M.load_sps_mapping(_sps_file(tmp_path, "V_library", rows))
    assert m == {"CTTGACTCGACAGTGAGAGC": "V545", "CGTATGGAGATGACGTCAGG": "V546"}


def test_d11_empty_after_filtering_still_refuses(tmp_path):
    """The existing, CORRECT behaviour: never run unscoped while claiming to be scoped."""
    with pytest.raises(ValueError, match="empty after filtering"):
        M.load_sps_mapping(_sps_file(tmp_path, "V_library"), ["V999"])


# ===========================================================================
# D7 -- matching unscoped is the DEFAULT and says nothing about it
# ===========================================================================
#
# `04_map_to_oligos.py --sps-mapping` defaults to None, and the else-branch just prints
# "Matching donors ... via core.matching.match_dataframe". Measured consequence of matching every
# row against every subpool at once: 313,570 of 4,458,480 named rows with a readable SPS (7.0%)
# took a name from a subpool their own read did not come from.
#
# The run-level state IS already recorded -- `prov.set_parameters(sps_scoping=bool(...))` -- so
# this is NOT fixed by adding a column (a per-row `sps_scoped=False` on every row of an unscoped
# run is a constant column, which names nothing). What was missing is a loud, quantified warning,
# and in particular the case where the input CARRIES an SPS column and the mapping was simply not
# supplied: there the information needed to scope was present and discarded.

STEP04 = _SRC / "snakemake/04_map_to_oligos.py"


def test_d7_unscoped_run_is_reported_loudly_and_quantified():
    step04 = _load(STEP04, "bdb_step04_d7")
    msg, stats = step04.unscoped_matching_report(n_rows=1000, has_sps_column=False)
    assert stats["sps_scoped"] is False
    assert stats["rows_matched_unscoped"] == 1000
    assert "7.0" in msg, "the measured cross-subpool mis-naming rate must be in the message"
    assert "WARNING" in msg.upper()


def test_d7_an_unused_sps_column_is_the_loudest_case():
    """The input could have been scoped and was not: that is a different, worse condition."""
    step04 = _load(STEP04, "bdb_step04_d7b")
    msg, stats = step04.unscoped_matching_report(n_rows=1000, has_sps_column=True)
    assert stats["sps_column_present_but_unused"] is True
    assert "--sps-mapping" in msg


def test_d7_scoped_run_reports_no_warning():
    """The guard must be silent when the run IS scoped, or it becomes noise to be ignored."""
    step04 = _load(STEP04, "bdb_step04_d7c")
    msg, stats = step04.unscoped_matching_report(n_rows=1000, has_sps_column=True, scoped=True)
    assert stats["sps_scoped"] is True
    assert stats["rows_matched_unscoped"] == 0
    assert msg == ""


def test_d7_require_sps_mapping_flag_exists_and_is_read():
    """A flag nothing reads is not a fix. --require-sps-mapping must turn the warning into an
    error, and the parser must expose it."""
    step04 = _load(STEP04, "bdb_step04_d7d")
    p = step04.build_parser()
    base = ["--input", "i.tsv", "--output", "o.tsv"]
    assert p.parse_args(base + ["--require-sps-mapping"]).require_sps_mapping is True
    assert p.parse_args(base).require_sps_mapping is False
    with pytest.raises(SystemExit):
        step04.enforce_sps_mapping_requirement(require=True, sps_mapping=None)
    # and it does NOT fire when a mapping was supplied
    step04.enforce_sps_mapping_requirement(require=True, sps_mapping="some/path.tsv")


# ===========================================================================
# D16 -- one trailing space destroys perfect-donor matching, and it reads as failed synthesis
# ===========================================================================
#
# REPRODUCED 2026-10-01. `match_sequence` upper-cases but never strips, and the vectorized fast
# path does `.str.upper()` only. Measured on the two-design fixture:
#
#   clean            perfect_donor=True   method=perfect
#   trailing space   perfect_donor=False  method=sliding_window     <-- the perfect route is GONE
#   trailing newline perfect_donor=False  method=sliding_window
#   leading space    perfect_donor=True   method=perfect            <-- survives the anchor trim
#
# So the defect is specifically TRAILING whitespace, and the anchor trim hides a leading one. The
# NN 241 attack measured the consequence on real data: perfect-donor matching 100% -> 0% and
# accuracy 68.6% -> 40.9%, which presents as a failed synthesis batch rather than as a parsing bug.


def test_d16_trailing_whitespace_used_to_destroy_the_perfect_route(tables):
    """Pins the mechanism: the stripped and unstripped donors must now match identically."""
    clean = _read(DONOR_A)
    for dirty in (clean + " ", " " + clean, clean + "\n", clean + "\t", "  " + clean + "  \r\n"):
        res = M.match_sequence(dirty, None, tables)
        assert res.perfect_donor is True, f"perfect route lost to whitespace: {dirty[-4:]!r}"
        assert res.match_method == "perfect"
        assert res.oligo_name == "design_A"


def test_d16_guide_whitespace_does_not_break_guide_matching(tables):
    """A padded guide must not read as an undesigned guide -- that would silently flip
    guide_status to guide_not_designed and suppress the D12 conflict guard."""
    res = M.match_sequence(_read(DONOR_A), " " + GUIDE_A + " ", tables)
    assert res.guide_status == M.GUIDE_MATCH
    res2 = M.match_sequence(_read(DONOR_A), GUIDE_B + "\n", tables)
    assert res2.guide_status == M.GUIDE_CONFLICT
    assert res2.is_unambiguous is False


def test_d16_whitespace_is_stripped_through_match_dataframe(tables):
    """Through the production entry point, including the vectorized fast path."""
    df = pd.DataFrame({"donor": [_read(DONOR_A), _read(DONOR_A) + " ", " " + _read(DONOR_A) + "\n"]})
    out = M.match_dataframe(df, tables)
    assert out["perfect_donor"].tolist() == [True, True, True]
    assert out["oligo_name"].tolist() == ["design_A"] * 3
    assert out["match_method"].tolist() == ["perfect"] * 3


def test_d16_donor_length_assertion_catches_a_systematic_shift(tables):
    """THE SECOND HALF of D16: a donor-length check, so a systematic length problem is REPORTED
    rather than presenting as a biological failure.

    A whitespace or truncation problem shifts the observed trimmed-donor lengths off the designed
    set. The report must say so and count it.
    """
    designed = {len(DONOR_A), len(DONOR_B)}
    ok = M.donor_length_report([129, 129, 129], designed)
    assert ok["off_design_lengths"] == 0
    assert ok["designed_lengths"] == sorted(designed)

    bad = M.donor_length_report([130, 130, 129], designed)
    assert bad["off_design_lengths"] == 2
    assert bad["off_design_fraction"] == pytest.approx(2 / 3)
    assert bad["length_histogram"][130] == 2


def test_d16_donor_length_report_is_wired_into_match_dataframe(tables, caplog):
    """A report nothing calls is not a fix. An all-off-design run must WARN."""
    import logging

    # every donor one base long -- nothing can match, and the lengths are off-design
    df = pd.DataFrame({"donor": [_read(DONOR_A) + "GGG"] * 3})
    with caplog.at_level(logging.WARNING, logger=M.logger.name):
        M.match_dataframe(df, tables)
    assert any("off-design" in r.message.lower() or "donor length" in r.message.lower()
               for r in caplog.records), [r.message for r in caplog.records]


def test_d11_the_SECOND_loader_in_step05_has_the_same_defects(tmp_path):
    """🔴 There are TWO `load_sps_mapping`. Fixing `core.matching`'s left the ASSIGNMENT path --
    `05_assign_sublibrary.py:99`, used by `assign_v_library_by_sps` -- with both original defects:
    it demands `library_id` (so the shipped `V_library` file raises) and it collapses duplicate
    SPS last-wins via `dict(zip(...))`. A fix to one of two definitions is not a fix.
    """
    step05 = _load(STEP05, "bdb_step05_d11")
    # (a) the shipped column name must be accepted
    m = step05.load_sps_mapping(_sps_file(tmp_path, "V_library", [
        ("CGTATGGAGATGACGTCAGG", "V546")]))
    assert m == {"CGTATGGAGATGACGTCAGG": "V546"}
    # (b) a duplicate SPS must NOT be resolved by file order
    with pytest.raises(ValueError, match="shared by several libraries|uniquely names"):
        step05.load_sps_mapping(_sps_file(tmp_path, "V_library"))


def test_d11_step05_still_reads_the_yaml_form(tmp_path):
    """05's loader also accepts a legacy YAML mapping; delegating the TSV path must not lose it."""
    import yaml

    step05 = _load(STEP05, "bdb_step05_d11_yaml")
    p = tmp_path / "m.yaml"
    p.write_text(yaml.safe_dump({"libraries": {"V546": {"sps": "cgtatggagatgacgtcagg"}}}))
    assert step05.load_sps_mapping(str(p)) == {"CGTATGGAGATGACGTCAGG": "V546"}


def test_c19_weighting_is_opt_in_and_the_default_is_unchanged():
    """🔴 THE HAZARD THE OPT-IN CLOSES. bdb stage 01 reads PAIRED R1/R2, which a collapsed FASTQ
    is not, and `bare_int` matches ANY all-digit read name -- an integer-renamed FASTQ
    (`@1, @2, @3, ...`) would be weighted by its RECORD INDEX. So acting on the header is opt-in;
    detection is always on."""
    step01 = _load(STEP01, "bdb_step01_c19_optin")
    p = step01.parse_args.__globals__["argparse"]
    assert "collapsed_counts" in step01.process_files.__code__.co_varnames
    # detection is unconditional; the ACTION is gated
    assert step01.read_weight_from_header("@seq_1_count=208") == (208, "count_kv")


def test_c19_integer_renamed_fastq_is_detected_as_a_rename():
    """A rename is strictly increasing and all-distinct; collapse counts repeat and are
    abundance-distributed. Only ever used to REFUSE, never to reinterpret."""
    step01 = _load(STEP01, "bdb_step01_c19_rename")
    assert step01.looks_like_an_integer_renamed_fastq(list(range(1, 101))) is True
    # a real collapse-count series repeats values
    assert step01.looks_like_an_integer_renamed_fastq(
        [208, 12, 5, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]) is False
    # too short to judge -> do not accuse
    assert step01.looks_like_an_integer_renamed_fastq([1, 2, 3]) is False


# ===========================================================================
# D15 -- the variant id is neither chromosome-qualified nor necessarily positional
# ===========================================================================
#
# MEASURED on the real shipped harmonized table (146,250 rows):
#   101,681 positional rows (`<pos>_<ref>_<alt>_SNP`), 22,252 distinct ids,
#       of which 50 APPEAR ON MORE THAN ONE CHROMOSOME  -> two different edits merge into one
#   26,101 NON-positional rows across 3,874 class labels -> `PTC` alone covers 3,054 rows
#   18,468 blank
#
# `individual_variants` is the variant identity used for unambiguity. Both failure modes make an
# AMBIGUOUS match report `is_unambiguous=True`:
#   * a bare position collides across chromosomes (50 ids measured);
#   * a CLASS LABEL is not an identity at all -- treating `PTC` as one collapses 3,054 genuinely
#     different premature-stop designs into a single "variant", so a match among any of them looks
#     unambiguous.
#
# `variant_coord` (col 28) is populated on 22,434 of the 26,101 class-label rows, so most of them
# have a real position available; and there IS a `chrom` column (col 12) to qualify with. The
# remainder fall through to the donor-sequence fallback that `build_lookup_tables` already has.


@pytest.mark.parametrize(
    "chrom,iv,vcoord,expect",
    [
        # positional AND chromosome-qualified. 🔴 D15 (Kevin, 2026-10-02): `variant_coord` is NOT
        # in the id. Keying on it split one real edit -- one donor, ~5 redundant guides whose
        # coords span up to 2,383 bp -- into several ids and manufactured false ambiguity.
        ("chrIV", "516707_T_C_SNP", "516707", "chrIV:516707_T_C_SNP"),
        # a class label is NOT rescued by a coordinate; it goes to the donor fallback. A
        # coordinate does not turn `PTC` into an identity.
        ("chrVII", "PTC", "392020", None),
        # positional id, no variant_coord -> still qualify with the chromosome
        ("chrII", "559546_G_A_SNP", None, "chrII:559546_G_A_SNP"),
        # a class label with NO coordinate names nothing positional -> MISSING
        ("chrIV", "PTC", None, None),
        ("chrIV", "synonymous_ctl", "", None),
        ("chrIV", "PDR5_aa_mutation_controls", None, None),
        # no chromosome -> not chromosome-qualified -> MISSING
        (None, "516707_T_C_SNP", "516707", None),
        ("", "516707_T_C_SNP", "516707", None),
        # blank
        ("chrIV", None, None, None),
        ("chrIV", "", "", None),
    ],
)
def test_d15_variant_id_must_be_chromosome_qualified_and_positional(chrom, iv, vcoord, expect):
    assert M.qualified_variant_id(chrom, iv, vcoord) == expect


def test_d15_the_cross_chromosome_collision_is_closed():
    """The measured defect: one positional id on two chromosomes must be TWO variants."""
    a = M.qualified_variant_id("chrIV", "68185_C_T_SNP", "68185")
    b = M.qualified_variant_id("chrVII", "68185_C_T_SNP", "68185")
    assert a is not None and b is not None and a != b


def test_d15_a_class_label_no_longer_collapses_thousands_of_edits(tmp_path):
    """`PTC` as a variant identity made every premature-stop design one variant.

    Two designs that install DIFFERENT edits but share the class label `PTC` must not be reported
    as one unambiguous variant when a donor matches both.
    """
    # 🔴 FIXTURE CORRECTED 2026-10-02 (NN 253 §S2). This previously gave both designs the SAME
    # donor, which does not occur: in the real harmonized table ZERO shared donors span two
    # chromosomes. With a shared donor the two designs are genuinely indistinguishable from the
    # read alone, so the old fixture was demanding that the matcher separate what the data cannot.
    # Under Kevin's D15 ruling a class label is routed to the DONOR fallback, so the donor is what
    # must distinguish them -- and in real data it does.
    donor_a = ("ACGT" * 33)[:129]
    donor_b = ("TTGCA" * 27)[:129]
    assert donor_a != donor_b
    oligo_df = pd.DataFrame(
        {
            "oligo_name": ["ptc_gene_A", "ptc_gene_B"],
            "guide": [GUIDE_A, GUIDE_B],
            "donor": [donor_a, donor_b],
            "chrom": ["chrIV", "chrVII"],
            "variant_coord": ["1000", "2000"],
            M.VARIANT_COLUMN: ["PTC", "PTC"],
        }
    )
    tables = M.build_lookup_tables(oligo_df)
    res_a = M.match_sequence(_read(donor_a), None, tables)
    res_b = M.match_sequence(_read(donor_b), None, tables)
    # the label `PTC` is shared, so if it were the identity these would be ONE variant
    assert res_a.num_candidate_oligos == 1 and res_b.num_candidate_oligos == 1
    assert res_a.oligo_name != res_b.oligo_name, (
        "two different PTC designs collapsed to one oligo -- the class label was used as an "
        "identity instead of falling back to the donor"
    )
    assert res_a.is_unambiguous is True and res_b.is_unambiguous is True


def test_d15_identical_edits_still_collapse_to_one_variant(tmp_path):
    """The DESIGNED behaviour must survive: redundant guides installing the SAME edit are one
    variant, so the match stays unambiguous. Without this the fix would destroy the feature
    variant-level unambiguity exists for."""
    shared_donor = ("ACGT" * 33)[:129]
    oligo_df = pd.DataFrame(
        {
            "oligo_name": ["guideA_same_edit", "guideB_same_edit"],
            "guide": [GUIDE_A, GUIDE_B],
            "donor": [shared_donor, shared_donor],
            "chrom": ["chrIV", "chrIV"],
            "variant_coord": ["1000", "1000"],
            M.VARIANT_COLUMN: ["1000_T_C_SNP", "1000_T_C_SNP"],
        }
    )
    tables = M.build_lookup_tables(oligo_df)
    res = M.match_sequence(_read(shared_donor), None, tables)
    assert res.num_candidate_oligos == 2
    assert res.num_candidate_variants == 1
    assert res.is_unambiguous is True


def test_d15_non_positional_rows_fall_back_to_the_donor_not_to_the_label():
    """A row whose identity cannot be made positional uses the DONOR sequence, the fallback the
    builder already had -- it must not keep the class label as an identity."""
    d1 = ("ACGT" * 33)[:129]
    d2 = ("TTGCA" * 26)[:129]
    oligo_df = pd.DataFrame(
        {
            "oligo_name": ["ctl_1", "ctl_2"],
            "guide": [GUIDE_A, GUIDE_B],
            "donor": [d1, d2],
            "chrom": ["chrIV", "chrIV"],
            "variant_coord": [None, None],
            M.VARIANT_COLUMN: ["PDR5_aa_mutation_controls", "PDR5_aa_mutation_controls"],
        }
    )
    tables = M.build_lookup_tables(oligo_df)
    assert tables.oligo_to_variant["ctl_1"] == d1
    assert tables.oligo_to_variant["ctl_2"] == d2
    assert tables.oligo_to_variant["ctl_1"] != tables.oligo_to_variant["ctl_2"]


# ===========================================================================
# C20 -- duplicate oligo names: attributes silently taken from whichever row came first
# ===========================================================================
#
# MEASURED on the real shipped harmonized table: 146,250 rows, 145,637 distinct `oligo_name`,
# **613 duplicated names, of which 0 conflict on `guide`** and 0 span multiple
# `guide_nuclease_PAM`. So the right behaviour is to DEDUPE benign duplicates with a count and
# REFUSE only conflicting ones: a blanket raise would break today's run on 613 harmless rows,
# while `drop_duplicates(subset=["oligo_name"])` -- what the code did -- is first-wins by file
# order and would silently pick one of two different guides if a conflict ever appeared.


def test_c20_benign_duplicates_are_deduped_and_counted():
    """613 of these exist in the real table today, all with identical attributes."""
    step04 = _load(STEP04, "bdb_step04_c20")
    h = pd.DataFrame({
        "oligo_name": ["o1", "o1", "o2"],
        "guide": [GUIDE_A, GUIDE_A, GUIDE_B],
        "lib": ["V545", "V545", "V546"],
    })
    m, stats = step04.resolve_oligo_attr(h, "oligo_name", ["guide", "lib"])
    assert m["guide"] == {"o1": GUIDE_A, "o2": GUIDE_B}
    assert stats["duplicated_names"] == 1
    assert stats["benign_duplicates"] == 1
    assert stats["conflicting_names"] == 0


def test_c20_a_conflicting_duplicate_is_refused_not_guessed():
    """Two different guides under one oligo_name: first-wins is a silent wrong attribution."""
    step04 = _load(STEP04, "bdb_step04_c20b")
    h = pd.DataFrame({
        "oligo_name": ["o1", "o1"],
        "guide": [GUIDE_A, GUIDE_B],
        "lib": ["V545", "V545"],
    })
    with pytest.raises(ValueError, match="conflict"):
        step04.resolve_oligo_attr(h, "oligo_name", ["guide", "lib"])


def test_c20_conflict_is_checked_on_EVERY_attribute_not_just_guide():
    """`oligo_pool` comes from the same map and the table merges the 2021 and 2024 pools, so a
    library conflict matters as much as a guide conflict."""
    step04 = _load(STEP04, "bdb_step04_c20c")
    h = pd.DataFrame({
        "oligo_name": ["o1", "o1"],
        "guide": [GUIDE_A, GUIDE_A],
        "lib": ["V545", "V629"],          # same guide, DIFFERENT library
    })
    with pytest.raises(ValueError, match="conflict"):
        step04.resolve_oligo_attr(h, "oligo_name", ["guide", "lib"])


def test_c20_lookup_table_dicts_are_not_silently_last_wins():
    """`guide_donor_to_oligo` is a dict(zip(...)), i.e. last-wins. A collision there must be
    counted and reported in the build stats, not resolved by row order."""
    shared_guide, shared_donor = GUIDE_A, ("ACGT" * 33)[:129]
    oligo_df = pd.DataFrame({
        "oligo_name": ["first_in_file", "second_in_file"],
        "guide": [shared_guide, shared_guide],
        "donor": [shared_donor, shared_donor],
        "chrom": ["chrIV", "chrIV"],
        "variant_coord": ["1000", "1000"],
        M.VARIANT_COLUMN: ["1000_T_C_SNP", "1000_T_C_SNP"],
    })
    tables = M.build_lookup_tables(oligo_df)
    stats = tables.get_stats()
    assert stats["guide_donor_key_collisions"] == 1, (
        "two oligos share one guide+donor key and the dict kept whichever came last, uncounted"
    )


# ===========================================================================
# D13 -- PAM-tuple agreement, carried as a FLAG (see 250_*.md §9 for why not a veto)
# ===========================================================================


def _pam_oligo_df():
    """Two oligos, same edit, DIFFERENT PAM tuple -- the designed redundant-guide case, which is
    97.5% of multi-oligo donor groups in the real library."""
    shared_donor = ("ACGT" * 33)[:129]
    return pd.DataFrame({
        "oligo_name": ["pam_AGTA", "pam_AGGA"],
        "guide": [GUIDE_A, GUIDE_B],
        "donor": [shared_donor, shared_donor],
        "chrom": ["chrI", "chrI"],
        "variant_coord": ["192337", "192337"],
        M.VARIANT_COLUMN: ["192337_T_TAA_SNP", "192337_T_TAA_SNP"],
        "PAM": ["AGTA", "AGGA"],
        "PAM_coord": ["192337", "192344"],
        "PAM_strand": ["-", "+"],
    }), shared_donor


def test_d13_pam_disagreement_is_flagged_but_does_not_veto_unambiguity():
    """🔶 THE DEVIATION, tested as implemented. The edit is unambiguous; the PAM attribution is not.

    Vetoing here would strip is_unambiguous from 97.5% of multi-oligo donor groups, because the
    library is DESIGNED with redundant guides at different PAM sites installing the same edit.
    """
    oligo_df, shared_donor = _pam_oligo_df()
    tables = M.build_lookup_tables(oligo_df)
    res = M.match_sequence(_read(shared_donor), None, tables)
    assert res.num_candidate_oligos == 2
    assert res.num_candidate_variants == 1
    assert res.is_unambiguous is True, "the EDIT is unambiguous and must stay so"
    assert res.candidates_agree_on_pam is False, "the PAM attribution is NOT determined"
    assert res.num_candidate_pams == 2


def test_d13_agreeing_pams_are_flagged_true():
    oligo_df, shared_donor = _pam_oligo_df()
    oligo_df.loc[:, ["PAM", "PAM_coord", "PAM_strand"]] = ["AGTA", "192337", "-"]
    tables = M.build_lookup_tables(oligo_df)
    res = M.match_sequence(_read(shared_donor), None, tables)
    assert res.candidates_agree_on_pam is True
    assert res.num_candidate_pams == 1


def test_d13_pam_flag_is_None_when_the_table_has_no_pam_columns(tables):
    """Unknown is a THIRD state and must not be reported as agreement."""
    res = M.match_sequence(_read(DONOR_A), None, tables)
    assert res.candidates_agree_on_pam is None
    assert res.num_candidate_pams == 0


def test_d13_pam_flag_reaches_the_dataframe():
    """A flag nothing can read is not a fix: it must be a column on the output."""
    oligo_df, shared_donor = _pam_oligo_df()
    tables = M.build_lookup_tables(oligo_df)
    out = M.match_dataframe(pd.DataFrame({"donor": [_read(shared_donor)]}), tables)
    assert "candidates_agree_on_pam" in out.columns
    assert bool(out.loc[0, "candidates_agree_on_pam"]) is False
    assert int(out.loc[0, "num_candidate_pams"]) == 2


# ===========================================================================
# D14 -- a run-level guide_not_designed fraction assertion
# ===========================================================================
#
# The NN 241 cross-library-contamination attack absorbed at 46.7% assigned with 47 rows flagged
# unambiguous, and nothing in the run said the guides did not belong to the library being matched
# against. `guide_not_designed` is already computed per row; what was missing is a run-level
# check on its FRACTION.


def test_d14_a_high_guide_not_designed_fraction_is_reported(tables, caplog):
    import logging

    # every guide is a real 20-mer but belongs to no design in this lookup
    df = pd.DataFrame({
        "donor": [_read(DONOR_A)] * 4,
        "guide": ["GGGGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCCCC",
                  "AAAAAAAAAAAAAAAAAAAA", GUIDE_A],
    })
    with caplog.at_level(logging.WARNING, logger=M.logger.name):
        out = M.match_dataframe(df, tables)
    assert (out["guide_status"] == M.GUIDE_NOT_DESIGNED).sum() == 3
    assert any("guide_not_designed" in r.message for r in caplog.records), \
        [r.message for r in caplog.records]


def test_d14_report_function_quantifies_and_has_a_documented_threshold():
    rep = M.guide_provenance_report({M.GUIDE_NOT_DESIGNED: 75, M.GUIDE_MATCH: 25})
    assert rep["guide_not_designed_fraction"] == 0.75
    assert rep["exceeds_threshold"] is True
    assert rep["threshold"] == M.GUIDE_NOT_DESIGNED_WARN_FRACTION
    ok = M.guide_provenance_report({M.GUIDE_NOT_DESIGNED: 1, M.GUIDE_MATCH: 99})
    assert ok["exceeds_threshold"] is False
    # rows with no guide at all are excluded from the denominator -- they are not evidence
    noguide = M.guide_provenance_report({M.GUIDE_NO_GUIDE: 1000, M.GUIDE_NOT_DESIGNED: 1})
    assert noguide["guide_bearing_rows"] == 1


# ===========================================================================
# D9 -- link_collapse.py:571 `sort_values(rank_col)` with no kind="mergesort" (SUSPECTED)
# ===========================================================================
#
# `merge_orphan_bc1` builds each bc1's per-link list "most abundant first" with
# `sub.sort_values(rank_col, ascending=False)`. numpy's default quicksort is NOT stable, so links
# TIED on rank_col can be ordered differently run to run -- and the pair-acceptance loop takes the
# FIRST compatible pair it finds, so a different order could in principle accept a different link.
#
# This was reported as SUSPECTED and NOT reproduced, and NN 240 found no depth-dependent
# mis-assignment anywhere (100% recovery at 1, 2, 5, 20 and 100 reads). It is tested here rather
# than patched on suspicion.
#
# 🔴 THE TRAP IN TESTING IT: numpy's quicksort falls back to INSERTION SORT on small partitions
# (<= 16 elements), which IS stable. A small fixture therefore looks deterministic whether or not
# the sort is stable, and would produce a false "NOT REPRODUCED". These fixtures use >16 tied
# links inside one bc1 so the real quicksort path is exercised.

from magestic.pipelines.bc1_donor_bc0.core import link_collapse as LC  # noqa: E402

_LC_KW = dict(rank_col="counts", sum_cols=["counts"], or_cols=[], preserve_cols=[])


def _orphan_fixture(n_tied=40, seed=0):
    """Two bc1 at Hamming 1, each with n_tied links ALL TIED on counts.

    Every link of the orphan is bc0-compatible with every link of the root, so which pair the
    acceptance loop finds FIRST is decided purely by the per-link ordering -- the thing the
    missing `kind="mergesort"` would make non-deterministic.
    """
    root = "A" * 20
    orphan = "A" * 19 + "C"          # Hamming 1 from root
    rows = []
    for i in range(n_tied):
        d = f"{'ACGT' * 10}{i:03d}"
        rows.append({"bc1": root, "donor": d, "bc0": f"BC{i:08d}", "counts": 7})
        rows.append({"bc1": orphan, "donor": d, "bc0": f"BC{i:08d}", "counts": 7})
    # one extra read so the root wins the component-root tie deterministically
    rows.append({"bc1": root, "donor": "ACGT" * 10 + "999", "bc0": "BC99999999", "counts": 9})
    rng = __import__("random").Random(seed)
    rng.shuffle(rows)
    return pd.DataFrame(rows)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 7, 13, 42, 99, 1234, 31337])
def test_d9_orphan_merge_is_permutation_invariant_over_many_tied_links(seed):
    """The property that matters. If this ever fails, D9 is real and the fix is a TOTAL ordering
    (secondary keys), not merely `kind="mergesort"` -- a stable sort only inherits input order."""
    base = _orphan_fixture(n_tied=40, seed=0)
    shuffled = _orphan_fixture(n_tied=40, seed=seed)

    a_df, a_stats, a_map = LC.merge_orphan_bc1(base.copy(), **_LC_KW)
    b_df, b_stats, b_map = LC.merge_orphan_bc1(shuffled.copy(), **_LC_KW)

    for k in ("accepted_pairs", "components_merged", "bc1_removed", "reads_moved",
              "rows_snapped_to_root_link", "bc1_out", "reads_out"):
        assert a_stats[k] == b_stats[k], f"{k} depends on input row order (seed={seed})"
    assert (dict(zip(a_map["bc1"], a_map["bc1_root"]))
            == dict(zip(b_map["bc1"], b_map["bc1_root"]))), "bc1->root map is order-dependent"


def test_d9_repeated_runs_on_one_input_are_identical():
    """Same input, many runs: an unstable sort can reorder equal keys between calls."""
    df = _orphan_fixture(n_tied=40, seed=5)
    outs = []
    for _ in range(8):
        _d, st, mp = LC.merge_orphan_bc1(df.copy(), **_LC_KW)
        outs.append((st["accepted_pairs"], st["reads_moved"], st["bc1_out"],
                     tuple(sorted(zip(mp["bc1"], mp["bc1_root"])))))
    assert len(set(outs)) == 1, "merge_orphan_bc1 is not deterministic across repeated runs"


def test_d9_no_reads_are_lost_whatever_the_order():
    """Whatever the tie-break does, the read accounting must balance -- no silent drop."""
    for seed in (0, 3, 42):
        df = _orphan_fixture(n_tied=40, seed=seed)
        _d, st, _m = LC.merge_orphan_bc1(df.copy(), **_LC_KW)
        assert st["reads_in"] == st["reads_out"], st


# ===========================================================================
# D10 -- unjustified, unexposed thresholds
# ===========================================================================
#
# Four magic numbers decided matching behaviour with no way to change or inspect them:
#   DEFAULT_WINDOW_SIZE = 60          sliding-window size
#   DEFAULT_ANCHOR_POS_SLACK = 10     how far either side of the canonical anchor position to scan
#   DEFAULT_ANCHOR_MAX_MISMATCH = 1   substitutions tolerated in the 11-bp anchor
#   BC1_LENGTH_TOLERANCE = 1          (stage 01; see its own note)
# plus the 30/10 fragment-key composition, now named in config (DONOR_FRAGMENT_LENGTH/BC0_LENGTH).
#
# Exposed on the CLI with the CURRENT values as defaults, so nothing changes unless asked. Each
# test passes a NON-default value and asserts BEHAVIOUR changes -- asserting that argparse
# accepts a flag proves only that argparse works.


def test_d10_matching_thresholds_are_exposed_with_current_values_as_defaults():
    step04 = _load(STEP04, "bdb_step04_d10")
    a = step04.build_parser().parse_args(["--input", "i.tsv", "--output", "o.tsv"])
    assert a.window_size == M.DEFAULT_WINDOW_SIZE == 60
    assert a.anchor_pos_slack == M.DEFAULT_ANCHOR_POS_SLACK == 10
    assert a.anchor_max_mismatch == M.DEFAULT_ANCHOR_MAX_MISMATCH == 1


def test_d10_window_size_actually_changes_matching():
    """A smaller window must admit a more truncated donor; a larger one must refuse it."""
    oligo_df = _oligo_df()
    short_obs = _read(DONOR_A[:40])          # 40 nt of variable region survives

    t60 = M.build_lookup_tables(oligo_df, window_size=60)
    assert M.match_sequence(short_obs, None, t60).oligo_name is None, \
        "a 40-nt donor cannot fill a 60-nt window"

    t30 = M.build_lookup_tables(oligo_df, window_size=30)
    assert M.match_sequence(short_obs, None, t30).oligo_name == "design_A"
    assert t30.window_size == 30 and t60.window_size == 60


def test_d10_anchor_slack_actually_changes_trimming():
    """With the anchor displaced beyond the slack, the trim must FAIL; widening it must rescue."""
    # anchor pushed 25 nt later than DEFAULT_ANCHOR_POS, and mutated by 1 so only the tolerant
    # (slack-scanning) pass can find it
    mutated = ANCHOR[:5] + ("A" if ANCHOR[5] != "A" else "C") + ANCHOR[6:]
    donor = "G" * (M.DEFAULT_ANCHOR_POS + 25) + mutated + DONOR_A

    tight = M.trim_to_variable_region(donor, pos_slack=10)
    wide = M.trim_to_variable_region(donor, pos_slack=40)
    assert tight == donor, "slack 10 should not reach an anchor displaced by 25"
    assert wide == DONOR_A, "slack 40 should find it"


def test_d10_fragment_key_composition_is_named_not_magic():
    """The 30/10 in `donor[-30:] + SPS + bc0` are now named constants, and they sum to 60."""
    assert CFG.DONOR_FRAGMENT_LENGTH == 30
    assert CFG.BC0_LENGTH == 10
    assert CFG.SPS_LENGTH == 20
    assert (CFG.DONOR_FRAGMENT_LENGTH + CFG.SPS_LENGTH + CFG.BC0_LENGTH
            == CFG.DONOR_BC0_FRAGMENT_LENGTH == 60)


def test_d10_flags_reach_the_SCOPED_build_closure_too():
    """🔴 The trap: 04's scoped path rebuilds a lookup per subpool inside a closure that called
    `load_and_build_lookup` WITHOUT these arguments, so a flag would work unscoped and be
    silently ignored the moment --sps-mapping was given."""
    import inspect

    step04 = _load(STEP04, "bdb_step04_d10b")
    src = inspect.getsource(step04.main)
    closure = src[src.index("def _build("):]
    closure = closure[: closure.index("out = match_dataframe_scoped")]
    assert "window_size=" in closure, (
        "the scoped per-subpool lookup ignores --window-size, so scoping silently reverts it"
    )


# ===========================================================================
# NN249-1 -- non-N20 bc1 calls must be REJECTED, not padded
# ===========================================================================
#
# NN 249, measured over 768 wells: **94 of 649 bc1 calls are not 20-mers, and 39.4% of those have
# NO supporting read against 0.4% for N20 calls -- OR 144, p = 4e-32.** That is the
# strongest-evidenced defect in the whole bc1 track, and it is the same phenomenon as C18 seen one
# layer later: a barcode whose boundary is wrong is not a barcode.
#
# A reusable, COUNTING guard, because the bc1 string is the join spine: a silent length filter
# anywhere here is how a denominator changes without anyone noticing.


def test_nn249_exact_n20_mask_separates_the_three_states():
    s = pd.Series(["A" * 20, "A" * 18, "A" * 21, None, "", "A" * 19 + "N", "a" * 20])
    mask, stats = B.exact_n20_mask(s)
    assert mask.tolist() == [True, False, False, False, False, False, True]
    assert stats["n_rows"] == 7
    assert stats["n_exact_n20"] == 2
    assert stats["n_wrong_length"] == 3        # 18, 21, and the empty string
    assert stats["n_missing"] == 1
    assert stats["n_non_acgt"] == 1            # the N-containing 20-mer
    # nothing is unaccounted for: every row lands in exactly one bucket
    assert (stats["n_exact_n20"] + stats["n_wrong_length"] + stats["n_missing"]
            + stats["n_non_acgt"]) == stats["n_rows"]


def test_nn249_mask_is_case_insensitive_but_not_whitespace_blind():
    """Lowercase is a real barcode; a padded one is the D16 whitespace defect and must be caught
    as a length problem rather than silently stripped here."""
    mask, stats = B.exact_n20_mask(pd.Series(["a" * 20, " " + "A" * 20, "A" * 20 + "\n"]))
    assert mask.tolist() == [True, False, False]
    assert stats["n_whitespace"] == 2


def test_nn249_rejection_is_counted_and_the_rate_is_reported():
    """The measured NN 249 rate (94/649 = 14.5%) must be computable from the stats."""
    s = pd.Series(["A" * 20] * 555 + ["A" * 18] * 94)
    mask, stats = B.exact_n20_mask(s)
    assert stats["n_rows"] == 649 and stats["n_exact_n20"] == 555
    assert stats["non_n20_fraction"] == pytest.approx(94 / 649)
    assert stats["exceeds_threshold"] is True       # 14.5% >> the documented floor


def test_nn249_a_clean_library_does_not_trip_the_guard():
    mask, stats = B.exact_n20_mask(pd.Series(["A" * 20] * 100))
    assert mask.all() and stats["non_n20_fraction"] == 0.0
    assert stats["exceeds_threshold"] is False


# ===========================================================================
# D17 / NN249-4 -- a provenance field that can CARRY a call, and never a constant
# ===========================================================================
#
# NN 249 explains why no field could have named the template from sequence alone: the
# genome-integrated and plasmid bc1 cassettes are THE SAME CONSTRUCT in opposite orientations,
# sequence-identical for 80 bp left and 2,421 bp right of the N20, so with 2x150 almost every read
# pair at the locus is genuinely uninformative. Unambiguous genomic support needs a read reaching
# reference position <=10189 or >=12711 -- i.e. the discrimination is MATE GEOMETRY, computed by
# NN 249, not sequence.
#
# The existing `bc1_from_wgs_source` is a CONSTANT (`G` x615) and names nothing. So this is
# plumbing: accept NN 249's per-bc1 call, join it, and ASSERT IT VARIES. When no call is supplied,
# emit NO column -- an always-NaN column is just as empty a claim as a constant one.


def test_d17_provenance_join_requires_variation():
    """🔴 A constant column names nothing. If every row gets the same call, refuse."""
    df = pd.DataFrame({"bc1": ["A" * 20, "C" * 20, "G" * 20]})
    allsame = pd.DataFrame({"bc1": ["A" * 20, "C" * 20, "G" * 20],
                            "bc1_from_wgs_source": ["G", "G", "G"]})
    with pytest.raises(ValueError, match="constant|varies|does not vary"):
        M.attach_bc1_provenance(df, allsame)


def test_d17_a_varying_call_is_joined_and_counted():
    df = pd.DataFrame({"bc1": ["A" * 20, "C" * 20, "G" * 20, "T" * 20]})
    calls = pd.DataFrame({
        "bc1": ["A" * 20, "C" * 20, "G" * 20],
        "bc1_from_wgs_source": ["genomic", "plasmid", "uninformative"],
    })
    out, stats = M.attach_bc1_provenance(df, calls, require_variation=True, return_stats=True)
    got = out["bc1_from_wgs_source"].tolist()
    assert got[:3] == ["genomic", "plasmid", "uninformative"]
    assert pd.isna(got[3]), "an uncalled bc1 must be null, not invented"
    assert stats["n_called"] == 3
    assert stats["n_uncalled"] == 1            # counted, not silently blank
    assert stats["distinct_calls"] == 3
    assert stats["call_counts"]["genomic"] == 1


def test_d17_no_call_supplied_emits_no_column():
    """Absent evidence must not be dressed as a column of NaN."""
    df = pd.DataFrame({"bc1": ["A" * 20]})
    out, stats = M.attach_bc1_provenance(df, None, return_stats=True)
    assert "bc1_from_wgs_source" not in out.columns
    assert stats["supplied"] is False


# ===========================================================================
# NN249-3 -- barcodes recurring across wells must be flagged
# ===========================================================================
#
# The decisive evidence for "mapping error, not contamination" in NN 249 was that 3 of 3 recurring
# barcodes mis-map to the SAME wrong gene in both clones carrying them -- contamination cannot be
# reproducible. A recurrence flag surfaces that class directly.


def test_nn249_recurrence_across_wells_is_flagged_and_counted():
    df = pd.DataFrame({
        "bc1": ["A" * 20, "A" * 20, "C" * 20, "G" * 20, "G" * 20, "G" * 20],
        "well": ["w1", "w2", "w3", "w4", "w4", "w5"],
    })
    out, stats = M.flag_bc1_recurrence(df, bc1_col="bc1", well_col="well")
    # A appears in 2 distinct wells; G in 2 distinct wells (w4 twice is ONE well); C in 1
    assert out["bc1_n_wells"].tolist() == [2, 2, 1, 2, 2, 2]
    assert out["bc1_recurs_across_wells"].tolist() == [True, True, False, True, True, True]
    assert stats["recurring_bc1"] == 2
    assert stats["distinct_bc1"] == 3
    # the flag VARIES on this input -- a flag that cannot vary is not a flag
    assert out["bc1_recurs_across_wells"].nunique() == 2


@pytest.mark.parametrize("vc", [
    None, "", 706445, 706445.0,               # pandas dtype inference gives a FLOAT
    "706445.0",                               # and the file itself stores it that way
    "706445", " 706445 ",
    706431, "706447.0",                       # a DIFFERENT coordinate for the same edit
    "706445.5",                               # genuinely fractional
    "chrII_extra",                            # non-numeric junk
])
def test_d15_the_id_is_stable_whatever_variant_coord_says(vc):
    """🔴 D15, Kevin's ruling 2026-10-02: `variant_coord` is accepted and IGNORED.

    This supersedes an earlier test that asserted the coordinate was *normalised into* the id
    (`chrII:706445:deletion`). Normalising it only papered over the dtype problem -- the harmonized
    table stores `706445.0`, so the id came out `chrII:706445.0:deletion`, stable only while pandas
    inferred the same dtype. Removing the coordinate from the id solves that outright AND fixes the
    larger defect it was hiding: one donor with ~5 redundant guides has coordinates up to 2,383 bp
    apart, so coordinate-keying SPLIT one real edit into several ids and reported false ambiguity
    (1,492 production rows, 708 `perfect_ambiguous`, downgraded on a single plate).

    So every value below -- including two genuinely different coordinates -- must give ONE id.
    """
    assert M.qualified_variant_id("chrII", "706445_G_A_SNP", vc) == "chrII:706445_G_A_SNP"


def test_d15_a_class_label_is_never_rescued_by_a_coordinate():
    """The other half of the same ruling: a coordinate does not make a label an identity."""
    for vc in (None, "", "392020", 392020.0):
        assert M.qualified_variant_id("chrVII", "PTC", vc) is None


def test_d13_pam_columns_are_requested_from_the_real_oligo_table():
    """🔴 BLOCKER FOUND BY REVIEW, not by my own tests. `COLUMNS_BC1_DONOR_BC0` does not load
    PAM/PAM_coord/PAM_strand and neither did my optional list, so in PRODUCTION
    `tables.oligo_to_pam` is {} and `match_dataframe` wrote `candidates_agree_on_pam=None` and
    `num_candidate_pams=0` on EVERY row -- the banned constant column, and exactly the mistake I
    had already caught for D15. My D13 tests passed only because the fixture hand-supplies the
    columns."""
    from magestic.utils import oligo_annotations as OA

    requested = set(OA.get_annotation_columns_for_bc1_donor_bc0())
    for c in ("PAM", "PAM_coord", "PAM_strand"):
        assert c in requested, f"{c} is never loaded, so oligo_to_pam is always empty"
    assert {"PAM", "PAM_coord", "PAM_strand"} <= set(OA.OPTIONAL_COLUMNS_BC1_DONOR_BC0), (
        "must be OPTIONAL: required would make usecols raise on every table predating them"
    )


def test_d13_no_pam_data_means_NO_COLUMN_rather_than_a_constant_one(tables):
    """A column that is None on every row names nothing. Emit it only when there is PAM data."""
    out = M.match_dataframe(pd.DataFrame({"donor": [_read(DONOR_A), _read(DONOR_B)]}), tables)
    assert "candidates_agree_on_pam" not in out.columns
    assert "num_candidate_pams" not in out.columns


def test_d13_with_pam_data_the_columns_appear_and_vary():
    """And when there IS data, the column must be present AND actually vary."""
    oligo_df, shared_donor = _pam_oligo_df()
    tables = M.build_lookup_tables(oligo_df)
    df = pd.DataFrame({"donor": [_read(shared_donor), _read(DONOR_B)]})
    # DONOR_B is not in this lookup -> no candidates -> None; the shared donor -> False
    out = M.match_dataframe(df, tables)
    assert "candidates_agree_on_pam" in out.columns
    assert out["candidates_agree_on_pam"].nunique(dropna=False) >= 2, (
        "the flag does not vary on an input where it should"
    )


# ===========================================================================
# NN249-2 -- a minimum screen-side read-support floor on bc1 -> oligo
# ===========================================================================
#
# NN 249: the 21 mis-assigned wells' barcodes are 13-20x more likely to be RARE in the screen
# reference -- OR 13.12 [2.64, 65.13], p = 0.0019 against NE-other (which matches on outcome
# class, removing the CE-selection objection). The four "rare in screen" cases had screen-side
# counts of 1, 2, 2 and 7: "the oligo call rests on a handful of reads".
#
# Implemented as a counted LABEL with the floor exposed, not a filter: this codebase's standing
# rule is LABEL, never filter (02_combine_and_collapse), and dropping these rows would change a
# denominator silently.


def test_nn249_support_floor_labels_and_counts():
    df = pd.DataFrame({
        "bc1": ["A" * 20, "C" * 20, "G" * 20, "T" * 20],
        "screen_counts": [1, 7, 8, 28228],
    })
    out, stats = M.label_low_screen_support(df, count_col="screen_counts")
    # the floor captures the measured rare cases (1, 2, 2, 7) and leaves 8+ alone
    assert out["bc1_low_screen_support"].tolist() == [True, True, False, False]
    assert stats["floor"] == M.SCREEN_SUPPORT_FLOOR == 8
    assert stats["n_below_floor"] == 2
    assert stats["fraction_below_floor"] == 0.5
    # the label VARIES on this input -- a label that cannot vary is not a label
    assert out["bc1_low_screen_support"].nunique() == 2


def test_nn249_support_floor_is_configurable_and_counts_missing_separately():
    df = pd.DataFrame({"bc1": ["A" * 20, "C" * 20], "screen_counts": [5, None]})
    out, stats = M.label_low_screen_support(df, count_col="screen_counts", floor=3)
    assert stats["floor"] == 3
    assert bool(out["bc1_low_screen_support"].iloc[0]) is False   # 5 >= 3
    assert stats["n_missing_support"] == 1                        # counted, not called low
    assert pd.isna(out["bc1_low_screen_support"].iloc[1])


def test_nn249_support_floor_emits_no_column_without_a_count_column():
    """No support data means no claim -- not a column of False."""
    df = pd.DataFrame({"bc1": ["A" * 20]})
    out, stats = M.label_low_screen_support(df, count_col="screen_counts")
    assert "bc1_low_screen_support" not in out.columns
    assert stats["supplied"] is False
