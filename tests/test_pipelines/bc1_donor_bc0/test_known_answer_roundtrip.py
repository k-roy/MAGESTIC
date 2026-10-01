"""Known-answer round trip for the bc1_donor_bc0 pipeline (NN 240, 2026-10-01).

Evidence here is RECOVERY AGAINST PLANTED TRUTH, not inspection. The synthetic template is
anchored in real data -- see `synth_bc1_library` for the measurements.

Findings that are DEFECTS are marked `xfail(strict=True)`: the shared suite stays green and
each test flips to a failure the moment the bug is fixed, which is the signal to delete the
marker. Nothing here refactors the pipeline.

Durable write-up: SGTC `planning/222_triumph_embl_library_design/240_bc1_known_answer_roundtrip.md`
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

import synth_bc1_library as S  # noqa: E402
import roundtrip_harness as H  # noqa: E402

from magestic.pipelines.bc1_donor_bc0.readers.data_loader import (  # noqa: E402
    _GDB_SUPPORT_COLUMNS, collapse_gdb_bridge, load_gdb_data)
from magestic.pipelines.bc1_donor_bc0.core.matching import (  # noqa: E402
    trim_to_variable_region)


# =================================================================================
# The template itself -- these pin the facts measured against real data, so a future
# change to the parser or the matcher that breaks the convention fails here.
# =================================================================================

REAL_PARSED_DONOR = (
    "ATCCCCTATTTCATAGGGCACTATCCATCCTTACCCACTCCTCCGTGGGTTTCTTCCACTCAGGTGCGCACCC"
    "TTACGTCAGAAGAAACGGGTTTCCTGTTTGAAGAGCTTTTTTTAGCACATTACCATAAGTAGGAATACAGATC"
    "TGACCTTTCCACGGAATCACGAGCCCAAAGGTGCTGAATCGTGATCACGACCCGCGGGTGAGCCCACACTCAC"
    "CGGAGTCTTAGTAAATCTG")
REAL_DESIGN_DONOR = (
    "TTTTTTTAGCACATTACCATAAGTAGGAATACAGATCTGACCTTTCCACGGAATCACGAGCCCAAAGGTGCTG"
    "AATCGTGATCACGACCCGCGGGTGAGCCCACACTCACCGGAGTCTTAGTAAATCTG")
REAL_SPS = "CTTGACTCGACAGTGAGAGC"
REAL_BC0 = "CCTCAAAGGT"
REAL_FRAGMENT = "CCCACACTCACCGGAGTCTTAGTAAATCTGCTTGACTCGACAGTGAGAGCCCTCAAAGGT"


def test_real_parsed_donor_trims_to_the_design_donor():
    """The design-donor convention, measured from a real yL406 row rather than inferred."""
    assert len(REAL_PARSED_DONOR) == 238
    assert len(REAL_DESIGN_DONOR) == 129
    assert REAL_PARSED_DONOR.find("GTTTGAAGAGC") == 98   # == DEFAULT_ANCHOR_POS
    assert trim_to_variable_region(REAL_PARSED_DONOR) == REAL_DESIGN_DONOR
    assert trim_to_variable_region(REAL_DESIGN_DONOR) == REAL_DESIGN_DONOR  # idempotent


def test_real_fragment_is_donor30_plus_sps_plus_bc0():
    assert REAL_DESIGN_DONOR[-30:] + REAL_SPS + REAL_BC0 == REAL_FRAGMENT
    assert len(REAL_FRAGMENT) == 60


def test_synthetic_amplicon_parses_back_to_planted_truth():
    """Stage 01 on a synthetic read pair returns the planted bc1/bc0/SPS/donor/fragment."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "p01", H.SNK / "01_parse_bc1_donor_bc0.py")
    p01 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p01)

    d = S.make_design(n_loci=1, variants_per_locus=1, guides_per_variant=1)
    S.plant_bc1s(d)
    bc1 = next(iter(d.links))
    o = d.by_name(d.links[bc1][0])
    r1, r2 = S.cut_reads(S.build_amplicon(bc1, o))
    cfg = p01.ReadConfig(r1_length=300, r2_length=300, orientation="standard")
    got = p01.process_read_pair(r1, r2, cfg)

    assert got["anchor_tier"] == "confirmed"
    assert got["bc1"] == bc1
    assert got["bc0"] == o.bc0
    assert got["SPS"] == o.sps
    assert len(got["donor"]) == 238
    assert trim_to_variable_region(got["donor"]) == o.donor
    assert got["donor_bc0_fragment"] == o.fragment
    assert p01.detect_orientation([r1] * 50, [r2] * 50) == "standard"


# =================================================================================
# DEFECT 1 -- collapse_gdb_bridge is blind to the real gdb count columns
# =================================================================================

REAL_GDB_COUNT_COLUMNS = ("total_counts", "total_fwd_counts", "total_rev_counts")


def test_gdb_support_columns_cover_the_real_schema():
    """REGRESSION GUARD. Was a DEFECT reproduced at 01:0x on 2026-10-01 and FIXED the same
    session (`total_counts` / `counts` / `grouped_guide_bc0_counts` added). Before the fix the
    collapse found no support column, skipped its sort, and 'most-supported wins' degraded to
    'first row in the file wins' on every real gdb table. Do not shrink this tuple."""
    assert any(c in _GDB_SUPPORT_COLUMNS for c in REAL_GDB_COUNT_COLUMNS), (
        f"_GDB_SUPPORT_COLUMNS={_GDB_SUPPORT_COLUMNS} matches none of the real gdb count "
        f"columns {REAL_GDB_COUNT_COLUMNS}")


def test_collapse_keeps_best_supported_on_real_schema(tmp_path):
    """REGRESSION GUARD, the end-to-end half of the test above. Exercises the REAL
    `total_counts` schema that a guide_donor_bc0_counts.tsv actually carries, which the
    sibling's 9 unit tests do not (their fixture column is named `count`)."""
    d = S.make_design(n_loci=2, variants_per_locus=1, guides_per_variant=1)
    S.plant_bc1s(d)
    victim, wrong = d.oligos[0], d.oligos[1].oligo_name
    gdb = tmp_path / "gdb"
    # the wrong oligo is written FIRST with support 1; the truth follows with 10/11/12
    S.write_gdb_dir(d, gdb, reads_per_fragment=3, count_col="total_counts",
                    ambiguous={victim.fragment: wrong})
    raw, _ = load_gdb_data(gdb, collapse_bridge=False)
    collapsed, stats = collapse_gdb_bridge(raw)
    kept = collapsed[collapsed["donor_bc0_fragment"] == victim.fragment]
    assert stats["ambiguous_keys"] == 1          # the ambiguity IS counted
    assert str(kept.iloc[0]["oligo_name"]) == victim.oligo_name


def test_collapse_does_keep_best_supported_with_the_count_column(tmp_path):
    """The control for the test above: rename the column to `count` and it works. This is
    what makes the failure above a schema-coverage bug and not a logic bug."""
    d = S.make_design(n_loci=2, variants_per_locus=1, guides_per_variant=1)
    S.plant_bc1s(d)
    victim, wrong = d.oligos[0], d.oligos[1].oligo_name
    gdb = tmp_path / "gdb"
    S.write_gdb_dir(d, gdb, reads_per_fragment=3, count_col="count",
                    ambiguous={victim.fragment: wrong})
    raw, _ = load_gdb_data(gdb, collapse_bridge=False)
    collapsed, stats = collapse_gdb_bridge(raw)
    kept = collapsed[collapsed["donor_bc0_fragment"] == victim.fragment]
    assert stats["ambiguous_keys"] == 1
    assert str(kept.iloc[0]["oligo_name"]) == victim.oligo_name


# =================================================================================
# DEFECT 2 -- the gdb fragment fallback key is 40 nt, stage 01's key is 60 nt
# =================================================================================

# FIXED 2026-10-01 by NN 250 / C17 -- this was `xfail(strict=True)` while the defect stood, and
# it XPASSed the moment the fallback was repointed at config.build_donor_bc0_fragment. The marker
# is removed rather than kept, so the test is now a live regression guard instead of a record of a
# defect: the ONE builder lives in config.py and both drifted fallbacks (data_loader.py:585 and
# :1019) call it. Original reason, for the record: load_gdb_data's donor_bc0_fragment FALLBACK was
# donor[-30:] + bc0 = 40 nt while stage 01 writes donor[-30:] + SPS + bc0 = 60 nt; collapse_gdb_bridge
# and merge.py key on donor_bc0_fragment, so a gdb table without that column bridged ZERO bc1s and
# the run read as a biological null. Further guards: test_nn250_bc1_fixes.py -k c17.
def test_gdb_fragment_fallback_matches_stage01_key_length(tmp_path):
    d = S.make_design(n_loci=2, variants_per_locus=1, guides_per_variant=1)
    S.plant_bc1s(d)
    gdb = tmp_path / "gdb"
    S.write_gdb_dir(d, gdb, reads_per_fragment=1, omit_fragment_column=True)
    raw, _ = load_gdb_data(gdb, collapse_bridge=False)
    frags = raw["donor_bc0_fragment"].dropna().astype(str)
    assert len(frags) > 0
    assert set(frags.str.len()) == {60}, (
        f"fallback fragment lengths {sorted(set(frags.str.len()))}, stage 01 writes 60; "
        f"example fallback {frags.iloc[0]!r} vs truth {d.oligos[0].fragment!r}")


# =================================================================================
# Behaviour that PASSES -- pinned so a regression is caught
# =================================================================================

@pytest.mark.slow
def test_clean_input_recovers_every_planted_link(tmp_path):
    """Headline: recovery and precision are both 100% on noise-free input, and no stage
    loses a bc1 without saying so. 48 planted bc1s / 48 designed oligos / 2 replicates."""
    d = S.make_design(n_loci=8, variants_per_locus=3, guides_per_variant=2)
    S.plant_bc1s(d)
    mols = [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=5) for b, n in d.links.items()]
    r = H.run_route_a(d, mols, tmp_path / "clean", n_samples=2)

    for stage, (planted, spurious, _rows) in r.survival.items():
        assert planted == len(d.links), f"{stage} lost planted bc1s: {planted}/{len(d.links)}"
        assert spurious == 0, f"{stage} invented {spurious} bc1s"

    sc = H.score_table(r.stage_frames["04"], d, "04")
    assert sc.counts["missing"] == 0
    assert sc.counts["wrong_variant"] == 0
    assert sc.counts["exact_oligo"] == len(d.links), sc.line()
    assert sc.recovery == 1.0 and sc.precision_oligo == 1.0


@pytest.mark.slow
def test_guide_bridge_is_load_bearing(tmp_path):
    """Non-vacuity control. Each variant has two guides sharing one donor, so without the
    gdb bridge the donor alone is 2-way ambiguous and precision_oligo must drop to ~50%.
    If this ever reads 100%, the round trip above has stopped testing the bridge."""
    d = S.make_design(n_loci=4, variants_per_locus=2, guides_per_variant=2)
    S.plant_bc1s(d)
    mols = [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=5) for b, n in d.links.items()]
    r = H.run_route_a(d, mols, tmp_path / "bridge", n_samples=2, run_05=False)

    with_bridge = H.score_table(r.stage_frames["04"], d, "with_bridge")
    out = tmp_path / "bridge" / "mapped_nobridge.tsv"
    H.run([H.PY, H.SNK / "04_map_to_oligos.py", "--input", r.collapsed, "--output", out],
          tmp_path / "bridge" / "common", tmp_path / "bridge" / "log_nobridge.txt")
    without = H.score_table(pd.read_csv(out, sep="\t", low_memory=False), d, "no_bridge")

    assert with_bridge.precision_oligo == 1.0
    assert without.precision_variant == 1.0          # the EDIT is still right
    assert without.precision_oligo < 0.75, without.line()


@pytest.mark.slow
def test_bc1_collision_is_not_silently_collapsed(tmp_path):
    """A bc1 carrying two clones -- donors one nt apart, inside --donor-edit-distance 5 --
    must stay two rows with both oligos named and bc1_purity < 1, not merge into one."""
    d = S.make_design(n_loci=4, variants_per_locus=2, guides_per_variant=1)
    S.plant_bc1s(d)
    o1, o2 = d.oligos[0], d.oligos[1]
    assert sum(x != y for x, y in zip(o1.donor, o2.donor)) == 1
    shared = next(iter(d.links))
    own2 = [k for k, v in d.links.items() if v[0] == o2.oligo_name][0]
    del d.links[own2]
    d.links[shared] = [o1.oligo_name, o2.oligo_name]
    mols = [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=5) for b, n in d.links.items()]
    mols.append(S.Molecule(bc1=shared, oligo=o2, depth=5))
    r = H.run_route_a(d, mols, tmp_path / "bc1col", n_samples=2, run_05=False)

    rows = r.stage_frames["04"]
    rows = rows[rows["bc1"] == shared]
    assert len(rows) == 2, f"the two clones collapsed into {len(rows)} row(s)"
    assert set(rows["oligo_name"].astype(str)) == {o1.oligo_name, o2.oligo_name}
    assert (rows["num_distinct_donors"] == 2).all()
    assert (rows["bc1_purity"] < 1.0).all()


@pytest.mark.slow
def test_bc0_collision_is_flagged_not_silent(tmp_path):
    """Two oligos sharing a bc0 AND a donor[-30:] have identical 60-nt fragments, so 03's
    40-bp linking key collides and one of them receives the other's guide. The donor still
    names the oligo correctly, and the disagreement must be RECORDED -- guide_link_n_guides
    > 1 and guide_status == 'conflict' -- all the way through to 05."""
    d = S.make_design(n_loci=4, variants_per_locus=2, guides_per_variant=1)
    S.plant_bc1s(d)
    a, b = d.oligos[0], d.oligos[1]
    assert a.locus == b.locus and a.donor[-30:] == b.donor[-30:]
    b.bc0 = a.bc0
    assert a.fragment == b.fragment
    mols = [S.Molecule(bc1=k, oligo=d.by_name(v[0]), depth=5) for k, v in d.links.items()]
    r = H.run_route_a(d, mols, tmp_path / "bc0col", n_samples=2)

    f = r.stage_frames["04"]
    sc = H.score_table(f, d, "04")
    assert sc.counts["wrong_variant"] == 0 and sc.recovery == 1.0, sc.line()
    hit = f[f["oligo_name"].isin([a.oligo_name, b.oligo_name])]
    assert (hit["guide_link_n_guides"] > 1).all(), "the key collision was not recorded"
    assert "conflict" in set(hit["guide_status"].astype(str)), "no guide conflict flagged"
    assert "guide_link_n_guides" in r.stage_frames["05"].columns, "flag lost by 05"


@pytest.mark.slow
@pytest.mark.parametrize("depth", [1, 2, 5, 20])
def test_recovery_is_depth_independent_in_route_a(tmp_path, depth):
    """Route A assigns correctly down to ONE read per bc1 per replicate.
    --min-pcr-replicates is a LABEL (support_tier), never a filter: 02's own help says
    'no row is dropped'. Any depth dependence in real data comes from replicate count,
    run_pipeline's --min-reads, or a downstream support_tier filter -- not from here."""
    d = S.make_design(n_loci=4, variants_per_locus=2, guides_per_variant=1)
    S.plant_bc1s(d)
    mols = [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=depth)
            for b, n in d.links.items()]
    r = H.run_route_a(d, mols, tmp_path / f"d{depth}", n_samples=2, run_05=False)
    sc = H.score_table(r.stage_frames["04"], d, f"depth{depth}")
    assert sc.recovery == 1.0 and sc.precision_oligo == 1.0, sc.line()


@pytest.mark.slow
def test_a_synthesis_error_never_yields_a_confidently_wrong_edit(tmp_path):
    """The important precision property. With 1-3 substitutions in the donor the matcher
    drops out of the `perfect` tier into sliding_window / ambiguous / unnamed, but the
    measured `wrong_variant` count stays 0 at every error level -- it fails loudly, never
    silently onto a different edit."""
    for nmm in (1, 2, 3):
        d = S.make_design(n_loci=4, variants_per_locus=2, guides_per_variant=2)
        S.plant_bc1s(d)
        mols = [S.Molecule(bc1=b, oligo=d.by_name(n[0]), depth=5, synth_mismatches=nmm)
                for b, n in d.links.items()]
        r = H.run_route_a(d, mols, tmp_path / f"mm{nmm}", n_samples=2, run_05=False)
        sc = H.score_table(r.stage_frames["04"], d, f"mm{nmm}")
        assert sc.counts["wrong_variant"] == 0, sc.line()
        assert sc.counts["ambiguous_excl_truth"] == 0, sc.line()
        assert not r.stage_frames["04"]["perfect_donor"].any(), (
            "a mismatched donor must not be reported as a perfect match")


@pytest.mark.slow
def test_empty_well_writes_the_real_schema_and_does_not_crash(tmp_path):
    d = S.make_design(n_loci=1, variants_per_locus=1, guides_per_variant=1)
    S.plant_bc1s(d)
    common = tmp_path / "common"
    S.write_common_dir(d, common)
    (tmp_path / "e_R1.fastq").write_text("")
    (tmp_path / "e_R2.fastq").write_text("")
    out = tmp_path / "empty.tsv"
    p = H.run([H.PY, H.SNK / "01_parse_bc1_donor_bc0.py", tmp_path / "e_R1.fastq",
               tmp_path / "e_R2.fastq", out], common, tmp_path / "log.txt", check=False)
    assert p.returncode == 0
    lines = out.read_text().splitlines()
    assert len(lines) == 1, "an empty well must write a header-only table, not nothing"
    assert lines[0].split("\t")[:6] == ["counts", "bc1", "donor", "bc0", "SPS",
                                        "donor_bc0_fragment"]


@pytest.mark.slow
def test_bc1_sequencing_error_phantoms_are_separable_by_support(tmp_path):
    """The actionable low-support result (NN 240 section 9b).

    A 2-nt sequencing error in a bc1 survives 02's --bc1-edit-distance 1 as a NEW lineage and
    is then confidently named by 04 (`is_unambiguous` True). It is NOT, however, unmarked:
    every phantom lands in `support_tier == 'provisional'` with `counts == 1`, while every real
    link is `confident`, because a per-read error does not recur in an independent PCR
    replicate. Filtering on support removes 100% of them at no cost to recovery.
    """
    d = S.make_design(n_loci=4, variants_per_locus=2, guides_per_variant=1)
    S.plant_bc1s(d)
    mols = []
    for b, n in d.links.items():
        o = d.by_name(n[0])
        mols.append(S.Molecule(bc1=b, oligo=o, depth=20))
        mols.append(S.Molecule(bc1=b, oligo=o, depth=2, bc1_seq_errors=2))
    r = H.run_route_a(d, mols, tmp_path / "phantom", n_samples=2, run_05=False)

    f = r.stage_frames["04"]
    planted = set(d.links)
    real = f[f["bc1"].astype(str).isin(planted)]
    phantom = f[~f["bc1"].astype(str).isin(planted)]

    assert len(phantom) > 0, "the 2-error construction produced no phantom lineages"
    assert set(real["support_tier"].astype(str)) == {"confident"}
    assert set(phantom["support_tier"].astype(str)) == {"provisional"}, (
        "phantoms are no longer separable by support_tier -- the filtering rule in NN 240 "
        "section 9b would stop working")
    assert real["counts"].min() > phantom["counts"].max(), "count ranges now overlap"
    # and the planted links themselves are still perfect
    sc = H.score_table(f, d, "phantom")
    assert sc.counts["exact_oligo"] == len(d.links) and sc.counts["wrong_variant"] == 0
