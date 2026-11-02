"""Regression tests for the five GDB root causes arm G reproduced at master (planning/288 §3), fixed in
planning/378. Each test is arm G's counter-fixture, adapted to import the package directly."""
import random

import pandas as pd
import pytest

pytest.importorskip("pysam")

from magestic.pipelines.guide_donor_bc0.core import oligo_matching as om  # noqa: E402


def _rand(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def _pool():
    """A, B share a guide; A, C share the central donor window; observed = A's donor with an error outside it."""
    s, e = om.MIDDLE_DONOR_START, om.MIDDLE_DONOR_END
    d_a = _rand(129, 1)
    d_c = _rand(s, 2) + d_a[s:e] + _rand(129 - e, 3)
    d_b = _rand(129, 4)
    g_ab, g_c = _rand(20, 5), _rand(20, 6)
    df = pd.DataFrame({"oligo_name": ["A", "B", "C"], "guide": [g_ab, g_ab, g_c], "donor": [d_a, d_b, d_c]})
    df["middle_donor"] = df["donor"].str.slice(s, e)
    observed = d_a[:5] + ("A" if d_a[5] != "A" else "C") + d_a[6:]
    return df, g_ab, observed


def test_c002_core_resolver_takes_the_guide_window_intersection():
    df, guide, donor = _pool()
    res = om.resolve_oligo_name(guide, donor, om.build_lookup_dictionaries(df))
    assert (res.oligo_name, res.match_type) == ("A", "intersection")


def test_c002_vectorized_counts_path_agrees_with_the_scalar_resolver(tmp_path):
    df, guide, donor = _pool()
    lk = om.build_lookup_dictionaries(df)
    pd.DataFrame({"library_ID": ["V1"], "total_counts": [7], "leader": ["x"], "guide": [guide], "donor": [donor],
                  "SPS": ["s"], "bc0": ["b"], "msd": ["m"], "donor_bc0_fragment": ["f"], "seq": ["q"]}).to_csv(
        tmp_path / "in.tsv", sep="\t", index=False)
    om.process_counts_file(tmp_path / "in.tsv", tmp_path / "out.tsv", lk)
    out = pd.read_csv(tmp_path / "out.tsv", sep="\t").iloc[0]
    assert (out["oligo_name"], out["match_type"]) == ("A", "intersection")


def test_c002_unique_guide_and_unique_donor_that_disagree_stay_unassigned():
    df, guide, _ = _pool()
    res = om.resolve_oligo_name(df.loc[2, "guide"], df.loc[1, "donor"], om.build_lookup_dictionaries(df))
    assert res.oligo_name is None and res.match_type == "conflicting"


def test_c003_identical_guide_donor_alias_does_not_depend_on_row_order():
    g, d = _rand(20, 7), _rand(129, 8)
    rows = pd.DataFrame({"oligo_name": ["ALIAS_Y", "ALIAS_X"], "guide": [g, g], "donor": [d, d]})
    rows["middle_donor"] = rows["donor"].str.slice(om.MIDDLE_DONOR_START, om.MIDDLE_DONOR_END)
    a = om.resolve_oligo_name(g, d, om.build_lookup_dictionaries(rows))
    b = om.resolve_oligo_name(g, d, om.build_lookup_dictionaries(rows.iloc[::-1].reset_index(drop=True)))
    assert a.oligo_name == b.oligo_name == "ALIAS_X"


def test_c009_dominant_oligo_is_read_support_not_name_order():
    counts = pd.DataFrame({"bc0": ["b1", "b1"], "oligo_name": ["A_high", "Z_low"], "count": [100, 1]})
    out = om.calculate_bc0_purity(counts).set_index("bc0")
    assert out.loc["b1", "dominant_oligo"] == "A_high" and out.loc["b1", "purity"] == pytest.approx(100 / 101)


def _fastq(path, n):
    with open(path, "w") as fh:
        for i in range(n):
            seq = _rand(240, 100 + i)
            fh.write(f"@seq_{i + 1}_count=1\n{seq}\n+\n{'I' * len(seq)}\n")


def test_resume_reparses_a_changed_fastq_at_the_same_path(tmp_path):
    from magestic.pipelines.guide_donor_bc0.core import sharded_combine as sc
    fq = tmp_path / "lib1_sample_1.fastq"
    _fastq(fq, 3)
    row = {"shard_id": "1", "fastq": str(fq), "library_ID": "lib1", "sample_number": "sample_1",
           "sequencing_date": "20240906", "inner_primers": "KR1882-KR1883", "counts_col": "fwd_0906",
           "direction": "fwd"}
    assert sc.parse_shard(row, tmp_path / "shards")["rows_out"] == 3
    assert sc.parse_shard(row, tmp_path / "shards")["skipped"]          # unchanged input: reused
    _fastq(fq, 11)
    second = sc.parse_shard(row, tmp_path / "shards")
    assert not second["skipped"] and second["rows_out"] == 11


def test_qc_reads_production_column_names_and_bins_match_labels():
    from magestic.pipelines.guide_donor_bc0.qc.pipeline_qc import (QCConfig, validate_bc0_purity,
                                                                    validate_library_representation)
    m = validate_bc0_purity(pd.DataFrame({"guide": ["g1", "g2"], "bc0": ["b1", "b2"],
                                          "guide_bc0_purity": [0.2, 0.3]}), QCConfig())
    assert m.mean_bc0_purity == pytest.approx(0.25)
    r = validate_library_representation(pd.DataFrame({"oligo_name": ["o1", "o2"], "counts": [0, 50]}), 2, QCConfig())
    assert r.n_oligos_covered and r.n_oligos_covered > 0
    rep = pd.DataFrame({"oligo_name": list("abcd"), "reads": [0, 1, 10, 100]})
    dist = validate_library_representation(rep, 4, QCConfig()).coverage_depth_distribution
    assert (dist.get("0"), dist.get("1-9"), dist.get("10-99"), dist.get("100-999")) == (1, 1, 1, 1)
