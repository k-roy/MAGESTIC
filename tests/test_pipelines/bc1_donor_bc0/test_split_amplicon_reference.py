"""Golden fixture for the split-amplicon bc1 reference route (magestic.pipelines.bc1_donor_bc0.split_amplicon).

Synthetic reads are built from the read geometry the route parses, so each assertion names a value a person can
check by hand. The real-data parity test against the NRD1 2025 table is planning/378 (it needs Sherlock scratch).
"""
import gzip
import random

import pandas as pd
import pytest

pytest.importorskip("pysam")
pytest.importorskip("rapidfuzz")

from magestic.pipelines.bc1_donor_bc0 import split_amplicon as sa  # noqa: E402

rng = random.Random(378)


def _seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


SPS = _seq(20)
MSD = _seq(12)
FWD, REV = _seq(20), _seq(20)
OLIGOS = {}
for name in ("oligo_001", "oligo_002", "oligo_003"):
    guide, donor, bc0 = _seq(20), _seq(129), _seq(10)
    OLIGOS[name] = {"guide": guide, "donor": donor, "bc0": bc0,
                    "oligo_seq": FWD + guide + sa.RE_SITE.lower() + donor + REV}
BC1 = {"A": _seq(20), "B": _seq(20), "C": _seq(20)}


def fragment(name):
    o = OLIGOS[name]
    return o["donor"][-30:] + SPS + o["bc0"]


def r1(bc1):
    return ("TCGA" + sa.FWD_BC1_PREFIX + bc1 + sa.SACAS9_SCEI_FRAGMENT + _seq(150))[:150]


def r2(name, mutate=False):
    o = OLIGOS[name]
    donor = o["donor"]
    if mutate:  # one substitution inside the 60-nt join key
        i = len(donor) - 5
        donor = donor[:i] + ("A" if donor[i] != "A" else "C") + donor[i + 1:]
    return sa.rev_comp(donor + SPS + o["bc0"] + MSD)[:150]


def merged(name):
    o = OLIGOS[name]
    return sa.rev_comp("TCC" + sa.RE_SITE + o["donor"] + SPS + o["bc0"] + _seq(9))


def write_pair(d, prefix, pairs):
    for k, idx in (("R1", 0), ("R2", 1)):
        with gzip.open(d / f"{prefix}_trimmed_{k}.fastq.gz", "wt") as fh:
            for i, p in enumerate(pairs):
                fh.write(f"@r{i}\n{p[idx]}\n+\n{'I' * len(p[idx])}\n")


def _make_project(tmp_path, extra=()):
    fq, coll = tmp_path / "fragment_fastq", tmp_path / "donor_bc0_fastq"
    fq.mkdir()
    coll.mkdir()
    s1 = ([(r1(BC1["A"]), r2("oligo_001"))] * 10 + [(r1(BC1["A"]), r2("oligo_001", mutate=True))]
          + [(r1(BC1["B"]), r2("oligo_002"))] * 5 + [(r1(BC1["B"]), r2("oligo_003"))] * 5
          + [(r1(BC1["C"]), r2("oligo_001"))] * 7)
    s2 = ([(r1(BC1["A"]), r2("oligo_001"))] * 8
          + [(r1(BC1["B"]), r2("oligo_002"))] * 5 + [(r1(BC1["B"]), r2("oligo_003"))] * 5)
    write_pair(fq, "20250815_plate_5B_sample_1", s1 + list(extra))
    write_pair(fq, "20250815_plate_5B_sample_2", s2 + list(extra))
    with open(coll / "20250815_plate_1C_sample_1_collapsed.fastq", "w") as fh:
        for n, (name, c) in enumerate((("oligo_001", 40), ("oligo_002", 30), ("oligo_003", 20))):
            s = merged(name)
            fh.write(f"@{c}\n{s}\n+\n{'I' * len(s)}\n")
    key = pd.DataFrame({
        "sequencing_date": [20250815] * 3,
        "description": ["bc1-donor-bc0", "bc1-donor-bc0", "donor-bc0"],
        "outer_primers": ["plate_5B", "plate_5B", "plate_1C"],
        "sample_number": ["sample_1", "sample_2", "sample_1"],
        "stage": ["plate wash"] * 3, "yL": ["yL435"] * 3, "polymerase": ["KAPA"] * 3,
    })
    key.to_csv(tmp_path / "key.tsv", sep="\t", index=False)
    pd.DataFrame({"oligo_name": list(OLIGOS), "oligo_seq": [o["oligo_seq"] for o in OLIGOS.values()]}).to_csv(
        tmp_path / "pool.tsv", sep="\t", index=False)
    pd.DataFrame({"library_ID": ["V625"] * 3, "guide": [o["guide"] for o in OLIGOS.values()],
                  "donor": [o["donor"] for o in OLIGOS.values()],
                  "donor_bc0_fragment": [fragment(n) for n in OLIGOS], "oligo_name": list(OLIGOS)}).to_csv(
        tmp_path / "gdb.tsv", sep="\t", index=False)
    pd.DataFrame({"oligo_name": list(OLIGOS), "common_ORF_name": ["NRD1"] * 3}).to_csv(
        tmp_path / "info.tsv", sep="\t", index=False)
    cfg = {"out_dir": str(tmp_path / "out"), "key": str(tmp_path / "key.tsv"), "fragment_fastq_dir": str(fq),
           "donor_bc0_fastq_dir": str(coll), "gdb_table": str(tmp_path / "gdb.tsv"),
           "oligo_pool": str(tmp_path / "pool.tsv"), "oligo_info": str(tmp_path / "info.tsv"),
           "table_name": "ref.tsv"}
    return cfg


@pytest.fixture
def project(tmp_path):
    return _make_project(tmp_path)


def _r2_error_variant_sorting_first(name):
    """R2 for `name` with one substitution in the join key that makes the 60-nt fragment sort BEFORE the
    designed one, so a stable tie order picks the error variant."""
    o = OLIGOS[name]
    d = o["donor"]
    i = next(k for k in range(len(d) - 30, len(d)) if d[k] != "A")
    bad = d[:i] + "A" + d[i + 1:]
    assert bad[-30:] + SPS + o["bc0"] < fragment(name)
    return sa.rev_comp(bad + SPS + o["bc0"] + MSD)[:150]


def test_stage_b_counts_and_geometry(project):
    sa.run(project, stages=["b"], n_jobs=1)
    b = pd.read_csv(f"{project['out_dir']}/bc1_donor_bc0_fragment_counts/20250815_plate_5B_sample_1"
                    f"{sa.FRAGMENT_SUFFIX}", sep="\t")
    assert list(b.columns) == sa.FRAGMENT_HEADER
    assert b["counts"].sum() == 28
    top = b.iloc[0]
    assert top["bc1"] == BC1["A"] and top["counts"] == 10
    assert top["donor_bc0_fragment"] == fragment("oligo_001")
    assert top["PS1"] == "TCGA" + sa.FWD_BC1_PREFIX and bool(top["SaCas9_SceI_fragment_found"])


def test_full_route_values(project):
    s = sa.run(project, n_jobs=1)
    t = pd.read_csv(f"{project['out_dir']}/data_tables/ref.tsv", sep="\t")
    a = t[t["bc1"] == BC1["A"]]
    assert len(a) == 1  # the 1-mismatch fragment collapsed into the 18-read one
    a = a.iloc[0]
    assert a["total_bc1_donor_bc0_counts"] == 19 and a["num_PCR_replicates"] == 2 and a["bc1_purity"] == 1.0
    assert a["db_oligo_name"] == "oligo_001" and a["gdb_oligo_name"] == "oligo_001" and bool(a["bc1_unique"])
    assert bool(a["unambiguous_oligo_name"]) and bool(a["db_perfect_donor"])
    assert BC1["C"] not in set(t["bc1"])  # one PCR replicate only
    b = t[t["bc1"] == BC1["B"]]
    assert sorted(b["db_oligo_name"]) == ["oligo_002", "oligo_003"]
    assert (b["bc1_purity"] == 0.5).all() and not b["bc1_unique"].any()
    assert (t["bc1_unique"] == t["multiple_donor_bc0_fragments_passing_filter"]).all()
    assert s["c"]["bc1_top_tied"] == 1 and s["h"]["bc1_out"] == 2


def test_top_fragment_ties_are_reported_and_deterministic(project, tmp_path):
    sa.run(project, stages=["b", "e"], n_jobs=1)
    out1, out2 = tmp_path / "c1", tmp_path / "c2"
    key = pd.read_csv(project["key"], sep="\t")
    frag = f"{project['out_dir']}/bc1_donor_bc0_fragment_counts"
    sa.combine_bc1_fragments(key, frag, str(out1), n_jobs=1)
    sa.combine_bc1_fragments(key, frag, str(out2), n_jobs=2)
    a = (out1 / "top_bc1_donor_bc0_fragments.tsv").read_bytes()
    assert a == (out2 / "top_bc1_donor_bc0_fragments.tsv").read_bytes()
    ties = pd.read_csv(out1 / "top_fragment_ties.tsv", sep="\t")
    assert set(ties["bc1"]) == {BC1["B"]} and len(ties) == 2
    assert ties["top_bc1_donor_fragment"].nunique() == 1


def test_stage_g_ambiguous_donor_takes_first_in_file_order(tmp_path):
    donor = _seq(129)
    pool = pd.DataFrame({"oligo_name": ["X_first", "X_second"],
                         "oligo_seq": [FWD + _seq(20) + sa.RE_SITE.lower() + donor + REV] * 2})
    pool.to_csv(tmp_path / "pool.tsv", sep="\t", index=False)
    frag = donor[-30:] + SPS + "ACGTACGTAC"
    pd.DataFrame({"RE_site_found": [True], "RE_site": [sa.RE_SITE], "donor": [donor], "SPS": [SPS],
                  "bc0": ["ACGTACGTAC"], "donor_bc0_fragment": [frag], "counts": [5]}).to_csv(
        tmp_path / "combined.tsv", sep="\t", index=False)
    sa.map_donor_bc0_to_oligos(str(tmp_path / "combined.tsv"), str(tmp_path / "pool.tsv"), str(tmp_path))
    m = pd.read_csv(tmp_path / "donor_bc0_to_oligo_mapping.tsv", sep="\t").iloc[0]
    assert m["db_oligo_name"] == "X_first" and not bool(m["unambiguous_oligo_name"])


def test_tie_order_is_validated(project):
    project["tie_order"] = "quick"
    with pytest.raises(ValueError):
        sa.run(project, stages=["c"], n_jobs=1)


@pytest.mark.parametrize("tie_order, expect", [("designed", "oligo_001"), ("stable", "")])
def test_designed_tie_order_keeps_the_designed_fragment_of_a_one_read_tie(tmp_path, tie_order, expect):
    # bc1 D: in each sample one read of oligo_001's fragment and one read of a 1-edit error variant. The two
    # collapse into one; which names the group is a 1-vs-1 tie. Stable order picks the variant (it sorts
    # first) and the bc1 loses its design; the designed order keeps oligo_001.
    d = _seq(20)
    cfg = _make_project(tmp_path, extra=[(r1(d), r2("oligo_001")), (r1(d), _r2_error_variant_sorting_first("oligo_001"))])
    cfg["tie_order"] = tie_order
    sa.run(cfg, n_jobs=1)
    t = pd.read_csv(f"{cfg['out_dir']}/data_tables/ref.tsv", sep="\t", keep_default_na=False, dtype=str)
    row = t[t["bc1"] == d]
    assert len(row) == 1 and row.iloc[0]["total_bc1_donor_bc0_counts"] == "4"
    assert row.iloc[0]["db_oligo_name"] == expect
