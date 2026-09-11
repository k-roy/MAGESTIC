"""gdb step 02 after the GDB-01 root fix: parse_sequence extracts G-ending Cas12a guides exactly, the parse
accounting closes, and the sharded parse+merge path reproduces process_library's table exactly."""
import json
import random
from pathlib import Path

import pandas as pd
import pytest

from magestic.pipelines.guide_donor_bc0.core import fastq_processing as fp
from magestic.pipelines.guide_donor_bc0.core import sharded_combine as sc
from magestic.utils.sequence_anchors import RE_SITE_PATTERNS

SP_SITE = dict(RE_SITE_PATTERNS)["SpCas9"]
LB_SITE = dict(RE_SITE_PATTERNS)["LbCas12a"]
MSD = fp.SEQ_CONST.msd
rng = random.Random(7)


def _seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _amplicon(guide, site, donor, sps, bc0, leader="GACTAC"):
    return leader + guide + site + donor + sps + bc0 + MSD


def test_parse_sequence_g_ending_cas12a_guide_is_exact_and_accounting_closes():
    stats = fp.new_parse_stats()
    guide = _seq(19) + "G"
    donor, sps, bc0 = _seq(117), _seq(20), _seq(10)
    read = _amplicon(guide, LB_SITE, donor, sps, bc0)
    leader, g, d, s, b, m, frag = fp.parse_sequence(read, stats=stats)
    assert (leader, g, d, s, b, m) == ("GACTAC", guide, donor, sps, bc0, MSD)
    assert frag == read[-12 - 60:-12]
    # an SpCas9 read, one with a 1-nt site error, one with no msd, one with no site
    fp.parse_sequence(_amplicon(_seq(20), SP_SITE, _seq(129), sps, bc0), stats=stats)
    fp.parse_sequence(_amplicon(_seq(20), "GTTTGAAGAGA", _seq(129), sps, bc0), stats=stats)
    fp.parse_sequence(_amplicon(_seq(20), SP_SITE, _seq(129), sps, bc0)[:-12] + "CCCCCCCCCCCC", stats=stats)
    fp.parse_sequence(_amplicon(_seq(20), "AAAAAAAAAAA", _seq(129), sps, bc0), stats=stats)
    assert stats["reads"] == 5
    assert stats["reads"] == stats["no_msd"] + stats["no_re_site"] + stats["bad_frame"] + stats["parsed"]
    assert stats["parsed"] == 3 and stats["no_msd"] == 1 and stats["no_re_site"] == 1
    assert stats["re_site_LbCas12a_exact"] == 1 and stats["re_site_SpCas9_exact"] == 1 and stats["re_site_SpCas9_mm1"] == 1
    assert fp.parse_parameters()["re_site_max_mismatches"] == 1


def test_bad_frame_is_counted():
    stats = fp.new_parse_stats()
    read = "AC" + LB_SITE + _seq(40) + _seq(20) + _seq(10) + MSD    # site too close to the 5' end for a 20-nt guide
    assert fp.parse_sequence(read, stats=stats) == ("NA",) * 7
    assert stats["bad_frame"] == 1


def _write_fastq(path, records):
    with open(path, "w") as fh:
        for count, seq in records:
            fh.write(f"@{count}\n{seq}\n+\n{'J' * len(seq)}\n")


@pytest.fixture
def tiny_project(tmp_path):
    """Two libraries, fwd + rev samples over two dates; the rev sample's reads are stored reverse-complemented
    and one read is shared between a fwd and a rev file so the merge has to sum across shards."""
    fq_dir = tmp_path / "collapsed_fastq"
    fq_dir.mkdir()
    sps_lb, sps_sp = _seq(20), _seq(20)
    lb_reads = [(50, _amplicon(_seq(19) + "G", LB_SITE, _seq(117), sps_lb, _seq(10))) for _ in range(5)]
    lb_reads.append((7, "ACGT" * 30))                                  # unparseable, must still be counted
    sp_reads = [(30, _amplicon(_seq(20), SP_SITE, _seq(129), sps_sp, _seq(10))) for _ in range(4)]
    shared = lb_reads[0][1]
    _write_fastq(fq_dir / "20210906_plate_5D_sample_1_collapsed.fastq", lb_reads)                      # V454 fwd
    _write_fastq(fq_dir / "20211012_plate_5D_sample_1_collapsed.fastq", lb_reads[2:4])                 # V454 fwd, date 2
    _write_fastq(fq_dir / "20210906_plate_5D_sample_2_collapsed.fastq",
                 [(9, fp.rev_comp(shared)), (3, fp.rev_comp(lb_reads[4][1]))])                         # V454 rev
    _write_fastq(fq_dir / "20210906_plate_5D_sample_3_collapsed.fastq", sp_reads)                      # V429 fwd
    keyfile = tmp_path / "key.tsv"
    keyfile.write_text(
        "sample_number\tlibrary_ID\tinner_primers\touter_primers\tsequencing_dates\tpolymerase\n"
        "sample_1\tV454\tKR1952-KR1590\tplate_5D\t20210906,20211012\tQ5\n"
        "sample_2\tV454\tKR1884-KR1953\tplate_5D\t20210906\tQ5\n"
        "sample_3\tV429\tKR1952-KR1590\tplate_5D\t20210906\tQ5\n"
        "sample_4\tempty\tKR1952-KR1590\tplate_5D\t20210906\tQ5\n")
    return tmp_path, keyfile, fq_dir, shared


def test_sharded_parse_and_merge_reproduce_process_library(tiny_project):
    tmp, keyfile, fq_dir, shared = tiny_project
    # reference: the in-process path
    lib2s, all_dates = fp.load_sample_keyfile(keyfile)
    ref_dir = tmp / "ref"
    ref_dir.mkdir()
    ref_stats = {lib: fp.process_library(lib, lib2s[lib], fq_dir, ref_dir, all_dates) for lib in ("V454", "V429")}
    assert ref_stats["V454"]["parse_stats"]["re_site_LbCas12a_exact"] == 9   # 5 + 2 + 2 parsed reads
    assert ref_stats["V454"]["parse_stats"]["no_msd"] == 1

    # sharded path
    manifest = sc.build_manifest(keyfile, fq_dir)
    assert [r["library_ID"] for r in manifest] == ["V454", "V454", "V454", "V429"]
    assert [r["counts_col"] for r in manifest] == ["fwd_0906", "fwd_1012", "rev_0906", "fwd_0906"]
    mpath = tmp / "manifest.tsv"
    sc.write_manifest(manifest, mpath)
    manifest = sc.read_manifest(mpath)
    shard_dir = tmp / "shards"
    for r in manifest:
        res = sc.parse_shard(r, shard_dir)
        assert res["skipped"] is False
        assert sc.parse_shard(r, shard_dir)["skipped"] is True        # idempotent on rc=0
        assert Path(res["tsv"] + ".provenance.json").exists()
    out_dir = tmp / "merged"
    for lib in ("V454", "V429"):
        res = sc.merge_library(lib, manifest, shard_dir, out_dir)
        assert res["skipped"] is False and res["reads_in"] == res["reads_out"]
        assert sc.merge_library(lib, manifest, shard_dir, out_dir)["skipped"] is True

    for lib in ("V454", "V429"):
        a = pd.read_csv(ref_dir / f"{lib}_guide_donor_bc0_counts.tsv", sep="\t")
        b = pd.read_csv(out_dir / f"{lib}_guide_donor_bc0_counts.tsv", sep="\t")
        assert list(a.columns) == list(b.columns)
        key = ["seq"]
        pd.testing.assert_frame_equal(a.sort_values(key).reset_index(drop=True),
                                      b.sort_values(key).reset_index(drop=True))
    v454 = pd.read_csv(out_dir / "V454_guide_donor_bc0_counts.tsv", sep="\t")
    # the read written to both a fwd and a rev file must be ONE row carrying both counts
    shared_row = v454[v454["seq"] == shared]
    assert len(shared_row) == 1
    assert int(shared_row["fwd_0906"].iloc[0]) == 50 and int(shared_row["rev_0906"].iloc[0]) == 9
    assert int(shared_row["total_counts"].iloc[0]) == 59
    assert ((v454["fwd_0906"] > 0) & (v454["rev_0906"] > 0)).sum() == 2   # `shared` and lb_reads[4]
    assert len(v454) == 6 and int(v454["total_counts"].sum()) == 50 * 5 + 7 + 50 * 2 + 9 + 3
    # G-ending Cas12a guides are exact 20-mers ending in G, never 'T' + 19.
    # (pandas reads the literal string "NA" back as NaN, hence dropna rather than != "NA".)
    guides = v454["guide"].dropna()
    assert len(guides) == 5
    assert (guides.str.len() == 20).all() and guides.str.endswith("G").all()
    # the decisive check: every extracted guide is the 20 nt that really precede the cloning site
    # in its own read (a 1-base shift, the GDB-01 signature, breaks this even though the shifted
    # string is still a 20-mer). The leader is 6 nt, so the guide is seq[6:26].
    parsed = v454.dropna(subset=["guide"])
    assert (parsed["seq"].str.slice(6, 26) == parsed["guide"]).all()
    assert (parsed["seq"].str.slice(26, 37) == LB_SITE).all()

    # provenance: inputs carry md5s, counts close, the policy ruling is quoted
    prov = json.load(open(str(out_dir / "V454_guide_donor_bc0_counts.tsv") + ".provenance.json"))
    assert prov["schema"].startswith("magestic-step-provenance")
    assert all("md5" in i for i in prov["inputs"]) and len(prov["inputs"]) == 3
    assert prov["counts"]["reads_in"] == prov["counts"]["reads_out"]
    assert prov["parameters"]["re_site_max_mismatches"] == 1
    assert any("hodge-podge" in d["quote"] for d in prov["decisions"])
    assert prov["outputs"][0]["rows"] == 6


def test_manifest_refuses_missing_fastq(tiny_project):
    tmp, keyfile, fq_dir, _shared = tiny_project
    (fq_dir / "20211012_plate_5D_sample_1_collapsed.fastq").unlink()
    with pytest.raises(FileNotFoundError):
        sc.build_manifest(keyfile, fq_dir)


def test_merge_refuses_shard_without_sentinel(tiny_project):
    tmp, keyfile, fq_dir, _shared = tiny_project
    manifest = sc.build_manifest(keyfile, fq_dir)
    with pytest.raises(RuntimeError):
        sc.merge_library("V429", manifest, tmp / "no_shards", tmp / "out")
