"""Fragment identity must survive the real Stage3-to-Stage4 cleaned handoff.

All alleles, positions and expected counts are declared before execution. No
private artifacts, audit-path imports, network or mocked samtools are used.
Historical parity intentionally retains its QNAME-only grouping contract.
"""
from collections import Counter
import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import shutil

import pysam
import pytest
import yaml

from magestic.pipelines.wgs.config import load_config
from magestic.pipelines.wgs.pipeline import step_3_sv_plasmid as stage3
from magestic.pipelines.wgs.pipeline import step_4_reconcile as stage4
from trace_crispr.wgs import sv_plasmid as sv


def sam_record(name="shared", rg="genome", *, template=False, flag=99, mapq=60,
               optional=()):
    pos, length, mate, tlen = (301, 10, 301, 10) if template else (296, 20, 351, 75)
    columns = [name, str(flag), "chrII", str(pos), str(mapq), f"{length}M",
               "=", str(mate), str(tlen), "A" * length, "D" * length]
    if rg is not None:
        columns.append("RG:Z:" + rg)
    return sv.parse_sam_line("\t".join(columns + list(optional)) + "\n")


def classify(records, *, mode="standard"):
    return sv.classify_donor_reads(records, "chrII", 301, 310, {}, dup_mode=mode)


@pytest.mark.parametrize("bad", [
    {"template": True}, {"mapq": 19}, {"mapq": 255}, {"flag": 99 | 8},
])
def test_discard_in_other_read_group_cannot_poison_valid_fragment(bad):
    result = classify([sam_record(rg="genome"), sam_record(rg="reject", **bad)])
    assert result.frags_genome == 1
    assert result.frags_plasmid == int(bad.get("template", False))
    assert result.variant_cov == 1


@pytest.mark.parametrize("rg_pair", [("lane1", "lane2"), (None, "lane2"), ("lane:1", "lane")])
def test_retained_fragment_denominator_distinguishes_read_groups(rg_pair):
    result = classify([sam_record(rg=rg) for rg in rg_pair])
    assert (result.frags_genome, result.frags_plasmid, result.variant_cov) == (2, 0, 2)


@pytest.mark.parametrize("rg", [None, "same_lane"])
def test_same_namespace_retains_pair_discard_precedence(rg):
    result = classify([sam_record(rg=rg), sam_record(rg=rg, template=True, flag=147)])
    assert (result.frags_genome, result.frags_plasmid, result.variant_cov) == (0, 1, 0)


@pytest.mark.parametrize("flag", [99 | 256, 99 | 2048])
def test_nonprimary_same_name_does_not_override_valid_primary(flag):
    result = classify([sam_record(), sam_record(template=True, flag=flag)])
    assert (result.frags_genome, result.frags_plasmid, result.variant_cov) == (1, 0, 1)


def test_unpaired_rejection_cannot_poison_paired_same_name_same_rg():
    # Stage4 treats an unpaired alignment as its own observation. HalfB may reject
    # that unsupported single read, but it must not reject the separate valid pair.
    result = classify([sam_record(), sam_record(flag=0)])
    assert (result.frags_genome, result.frags_plasmid, result.variant_cov) == (1, 0, 1)


@pytest.mark.parametrize("optional", [
    ("RG:Z:second",), ("RG:i:1",), ("RG:Z:",),
])
def test_ambiguous_or_malformed_rg_namespace_is_rejected(optional):
    # No new SamRecord field is assumed here: the public SAM parser is the API.
    # The first case duplicates a valid RG; the others declare invalid types/IDs.
    rg = "first" if optional[0] == "RG:Z:second" else None
    with pytest.raises(ValueError, match="[Rr]ead.group|RG"):
        sam_record(rg=rg, optional=optional)


def test_parity_explicitly_preserves_historical_qname_namespace():
    result = classify([sam_record(rg="genome"), sam_record(rg="template", template=True)], mode="parity")
    assert result.unique_genome == result.unique_plasmid == {"shared"}
    assert result.included == result.discarded == {"shared"}
    assert (result.frags_genome, result.frags_plasmid, result.variant_cov) == (1, 1, 1)


BC = "ACGTCAGTTCGATGACCTGA"
GLEFT, GRIGHT = "ACCCTGTTAT", "GAGCATGACC"
PLEFT, PRIGHT = "GGTCATGCTC", "ATAACAGGGT"


def make_alignment(header, *, name, contig, start, sequence, cigar, flag=99,
                   mate_start=180, tlen=160, rg=None):
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.query_sequence = sequence
    read.flag = flag
    read.reference_id = header.get_tid(contig)
    read.reference_start = start
    read.mapping_quality = 60
    read.cigarstring = cigar
    read.next_reference_id = header.get_tid(contig) if flag & 1 else -1
    read.next_reference_start = mate_start if flag & 1 else -1
    read.template_length = tlen if flag & 1 else 0
    read.query_qualities = pysam.qualitystring_to_array("D" * len(sequence))
    if rg is not None:
        read.set_tag("RG", rg)
    return read


def build_case(root, *, collision=True, rejected_unpaired=False, extra_off_target=False):
    """Four retained genomic pairs and four discarded fragments/single reads."""
    samtools = shutil.which("samtools")
    assert samtools, "the real Stage3-to-Stage4 regression requires samtools"
    root.mkdir()
    refs = {
        "chrII": "A" * 90 + GLEFT + "N" * 20 + GRIGHT + "A" * 370,
        "plasmid_demo": "A" * 90 + PLEFT + "N" * 20 + PRIGHT + "A" * 370,
    }
    fasta = root / "ref.fa"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in refs.items()))
    pysam.faidx(str(fasta))
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": len(seq), "M5": hashlib.md5(seq.encode()).hexdigest()}
               for name, seq in refs.items()],
        "RG": [{"ID": name, "SM": "well_A1"} for name in ("genome", "template", "off_target")],
        "PG": [{"ID": "fixture_generator", "PN": "synthetic_fixture", "CL": f"synthetic fixture reference {fasta}"}],
    })
    reverse_bc = BC.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    reads = []
    for contig, left, barcode, prefix in (
        ("chrII", GLEFT, reverse_bc, "g"), ("plasmid_demo", PLEFT, BC, "p"),
    ):
        reads.extend([
            make_alignment(header, name=prefix + "1", contig=contig, start=90,
                           sequence=left + barcode, cigar="10M20S"),
            make_alignment(header, name=prefix + "2", contig=contig, start=89,
                           sequence="A" + left + barcode, cigar="11M20S", tlen=161),
        ])
    for i in range(4):
        genome_name = f"q{i}" if collision else f"genome_q{i}"
        reads.extend([
            make_alignment(header, name=genome_name, rg="genome", contig="chrII", start=295,
                           sequence="A" * 5 + "C" + "A" * 14, cigar="20M", flag=99, mate_start=350, tlen=75),
            make_alignment(header, name=genome_name, rg="genome", contig="chrII", start=350,
                           sequence="A" * 20, cigar="20M", flag=147, mate_start=295, tlen=-75),
        ])
        if rejected_unpaired:
            reads.append(make_alignment(header, name=f"q{i}", rg="genome", contig="chrII", start=300,
                                        sequence="C" + "A" * 9, cigar="10M", flag=0))
        else:
            for flag, tlen in ((99, 10), (147, -10)):
                reads.append(make_alignment(header, name=f"q{i}", rg="template", contig="chrII", start=300,
                                            sequence="C" + "A" * 9, cigar="10M", flag=flag, mate_start=300, tlen=tlen))
    if extra_off_target:
        reads.append(make_alignment(header, name="q0", rg="off_target", contig="chrII", start=420,
                                    sequence="A" * 20, cigar="20M", flag=0))
    bam = root / "input.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
        for read in sorted(reads, key=lambda record: (record.reference_id, record.reference_start)):
            handle.write(read)
    pysam.index(str(bam))
    manifest = root / "manifest.tsv"
    manifest.write_text(f"sample\tbam\treference\nwell_A1\t{bam}\t{fasta}\n")
    windows = root / "windows.tsv"
    windows.write_text('oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\tintended_variants\n'
                       'oligo_1\tchrII\t301\t310\t[{"pos":301,"ref":"A","alt":"C"}]\n')
    assignments = root / "assignments.tsv"
    assignments.write_text("sample\toligo_name\nwell_A1\toligo_1\n")
    bridge = root / "bridge.tsv"
    bridge.write_text(f"bc1\toligo_name\n{BC}\toligo_1\n")
    blacklist = root / "blacklist.tsv"
    blacklist.write_text("plasmid_demo\t1\t500\n")
    config_path = root / "config.yaml"
    config_path.write_text(yaml.safe_dump({
        "cohort_name": "synthetic_namespace", "sequencing_date": "20260920",
        "paths": {"project_root": str(root), "manifest": str(manifest)},
        "locus": {"barcode_locus": "chrII:101-120", "barcode_left_sequence": GLEFT,
                  "barcode_right_sequence": GRIGHT, "plasmid_bc_locus": "plasmid_demo:101-120",
                  "plasmid_bc_left_sequence": PLEFT, "plasmid_bc_right_sequence": PRIGHT},
        "bc1": {"use_plasmid": True, "bridge_table": str(bridge)},
        "reconcile": {"oligo_windows": str(windows), "oligo_assignments": str(assignments),
                      "read_source": "donor_cleaned"},
        "sv_plasmid": {"blacklist": str(blacklist), "blacklist_reference": str(fasta), "samtools": samtools},
    }))
    return load_config(config_path), bam, len(reads)


def only_row(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["sample"] == "well_A1"
    return rows[0]


@pytest.mark.parametrize("collision,rejected_unpaired", [(False, False), (True, False), (True, True)])
def test_public_cleaned_handoff_preserves_four_genomic_fragments(tmp_path, collision, rejected_unpaired):
    config, bam_path, input_count = build_case(tmp_path / "case", collision=collision,
                                              rejected_unpaired=rejected_unpaired)
    assert input_count == (16 if rejected_unpaired else 20)
    with pysam.AlignmentFile(str(bam_path)) as handle:
        assert sum(1 for _ in handle) == input_count
    # Exact same public consumer, changing only declared input policy. Raw control
    # establishes eight callable fragments before any attempted donor cleaning.
    raw_config = replace(config, reconcile=replace(config.reconcile, read_source="raw"))
    raw = only_row(stage4.run_step_4(raw_config))
    raw_evidence = json.loads(raw["allele_counts"])
    assert len(raw_evidence) == 1
    assert (raw_evidence[0]["alt_count"], raw_evidence[0]["ref_count"], raw_evidence[0]["depth"]) == (8, 0, 8)
    assert "RAW_DONOR_READS" in raw["flags"].split(",")

    cleaned = only_row(stage4.run_step_4(config))
    evidence = json.loads(cleaned["allele_counts"])
    producer = only_row(stage3.shard_path(config, "well_A1"))
    cleaned_path = stage3.cleaned_bam_path(config, "well_A1")
    with pysam.AlignmentFile(str(cleaned_path)) as handle:
        assert handle.check_index()
        kept = list(handle)
        indexed = list(handle.fetch("chrII", 294, 370))
    genome_name = (lambda i: f"q{i}") if collision else (lambda i: f"genome_q{i}")
    expected_mates = Counter(("genome", genome_name(i), flag, start)
                             for i in range(4) for flag, start in ((99, 295), (147, 350)))
    actual_mates = Counter((read.get_tag("RG"), read.query_name, read.flag, read.reference_start)
                          for read in indexed)
    assert actual_mates == expected_mates
    assert len(kept) == 12  # Eight genome-mate alignments plus four raw barcode alignments.
    assert {read.query_name for read in kept if not read.has_tag("RG")} == {"g1", "g2", "p1", "p2"}
    assert producer["half_b_status"] == "run"
    assert [producer[key] for key in ("read_count_mid_point", "frags_genome", "frags_plasmid", "frags_chimeric")] == ["4", "4", "0" if rejected_unpaired else "4", "0"]
    assert len(evidence) == 1
    assert (evidence[0]["alt_count"], evidence[0]["ref_count"], evidence[0]["depth"]) == (4, 0, 4)
    assert (cleaned["label"], cleaned["n_nocov"], cleaned["n_genome_bc"], cleaned["n_plasmid_bc"]) == ("SINGLE_EDIT", "0", "1", "1")
    assert "DONOR_CLEANED_READS" in cleaned["flags"].split(",")


def test_full_stream_retains_unrelated_offtarget_collision(tmp_path):
    config, _, input_count = build_case(tmp_path / "case", collision=False, extra_off_target=True)
    assert input_count == 21
    stage3.run_step_3(config)
    with pysam.AlignmentFile(str(stage3.cleaned_bam_path(config, "well_A1"))) as handle:
        assert handle.check_index()
        off_target = list(handle.fetch("chrII", 420, 440))
    assert [(read.get_tag("RG"), read.query_name, read.flag, read.reference_start) for read in off_target] == [("off_target", "q0", 0, 420)]
