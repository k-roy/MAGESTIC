"""An externally supplied per-sample VCF is a Stage2 input, never a Stage2 product.

With `donor_precall.plate_vcfs: []` the `.fil.vcf` files already exist (the route the
empirical v11 harness uses). Stage2 must hash them as consumed inputs, must not
write or demand a producer receipt beside them, and must still reuse a completed
sample. The external files here are made by the real plate adapter before Stage2
runs, so their content is exactly what the plate route would consume. Outcome
expectations are the same independently specified pair as the provenance contract:
sample_B installs the designed A>G; sample_A stays unknown on a 0/. genotype.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pysam
import pytest
import yaml

from magestic.pipelines.wgs.config import load_config
from magestic.pipelines.wgs.core.vcf_adapter import split_plate_vcf
from magestic.pipelines.wgs.pipeline import step_2_donor_precall as stage

CHROM = "integrated_plasmid_oligo_1"
WT, DESIGN = "ACGTACGTAA", "ACGTACGTAG"
REFERENCE = "N" * 10 + WT + "N" * 10
SAMPLES = ("sample_B", "sample_A")


def receipt(path):
    return Path(str(path) + ".provenance.json")


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


@pytest.fixture
def external(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(f">{CHROM}\n{REFERENCE}\n")
    pysam.faidx(str(reference))
    rows = []
    for sample in SAMPLES:
        bam = tmp_path / f"{sample}.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": CHROM, "LN": len(REFERENCE), "M5": hashlib.md5(REFERENCE.encode()).hexdigest()}],
                  "RG": [{"ID": sample, "SM": sample}],
                  "PG": [{"ID": "synthetic_fixture", "PN": "synthetic_fixture",
                          "CL": "generate known empty BAM with exact SQ M5"}]}
        with pysam.AlignmentFile(str(bam), "wb", header=header):
            pass
        pysam.index(str(bam))
        rows.append(f"{sample}\t{bam}\t{reference}\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("sample\tbam\treference\n" + "".join(rows))
    donors = tmp_path / "donors.tsv"
    donors.write_text("oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\tdonor_strand\treference\n"
                      f"oligo\t{CHROM}\t11\t20\t{WT}\t{DESIGN}\t+\t{reference}\n")
    assignments = tmp_path / "assignments.tsv"
    assignments.write_text("sample\toligo_name\n" + "".join(f"{s}\toligo\n" for s in SAMPLES))
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump({
        "cohort_name": "stage2_external_vcf", "sequencing_date": "20260921",
        "paths": {"project_root": str(tmp_path), "manifest": str(manifest)},
        "donor_precall": {"donor_table": str(donors), "plate_vcfs": [],
                          "synthetic_donor_source": "design_donor"},
        "reconcile": {"oligo_assignments": str(assignments)}}))
    config = load_config(config_path)
    workdir = config.paths.donor_precall_workdir
    # Produce the external per-sample VCFs with the real adapter OUTSIDE Stage2, from a
    # plate file that the Stage2 config never names.
    plate = tmp_path / "outside_plate.vcf"
    plate.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\n"
                     f"{CHROM}\t20\t.\tA\tC,G\t50\tPASS\t.\tAD:GT:DP\t9,0,0:0/.:9\t0,0,9:2/2:9\n")
    split_plate_vcf(plate, workdir, samples=list(SAMPLES), ad_slot=config.donor_precall.ad_format_slot_out)
    for sample in SAMPLES:
        gvcf = workdir / "gvcf" / f"{sample}.g.vcf"
        gvcf.parent.mkdir(parents=True, exist_ok=True)
        text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n"
        if sample == "sample_A":
            text += f"{CHROM}\t11\t.\tA\t<NON_REF>\t.\tPASS\tEND=20\tGT:MIN_DP\t0:9\n"
        gvcf.write_text(text)
    vcfs = {s: workdir / "vcf" / f"{s}.fil.vcf" for s in SAMPLES}
    for path in vcfs.values():
        # Anything beside an external input is not Stage2's; start from none.
        assert not receipt(path).exists()
    return config, vcfs


def outcomes(path):
    with Path(path).open() as handle:
        rows = list(csv.DictReader((l for l in handle if not l.startswith("# synthetic_error_axis=")), delimiter="\t"))
    return [(r["#sample_name"], r["editing_outcome"], r["designed_variants"]) for r in rows]


EXPECTED = [("sample_B", "designed_variant_installed", f"{CHROM}_20_A_to_G"),
            ("sample_A", "unknown;ambiguous_genotype", "")]


def test_external_vcf_is_consumed_unchanged_and_never_claimed_as_a_product(external):
    config, vcfs = external
    before = {s: p.read_bytes() for s, p in vcfs.items()}
    assert outcomes(stage.run_step_2(config)) == EXPECTED
    for sample, path in vcfs.items():
        assert path.read_bytes() == before[sample], "Stage2 must not rewrite an external input"
        assert not receipt(path).exists(), "Stage2 must not write a producer receipt for an input it did not produce"


def test_external_vcf_is_hashed_as_an_input_of_sample_and_cohort_receipts(external):
    config, vcfs = external
    primary = stage.run_step_2(config)
    cohort = json.loads(receipt(primary).read_text())
    cohort_inputs = {Path(r["path"]).resolve(): r for r in cohort["inputs"]}
    for sample, vcf in vcfs.items():
        shard = json.loads(receipt(stage.shard_path(config, sample)).read_text())
        shard_inputs = {Path(r["path"]).resolve(): r for r in shard["inputs"]}
        for inputs in (shard_inputs, cohort_inputs):
            row = inputs[vcf.resolve()]
            assert row.get("role") and row["md5"] == md5(vcf) and row["bytes"] == vcf.stat().st_size


def test_completed_external_route_is_reused_without_a_split_receipt(external, monkeypatch):
    config, _ = external
    stage.run_step_2(config)
    calls = []
    original = stage.donor_outcome.precall
    def counted(*args, **kwargs):
        calls.extend(sorted(kwargs["samples"]))
        return original(*args, **kwargs)
    monkeypatch.setattr(stage.donor_outcome, "precall", counted)
    assert outcomes(stage.run_step_2(config)) == EXPECTED
    assert calls == [], "a receipt Stage2 never owns must not be demanded for completion"


def test_changed_external_vcf_recomputes_only_its_sample(external, monkeypatch):
    # Positive control for the existing input-identity cache.
    config, vcfs = external
    stage.run_step_2(config)
    vcfs["sample_A"].write_text(vcfs["sample_A"].read_text())  # new identity, same content
    calls = []
    original = stage.donor_outcome.precall
    def counted(*args, **kwargs):
        calls.extend(sorted(kwargs["samples"]))
        return original(*args, **kwargs)
    monkeypatch.setattr(stage.donor_outcome, "precall", counted)
    assert outcomes(stage.run_step_2(config)) == EXPECTED
    assert calls == ["sample_A"]
