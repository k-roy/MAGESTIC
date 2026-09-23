"""Portable before/fix contracts for actual Stage2 artifact provenance.

Scientific expectations are independently specified: manifest B,A; plate A,B;
AD-first multiallelic ALT C,G; B selects G, A remains uncertain despite coverage.
No planning fixture, external caller, network or real biological data is used.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest
import yaml

from magestic.pipelines.wgs.config import load_config
from magestic.pipelines.wgs.pipeline import step_2_donor_precall as stage
from magestic.utils import provenance as shared_provenance


CHROM = "integrated_plasmid_oligo_1"
WT = "ACGTACGTAA"
DESIGN = "ACGTACGTAG"
REFERENCE = "N"*10 + WT + "N"*10
SAMPLES = ("sample_B", "sample_A")


def receipt(path):
    return Path(str(path)+".provenance.json")


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def outcome_rows(path):
    # A real TRACE header begins '#sample_name'; only the axis stamp is skipped.
    with Path(path).open() as handle:
        return list(csv.DictReader((line for line in handle if not line.startswith("# synthetic_error_axis=")), delimiter="\t"))


@pytest.fixture
def cohort(tmp_path):
    reference = tmp_path/"reference.fa"
    reference.write_text(f">{CHROM}\n{REFERENCE}\n")
    pysam.faidx(str(reference))
    bams = {}
    for sample in SAMPLES:
        bam = tmp_path/f"{sample}.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": CHROM, "LN": len(REFERENCE),
                          "M5": hashlib.md5(REFERENCE.encode()).hexdigest()}],
                  "RG": [{"ID": sample, "SM": sample}],
                  "PG": [{"ID": "synthetic_fixture", "PN": "synthetic_fixture",
                          "CL": "generate known empty BAM with exact SQ M5"}]}
        with pysam.AlignmentFile(str(bam), "wb", header=header):
            pass
        pysam.index(str(bam))
        bams[sample] = bam
    manifest = tmp_path/"manifest.tsv"
    manifest.write_text("sample\tbam\treference\n"+"".join(
        f"{sample}\t{bams[sample]}\t{reference}\n" for sample in SAMPLES))
    donors = tmp_path/"donors.tsv"
    donors.write_text("oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\tdonor_strand\treference\n"
                      f"oligo\t{CHROM}\t11\t20\t{WT}\t{DESIGN}\t+\t{reference}\n")
    assignments = tmp_path/"assignments.tsv"
    assignments.write_text("sample\toligo_name\n"+"".join(f"{sample}\toligo\n" for sample in SAMPLES))
    plate = tmp_path/"plate.vcf"
    plate.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A\tsample_B\n"
                     f"{CHROM}\t20\t.\tA\tC,G\t50\tPASS\t.\tAD:GT:DP\t9,0,0:0/.:9\t0,0,9:2/2:9\n")
    config_path = tmp_path/"config.yaml"
    config_path.write_text(yaml.safe_dump({
        "cohort_name": "stage2_provenance_contract", "sequencing_date": "20260920",
        "paths": {"project_root": str(tmp_path), "manifest": str(manifest)},
        "donor_precall": {"donor_table": str(donors), "plate_vcfs": [str(plate)],
                          "synthetic_donor_source": "design_donor"},
        "reconcile": {"oligo_assignments": str(assignments)}}))
    config = load_config(config_path)
    gvcfs = {}
    for sample in SAMPLES:
        path = config.paths.donor_precall_workdir/"gvcf"/f"{sample}.g.vcf"
        path.parent.mkdir(parents=True, exist_ok=True)
        text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sample+"\n"
        if sample == "sample_A":
            text += f"{CHROM}\t11\t.\tA\t<NON_REF>\t.\tPASS\tEND=20\tGT:MIN_DP\t0:9\n"
        path.write_text(text)
        gvcfs[sample] = path
    return SimpleNamespace(config=config, manifest=manifest, donors=donors,
                           assignments=assignments, plate=plate, reference=reference,
                           bams=bams, gvcfs=gvcfs)


def artifact(case, sample, kind):
    if kind == "shard":
        return stage.shard_path(case.config, sample)
    if kind == "sheet":
        return stage.samplesheet_path(case.config, sample)
    if kind == "split":
        return case.config.paths.donor_precall_workdir/"vcf"/f"{sample}.fil.vcf"
    raise AssertionError(kind)


def assert_receipt(path):
    path = Path(path)
    metadata = receipt(path)
    assert metadata.is_file(), f"producer receipt missing for {path}"
    record = json.loads(metadata.read_text())
    assert record["schema"] == "magestic-step-provenance/1"
    assert record["pipeline"] == "wgs"
    assert Path(record["output"]).resolve() == path.resolve()
    own = [row for row in record["outputs"] if Path(row["path"]).resolve() == path.resolve()]
    assert len(own) == 1
    assert own[0]["md5"] == md5(path)
    assert own[0]["bytes"] == path.stat().st_size
    assert record["code"]["commit"]
    assert record["counts"], "a receipt must state the relevant output denominator"
    return record


def assert_inputs(record, required):
    observed = {Path(row["path"]).resolve(): row for row in record["inputs"]}
    required = {Path(path).resolve() for path in required}
    assert required <= set(observed), f"missing source input roles: {sorted(str(p) for p in required-set(observed))}"
    for path in required:
        row = observed[path]
        assert row.get("role"), f"input has no role: {path}"
        assert row["md5"] == md5(path)
        assert row["bytes"] == path.stat().st_size


def assert_science(path):
    rows = outcome_rows(path)
    assert [row["#sample_name"] for row in rows] == list(SAMPLES)
    assert rows[0]["editing_outcome"] == "designed_variant_installed"
    assert rows[0]["designed_variants"] == f"{CHROM}_20_A_to_G"
    assert rows[1]["editing_outcome"] == "unknown;ambiguous_genotype"
    assert rows[1]["designed_variants"] == ""
    assert json.loads(rows[1]["ambiguous_genotypes"]) == [{"chrom": CHROM, "pos": 20, "GT": "0/."}]
    assert all(row["synthetic_errors"] == row["spontaneous_variants"] == "" for row in rows)


def test_existing_coupled_scientific_contract_is_preserved(cohort):
    assert_science(stage.run_step_2(cohort.config))


def test_every_stage2_owned_data_artifact_has_its_own_receipt(cohort):
    primary = stage.run_step_2(cohort.config)
    owned = [primary]+[artifact(cohort, sample, kind) for sample in SAMPLES for kind in ("shard", "sheet", "split")]
    assert len(owned) == 7
    missing = [str(path) for path in owned if not receipt(path).is_file()]
    assert missing == [], f"{len(missing)}/7 owned artifacts lack receipts: {missing}"
    for path in owned:
        assert_receipt(path)


def test_cohort_receipt_has_source_closure_and_contributing_shards(cohort):
    primary = stage.run_step_2(cohort.config)
    required = {cohort.manifest, cohort.donors, cohort.assignments, cohort.plate, cohort.reference,
                *cohort.gvcfs.values(),
                *(artifact(cohort, sample, "split") for sample in SAMPLES),
                *(artifact(cohort, sample, "shard") for sample in SAMPLES)}
    assert len(required) == 11
    assert_inputs(assert_receipt(primary), required)


@pytest.mark.parametrize("kind", ["split", "sheet", "shard"])
def test_each_producer_receipt_records_its_actual_input_roles(cohort, kind):
    stage.run_step_2(cohort.config)
    sample = "sample_A"
    required = {
        "split": {cohort.plate},
        "sheet": {cohort.manifest, cohort.donors, cohort.assignments},
        "shard": {artifact(cohort, sample, "sheet"), artifact(cohort, sample, "split"), cohort.gvcfs[sample], cohort.reference},
    }[kind]
    assert_inputs(assert_receipt(artifact(cohort, sample, kind)), required)


@pytest.mark.parametrize("kind", ["split", "sheet", "shard"])
@pytest.mark.parametrize("damage", ["missing", "foreign_content"])
def test_receipt_damage_invalidates_only_affected_sample_completion(cohort, monkeypatch, kind, damage):
    stage.run_step_2(cohort.config)
    target = artifact(cohort, "sample_A", kind)
    if damage == "missing":
        receipt(target).unlink(missing_ok=True)
    else:
        receipt(target).write_text(json.dumps({"schema": "magestic-step-provenance/1", "pipeline": "wgs",
                                              "output": "foreign_artifact.tsv", "outputs": []}))
    calls = []
    original = stage.donor_outcome.precall
    def counted(*args, **kwargs):
        calls.extend(sorted(kwargs["samples"]))
        return original(*args, **kwargs)
    monkeypatch.setattr(stage.donor_outcome, "precall", counted)
    primary = stage.run_step_2(cohort.config)
    assert calls == ["sample_A"], "receipt damage must not be skipped, repaired post-hoc, or invalidate unaffected sample_B"
    assert_receipt(target)
    assert_science(primary)


def test_receipt_publication_failure_cannot_commit_sample_success(cohort, monkeypatch):
    target = artifact(cohort, "sample_A", "shard")
    original = shared_provenance.write_json_atomic
    def reject_target(record, destination):
        if Path(destination).resolve() == receipt(target).resolve():
            raise OSError("injected sample receipt publication failure")
        return original(record, destination)
    with monkeypatch.context() as patch:
        patch.setattr(shared_provenance, "write_json_atomic", reject_target)
        with pytest.raises(OSError, match="sample receipt publication failure"):
            stage.run_step_2(cohort.config)
    assert not Path(str(target)+".cache.json").exists()
    assert not stage.table_path(cohort.config).exists()
    primary = stage.run_step_2(cohort.config)
    assert_science(primary)
    assert_receipt(target)


def test_neutral_donor_file_change_refreshes_input_identity_without_changing_calls(cohort):
    primary = stage.run_step_2(cohort.config)
    original_data, original_hash = primary.read_bytes(), md5(cohort.donors)
    # A blank physical line changes the source artifact without changing a donor.
    cohort.donors.write_text(cohort.donors.read_text()+"\n")
    assert md5(cohort.donors) != original_hash
    primary = stage.run_step_2(cohort.config)
    assert primary.read_bytes() == original_data
    assert_inputs(assert_receipt(primary), {cohort.donors})
    for sample in SAMPLES:
        assert_inputs(assert_receipt(artifact(cohort, sample, "sheet")), {cohort.donors})
