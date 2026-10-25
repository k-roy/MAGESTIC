"""A rescued minority allele reports as a mixture, never as unedited (planning/268 §2e).

Kevin, 2026-10-24, agreed option 2. A haploid caller writes a 10-30% allele as `GT 0`, so
the opt-in allele-fraction rescue is the only thing that can see it. The pinned TRACE
selects alleles by GT, so a rescued record left at `GT 0` selects no ALT, and the well
read `designed_variant_unedited`. The rescue now writes the heterozygous genotype the
read counts show, and keeps the original in `MAG_RESCUED_GT`. TRACE then refuses the call
as `unknown;ambiguous_genotype`, and step_5 reports a well whose ambiguous donor-window
genotype is itself a `MAG_MINOR` allele as `mixture_signal`. Nothing here asserts an edit.
"""
import json
from pathlib import Path

from trace_crispr.wgs.donor_outcome import DonorSpec, precall_sample

from magestic.pipelines.wgs.core.vcf_adapter import (
    MAG_MINOR_TAG,
    MAG_RESCUED_GT_TAG,
    RESCUED_GT,
    split_plate_vcf,
)
from magestic.pipelines.wgs.pipeline import step_2_donor_precall as s2
from magestic.pipelines.wgs.pipeline import step_4_reconcile as s4
from magestic.pipelines.wgs.pipeline import step_5_join_outcomes as s5

from . import test_255_blocker_fixes as t255

WT = "ACGTACGTACGTACGTACGTA"  # chrXV:1000-1020; G at 1010
DESIGN = WT[:10] + "T" + WT[11:]  # the design is G>T at 1010
CHROM_SEQ = "N" * 999 + WT + "N" * 200


def _plate(path: Path, gt: str, ad: str) -> Path:
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chrXV>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        f"chrXV\t1010\t.\tG\tT\t50\t.\t.\tGT:AD\t{gt}:{ad}\n")
    return path


def _gvcf(path: Path) -> Path:
    rows = [f"chrXV\t{p}\t.\t{CHROM_SEQ[p - 1]}\t<NON_REF>\t.\t.\t.\tGT:AD\t0/0:30,0"
            for p in range(1000, 1021)]
    path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                    "\tFORMAT\ts1\n" + "\n".join(rows) + "\n")
    return path


def _record(vcf: Path) -> list[str]:
    data = [ln.split("\t") for ln in vcf.read_text().splitlines() if not ln.startswith("#")]
    return data[0] if data else []


def test_the_rescue_rewrites_gt_to_heterozygous_and_keeps_the_original(tmp_path):
    plate = _plate(tmp_path / "p.vcf", "0", "27,3")  # a 10% allele, haploid-called
    off = split_plate_vcf(plate, tmp_path / "off")
    assert _record(Path(off.vcf_paths["s1"])) == []  # GT 0, so dropped: the gap itself
    on = split_plate_vcf(plate, tmp_path / "on", min_allele_fraction=0.05)
    rec = _record(Path(on.vcf_paths["s1"]))
    fmt = dict(zip(rec[8].split(":"), rec[9].split(":"), strict=True))
    assert fmt["GT"] == RESCUED_GT == "0/1"
    info = rec[7].split(";")
    assert f"{MAG_RESCUED_GT_TAG}=0" in info and MAG_MINOR_TAG in info
    assert on.minor_allele_variants["s1"] == ("chrXV_1010_G_to_T",)


def test_a_rescued_minority_allele_is_ambiguous_not_unedited(tmp_path):
    """BEFORE planning/268 §2e this was `designed_variant_unedited`: a false 'no edit'."""
    res = split_plate_vcf(_plate(tmp_path / "p.vcf", "0", "21,9"), tmp_path / "wd",
                          min_allele_fraction=0.05)
    spec = DonorSpec("s1", "chrXV", 1000, 1020, WT, DESIGN, DESIGN)
    out = precall_sample(spec, Path(res.vcf_paths["s1"]), _gvcf(tmp_path / "g.vcf"), CHROM_SEQ)
    assert out.editing_outcome == "unknown;ambiguous_genotype"
    assert {"chrom": "chrXV", "pos": 1010, "GT": "0/1"} in out.ambiguous_genotypes


def _row(outcome, minor, ambiguous):
    return {"editing_outcome": outcome, "designed_variants": "", "minor_allele_variants": minor,
            "ambiguous_genotypes": json.dumps(ambiguous)}


def test_the_mixture_rule_needs_the_minor_allele_AT_the_ambiguous_position():
    at = [{"chrom": "chrXV", "pos": 1010, "GT": "0/1"}]
    hit = s5._minor_allele_at_ambiguous_genotype
    assert hit(_row("unknown;ambiguous_genotype", "chrXV_1010_G_to_T", at))
    # an underscore-bearing contig parses from the right
    contig = "integrated_plasmid_o1"
    assert hit(_row("unknown;ambiguous_genotype", f"{contig}_1010_G_to_T",
                    [{"chrom": contig, "pos": 1010, "GT": "0/1"}]))
    # NEGATIVE controls: a minor allele elsewhere, an unresolved ALT, a called outcome
    assert not hit(_row("unknown;ambiguous_genotype", "chrIV_9999_C_to_T", at))
    assert not hit(_row("unknown;ambiguous_genotype", "chrXV_1010_G_to_T",
                        [{"chrom": "chrXV", "pos": 1010, "GT": "1", "reason": "unresolved_allele"}]))
    assert not hit(_row("designed_variant_unedited", "chrXV_1010_G_to_T", at))
    assert not hit(_row("unknown;ambiguous_genotype", "", at))


def test_step_5_reports_it_as_mixture_signal_with_a_reason(tmp_path):
    config = t255._config(tmp_path, samples=("s1", "s2"))
    t255._stage_4_audited_table(s4.table_path(config))
    cols = (["sample_name", "editing_outcome", "designed_variants", "synthetic_errors",
             "spontaneous_variants", "ambiguous_genotypes"] + list(s2.OUTCOME_REASON_COLUMNS))
    amb = json.dumps([{"chrom": "chrXV", "pos": 1010, "GT": "0/1"}])
    rows = [
        ["s1", "unknown;ambiguous_genotype", "", "", "", amb,
         "false", "uncallable", "10", "10", "0", "chrXV_1010_G_to_T"],
        ["s2", "unknown;ambiguous_genotype", "", "", "", amb,
         "false", "uncallable", "10", "10", "0", ""],
    ]
    s2.table_path(config).write_text("# synthetic_error_axis=disabled\n#" + "\t".join(cols) + "\n"
                                     + "\n".join("\t".join(r) for r in rows) + "\n")
    written = s5.run_step_5(config)
    out = [ln.split("\t") for ln in written[0].read_text().splitlines()]
    header, data = out[0], out[1:]
    state = header.index("precall_state")
    reason = header.index("precall_state_reason")
    assert (data[0][state], data[0][reason]) == ("mixture_signal",
                                                 s5.PRECALL_STATE_MINOR_ALLELE_AMBIGUOUS)
    assert (data[1][state], data[1][reason]) == ("undetermined", s5.PRECALL_STATE_FROM_OUTCOME)
