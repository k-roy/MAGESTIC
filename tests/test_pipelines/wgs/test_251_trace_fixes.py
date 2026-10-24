"""NN 251 — one reproduction + one fix assertion per defect the WGS trace audits found.

METHOD (non-negotiable, and the reason this file is structured in pairs): every defect
gets a ``test_repro_*`` that drives an **existing** entry point and asserts the CORRECT
result, so it fails on wrong *output* rather than on a missing import; and a
``test_fixed_*`` that asserts the adapter's behaviour. A defect with no reproduction is
recorded as ``NOT REPRODUCED`` in `planning/222_triumph_embl_library_design/251_wgs_trace_fixes.md`
and is not "fixed".

Inputs: `222b_audit_variant_processing.md`, `234_magestic_bug_fixes.md`,
`243_wgs_trace_adversarial.md` (+ its runnable repros). TRACE is pinned at `7c4c10f1`
and is NOT edited; each defect that lives in TRACE is closed magestic-side by an adapter
and the TRACE-side change is specified in the planning doc.

Everything here is synthetic. No Oak/Sherlock path appears in this file.

Author: NN 251 (wgs-fixes)
"""

from __future__ import annotations

import math
import random
from pathlib import Path

import pysam
import pytest

from magestic.pipelines.wgs.core import _trace

TR = _trace.reconcile()
DO = _trace.donor_outcome() if hasattr(_trace, "donor_outcome") else None

# The five composite contigs the WGS reference-mapping policy mandates. Four contain
# an underscore, which is defect 3.
POLICY_CONTIGS = (
    "chrXV",
    "unplaced_contig_34",
    "blank_MAGESTIC3_pF835_guide_donor_plasmid",
    "native_2micron",
    "integrated_plasmid_o1",
)

CONTIG = "chrXV"
UNDERSCORE_CONTIG = "blank_MAGESTIC3_pF835_guide_donor_plasmid"


def defect_persists(defect: str, *, closed_by: str):
    """Mark a REPRODUCTION of a defect that still lives inside pinned TRACE.

    `strict=True` is the whole point. The test keeps running, so the defect stays
    documented by executable evidence rather than by prose; it reports `xfail`, so
    the suite is honest about it being a known defect rather than a regression; and
    if TRACE is ever upgraded such that the defect disappears, the test **XPASSes
    and that is a FAILURE**, which forces this file and the planning doc to be
    updated instead of silently drifting.

    🔴 `raises=AssertionError` is as load-bearing as `strict`. Without it, a repro that
    died with `ImportError`, `AttributeError` or `TypeError` would still count as
    "defect reproduced" — so a renamed TRACE symbol would masquerade as a persisting
    defect, which is precisely the fail-on-wrong-output rule this task was given.

    `closed_by` names the magestic-side adapter that closes it for this pipeline.
    """
    return pytest.mark.xfail(
        strict=True,
        raises=AssertionError,
        reason=(
            f"defect {defect} is in TRACE @7c4c10f1 (pinned, not edited); closed "
            f"magestic-side by {closed_by}. TRACE-side fix specified in that module's "
            "docstring and in planning/222_triumph_embl_library_design/251_wgs_trace_fixes.md"
        ),
    )


# --------------------------------------------------------------------------- #
# Synthetic fixture builders (adapted from 243_wgs_trace_adversarial_repro)     #
# --------------------------------------------------------------------------- #


def _seq(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _other(b: str) -> str:
    return {"A": "G", "G": "A", "C": "T", "T": "C", "N": "A"}[b.upper()]


def write_fasta(path: Path, contigs: dict[str, str]) -> Path:
    with path.open("w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")
    pysam.faidx(str(path))
    return path


def write_bam(
    path: Path,
    truth: str,
    *,
    contig: str = CONTIG,
    spans=((3001, 3120),),
    edits: dict[int, str] | None = None,
    deletion: tuple[int, int] | None = None,
    skip_positions: tuple[int, ...] = (),
    read_len: int = 60,
    depth: int = 14,
) -> Path:
    """A coordinate-sorted, indexed BAM of `read_len` reads tiling each span.

    `skip_positions` drops any read that would cover one of those positions, which is
    how the A10/A12 "window maximum hides zero depth at the edit" case is built.
    `deletion` applies a CIGAR deletion `(start, length)` to every covering read.
    """
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": len(truth)}]}
    rng = random.Random(11)
    recs = []
    n = 0
    for ws, we in spans:
        starts = [rng.randint(max(1, ws - read_len + 1), we) for _ in range(depth * 4)]
        for k, st in enumerate(starts):
            st = max(1, min(st, len(truth) - read_len + 1))
            if any(st <= p < st + read_len for p in skip_positions):
                continue
            s = list(truth[st - 1 : st - 1 + read_len])
            for pos, alt in (edits or {}).items():
                if st <= pos < st + read_len:
                    s[pos - st] = alt
            a = pysam.AlignedSegment()
            a.query_name = f"r{n}"
            n += 1
            a.flag = 16 if k % 2 else 0
            a.reference_id = 0
            a.reference_start = st - 1
            a.mapping_quality = 60
            if deletion is not None and st <= deletion[0] < st + read_len - deletion[1]:
                dstart, dlen = deletion
                keep = "".join(s)
                left = dstart - st
                a.query_sequence = (keep[:left] + keep[left + dlen :]).replace("N", "A")
                a.cigartuples = [(0, left), (2, dlen), (0, read_len - dlen - left)]
                a.query_qualities = pysam.qualitystring_to_array("I" * (read_len - dlen))
            else:
                a.query_sequence = "".join(s).replace("N", "A")
                a.cigartuples = [(0, read_len)]
                a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            recs.append(a)
    recs.sort(key=lambda r: r.reference_start)
    tmp = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(tmp), "wb", header=header) as out:
        for r in recs:
            out.write(r)
    pysam.sort("-o", str(path), str(tmp))
    pysam.index(str(path))
    tmp.unlink()
    return path


def write_window_table(path: Path, rows, header=None) -> Path:
    header = header or ("oligo_name", "chrom", "donor_start_coord", "donor_end_coord")
    with path.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path


def write_vcf(path: Path, records, *, contig=CONTIG) -> Path:
    """A single-sample plain-text `.fil.vcf` as TRACE's `_iter_vcf_records` reads it."""
    with path.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={contig}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for chrom, pos, ref, alt in records:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT:AD\t1/1:0,30\n")
    return path


def write_gvcf(path: Path, blocks, *, contig=CONTIG, ref: str | None = None,
               min_dp: bool = True) -> Path:
    """gVCF reference blocks as `(pos, end_or_None, ref_depth)`.

    `end_or_None` writes `INFO END=<end>`, which is exactly how GATK/`bcftools -g`
    emit a *compressed* reference block — the shape defect 5 is about.

    `ref` (the contig sequence) supplies each block's REF base. planning/268: the pinned
    TRACE validates every gVCF REF against the reference, and the former constant `A`
    was a fixture no caller could write. An END block also carries `MIN_DP` unless
    `min_dp=False`, as `bcftools --gvcf` and GATK write it: the pinned TRACE counts a
    block's interior only from `MIN_DP`, its lower bound over the interval.
    """
    with path.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Block end">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for pos, end, depth in blocks:
            info = f"END={end}" if end is not None else "."
            base = ref[pos - 1] if ref else "A"
            fmt, val = ("GT:AD:MIN_DP", f"0/0:{depth},0:{depth}") if end is not None and min_dp \
                else ("GT:AD", f"0/0:{depth},0")
            fh.write(f"{contig}\t{pos}\t.\t{base}\t<NON_REF>\t.\t.\t{info}\t{fmt}\t{val}\n")
    return path


def qwin(dstart: int, dend: int, size: int = 10) -> tuple[int, int]:
    """TRACE's own quantification-window formula, restated so a test can target it."""
    start = math.ceil(dend - (dend - dstart) / 2 - size / 2)
    return start, start + size - 1


# --------------------------------------------------------------------------- #
# Defect 1 — diploid default on haploid yeast; no allele-fraction threshold     #
# --------------------------------------------------------------------------- #


class TestDefect1AlleleFraction:
    """`call -mv` at the diploid default emits `GT 0/0` for a 33% minor allele, so
    `carries_alt()` drops the record and Stage 2 writes a bare `unedited`."""

    @staticmethod
    def _plate(tmp_path: Path) -> Path:
        p = tmp_path / "plate.vcf"
        with p.open("w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
            # S1: minor allele, 10/30 = 33%, called 0/0 by the DIPLOID caller.
            # S2: clean homozygous alt.
            fh.write(
                "chrXV\t3060\t.\tA\tG\t50\t.\t.\tGT:PL:DP:AD"
                "\t0/0:0,30,100:30:20,10\t1/1:100,30,0:30:0,30\n"
            )
            # A record with no AD at all, so the no-AD denominator is exercised.
            fh.write(
                "chrXV\t3070\t.\tC\tT\t50\t.\t.\tGT:PL:DP:AD"
                "\t0/0:0,30,100:30:.\t1/1:100,30,0:30:0,30\n"
            )
        return p

    @defect_persists("1 (diploid GT drops a haploid minor allele)",
                     closed_by="vcf_adapter.split_plate_vcf(min_allele_fraction=...)")
    def test_repro_minor_allele_is_silently_dropped(self, tmp_path):
        """REPRODUCTION: a 33% alt in S1 produces ZERO records for S1."""
        from magestic.pipelines.wgs.core.vcf_adapter import split_plate_vcf

        res = split_plate_vcf(self._plate(tmp_path), tmp_path / "wd")
        assert res.n_records_in == 2
        # S2 is kept by GT; S1 is dropped entirely despite carrying 10 alt reads.
        assert res.n_records_out["S2"] == 2
        assert res.n_records_out["S1"] == 1, (
            "a 33% minor allele with GT 0/0 is dropped, so Stage 2 sees no variant "
            "and emits a bare `unedited` — the silent conversion of a real partial "
            "edit into a no-edit"
        )

    def test_fixed_allele_fraction_rescues_and_flags(self, tmp_path):
        from magestic.pipelines.wgs.core.vcf_adapter import split_plate_vcf

        res = split_plate_vcf(
            self._plate(tmp_path), tmp_path / "wd", min_allele_fraction=0.10
        )
        assert res.n_records_out["S1"] == 1
        assert res.n_rescued_by_af == 1
        assert res.n_no_ad == 1, "the AD-less record must be counted, never silently kept"
        # The rescued allele is 33%, far below a full call: it must be FLAGGED, because
        # TRACE reads only columns 0/1/3/4 and would otherwise treat it as 100%.
        assert res.minor_allele_variants["S1"] == ("chrXV_3060_A_to_G",)
        text = (tmp_path / "wd" / "vcf" / "S1.fil.vcf").read_text()
        assert "MAG_AF=0.333" in text
        assert "MAG_MINOR" in text

    def test_fixed_default_is_unchanged(self, tmp_path):
        """Parity: without `min_allele_fraction` nothing moves (GT-only, as today)."""
        from magestic.pipelines.wgs.core.vcf_adapter import split_plate_vcf

        res = split_plate_vcf(self._plate(tmp_path), tmp_path / "wd")
        assert res.n_records_out["S1"] == 0
        assert res.n_rescued_by_af == 0

    def test_fixed_ploidy_one_is_owned_and_stated(self):
        """The mpileup/call command is owned here, at `--ploidy 1`, not left implicit."""
        from magestic.pipelines.wgs.core import vcf_adapter as va

        assert "--ploidy 1" in va.MPILEUP_CALL_COMMAND
        assert va.CALLER_PLOIDY == 1
        # And the limitation is stated in the module, not only in a planning doc.
        #
        # 🔴 UPDATED by NN 255: this asserted the word "untested", which was true only
        # while nobody had run the binary. NN 252 F13 ran it (bcftools 1.21, exit 0),
        # so the caveat now carries MEASURED floors and the test pins those.
        assert "MEASURED against bcftools 1.21" in va.PLOIDY_CAVEAT
        assert "private" in va.PLOIDY_CAVEAT.lower()


# --------------------------------------------------------------------------- #
# Defect 2 — the blacklist is silently inert on a contig-name mismatch          #
# --------------------------------------------------------------------------- #


class TestDefect2Blacklist:
    @staticmethod
    def _blacklist(tmp_path: Path) -> Path:
        p = tmp_path / "REDI_reads_mapped_regions.txt"
        p.write_text(
            "chr15\t100\t200\tCUP1\n"
            "chr1\t300\t400\tmarker\n"
            "yT177_after_bc1_integration 500 600\n"  # space separated, not tabs
        )
        return p

    def test_repro_all_names_absent_parses_clean(self, tmp_path):
        """REPRODUCTION: TRACE's loader accepts a 100%-inert blacklist without a word."""
        sv = _trace.sv_plasmid()
        # planning/268: the pinned TRACE refuses the space-separated row in its standard
        # mode, which closes that half of the defect upstream. The Perl-port parity mode
        # still parses it, and two well-formed rows naming absent contigs parse in either
        # mode, so the magestic guard below is still the only refusal of an inert list.
        with pytest.raises(ValueError, match="expected chrom<TAB>start<TAB>end"):
            sv.load_blacklist(str(self._blacklist(tmp_path)))
        bl = sv.load_blacklist(str(self._blacklist(tmp_path)), mode="parity")
        reference_contigs = {"chrXV": 10000, "chrI": 5000}
        present = [c for c in bl if c in reference_contigs]
        assert len(bl) == 3  # parsed happily
        assert present == [], "every blacklist name is absent from the reference"
        # Nothing in TRACE's own path turns that into a failure. That is the bug: the
        # run looks filtered and is not. The guard must raise its OWN error class --
        # a bare `pytest.raises(Exception)` would be satisfied by `ModuleNotFoundError`
        # and so would pass before the fix existed, proving nothing.
        from magestic.pipelines.wgs.core.blacklist_adapter import (
            BlacklistInapplicableError,
            require_blacklist_applicable,
        )

        with pytest.raises(BlacklistInapplicableError):
            require_blacklist_applicable(
                str(self._blacklist(tmp_path)), reference_contigs
            )

    def test_fixed_audit_counts_every_reject(self, tmp_path):
        from magestic.pipelines.wgs.core.blacklist_adapter import audit_blacklist

        audit = audit_blacklist(
            str(self._blacklist(tmp_path)), {"chrXV": 10000, "chrI": 5000}
        )
        assert audit.n_rows == 3
        assert audit.n_absent_contigs == 3
        assert audit.n_present_contigs == 0
        assert audit.n_space_separated == 1
        assert audit.n_malformed == 1  # the space-separated row has no usable interval
        assert audit.fraction_applicable == 0.0
        assert "chr15" in audit.absent_contigs

    def test_fixed_hard_check_fails_loudly(self, tmp_path):
        from magestic.pipelines.wgs.core.blacklist_adapter import (
            BlacklistInapplicableError,
            require_blacklist_applicable,
        )

        with pytest.raises(BlacklistInapplicableError) as exc:
            require_blacklist_applicable(
                str(self._blacklist(tmp_path)), {"chrXV": 10000, "chrI": 5000}
            )
        assert "0 of 3" in str(exc.value) or "0/3" in str(exc.value)

    def test_fixed_out_of_bounds_interval_is_reported(self, tmp_path):
        """A same-name blacklist in the WRONG coordinate frame (the PRE-coords copy)."""
        from magestic.pipelines.wgs.core.blacklist_adapter import audit_blacklist

        p = tmp_path / "bl2.txt"
        p.write_text("chrXV\t100\t200\tok\nchrXV\t99000\t99500\tbeyond_the_contig\n")
        audit = audit_blacklist(str(p), {"chrXV": 10000})
        assert audit.n_out_of_bounds == 1
        assert audit.n_present_contigs == 1

    def test_fixed_a_matching_blacklist_passes(self, tmp_path):
        from magestic.pipelines.wgs.core.blacklist_adapter import (
            require_blacklist_applicable,
        )

        p = tmp_path / "bl3.txt"
        p.write_text("chrXV\t100\t200\tok\n")
        audit = require_blacklist_applicable(str(p), {"chrXV": 10000})
        assert audit.fraction_applicable == 1.0


# --------------------------------------------------------------------------- #
# Defect 3 — split("_") on contig names containing underscores                  #
# --------------------------------------------------------------------------- #


class TestDefect3UnderscoreContigs:
    @staticmethod
    def _case(tmp_path: Path, contig: str):
        ref = _seq(5, 600)
        dstart, dend = 101, 160
        wt = ref[dstart - 1 : dend]
        local = 30
        alt = _other(wt[local])
        design = wt[:local] + alt + wt[local + 1 :]
        fasta = write_fasta(tmp_path / f"{contig}.fa", {contig: ref})
        vcf = write_vcf(tmp_path / "s.fil.vcf", [(contig, dstart + local, wt[local], alt)],
                        contig=contig)
        qs, qe = qwin(dstart, dend)
        gvcf = write_gvcf(
            tmp_path / "s.g.vcf", [(p, None, 30) for p in range(qs, qe + 1)], contig=contig,
            ref=ref,
        )
        spec = TR and None
        return ref, dstart, dend, wt, design, fasta, vcf, gvcf, contig

    @defect_persists("3 (positional split('_') on contig names)",
                     closed_by="donor_precall.precall_sample_aliased + ContigAlias")
    def test_repro_precall_raises_on_every_policy_contig(self, tmp_path):
        """REPRODUCTION: all four underscore-bearing mandated contigs crash Stage 2."""
        do = _trace.donor_outcome()
        broken = []
        for contig in POLICY_CONTIGS:
            d = tmp_path / contig.replace("/", "_")
            d.mkdir()
            ref, dstart, dend, wt, design, fasta, vcf, gvcf, _ = self._case(d, contig)
            spec = do.DonorSpec(
                sample="S1", chrom=contig, donor_start=dstart, donor_end=dend,
                donor_wt=wt, donor_design=design, donor_step2=design,
            )
            try:
                do.precall_sample(spec, vcf, gvcf, ref)
            except ValueError:
                broken.append(contig)
        assert broken == [], (
            f"{len(broken)} of {len(POLICY_CONTIGS)} mandated composite contigs crash "
            f"`precall_sample`'s positional split('_'): {broken}"
        )

    def test_fixed_underscore_contig_is_called_with_real_names_out(self, tmp_path):
        """The adapter aliases on the way in and DECODES on the way out."""
        from magestic.pipelines.wgs.core.donor_precall import precall_sample_aliased

        do = _trace.donor_outcome()
        ref, dstart, dend, wt, design, fasta, vcf, gvcf, contig = self._case(
            tmp_path, UNDERSCORE_CONTIG
        )
        spec = do.DonorSpec(
            sample="S1", chrom=contig, donor_start=dstart, donor_end=dend,
            donor_wt=wt, donor_design=design, donor_step2=design,
        )
        res = precall_sample_aliased(spec, vcf, gvcf, ref, contigs=POLICY_CONTIGS)
        assert res.outcome.editing_outcome.startswith("designed_variant_installed")
        # 🔴 The shipped table must carry the REAL contig name, never `ctg0001`.
        assert len(res.outcome.designed_variants) == 1
        shipped = res.outcome.designed_variants[0]
        assert shipped.startswith(UNDERSCORE_CONTIG + "_"), shipped
        assert "ctg" not in shipped

    def test_fixed_no_aliased_name_reaches_a_shipped_table(self, tmp_path):
        """Every field of every shipped row round-trips to a real contig name."""
        from magestic.pipelines.wgs.core.contig_alias import ALIAS_PREFIX
        from magestic.pipelines.wgs.core.donor_precall import precall_sample_aliased

        do = _trace.donor_outcome()
        for contig in POLICY_CONTIGS:
            d = tmp_path / ("c_" + contig)
            d.mkdir()
            ref, dstart, dend, wt, design, fasta, vcf, gvcf, _ = self._case(d, contig)
            spec = do.DonorSpec(
                sample="S1", chrom=contig, donor_start=dstart, donor_end=dend,
                donor_wt=wt, donor_design=design, donor_step2=design,
            )
            res = precall_sample_aliased(spec, vcf, gvcf, ref, contigs=POLICY_CONTIGS)
            row = do.outcome_to_row(res.outcome) if hasattr(do, "outcome_to_row") else [
                res.outcome.sample, res.outcome.editing_outcome,
                ";".join(res.outcome.designed_variants),
                ";".join(res.outcome.synthetic_errors),
                ";".join(res.outcome.spontaneous_variants),
            ]
            joined = "\t".join(str(x) for x in row)
            assert ALIAS_PREFIX not in joined, f"aliased name shipped for {contig}: {joined}"
            assert contig in joined

    def test_fixed_decode_variant_id_handles_lc_prefix(self):
        from magestic.pipelines.wgs.core.contig_alias import ContigAlias

        alias = ContigAlias.from_names(POLICY_CONTIGS)
        enc = alias.encode(UNDERSCORE_CONTIG)
        assert alias.decode_variant_id(f"(LC){enc}_131_A_to_G") == (
            f"(LC){UNDERSCORE_CONTIG}_131_A_to_G"
        )

    def test_fixed_aliasing_a_vcf_twice_is_idempotent(self, tmp_path):
        """`step_2` rewrites `{workdir}/vcf/{sample}.fil.vcf` IN PLACE, so a re-run
        after a partial failure aliases an already-aliased file.

        This was flagged as a self-identified regression risk (251 section 7 item 2)
        and is closed here rather than argued: `encode` is identity on an
        already-safe name, so a second pass is a no-op and the file is byte-identical.
        """
        from magestic.pipelines.wgs.core.contig_alias import ContigAlias
        from magestic.pipelines.wgs.core.donor_precall import alias_vcf_contigs

        alias = ContigAlias.from_names(POLICY_CONTIGS)
        src = write_vcf(
            tmp_path / "s.fil.vcf",
            [(UNDERSCORE_CONTIG, 131, "A", "G"), (CONTIG, 200, "C", "T")],
            contig=UNDERSCORE_CONTIG,
        )
        once = tmp_path / "once.vcf"
        twice = tmp_path / "twice.vcf"
        assert alias_vcf_contigs(src, once, alias) == 2
        assert alias_vcf_contigs(once, twice, alias) == 2
        assert once.read_text() == twice.read_text()
        # ...and the aliased names really are underscore-free on the data rows.
        for line in once.read_text().splitlines():
            if not line.startswith("#"):
                assert "_" not in line.split("\t")[0]
        # ...and a third pass still decodes back to the real name.
        got = {
            alias.decode(line.split("\t")[0])
            for line in twice.read_text().splitlines()
            if not line.startswith("#")
        }
        assert got == {UNDERSCORE_CONTIG, CONTIG}


# --------------------------------------------------------------------------- #
# Defects 4 + 5 — `unknown` with no reason; `flag_covered` ignores END          #
# --------------------------------------------------------------------------- #


class TestDefects4And5CoverageAndReason:
    @staticmethod
    def _unedited_case(tmp_path: Path):
        ref = _seq(9, 600)
        dstart, dend = 101, 160
        wt = ref[dstart - 1 : dend]
        design = wt[:30] + _other(wt[30]) + wt[31:]
        vcf = write_vcf(tmp_path / "s.fil.vcf", [])  # no called variants: a true WT cell
        return ref, dstart, dend, wt, design, vcf

    def test_repro_compressed_gvcf_block_yields_unknown(self, tmp_path):
        """REPRODUCTION (defect 5): ONE reference block spanning the whole
        quantification window still gives `unknown`, because only POS is compared.

        CLOSED UPSTREAM (planning/268): the pinned TRACE reads `END` with `MIN_DP`, so the
        correct result is asserted outright and `defect_persists` is removed. A block with
        no `MIN_DP` stays uncovered there by design (see the reason test below)."""
        do = _trace.donor_outcome()
        ref, dstart, dend, wt, design, vcf = self._unedited_case(tmp_path)
        qs, qe = qwin(dstart, dend)
        assert qs > dstart, "the block must start before the window for this to bite"
        gvcf = write_gvcf(tmp_path / "s.g.vcf", [(dstart, dend, 30)], ref=ref)  # END spans qwin
        spec = do.DonorSpec(
            sample="S1", chrom=CONTIG, donor_start=dstart, donor_end=dend,
            donor_wt=wt, donor_design=design, donor_step2=design,
        )
        out = do.precall_sample(spec, vcf, gvcf, ref)
        assert out.editing_outcome == "designed_variant_unedited", (
            "a normal COMPRESSED gVCF whose reference block covers the whole "
            f"quantification window [{qs},{qe}] still reports {out.editing_outcome!r}"
        )

    def test_fixed_the_compressed_block_counts_without_expansion(self, tmp_path):
        """WITHDRAWN (planning/268): "`expand_reference_blocks` makes the block count".
        The expander copies the block's first REF into every position, which the pinned
        TRACE refuses (`gVCF REF does not match reference`), so step_2 no longer calls it.
        The defect-5 claim it served holds on the raw block, and is asserted here."""
        do = _trace.donor_outcome()
        ref, dstart, dend, wt, design, vcf = self._unedited_case(tmp_path)
        gvcf = write_gvcf(tmp_path / "s.g.vcf", [(dstart, dend, 30)], ref=ref)
        spec = do.DonorSpec(
            sample="S1", chrom=CONTIG, donor_start=dstart, donor_end=dend,
            donor_wt=wt, donor_design=design, donor_step2=design,
        )
        out = do.precall_sample(spec, vcf, gvcf, ref)
        assert out.editing_outcome == "designed_variant_unedited"
        from magestic.pipelines.wgs.core.gvcf_adapter import expand_reference_blocks

        expanded = tmp_path / "s.expanded.g.vcf"
        expand_reference_blocks(gvcf, expanded)
        with pytest.raises(ValueError, match="gVCF REF does not match reference"):
            do.precall_sample(spec, vcf, expanded, ref)

    def test_fixed_expansion_is_a_noop_on_an_uncompressed_gvcf(self, tmp_path):
        """Parity: a per-position gVCF passes through byte-equivalently."""
        from magestic.pipelines.wgs.core.gvcf_adapter import expand_reference_blocks

        gvcf = write_gvcf(tmp_path / "s.g.vcf", [(p, None, 30) for p in range(101, 111)])
        out = tmp_path / "out.g.vcf"
        res = expand_reference_blocks(gvcf, out)
        assert res.n_blocks_expanded == 0
        assert res.n_records_in == res.n_records_out == 10
        assert out.read_text() == gvcf.read_text()

    @defect_persists("4 (`unknown` records no reason)",
                     closed_by="outcome_reason.explain_outcome")
    def test_repro_unknown_carries_no_reason(self, tmp_path):
        """REPRODUCTION (defect 4): `unknown` is a bare string. The single most
        load-bearing distinction in the analysis — uncallable vs unedited — is lost."""
        do = _trace.donor_outcome()
        ref, dstart, dend, wt, design, vcf = self._unedited_case(tmp_path)
        gvcf = write_gvcf(tmp_path / "s.g.vcf", [(dstart, None, 30)], ref=ref)  # 1 position only
        spec = do.DonorSpec(
            sample="S1", chrom=CONTIG, donor_start=dstart, donor_end=dend,
            donor_wt=wt, donor_design=design, donor_step2=design,
        )
        out = do.precall_sample(spec, vcf, gvcf, ref)
        assert out.editing_outcome == "unknown"
        assert getattr(out, "reason", None), (
            "`unknown` records no reason, so an uncallable cell is indistinguishable "
            "from a called one in every downstream table"
        )

    def test_fixed_reason_is_explicit_and_machine_readable(self, tmp_path):
        from magestic.pipelines.wgs.core.outcome_reason import (
            REASON_GVCF_ABSENT,
            REASON_QWIN_UNCOVERED,
            REASON_UNEDITED_COVERED,
            explain_outcome,
        )

        ref, dstart, dend, wt, design, vcf = self._unedited_case(tmp_path)
        gvcf = write_gvcf(tmp_path / "s.g.vcf", [(dstart, None, 30)], ref=ref)
        ex = explain_outcome(
            "unknown", gvcf_path=gvcf, chrom=CONTIG, donor_start=dstart, donor_end=dend
        )
        assert ex.reason == REASON_QWIN_UNCOVERED
        assert ex.callable_ is False
        assert ex.n_window_positions_covered == 0
        assert ex.n_window_positions == 10
        # ...and the same function says so when the gVCF does not exist at all.
        ex2 = explain_outcome(
            "unknown", gvcf_path=tmp_path / "nope.g.vcf", chrom=CONTIG,
            donor_start=dstart, donor_end=dend,
        )
        assert ex2.reason == REASON_GVCF_ABSENT
        assert ex2.callable_ is False
        # ...and an `unedited` is positively marked CALLABLE, so the two never collide.
        gvcf2 = write_gvcf(tmp_path / "ok.g.vcf", [(dstart, dend, 30)], ref=ref)
        ex3 = explain_outcome(
            "designed_variant_unedited", gvcf_path=gvcf2, chrom=CONTIG,
            donor_start=dstart, donor_end=dend,
        )
        assert ex3.reason == REASON_UNEDITED_COVERED
        assert ex3.callable_ is True

    def test_fixed_reason_agrees_with_traces_own_flag(self, tmp_path):
        """A reason that contradicts TRACE's flag would be a NEW silent defect.

        🔴 This is only a real test because `explain_outcome` derives its reason from
        the END-AWARE COVERAGE, not from the flag string. Were it derived from the
        flag it could not disagree, and this test would be a tautology.
        `test_fixed_reason_names_the_unexpanded_block_case` pins the one case where
        the two genuinely diverge, which is what proves the derivation is independent.
        """
        from magestic.pipelines.wgs.core.outcome_reason import explain_outcome

        do = _trace.donor_outcome()
        ref, dstart, dend, wt, design, vcf = self._unedited_case(tmp_path)
        for blocks, expect_callable in (
            ([(dstart, dend, 30)], True),
            ([(dstart, None, 30)], False),
            ([], False),
        ):
            # planning/268: read as written; the pinned TRACE honours END itself.
            ge = write_gvcf(tmp_path / f"g{len(blocks)}{expect_callable}.g.vcf", blocks, ref=ref)
            spec = do.DonorSpec(
                sample="S1", chrom=CONTIG, donor_start=dstart, donor_end=dend,
                donor_wt=wt, donor_design=design, donor_step2=design,
            )
            flag = do.precall_sample(spec, vcf, ge, ref).editing_outcome
            ex = explain_outcome(
                flag, gvcf_path=ge, chrom=CONTIG, donor_start=dstart, donor_end=dend
            )
            assert ex.callable_ is expect_callable
            assert (flag == "unknown") is (not expect_callable)
            # The numbers must support the reason, not merely sit beside it.
            assert (ex.n_window_positions_covered >= ex.n_window_positions) is (
                expect_callable
            )

    def test_fixed_a_block_without_min_dp_is_uncovered_on_both_sides(self, tmp_path):
        """WITHDRAWN (planning/268): "TRACE says `unknown` on a covered compressed block,
        so the reason must say the blocks were not expanded". The pinned TRACE reads
        `END`; it counts a block's interior only from `MIN_DP`, and
        `gvcf_adapter.window_reference_coverage` now follows that rule. So a block
        without `MIN_DP` is uncovered on BOTH sides and the reason says so, without
        arguing with its own numbers; `REASON_BLOCKS_NOT_EXPANDED` is no longer emitted.
        """
        from magestic.pipelines.wgs.core.outcome_reason import (
            REASON_BLOCKS_NOT_EXPANDED,
            REASON_QWIN_UNCOVERED,
            explain_outcome,
        )

        do = _trace.donor_outcome()
        ref, dstart, dend, wt, design, vcf = self._unedited_case(tmp_path)
        raw = write_gvcf(tmp_path / "raw.g.vcf", [(dstart, dend, 30)], ref=ref, min_dp=False)
        spec = do.DonorSpec(
            sample="S1", chrom=CONTIG, donor_start=dstart, donor_end=dend,
            donor_wt=wt, donor_design=design, donor_step2=design,
        )
        flag = do.precall_sample(spec, vcf, raw, ref).editing_outcome
        assert flag == "unknown"
        ex = explain_outcome(
            flag, gvcf_path=raw, chrom=CONTIG, donor_start=dstart, donor_end=dend
        )
        assert ex.n_window_positions_covered == 0
        assert ex.reason == REASON_QWIN_UNCOVERED != REASON_BLOCKS_NOT_EXPANDED
        assert ex.callable_ is False

    def test_fixed_the_reason_reaches_the_shipped_shard_not_only_the_log(self, tmp_path):
        """🔴 Defect 4 is only fixed if the reason is in the TABLE.

        A reason that reached only the run log would leave every downstream reader
        unable to tell "uncallable" from "unedited" — so `step_2` appends the columns
        to TRACE's own 5-column output. Appended, never inserted, because `step_5`
        looks `editing_outcome` up by NAME (`step_5_join_outcomes.py:139`).
        """
        from magestic.pipelines.wgs.pipeline.step_2_donor_precall import (
            OUTCOME_REASON_COLUMNS,
            _append_columns,
        )

        t = tmp_path / "out.tsv"
        t.write_text(
            "#sample_name\tediting_outcome\tdesigned_variants\tsynthetic_errors"
            "\tspontaneous_variants\n"
            "S1\tunknown\t\t\t\n"
        )
        _append_columns(t, OUTCOME_REASON_COLUMNS,
                        ["false", "uncallable_gvcf_absent", "0", "10", "2",
                         "chrXV_3060_A_to_G"])
        header, row = t.read_text().splitlines()
        assert header.split("\t")[:5] == [
            "#sample_name", "editing_outcome", "designed_variants",
            "synthetic_errors", "spontaneous_variants",
        ]
        cols = header.split("\t")
        assert cols[cols.index("callable")] == "callable"
        assert row.split("\t")[cols.index("outcome_reason")] == "uncallable_gvcf_absent"
        assert row.split("\t")[cols.index("callable")] == "false"
        assert row.split("\t")[1] == "unknown"  # the flag did not move
        assert len(row.split("\t")) == len(cols)
        # Idempotent: an array task re-running must not append twice.
        before = t.read_text()
        _append_columns(t, OUTCOME_REASON_COLUMNS, ["x"] * len(OUTCOME_REASON_COLUMNS))
        assert t.read_text() == before

    def test_fixed_step_2_wires_the_block_expansion_and_the_reason(self):
        """The adapters must be ON THE PIPELINE PATH, which is defect 14's own lesson.

        An adapter that exists and is called by nothing is exactly the state NN 251
        was sent to fix; asserting the call sites keeps this change set from repeating
        it.
        """
        from magestic.pipelines.wgs.pipeline import step_2_donor_precall as s2

        src = Path(s2.__file__).read_text()
        # WITHDRAWN (planning/268): "step_2 expands gVCF blocks" (defect 5 is closed in
        # the pinned TRACE, and the expander's copied REF is refused by it).
        assert "expand_reference_blocks(" not in src
        assert "explain_outcome(" in src  # defect 4
        assert "audit_vcf_window(" in src  # defect 8
        assert "n_multiallelic" in src  # defect 14
        assert "min_allele_fraction" in src  # defect 1
        # 🔴 Defect 1's second half: the RESCUED minor-allele ids must reach the shard.
        # TRACE reads only VCF columns 0/1/3/4, so a rescued 10% allele is treated as a
        # 100% call; without this column a mixture ships as `designed_variant_installed`.
        assert "minor_allele_variants" in src
        # planning/268: the ids are read back from the split VCF's MAG_MINOR tags.
        assert "minor_alleles=_minor_alleles_from_vcf(" in src
        from magestic.pipelines.wgs.pipeline.step_2_donor_precall import (
            OUTCOME_REASON_COLUMNS,
        )

        assert "minor_allele_variants" in OUTCOME_REASON_COLUMNS

        # Defect 17, checked on the AST rather than the text: the module's own
        # docstring quotes the broken `ad_index=` call to explain the bug, so a
        # substring search would match the explanation instead of the code.
        import ast

        calls = [
            n
            for n in ast.walk(ast.parse(src))
            if isinstance(n, ast.Call)
            and getattr(n.func, "id", getattr(n.func, "attr", None)) == "split_plate_vcf"
        ]
        assert len(calls) == 1, f"expected one split_plate_vcf call, found {len(calls)}"
        kwargs = {k.arg for k in calls[0].keywords}
        assert "ad_slot" in kwargs, "the parameter is ad_slot, not ad_index"
        assert "ad_index" not in kwargs
        slot = next(k.value for k in calls[0].keywords if k.arg == "ad_slot")
        # ...and it must be the DESTINATION slot (2), not the source slot (4).
        assert ast.unparse(slot).endswith("ad_format_slot_out"), ast.unparse(slot)

    def test_fixed_step_4_wires_the_depth_gate_and_the_frame_check(self):
        from magestic.pipelines.wgs.pipeline import step_4_reconcile as s4

        src = Path(s4.__file__).read_text()
        assert "verify_window_frame(" in src  # defect 10
        assert "load_windows_audited(" in src  # defect 12
        assert "incomparable_contigs(" in src  # defect 14
        assert "require_reference_consistency(" in src  # defect 14
        assert "designed_positions_by_oligo(" in src  # defect 11, for undeclared windows
        # WITHDRAWN (planning/268): the `depth_gate` switch (the read-evidence route is
        # always gated now), the `bc_confidence`-renamed writer (15) and
        # `unwindowed_oligos` (12): `build_edit_sites` raises for an oligo with no window,
        # so none can be dropped in silence.
        assert "depth_gate" not in src
        assert "write_reconciliations(" in src

    def test_fixed_a_sample_with_no_gvcf_still_gets_a_row(self, tmp_path):
        """TRACE drops such a sample with no row at all — the worst `uncallable`."""
        from magestic.pipelines.wgs.core.outcome_reason import rows_for_missing_inputs

        rows = rows_for_missing_inputs(
            ["S1", "S2", "S3"], processed=["S1"], workdir=tmp_path
        )
        assert [r.sample for r in rows] == ["S2", "S3"]
        assert all(r.callable_ is False and r.reason for r in rows)


# --------------------------------------------------------------------------- #
# Defects 6 + 11 — one number for SNV/indel/SV; min_cov on the window MAXIMUM   #
# --------------------------------------------------------------------------- #


class TestDefects6And11EditProfile:
    WIN = (3001, 3060)

    @staticmethod
    def _genome(tmp_path: Path):
        ref = _seq(21, 4000)
        return ref, write_fasta(tmp_path / "g.fa", {CONTIG: ref})

    @defect_persists("6 (SNV/indel/SV collapse to one number) in TRACE's LEGACY "
                     "window_edit_vaf",
                     closed_by="TRACE score_edit_site (allele-specific evidence), which "
                               "reconcile_driver.score_sites calls since planning/268")
    def test_repro_trace_collapses_a_deletion_and_a_substitution(self, tmp_path):
        """REPRODUCTION (defect 6), driving **TRACE's own functions** directly.

        A 6 bp deletion and a 1 bp substitution give the SAME `(vaf, cov, status)` from
        `window_edit_vaf` + `classify_site`, and nothing in TRACE's `SiteCall` records
        which it was — so an SV is indistinguishable from a substitution.

        This targets TRACE deliberately. `score_sites` is MAGESTIC code and is now
        fixed (see `test_fixed_score_sites_records_the_variant_class`), so a repro
        driving it would no longer be reproducing anything.
        """
        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        got = {}
        for name, bam_path in (
            ("del", write_bam(tmp_path / "del.bam", ref, spans=(self.WIN,),
                              deletion=(3030, 6))),
            ("snv", write_bam(tmp_path / "snv.bam", ref, spans=(self.WIN,),
                              edits={3030: _other(ref[3029])})),
        ):
            with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
                str(fasta)
            ) as fa:
                vaf, cov = TR.window_edit_vaf(bam, fa, CONTIG, ws, we)
                call = TR.SiteCall(TR.EditSite(CONTIG, ws, we, "O1"), vaf, cov,
                                   TR.classify_site(vaf, cov))
            got[name] = call
        assert got["del"].status == got["snv"].status == TR.SITE_FULL
        assert hasattr(got["del"], "variant_class"), (
            "TRACE's SiteCall records no variant type, so SNV, indel and SV all "
            "collapse to SINGLE_EDIT and an SV reads as a substitution"
        )

    @staticmethod
    def _allele_table(path: Path, ws: int, we: int, ref: str, *, pos: int = 3030) -> Path:
        """One window declaring one intended SNV at `pos` (planning/268: the pinned TRACE
        scores only declared `intended_variants`; an undeclared window is `unassessed`)."""
        import json

        allele = json.dumps([{"pos": pos, "ref": ref[pos - 1], "alt": _other(ref[pos - 1])}])
        return write_window_table(
            path, [("O1", CONTIG, ws, we, allele)],
            header=("oligo_name", "chrom", "donor_start_coord", "donor_end_coord",
                    "intended_variants"),
        )

    def test_fixed_score_sites_records_the_variant_class(self, tmp_path):
        """🔴 Defect 6 on the path Stage 4 ACTUALLY TAKES: a deletion must never come back
        as a substitution.

        planning/268: `score_sites` has two routes and both are pinned. A window with NO
        declared allele is classified from read evidence by `core.edit_profile` and
        records `variant_class` (Kevin 2026-10-24: an intended edit is not required to
        assign confidence). A window that declares the SNV is scored by TRACE per allele:
        FULL on the SNV's own fragments, never FULL from a 6 bp deletion over that base.
        WITHDRAWN: the `depth_gate=False` parity switch, removed with that config key.
        """
        from magestic.pipelines.wgs.core.edit_profile import (
            VARIANT_DELETION,
            VARIANT_SUBSTITUTION,
        )
        from magestic.pipelines.wgs.core.reconcile_driver import (
            FLAG_WINDOW_EVIDENCE,
            load_oligo_windows,
            score_sites,
        )

        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        bare = load_oligo_windows(
            str(write_window_table(tmp_path / "w.tsv", [("O1", CONTIG, ws, we)]))
        )
        declared = load_oligo_windows(str(self._allele_table(tmp_path / "a.tsv", ws, we, ref)))
        snv = write_bam(tmp_path / "s.bam", ref, spans=(self.WIN,),
                        edits={3030: _other(ref[3029])})
        dele = write_bam(tmp_path / "d.bam", ref, spans=(self.WIN,), deletion=(3030, 6))
        got = {}
        for name, bam_path, expect in (("del", dele, VARIANT_DELETION),
                                       ("snv", snv, VARIANT_SUBSTITUTION)):
            with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
                str(fasta)
            ) as fa:
                window = score_sites(bam, fa, ["O1"], bare)[0]
                got[name] = score_sites(bam, fa, ["O1"], declared)[0]
            assert window.variant_class == expect, (name, window.variant_class)
            assert window.reason, "no status without a reason"
            assert FLAG_WINDOW_EVIDENCE in window.flags
            assert window.trace_status == TR.SITE_UNASSESSED
        assert got["snv"].status == TR.SITE_FULL
        assert got["snv"].allele_evidence[0].alt_count == got["snv"].allele_evidence[0].depth > 0
        assert got["del"].status != TR.SITE_FULL
        assert all(e.alt_count == 0 for e in got["del"].allele_evidence)

    def test_fixed_profile_classifies_the_variant_type(self, tmp_path):
        from magestic.pipelines.wgs.core.edit_profile import (
            VARIANT_DELETION,
            VARIANT_SUBSTITUTION,
            window_edit_profile,
        )

        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        cases = {
            "del": (write_bam(tmp_path / "d.bam", ref, spans=(self.WIN,),
                              deletion=(3030, 6)), VARIANT_DELETION),
            "snv": (write_bam(tmp_path / "s.bam", ref, spans=(self.WIN,),
                              edits={3030: _other(ref[3029])}), VARIANT_SUBSTITUTION),
        }
        for name, (bam_path, expect) in cases.items():
            with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
                str(fasta)
            ) as fa:
                prof = window_edit_profile(bam, fa, CONTIG, ws, we)
            assert prof.variant_class == expect, f"{name}: {prof.variant_class}"
        # 🔴 the hard requirement: a deletion is NEVER reported as a substitution.
        assert cases["del"][1] != VARIANT_SUBSTITUTION

    def test_fixed_profile_is_vaf_parity_with_trace(self, tmp_path):
        """The headline VAF/cov provably does not move."""
        from magestic.pipelines.wgs.core.edit_profile import window_edit_profile

        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        for bam_path in (
            write_bam(tmp_path / "p1.bam", ref, spans=(self.WIN,),
                      edits={3030: _other(ref[3029])}),
            write_bam(tmp_path / "p2.bam", ref, spans=(self.WIN,)),
            write_bam(tmp_path / "p3.bam", ref, spans=(self.WIN,), deletion=(3030, 6)),
        ):
            with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
                str(fasta)
            ) as fa:
                prof = window_edit_profile(bam, fa, CONTIG, ws, we)
            with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
                str(fasta)
            ) as fa:
                vaf, cov = TR.window_edit_vaf(bam, fa, CONTIG, ws, we)
            assert (round(prof.max_vaf, 9), prof.cov_at_max) == (round(vaf, 9), cov)

    def test_fixed_profile_is_vaf_parity_with_trace_on_an_insertion(self, tmp_path):
        """🔴 The one place the re-walk could diverge from TRACE, closed by measurement.

        TRACE increments `nonref` under
        `query_sequence[qpos] != refbase or pr.indel > 0` — a single `or`, so a read
        whose anchor base BOTH mismatches and opens an insertion counts ONCE. The
        profile splits that counter into `n_substitution` / `n_insertion` to type the
        variant (defect 6), which risks counting such a read twice and inflating the
        VAF above TRACE's. `elif` ordering prevents it; this asserts it on a BAM with
        a real insertion CIGAR rather than arguing it (251 section 7 item 4).
        """
        from magestic.pipelines.wgs.core.edit_profile import (
            VARIANT_INSERTION,
            window_edit_profile,
        )

        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        # Build reads carrying a 3 bp insertion at 3030, with a MISMATCHED anchor.
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": CONTIG, "LN": len(ref)}]}
        recs = []
        for i in range(20):
            st = 3001 + (i % 10)
            left = 3030 - st + 1
            seg = list(ref[st - 1 : st - 1 + 50])
            seg[left - 1] = _other(seg[left - 1])  # anchor mismatches AND inserts
            a = pysam.AlignedSegment()
            a.query_name = f"ins{i}"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = st - 1
            a.mapping_quality = 60
            a.query_sequence = ("".join(seg[:left]) + "GGG" + "".join(seg[left:]))
            a.cigartuples = [(0, left), (1, 3), (0, 50 - left)]
            a.query_qualities = pysam.qualitystring_to_array("I" * (50 + 3))
            recs.append(a)
        unsorted = tmp_path / "ins.unsorted.bam"
        with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
            for r in sorted(recs, key=lambda r: r.reference_start):
                out.write(r)
        bam_path = tmp_path / "ins.bam"
        pysam.sort("-o", str(bam_path), str(unsorted))
        pysam.index(str(bam_path))

        with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            prof = window_edit_profile(bam, fa, CONTIG, ws, we)
        with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            vaf, cov = TR.window_edit_vaf(bam, fa, CONTIG, ws, we)
        assert (round(prof.max_vaf, 9), prof.cov_at_max) == (round(vaf, 9), cov)
        # No position may count a read twice: non-ref can never exceed the depth.
        for p in prof.positions:
            assert p.n_nonref <= p.depth, (p.pos, p.n_nonref, p.depth)
        assert prof.variant_class_at(3030) == VARIANT_INSERTION

    @defect_persists("11 (min_cov gated on the window MAXIMUM depth) in TRACE's LEGACY "
                     "window_edit_vaf",
                     closed_by="TRACE score_edit_site (depth AT each intended allele), which "
                               "reconcile_driver.score_sites calls since planning/268")
    def test_repro_trace_reports_a_confident_no_edit_over_a_hole(self, tmp_path):
        """REPRODUCTION (defect 11 / A10 / A12), driving **TRACE's own functions**.

        The window MAXIMUM depth clears `min_cov` while the designed base sits at
        depth 0, so `classify_site` returns `SITE_WT`, `reconcile_cell` reports
        `NO_EDIT` with `n_nocov=0`, and confidence is never downgraded.
        """
        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        hole = tuple(range(3028, 3034))
        bam_path = write_bam(tmp_path / "hole.bam", ref, spans=(self.WIN,),
                             skip_positions=hole)
        with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            vaf, cov = TR.window_edit_vaf(bam, fa, CONTIG, ws, we)
        status = TR.classify_site(vaf, cov)
        cell = TR.reconcile_cell(
            "S1", [TR.SiteCall(TR.EditSite(CONTIG, ws, we, "O1"), vaf, cov, status)],
            n_genome_bc=1, bc_confidence="high",
        )
        assert cov >= 4  # the window maximum clears min_cov
        assert cell.label == TR.NO_EDIT and cell.n_nocov == 0
        assert status == TR.SITE_NOCOV, (
            f"depth 0 across {hole} yet status={status} cov={cov} "
            f"label={cell.label} n_nocov={cell.n_nocov} confidence={cell.confidence} "
            "— a confident NO_EDIT over a hole"
        )

    def test_fixed_score_sites_gates_depth_at_the_hole(self, tmp_path):
        """🔴 Defect 11 on the path Stage 4 ACTUALLY TAKES.

        Both `score_sites` routes must turn that confident `NO_EDIT` into `LOW_CONF`
        (planning/268). Undeclared window: gated at the window minimum, with the
        reclassification recorded (TRACE alone would call it `unassessed`). Declared
        allele: TRACE counts fragments AT the allele, so the hole is `nocov`. The window
        route's headline `(vaf, cov)` stays TRACE's `window_edit_vaf`, so no number moves.
        """
        from magestic.pipelines.wgs.core.reconcile_driver import (
            load_oligo_windows,
            score_sites,
        )

        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        bam_path = write_bam(tmp_path / "hole.bam", ref, spans=(self.WIN,),
                             skip_positions=tuple(range(3028, 3034)))
        bare = load_oligo_windows(
            str(write_window_table(tmp_path / "w.tsv", [("O1", CONTIG, ws, we)]))
        )
        declared = load_oligo_windows(str(self._allele_table(tmp_path / "a.tsv", ws, we, ref)))
        with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            calls = score_sites(bam, fa, ["O1"], bare)
            allele = score_sites(bam, fa, ["O1"], declared)
            vaf, cov = TR.window_edit_vaf(bam, fa, CONTIG, ws, we)
        cell = TR.reconcile_cell("S1", calls, n_genome_bc=1, bc_confidence="high")
        assert calls[0].status == TR.SITE_NOCOV
        assert calls[0].reason == "window_min_depth_below_min_cov"
        assert calls[0].min_depth == 0
        assert calls[0].n_positions_uncovered >= 6
        assert calls[0].trace_status == TR.SITE_UNASSESSED
        assert cell.label == TR.LOW_CONF and cell.n_nocov == 1
        assert (round(calls[0].vaf, 9), calls[0].cov) == (round(vaf, 9), cov)
        assert allele[0].status == TR.SITE_NOCOV
        assert allele[0].allele_evidence[0].depth == 0
        assert "UNCOVERED_ALLELES" in allele[0].flags
        other = TR.reconcile_cell("S1", allele, n_genome_bc=1, bc_confidence="high")
        assert other.label == TR.LOW_CONF and other.confidence == "low"

    def test_fixed_designed_position_gate_reaches_score_sites(self, tmp_path):
        """With `wt_seq` + `design_donor` in the table, the gate narrows to the edit;
        with the allele declared, TRACE gates at it too (planning/268)."""
        from magestic.pipelines.wgs.core.reconcile_driver import (
            load_oligo_windows,
            score_sites,
        )
        from magestic.pipelines.wgs.core.window_table import (
            designed_positions_by_oligo,
            load_windows_audited,
            verify_window_frame,
        )

        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        wt_seq = ref[ws - 1 : we]
        design = wt_seq[:29] + _other(wt_seq[29]) + wt_seq[30:]
        tab = write_window_table(
            tmp_path / "w.tsv",
            [("O1", CONTIG, ws, we, wt_seq, design)],
            header=("oligo_name", "chrom", "donor_start_coord", "donor_end_coord",
                    "wt_seq", "design_donor"),
        )
        audit = load_windows_audited(str(tab))
        # The orientation comes from the frame check, exactly as `step_4` does it —
        # it is never assumed, because a plus-strand assumption on a minus-strand
        # design gates the mirror base.
        with pysam.FastaFile(str(fasta)) as fa:
            oriented = verify_window_frame(fa, CONTIG, ws, we, wt_seq).oriented
        dpos = designed_positions_by_oligo(audit, orientations={"O1": oriented})
        assert dpos == {"O1": (3030,)}
        # A hole at the DESIGNED base only.
        bam_path = write_bam(tmp_path / "h.bam", ref, spans=(self.WIN,),
                             skip_positions=(3030,))
        declared = load_oligo_windows(str(self._allele_table(tmp_path / "a.tsv", ws, we, ref)))
        with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            calls = score_sites(
                bam, fa, ["O1"], load_oligo_windows(str(tab)), designed_positions=dpos
            )
            allele = score_sites(bam, fa, ["O1"], declared)
        assert calls[0].status == TR.SITE_NOCOV
        assert calls[0].reason == "designed_position_depth_below_min_cov"
        assert calls[0].depth_at_designed == 0
        assert allele[0].status == TR.SITE_NOCOV
        assert allele[0].allele_evidence[0].depth == 0

    def test_fixed_min_depth_gates_the_call(self, tmp_path):
        from magestic.pipelines.wgs.core.edit_profile import (
            classify_window,
            window_edit_profile,
        )

        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        bam_path = write_bam(tmp_path / "hole.bam", ref, spans=(self.WIN,),
                             skip_positions=tuple(range(3028, 3034)))
        with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            prof = window_edit_profile(bam, fa, CONTIG, ws, we)
        assert prof.min_depth == 0
        assert prof.n_positions_uncovered >= 6
        call = classify_window(prof)
        assert call.status == TR.SITE_NOCOV
        assert call.reason == "window_min_depth_below_min_cov"

    def test_fixed_designed_position_depth_gates_when_known(self, tmp_path):
        """With `wt_seq`+design (defect 10) the gate is at the designed base itself."""
        from magestic.pipelines.wgs.core.edit_profile import (
            classify_window,
            window_edit_profile,
        )

        ref, fasta = self._genome(tmp_path)
        ws, we = self.WIN
        # A hole at 3030 only; the rest of the window is deep.
        bam_path = write_bam(tmp_path / "h2.bam", ref, spans=(self.WIN,),
                             skip_positions=(3030,))
        with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            prof = window_edit_profile(bam, fa, CONTIG, ws, we)
        call = classify_window(prof, designed_positions=(3030,))
        assert call.status == TR.SITE_NOCOV
        assert call.reason == "designed_position_depth_below_min_cov"
        assert call.depth_at_designed == 0


# --------------------------------------------------------------------------- #
# Defect 7 — no inversion/duplication detector                                  #
# --------------------------------------------------------------------------- #


class TestDefect7SvRefusal:
    def test_repro_no_detector_exists(self):
        """REPRODUCTION: TRACE has no inversion or tandem-duplication detector at all.

        `decode_duplicate` is a SAM dup-FLAG decoder, not an SV detector, so it is
        excluded by name — the search is for a *detect*/*call*/*find* entry point.
        """
        sv = _trace.sv_plasmid()
        detectors = [
            n
            for n in dir(sv)
            if ("invers" in n.lower() or "tandem" in n.lower()
                or ("duplicat" in n.lower() and n != "decode_duplicate"))
            and any(v in n.lower() for v in ("detect", "call", "find", "count", "score"))
        ]
        assert detectors == [], f"unexpected detector(s) appeared: {detectors}"
        # The translocation / large-deletion helpers DO exist but nothing calls them,
        # so there is no path from a BAM to an SV rate of any kind.
        assert hasattr(sv, "mate_pair_large_deletion")
        assert hasattr(sv, "within_read_large_deletions")
        from magestic.pipelines.wgs import pipeline

        text = "\n".join(
            p.read_text() for p in sorted(Path(pipeline.__file__).parent.glob("*.py"))
        )
        assert "mate_pair_large_deletion" not in text
        assert "within_read_large_deletions" not in text

    def test_fixed_pipeline_refuses_to_emit_an_sv_rate(self):
        """DELIBERATE CHOICE: refuse, rather than ship a half-built detector."""
        from magestic.pipelines.wgs.core.sv_status import (
            SV_UNMEASURED,
            SvRateRefused,
            require_no_sv_rate,
            sv_status,
        )

        assert sv_status() == SV_UNMEASURED
        assert "UNMEASURED" in SV_UNMEASURED
        require_no_sv_rate({"upstream_count": 3, "cells": 10})  # fine
        for bad in ("sv_rate", "inversion_rate", "duplication_rate", "sv_frequency",
                    "pct_sv", "structural_variant_rate"):
            with pytest.raises(SvRateRefused):
                require_no_sv_rate({bad: 0.0})

    def test_fixed_step_3_ships_an_explicit_sv_status_column(self):
        from magestic.pipelines.wgs.core.sv_status import SV_UNMEASURED
        from magestic.pipelines.wgs.pipeline import step_3_sv_plasmid as s3

        assert "sv_status" in s3.SHARD_COLUMNS
        assert s3.SV_STATUS_VALUE == SV_UNMEASURED

    def test_fixed_a_stale_shard_is_refused_on_resume(self, tmp_path):
        """🔴 Adding columns to a LIVE stage's shard breaks idempotent resume.

        Stage 3 is the live stage, so 10-column shards written before NN 251 may
        already exist on Oak. `_run_one_sample` skips an existing shard and
        `collect_shards` only catches drift BETWEEN shards, so a cohort of uniformly
        stale shards would be read with the NEW indices and `half_b_status` taken
        from `frags_genome`. That is a silent column shift, not a missing column.
        """
        from magestic.pipelines.wgs.pipeline.step_3_sv_plasmid import (
            SHARD_COLUMNS,
            _require_current_shard_schema,
        )
        from magestic.pipelines.wgs.pipeline._common import StepIncompleteError

        old = ("sample", "upstream_count", "downstream_count", "total_count",
               "n_integrated_contigs", "half_b_status", "read_count_mid_point",
               "frags_genome", "frags_plasmid", "frags_chimeric")
        # planning/268 adds the 205 lane's `half_a_status`: 13 columns.
        assert len(old) == 10 and len(SHARD_COLUMNS) == 13
        stale = tmp_path / "S1.sv_plasmid.tsv"
        stale.write_text("\t".join(old) + "\n" + "\t".join(["x"] * 10) + "\n")
        with pytest.raises(StepIncompleteError) as exc:
            _require_current_shard_schema(stale, "S1")
        assert "DIFFERENT" in str(exc.value) and "Delete this shard" in str(exc.value)
        # A current shard resumes cleanly.
        fresh = tmp_path / "S2.sv_plasmid.tsv"
        fresh.write_text(
            "\t".join(SHARD_COLUMNS) + "\n" + "\t".join(["x"] * len(SHARD_COLUMNS)) + "\n"
        )
        _require_current_shard_schema(fresh, "S2")

    def test_fixed_no_step_emits_a_numeric_sv_rate(self):
        """Source-text guard: no rate for a quantity with no detector, anywhere."""
        import re

        from magestic.pipelines.wgs import pipeline

        root = Path(pipeline.__file__).parent
        bad = re.compile(r"\"(sv|inversion|duplication|structural)\w*_(rate|freq\w*)\"")
        hits = [
            f"{p.name}:{i}"
            for p in sorted(root.glob("*.py"))
            for i, line in enumerate(p.read_text().splitlines(), 1)
            if bad.search(line)
        ]
        assert hits == [], hits


# --------------------------------------------------------------------------- #
# Defect 8 — variants outside the donor window are discarded invisibly          #
# --------------------------------------------------------------------------- #


class TestDefect8DistalDiscards:
    def test_repro_distal_variants_vanish(self, tmp_path):
        """REPRODUCTION: a variant 240 bp away is dropped with no count."""
        do = _trace.donor_outcome()
        ref = _seq(31, 900)
        dstart, dend = 101, 160
        wt = ref[dstart - 1 : dend]
        design = wt[:30] + _other(wt[30]) + wt[31:]
        vcf = write_vcf(
            tmp_path / "s.fil.vcf",
            [
                (CONTIG, dstart + 30, wt[30], _other(wt[30])),  # in window, designed
                (CONTIG, 400, ref[399], _other(ref[399])),  # distal, invisible
                ("otherchrom", 50, "A", "G"),  # off-contig, invisible
            ],
        )
        qs, qe = qwin(dstart, dend)
        gvcf = write_gvcf(tmp_path / "s.g.vcf", [(p, None, 30) for p in range(qs, qe + 1)], ref=ref)
        spec = do.DonorSpec(
            sample="S1", chrom=CONTIG, donor_start=dstart, donor_end=dend,
            donor_wt=wt, donor_design=design, donor_step2=design,
        )
        out = do.precall_sample(spec, vcf, gvcf, ref)
        total = len(out.designed_variants) + len(out.synthetic_errors) + len(
            out.spontaneous_variants
        )
        # 3 variants in the VCF, 1 in the outcome, and nothing anywhere says 2 vanished.
        assert total == 1  # only the in-window one survives
        assert sum(1 for _ in open(vcf) if not _.startswith("#")) == 3
        from magestic.pipelines.wgs.core.donor_precall import audit_vcf_window

        audit = audit_vcf_window(vcf, chrom=CONTIG, donor_start=dstart, donor_end=dend)
        assert audit.n_discarded == 2
        assert audit.n_distal_same_contig == 1
        assert audit.n_other_contig == 1
        assert f"{CONTIG}_400" in ";".join(audit.discarded_ids)

    def test_fixed_window_is_configurable_and_the_count_follows(self, tmp_path):
        from magestic.pipelines.wgs.core.donor_precall import audit_vcf_window

        ref = _seq(31, 900)
        dstart, dend = 101, 160
        vcf = write_vcf(
            tmp_path / "s.fil.vcf",
            [(CONTIG, 400, ref[399], _other(ref[399]))],
        )
        tight = audit_vcf_window(vcf, chrom=CONTIG, donor_start=dstart, donor_end=dend)
        assert tight.n_discarded == 1
        wide = audit_vcf_window(
            vcf, chrom=CONTIG, donor_start=dstart, donor_end=dend, flank=500
        )
        assert wide.n_discarded == 0
        assert wide.n_in_flank == 1


# --------------------------------------------------------------------------- #
# Defect 10 — Stage 4 is never given the design's wt_seq (the allele inversion)  #
# --------------------------------------------------------------------------- #


class TestDefect10FrameVerification:
    WIN = (3001, 3060)

    @defect_persists("10 / NN243 A11 (Stage 4 never sees the design's wt_seq)",
                     closed_by="window_table.verify_window_frame")
    def test_repro_an_inverted_reference_frame_is_undetectable(self, tmp_path):
        """REPRODUCTION (A11): reads that are WT against the TRUE genome score as a
        full edit against a reference that already carries the design's ALT — and
        nothing in Stage 4 can tell, because it never sees `wt_seq`."""
        from magestic.pipelines.wgs.core.reconcile_driver import (
            load_oligo_windows,
            score_sites,
        )

        truth = _seq(41, 4000)
        pos = 3030
        wrong = truth[: pos - 1] + _other(truth[pos - 1]) + truth[pos:]
        fasta = write_fasta(tmp_path / "wrong.fa", {CONTIG: wrong})
        ws, we = self.WIN
        bam_path = write_bam(tmp_path / "wt.bam", truth, spans=(self.WIN,))
        wt_tab = write_window_table(tmp_path / "w.tsv", [("O1", CONTIG, ws, we)])
        ow = load_oligo_windows(str(wt_tab))
        with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            calls = score_sites(bam, fa, ["O1"], ow)
            cell = TR.reconcile_cell("S1", calls, n_genome_bc=1, bc_confidence="high")
        # Unedited reads, yet a confident SINGLE_EDIT: the inversion.
        assert cell.label == TR.SINGLE_EDIT and calls[0].vaf > 0.9
        assert "FRAME" in " ".join(cell.flags) or False, (
            f"unedited reads scored {cell.label} vaf={calls[0].vaf:.2f} against a "
            "reference carrying the design's ALT, with no flag of any kind"
        )

    def test_fixed_wt_seq_detects_the_inverted_frame(self, tmp_path):
        from magestic.pipelines.wgs.core.window_table import (
            FrameMismatchError,
            verify_window_frame,
        )

        truth = _seq(41, 4000)
        pos = 3030
        wrong = truth[: pos - 1] + _other(truth[pos - 1]) + truth[pos:]
        fasta = write_fasta(tmp_path / "wrong.fa", {CONTIG: wrong})
        ws, we = self.WIN
        wt_seq = truth[ws - 1 : we]
        with pysam.FastaFile(str(fasta)) as fa:
            with pytest.raises(FrameMismatchError) as exc:
                verify_window_frame(fa, CONTIG, ws, we, wt_seq)
            assert "1 mismatch" in str(exc.value) or "mismatch" in str(exc.value)
            # ...and the correct reference verifies, in either orientation.
            ok = write_fasta(tmp_path / "right.fa", {CONTIG: truth})
        with pysam.FastaFile(str(ok)) as fa:
            assert verify_window_frame(fa, CONTIG, ws, we, wt_seq).oriented == "plus"
            rc = _trace.donor_outcome().reverse_complement(wt_seq)
            assert verify_window_frame(fa, CONTIG, ws, we, rc).oriented == "minus"

    def test_fixed_missing_wt_seq_is_explicit_not_silent(self, tmp_path):
        """If `wt_seq` is absent the guard must SAY it did not run."""
        from magestic.pipelines.wgs.core.window_table import load_windows_audited

        ws, we = self.WIN
        p = write_window_table(tmp_path / "w.tsv", [("O1", CONTIG, ws, we)])
        audit = load_windows_audited(str(p))
        assert audit.has_wt_seq is False
        assert audit.frame_status == "frame_unverified_no_wt_seq"

    def test_repro_a_minus_strand_design_gates_the_mirror_base(self, tmp_path):
        """🔴 A defect NN 251 ITSELF introduced, caught before shipping.

        `verify_window_frame` accepts a MINUS-oriented `wt_seq` (the design tables
        legitimately store either orientation). But `designed_positions` maps string
        index `i` to genomic `window_start + i`, which is only true on the plus
        strand: on the minus strand `wt[i]` is genomic `end - i`.

        So for a minus-strand donor the depth gate and the VAF read land on the
        MIRROR base. A hole at the real edit passes as covered, and a real edit reads
        as WT — the exact A11 class of failure this task exists to close, reintroduced
        by the fix for defect 11.
        """
        from magestic.pipelines.wgs.core.window_table import (
            designed_positions_by_oligo,
            load_windows_audited,
        )

        do = _trace.donor_outcome()
        ref = _seq(77, 4000)
        ws, we = 3001, 3060
        plus_wt = ref[ws - 1 : we]
        true_edit = 3030  # offset 29 from ws
        plus_design = plus_wt[:29] + _other(plus_wt[29]) + plus_wt[30:]
        # Store the SAME design minus-oriented, as a real table may.
        minus_wt = do.reverse_complement(plus_wt)
        minus_design = do.reverse_complement(plus_design)
        tab = write_window_table(
            tmp_path / "w.tsv",
            [("O1", CONTIG, ws, we, minus_wt, minus_design)],
            header=("oligo_name", "chrom", "donor_start_coord", "donor_end_coord",
                    "wt_seq", "design_donor"),
        )
        audit = load_windows_audited(str(tab))
        got = designed_positions_by_oligo(audit, orientations={"O1": "minus"})
        assert got == {"O1": (true_edit,)}, (
            f"minus-oriented design put the designed position at {got} instead of "
            f"{true_edit}; the mirror base is {ws + we - true_edit}, so the depth gate "
            "and the VAF read would land on the wrong base entirely"
        )

    def test_fixed_unverified_orientation_omits_the_oligo(self, tmp_path):
        """Never GUESS the orientation — fall back to the window minimum instead."""
        from magestic.pipelines.wgs.core.window_table import (
            designed_positions_by_oligo,
            load_windows_audited,
        )

        ref = _seq(78, 4000)
        ws, we = 3001, 3060
        wt = ref[ws - 1 : we]
        design = wt[:29] + _other(wt[29]) + wt[30:]
        tab = write_window_table(
            tmp_path / "w.tsv",
            [("O1", CONTIG, ws, we, wt, design)],
            header=("oligo_name", "chrom", "donor_start_coord", "donor_end_coord",
                    "wt_seq", "design_donor"),
        )
        audit = load_windows_audited(str(tab))
        # No orientation supplied at all -> omitted, so `score_sites` gates on the
        # window minimum rather than on a position derived from an assumed strand.
        assert designed_positions_by_oligo(audit) == {}
        # An oligo absent from a supplied mapping -> also omitted.
        assert designed_positions_by_oligo(audit, orientations={"O_other": "plus"}) == {}
        # Verified plus -> used.
        assert designed_positions_by_oligo(audit, orientations={"O1": "plus"}) == {
            "O1": (3030,)
        }

    def test_fixed_designed_positions_derive_from_wt_vs_design(self):
        from magestic.pipelines.wgs.core.window_table import designed_positions

        wt = "AAAACCCCGGGG"
        design = "AAAACGCCGGGG"
        assert designed_positions(3001, wt, design) == (3006,)


# --------------------------------------------------------------------------- #
# Defect 12 — a dropped donor window is a confident undercount                   #
# --------------------------------------------------------------------------- #


class TestDefect12WindowTableDenominators:
    ROWS = [
        ("O1", CONTIG, 3001, 3060),
        ("O2", CONTIG, "3001.0", 3060),  # float-formatted -> .isdigit() False
        ("O3", CONTIG, "", 3060),  # empty (50.2% of the real table)
        ("O4", CONTIG, "3,001", 3060),  # thousands separator
        ("O1", CONTIG, 9001, 9060),  # duplicate oligo: FIRST wins, silently
    ]

    def test_repro_rows_vanish_with_no_denominator(self, tmp_path):
        """REPRODUCTION: 5 rows in, 1 window out, nothing reported.

        CLOSED UPSTREAM (planning/268): the pinned TRACE's loader REFUSES the table at the
        first malformed row, so no row can vanish in silence. `defect_persists` removed;
        `load_windows_audited` (below) still supplies the denominators for the refusal.
        """
        from magestic.pipelines.wgs.core.reconcile_driver import load_oligo_windows

        p = write_window_table(tmp_path / "w.tsv", self.ROWS)
        with pytest.raises(ValueError, match="w.tsv:3: invalid window"):
            load_oligo_windows(str(p))

    def test_fixed_audited_loader_returns_every_denominator(self, tmp_path):
        from magestic.pipelines.wgs.core.window_table import load_windows_audited

        p = write_window_table(tmp_path / "w.tsv", self.ROWS)
        a = load_windows_audited(str(p))
        assert a.n_rows_read == 5
        assert a.n_windows == 1
        assert a.n_dropped_bad_coords == 3
        assert a.n_duplicate_oligos == 1
        assert sorted(a.dropped_oligos) == ["O2", "O3", "O4"]
        assert a.windows == {"O1": (CONTIG, 3001, 3060)}

    def test_fixed_strict_mode_refuses(self, tmp_path):
        from magestic.pipelines.wgs.core.window_table import (
            WindowTableError,
            load_windows_audited,
        )

        p = write_window_table(tmp_path / "w.tsv", self.ROWS)
        with pytest.raises(WindowTableError) as exc:
            load_windows_audited(str(p), strict=True)
        msg = str(exc.value)
        assert "3" in msg and "5" in msg

    def test_fixed_a_clean_table_passes_strict(self, tmp_path):
        from magestic.pipelines.wgs.core.window_table import load_windows_audited

        p = write_window_table(tmp_path / "w.tsv", [("O1", CONTIG, 3001, 3060)])
        a = load_windows_audited(str(p), strict=True)
        assert (a.n_rows_read, a.n_windows, a.n_dropped_bad_coords) == (1, 1, 0)

    def test_fixed_cr_line_endings_do_not_drop_a_row(self, tmp_path):
        """`invisible-line-endings-fake-clean-negatives`: a bare CR must not drop a row."""
        from magestic.pipelines.wgs.core.window_table import load_windows_audited

        p = tmp_path / "cr.tsv"
        p.write_bytes(
            b"oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\n"
            b"O1\tchrXV\t3001\t3060\r\n"
        )
        a = load_windows_audited(str(p))
        assert a.n_windows == 1 and a.n_dropped_bad_coords == 0

    def test_fixed_an_unwindowed_oligo_counts_as_nocov_not_absent(self, tmp_path):
        from magestic.pipelines.wgs.core.window_table import unwindowed_oligos

        p = write_window_table(tmp_path / "w.tsv", [("O1", CONTIG, 3001, 3060)])
        from magestic.pipelines.wgs.core.window_table import load_windows_audited

        a = load_windows_audited(str(p))
        missing = unwindowed_oligos(["O1", "O2", "O3"], a)
        assert missing == ("O2", "O3")


# --------------------------------------------------------------------------- #
# Defect 13 — Stage 3 reports 0 rather than "not measured"                       #
# --------------------------------------------------------------------------- #


class TestDefect13StageThreeNotMeasured:
    def test_repro_zero_is_indistinguishable_from_unmeasured(self, tmp_path):
        """REPRODUCTION: a BAM with no `integrated_plasmid_` contig returns 0/0, the
        same numbers a measured-and-negative sample returns.

        CLOSED UPSTREAM (planning/268): the pinned TRACE now REFUSES to measure a BAM with
        no integrated contig, so it cannot return that 0/0. The magestic-side guard
        below is kept, and reports the unmeasured state without making the call."""
        sv = _trace.sv_plasmid()
        ref = _seq(51, 2000)
        bam = write_bam(tmp_path / "noplasmid.bam", ref, spans=((101, 400),))
        with pysam.AlignmentFile(str(bam), "rb") as b:
            assert [r for r in b.references if r.startswith("integrated_plasmid_")] == []
        with pytest.raises(ValueError, match="integration is unmeasured"):
            sv.detect_plasmid_integration(str(bam))
        from magestic.pipelines.wgs.core.sv_status import integration_measurement

        m = integration_measurement(str(bam))
        assert m.measured is False
        assert m.total is None
        assert m.status == "not_measured_no_integrated_contig"

    def test_fixed_a_measured_sample_is_marked_measured(self, tmp_path):
        from magestic.pipelines.wgs.core.sv_status import integration_measurement

        # planning/268: the pinned TRACE requires a geometrically possible contig, i.e.
        # one longer than both junction arms (2 x 1129), so the fixture is 3000 bp.
        ref = _seq(52, 3000)
        bam = write_bam(
            tmp_path / "withplasmid.bam", ref, contig="integrated_plasmid_o1",
            spans=((101, 400),),
        )
        m = integration_measurement(str(bam))
        assert m.measured is True
        assert m.total == 0  # measured AND zero — a real negative
        assert m.status == "measured"

    def test_fixed_a_cohort_rate_refuses_unmeasured_samples(self):
        from magestic.pipelines.wgs.core.sv_status import (
            IntegrationNotMeasured,
            cohort_integration_rate,
        )

        assert cohort_integration_rate([("a", 3), ("b", 0)]) == pytest.approx(0.5)
        with pytest.raises(IntegrationNotMeasured):
            cohort_integration_rate([("a", 3), ("b", None)])


# --------------------------------------------------------------------------- #
# Defect 14 — four written guards are called by nothing                          #
# --------------------------------------------------------------------------- #


class TestDefect14GuardsWired:
    GUARDS = {
        "require_reference_consistency": "manifest/BAM reference mismatch",
        "incomparable_contigs": "a contig length that differs between references",
        "ContigAlias": "underscore contigs (GATE 3)",
        "n_multiallelic": "A7's multi-allelic inversion",
    }

    @staticmethod
    def _pipeline_text() -> str:
        from magestic.pipelines.wgs import pipeline

        root = Path(pipeline.__file__).parent
        return "\n".join(p.read_text() for p in sorted(root.glob("*.py")))

    @pytest.mark.parametrize("guard", sorted(GUARDS))
    def test_repro_and_fixed_guard_has_a_call_site_in_pipeline(self, guard):
        """REPRODUCTION + FIX in one: all four start with zero call sites under
        `pipeline/`. A guard nobody calls is documentation, not a guard."""
        text = self._pipeline_text()
        assert guard in text, f"{guard} is written and documented but called by nothing"

    def test_fixed_a_real_mismatch_is_reported_before_an_unobserved_one(self):
        """🔴 In a MIXED cohort the mismatch must win.

        Raising `ReferenceUnobservedError` first would hide a genuine wrong-reference
        alignment behind a "could not check" message. Mismatch means the data is
        wrong; unobserved only means the check could not run.
        """
        from magestic.pipelines.wgs.core.manifest import (
            ReferenceMismatchError,
            WGSSample,
            require_reference_consistency,
        )

        good = WGSSample(sample="S1", bam=Path("/x/S1.bam"), reference=Path("/r/ref.fa"))
        bad = WGSSample(sample="S2", bam=Path("/x/S2.bam"), reference=Path("/r/ref.fa"))
        with pytest.raises(ReferenceMismatchError) as exc:
            require_reference_consistency(
                [good, bad], {"S2": "/r/OTHER.fa"}  # S1 unobserved, S2 mismatched
            )
        assert "OTHER.fa" in str(exc.value)
        assert "could not be observed" in str(exc.value)  # mentioned, not raised

    def test_fixed_disjoint_contig_names_are_not_a_clean_pass(self):
        """🔴 An EMPTY intersection must be reported, not silently cleared.

        `incomparable_contigs` compares only shared names, which is right (a composite
        reference legitimately ADDS contigs). But with no shared name the loop finds
        nothing and would return "all clear" for the worst case of all — the
        `chr15`/`chrXV` naming mismatch that made the blacklist 100% inert.
        """
        from magestic.pipelines.wgs.core.manifest import incomparable_contigs

        assert incomparable_contigs({"chrXV": 1091291}, {"chrXV": 1091291}) == []
        # A composite that adds contigs is fine.
        assert incomparable_contigs(
            {"chrXV": 10, "native_2micron": 6318}, {"chrXV": 10}
        ) == []
        # A shared name with a different length is never benign (the CUP1 case).
        assert len(incomparable_contigs({"chrVIII": 10}, {"chrVIII": 12})) == 1
        # ...and a fully disjoint pair is a problem, not a pass.
        problems = incomparable_contigs({"chr15": 10}, {"chrXV": 10})
        assert len(problems) == 1
        assert "no contig name is shared" in problems[0]
        # Two empty mappings are not a naming problem, just nothing to compare.
        assert incomparable_contigs({}, {}) == []

    def test_fixed_missing_pg_line_is_its_own_error(self):
        """A BAM with no `@PG ID:bwa` must not read as a mismatch."""
        from magestic.pipelines.wgs.core.manifest import (
            ReferenceUnobservedError,
            require_reference_consistency,
        )
        from magestic.pipelines.wgs.core.manifest import WGSSample

        s = WGSSample(sample="S1", bam=Path("/x/S1.bam"), reference=Path("/r/ref.fa"))
        with pytest.raises(ReferenceUnobservedError):
            require_reference_consistency([s], {})


# --------------------------------------------------------------------------- #
# Defect 15 — `confidence` never reflects edit evidence                          #
# --------------------------------------------------------------------------- #


class TestDefect15EditEvidence:
    @defect_persists("15 (`confidence` never reflects edit evidence)",
                     closed_by="emit.RECONCILE_COLUMNS + edit_profile.edit_evidence")
    def test_repro_confidence_is_independent_of_the_edit_evidence(self):
        """REPRODUCTION: a 1-read site and a 400-read site both get `high`."""
        thin = TR.SiteCall(TR.EditSite(CONTIG, 3001, 3060, "O1"), 1.0, 4, TR.SITE_FULL)
        thick = TR.SiteCall(TR.EditSite(CONTIG, 3001, 3060, "O1"), 1.0, 400, TR.SITE_FULL)
        a = TR.reconcile_cell("S1", [thin], n_genome_bc=1, bc_confidence="high")
        b = TR.reconcile_cell("S2", [thick], n_genome_bc=1, bc_confidence="high")
        assert a.confidence == b.confidence == "high"
        assert "confidence" in TR.EMIT_COLUMNS
        assert "edit_evidence" in TR.EMIT_COLUMNS, (
            "the column a reader would take for edit confidence is purely the bc1 "
            "confidence; nothing records the edit evidence"
        )

    def test_fixed_edit_evidence_is_its_own_explicit_field(self):
        from magestic.pipelines.wgs.core.edit_profile import edit_evidence

        thin = TR.SiteCall(TR.EditSite(CONTIG, 3001, 3060, "O1"), 1.0, 4, TR.SITE_FULL)
        thick = TR.SiteCall(TR.EditSite(CONTIG, 3001, 3060, "O1"), 1.0, 400, TR.SITE_FULL)
        nocov = TR.SiteCall(TR.EditSite(CONTIG, 3001, 3060, "O1"), 0.0, 1, TR.SITE_NOCOV)
        assert edit_evidence([thin]) == "thin"
        assert edit_evidence([thick]) == "strong"
        assert edit_evidence([nocov]) == "none"
        assert edit_evidence([]) == "no_sites_scored"

    def test_fixed_shipped_row_carries_it_and_names_bc_confidence(self):
        from magestic.pipelines.wgs.core.emit import RECONCILE_COLUMNS

        assert "edit_evidence" in RECONCILE_COLUMNS
        assert "bc_confidence" in RECONCILE_COLUMNS
        # the ambiguous bare name must not be the only confidence column shipped
        assert RECONCILE_COLUMNS.index("bc_confidence") >= 0


# --------------------------------------------------------------------------- #
# Defect 16 — nothing in the repo executes the caller on a BAM                   #
# --------------------------------------------------------------------------- #


class TestDefect16EndToEndOnABam:
    WIN = (3001, 3060)

    def test_end_to_end_clean_edit_and_clean_wt(self, tmp_path):
        """The one test that would have caught A1, A10 and A1c: the caller, on a BAM.

        planning/268: run twice, once with the intended allele declared (TRACE's
        allele-specific route) and once without it (the read-evidence route, because an
        intended edit is not required to assign confidence, Kevin 2026-10-24). Both must
        give the same labels.
        """
        import json

        from magestic.pipelines.wgs.core.reconcile_driver import (
            load_oligo_windows,
            score_sites,
        )

        ref = _seq(61, 4000)
        fasta = write_fasta(tmp_path / "g.fa", {CONTIG: ref})
        ws, we = self.WIN
        allele = json.dumps([{"pos": 3030, "ref": ref[3029], "alt": _other(ref[3029])}])
        tab = write_window_table(
            tmp_path / "w.tsv", [("O1", CONTIG, ws, we, allele)],
            header=("oligo_name", "chrom", "donor_start_coord", "donor_end_coord",
                    "intended_variants"),
        )
        ow = load_oligo_windows(str(tab))
        bare = load_oligo_windows(
            str(write_window_table(tmp_path / "bare.tsv", [("O1", CONTIG, ws, we)]))
        )
        edited = write_bam(tmp_path / "e.bam", ref, spans=(self.WIN,),
                           edits={3030: _other(ref[3029])})
        wt = write_bam(tmp_path / "w.bam", ref, spans=(self.WIN,))
        got = {}
        for name, p in (("edited", edited), ("wt", wt)):
            with pysam.AlignmentFile(str(p), "rb") as bam, pysam.FastaFile(
                str(fasta)
            ) as fa:
                calls = score_sites(bam, fa, ["O1"], ow)
                got[name] = TR.reconcile_cell(
                    name, calls, n_genome_bc=1, bc_confidence="high"
                ).label
                undeclared = score_sites(bam, fa, ["O1"], bare)
                assert TR.reconcile_cell(
                    name, undeclared, n_genome_bc=1, bc_confidence="high"
                ).label == got[name]
        assert got == {"edited": TR.SINGLE_EDIT, "wt": TR.NO_EDIT}

    def test_end_to_end_audited_path_agrees_and_adds_the_guards(self, tmp_path):
        from magestic.pipelines.wgs.core.edit_profile import (
            classify_window,
            window_edit_profile,
        )
        from magestic.pipelines.wgs.core.window_table import (
            designed_positions,
            verify_window_frame,
        )

        ref = _seq(61, 4000)
        fasta = write_fasta(tmp_path / "g.fa", {CONTIG: ref})
        ws, we = self.WIN
        wt_seq = ref[ws - 1 : we]
        design = wt_seq[:29] + _other(wt_seq[29]) + wt_seq[30:]
        dpos = designed_positions(ws, wt_seq, design)
        assert dpos == (3030,)
        edited = write_bam(tmp_path / "e.bam", ref, spans=(self.WIN,),
                           edits={3030: _other(ref[3029])})
        with pysam.AlignmentFile(str(edited), "rb") as bam, pysam.FastaFile(
            str(fasta)
        ) as fa:
            assert verify_window_frame(fa, CONTIG, ws, we, wt_seq).oriented == "plus"
            prof = window_edit_profile(bam, fa, CONTIG, ws, we)
            call = classify_window(prof, designed_positions=dpos)
        assert call.status == TR.SITE_FULL
        assert call.variant_class == "SUBSTITUTION"
        assert call.depth_at_designed >= 4
        assert call.reason == "designed_position_full_vaf"
