"""NN 255 — the blockers NN 252's adversarial panel found, each with its failing-first test.

Every test here FAILED before its fix and the planning doc
(`planning/222_triumph_embl_library_design/255_wgs_blocker_fixes.md`) records the
measured before/after. Ordering mirrors the doc: by whether a wrong answer reaches a
human.

Two house rules this file is deliberate about:

* **A null result needs a positive control.** No test asserts "0 violations" without a
  known-real positive detected by the same code path in the same test.
* **A detector returning 0 may mean it was disabled.** Every guard is exercised on the
  DEFAULT code path, not only with its switch turned on — NN 252's F11b was exactly a
  guard that fired only on a non-default branch.
"""

from __future__ import annotations

import builtins
import importlib
import shutil
import subprocess
import symtable
from pathlib import Path

import pytest

from .conftest import WGS_PKG, import_wgs

# --------------------------------------------------------------------------- #
# F0 — an unbound name on the live path (found while enumerating F1)           #
# --------------------------------------------------------------------------- #


def _wgs_modules(wgs_src_dir: Path) -> list[str]:
    """Every importable module in the package, as dotted suffixes."""
    out = []
    for path in sorted(wgs_src_dir.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        rel = path.relative_to(wgs_src_dir).with_suffix("")
        parts = [p for p in rel.parts if p != "__init__"]
        out.append(".".join(parts))
    return out


# WITHDRAWN (planning/268): "step_2 binds `expand_reference_blocks`". The merge retired the
# gVCF block expansion from step_2: the pinned TRACE reads `INFO/END` itself and validates
# every gVCF REF, which the expander's copied block REF failed. The unbound-name class this
# test was written for is still guarded package-wide by the test below.


def test_f0_no_unbound_global_name_anywhere_in_the_package(wgs_src_dir):
    """The sweep that WOULD have caught F0, kept as a tripwire for the next one.

    `symtable` gives the names each scope reads as globals; anything not bound at
    module scope and not a builtin is an `F821` waiting for the right input. Two were
    present before NN 255 (`expand_reference_blocks` in `step_2`, `Mapping` in
    `core/window_table.py`), and the second was invisible only because
    `from __future__ import annotations` defers annotation evaluation.

    POSITIVE CONTROL in the same assertion style: `_known_bad` is compiled through the
    identical code path and must be reported, so a green result here cannot mean the
    scan silently found nothing.
    """

    def unbound(src: str, name: str, module_globals: set[str]) -> list[str]:
        table = symtable.symtable(src, name, "exec")
        top = set(table.get_identifiers())
        bad: list[str] = []
        stack = list(table.get_children())
        while stack:
            scope = stack.pop()
            stack.extend(scope.get_children())
            for sym in scope.get_symbols():
                if not sym.is_global() or sym.is_assigned():
                    continue
                n = sym.get_name()
                if n in top or n in module_globals or hasattr(builtins, n):
                    continue
                bad.append(f"{name}:{n}")
        return bad

    assert unbound(
        "def f():\n    return _definitely_not_defined_anywhere\n", "_known_bad", set()
    ) == ["_known_bad:_definitely_not_defined_anywhere"], (
        "the positive control was not detected, so this scan proves nothing"
    )

    offenders: list[str] = []
    for suffix in _wgs_modules(wgs_src_dir):
        try:
            mod = importlib.import_module(f"{WGS_PKG}.{suffix}" if suffix else WGS_PKG)
        except Exception:  # pragma: no cover - an unimportable module is another test
            continue
        path = Path(mod.__file__)
        offenders += unbound(path.read_text(), suffix or WGS_PKG, set(dir(mod)))
    assert offenders == [], f"unbound global name(s) on a live path: {offenders}"


# --------------------------------------------------------------------------- #
# F1 — Stage 5 could not run at all                                            #
# --------------------------------------------------------------------------- #


def _config(tmp_path, samples=()):
    """`samples` become manifest rows: since planning/268, step_5 refuses a join whose
    sample set differs from the manifest (the 205 lane's manifest-scoped join)."""
    cfg = import_wgs("config")
    root = tmp_path / "cohort"
    (root / "processed_data").mkdir(parents=True)
    (root / "keyfiles").mkdir(parents=True)
    (root / "intermediate_data").mkdir(parents=True)
    manifest = root / "keyfiles" / "manifest.tsv"
    manifest.write_text("sample\tbam\treference\n" + "".join(
        f"{s}\t{root}/{s}.bam\t{root}/ref.fa\n" for s in samples))
    return cfg.WGSConfig(
        cohort_name="nn255",
        sequencing_date="20260728",
        paths=cfg.WGSPaths(
            project_root=root,
            manifest=manifest,
            intermediate=root / "intermediate_data",
            processed=root / "processed_data",
            keyfiles=root / "keyfiles",
        ),
    )


def _stage_4_audited_table(path: Path, samples=("s1", "s2")) -> list[str]:
    """A Stage-4 table under the header `step_4` actually writes.

    Built from `emit.RECONCILE_COLUMNS` — i.e. from TRACE's own list with NN 251's
    rename applied — rather than from a literal header, so this fixture tracks the
    schema instead of freezing a copy of it.
    """
    emit = import_wgs("core.emit")
    columns = list(emit.RECONCILE_COLUMNS)
    labels = {"s1": "SINGLE_EDIT", "s2": "NO_EDIT"}
    rows = []
    for s in samples:
        row = []
        for c in columns:
            row.append(
                s if c == "sample"
                else labels[s] if c == "label"
                else "high" if c in ("bc_confidence", "confidence")
                else "moderate" if c == "edit_evidence"
                else "0"
            )
        rows.append("\t".join(row))
    path.write_text("\t".join(columns) + "\n" + "\n".join(rows) + "\n")
    return columns


def _stage_2_table(path: Path, samples=("s1", "s2"), *, minor=("", "")) -> None:
    """A Stage-2 table with TRACE's `#`-prefixed header plus step_2's appended columns."""
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    base = ["sample_name", "editing_outcome", "designed_variants",
            "synthetic_errors", "spontaneous_variants"]
    cols = base + list(step_2.OUTCOME_REASON_COLUMNS)
    outcomes = {"s1": "designed_variant_installed", "s2": "designed_variant_unedited"}
    lines = ["# synthetic_error_axis=disabled", "#" + "\t".join(cols)]
    for s, m in zip(samples, minor):
        lines.append("\t".join([
            s, outcomes[s], "chrXV_130_G_to_A", "", "",
            "true", "callable_covered", "10", "10", "0", m,
        ]))
    path.write_text("\n".join(lines) + "\n")


def test_f1_run_step_5_survives_the_renamed_confidence_column(tmp_path, monkeypatch):
    """🔴 HARD BLOCKER. `step_5:149` read `record["confidence"]`; step_4 writes
    `bc_confidence`, so Stage 5 raised `KeyError('confidence')` and the Stage-2 <->
    Stage-4 agreement matrix — the artifact the TRIUMPH rates are read off — could not
    be produced at all.

    Driven through `run_step_5`, not through a re-execution of its join body: NN 252
    reproduced the defect by re-running the body, and NN 251 missed it by grepping.
    """
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    step_4 = import_wgs("pipeline.step_4_reconcile")
    step_5 = import_wgs("pipeline.step_5_join_outcomes")
    config = _config(tmp_path, samples=("s1", "s2"))

    columns = _stage_4_audited_table(step_4.table_path(config))
    assert "bc_confidence" in columns and "confidence" not in columns, (
        "fixture is not exercising the rename"
    )
    _stage_2_table(step_2.table_path(config))

    written = step_5.run_step_5(config)
    rows = [ln.split("\t") for ln in written[0].read_text().splitlines()]
    header, data = rows[0], rows[1:]
    assert len(data) == 2
    assert header[header.index("reconcile_confidence")] == "reconcile_confidence"
    assert data[0][header.index("reconcile_confidence")] == "high"


def test_f1_legacy_unrenamed_confidence_still_joins(tmp_path):
    """The transition, not a flag day: a table written by the OLD writer still joins."""
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    step_4 = import_wgs("pipeline.step_4_reconcile")
    step_5 = import_wgs("pipeline.step_5_join_outcomes")
    config = _config(tmp_path, samples=("s1", "s2"))
    _stage_4_audited_table(step_4.table_path(config))
    legacy = step_4.table_path(config)
    legacy.write_text(legacy.read_text().replace("bc_confidence", "confidence", 1))
    _stage_2_table(step_2.table_path(config))
    written = step_5.run_step_5(config)
    assert len(written[0].read_text().splitlines()) == 3


def test_f1_a_missing_confidence_column_FAILS_LOUDLY_not_as_an_empty_string(tmp_path):
    """Neither name present must raise, naming both candidates.

    A `.get(name, "")` would have published a column of empty strings, which is the
    silent-wrong-answer class rather than a blocked deliverable.
    """
    common = import_wgs("pipeline._common")
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    step_4 = import_wgs("pipeline.step_4_reconcile")
    step_5 = import_wgs("pipeline.step_5_join_outcomes")
    config = _config(tmp_path, samples=("s1", "s2"))
    _stage_4_audited_table(step_4.table_path(config))
    p = step_4.table_path(config)
    p.write_text(p.read_text().replace("bc_confidence", "something_else", 1))
    _stage_2_table(step_2.table_path(config))
    with pytest.raises(common.StepIncompleteError) as exc:
        step_5.run_step_5(config)
    assert "bc_confidence" in str(exc.value) and "confidence" in str(exc.value)


def test_f1_every_renamed_stage_4_column_has_an_alias_entry():
    """The structural tripwire for the NEXT rename, not just this one.

    `emit.RECONCILE_COLUMNS` is derived from TRACE's `EMIT_COLUMNS`; any name that
    differs between the two is a rename, and a rename with no `FIELD_ALIASES` entry is
    an F1 waiting to happen. POSITIVE CONTROL: the assertion is proven non-vacuous by
    requiring that at least one rename exists and is covered.
    """
    emit = import_wgs("core.emit")
    trace = import_wgs("core._trace").reconcile()
    step_5 = import_wgs("pipeline.step_5_join_outcomes")
    audited = set(emit.RECONCILE_COLUMNS)
    renamed = [c for c in trace.EMIT_COLUMNS if c not in audited]
    assert renamed, "no rename detected — this guard would be vacuous"
    covered = {n for names in step_5.FIELD_ALIASES.values() for n in names}
    missing = [c for c in renamed if c not in covered]
    assert missing == [], (
        f"Stage-4 column(s) {missing} were renamed but no step_5 FIELD_ALIASES entry "
        "accepts the old name; that is exactly the F1 KeyError"
    )


# --------------------------------------------------------------------------- #
# F11b — defect 1's fix was inert on the input defect 1 describes              #
# --------------------------------------------------------------------------- #

def _working_bcftools() -> str | None:
    """The first candidate that RUNS, not the first that exists.

    🔴 A toolchain guard must execute the binary. `shutil.which("bcftools")` resolves
    `~/miniconda3/bin/bcftools` on this machine, which exits **134** (`libdeflate.so`
    not found) — NN 252 recorded that this broken copy is probably what produced the
    original "bcftools is absent from this toolchain" claim, and the same class of trap
    is documented for SCG's `module load samtools/1.21`. Existence is not availability.
    """
    for cand in (
        "/Users/kevinroy/miniconda3/envs/bcftools/bin/bcftools",
        shutil.which("bcftools"),
    ):
        if not cand or not Path(cand).exists():
            continue
        try:
            probe = subprocess.run(
                [cand, "--version"], capture_output=True, text=True, timeout=30
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if probe.returncode == 0 and "bcftools" in probe.stdout.lower():
            return cand
    return None


BCFTOOLS = _working_bcftools()
_needs_bcftools = pytest.mark.skipif(
    BCFTOOLS is None, reason="no working bcftools on this toolchain"
)

_REF = ("ACGT" * 60)[:200]
_EDIT_POS0 = 99  # 1-based 100


def _bam(td: Path, name: str, *, n_total: int, n_alt: int) -> tuple[Path, Path]:
    """One flat-coverage BAM, `n_alt` of `n_total` reads carrying the alt base."""
    import pysam

    fa = td / "ref.fa"
    if not fa.exists():
        fa.write_text(f">chrXV\n{_REF}\n")
        pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrXV", "LN": len(_REF)}],
              "RG": [{"ID": name, "SM": name}]}
    bam = td / f"{name}.bam"
    alt = "G" if _REF[_EDIT_POS0] != "G" else "T"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for i in range(n_total):
            a = pysam.AlignedSegment()
            a.query_name = f"{name}_{i}"
            start, length = 50, 100
            seq = list(_REF[start:start + length])
            if i < n_alt:
                seq[_EDIT_POS0 - start] = alt
            a.query_sequence = "".join(seq)
            a.flag, a.reference_id, a.reference_start = 0, 0, start
            a.mapping_quality = 60
            a.cigartuples = [(0, length)]
            a.query_qualities = pysam.qualitystring_to_array("I" * length)
            a.set_tag("RG", name)
            out.write(a)
    pysam.index(str(bam))
    return fa, bam


def _call(ref: Path, bams: list[Path], flags: str, out: Path) -> None:
    cmd = (
        f"{BCFTOOLS} mpileup -q 20 -Q 20 -d 500 -a FORMAT/AD,FORMAT/DP -f {ref} "
        + " ".join(str(b) for b in bams)
        + f" | {BCFTOOLS} call {flags} -Ov -o {out}"
    )
    subprocess.run(["bash", "-c", cmd], check=True, capture_output=True, text=True)


def _synthetic_plate(path: Path) -> Path:
    """A diploid-default plate VCF: a 10% `GT 0/1`, a 33% `GT 0/1`, a clean `1/1`.

    This is the shape NN 252 measured on the real binary (`call -mv` emits a 10% allele
    as `0/1`), written by hand so the assertion does not need bcftools.
    """
    with path.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
                 "\tsHOM\ts33\ts10\n")
        fh.write(
            "chrXV\t100\t.\tT\tG\t50\t.\tDP=120\tGT:PL:DP:AD"
            "\t1/1:219,90,0:40:0,40"      # clean full call  -> NEGATIVE control
            "\t0/1:142,0,184:40:27,13"    # 33% minor allele -> must be flagged
            "\t0/1:16,0,159:40:36,4\n"    # 10% minor allele -> must be flagged
        )
    return path


def test_f11b_a_carries_alt_minor_allele_is_flagged_on_the_DEFAULT_path(tmp_path):
    """🔴 THE ONE THAT MATTERS. `MAG_AF`/`MAG_MINOR` lived inside `if not carries_alt`.

    On a diploid-default plate VCF — what `vcf_adapter`'s own defect-1 docstring says
    this cohort has — a 10% allele is `GT 0/1`, so `carries_alt` is True and the record
    was passed through untouched: measured `n_rescued_by_af=0`,
    `minor_allele_variants={}`, no tag. The over-call defect 1 exists to catch was
    completely unflagged on exactly the input defect 1 describes.

    Run with `min_allele_fraction=None`, i.e. **the default code path**, because a guard
    that fires only with its switch on is the F11b failure itself.

    POSITIVE and NEGATIVE control in the same table: the 10% and 33% records must be
    flagged, and the clean `1/1` record must be written byte-identically with no tag —
    so a green result cannot mean "the tag is written to everything".
    """
    from magestic.pipelines.wgs.core.vcf_adapter import split_plate_vcf

    res = split_plate_vcf(_synthetic_plate(tmp_path / "plate.vcf"), tmp_path / "wd")

    assert res.n_af_computed == 3, "denominator: AD was usable for all three samples"
    assert res.n_rescued_by_af == 0, "no record was rescued; all three carry the alt"
    assert res.n_flagged_minor == 2
    assert res.n_flagged_minor_carrying_alt == 2, (
        "both minority alleles already named the alt in GT, which is precisely the "
        "class the old `if not carries_alt` branch could never reach"
    )
    assert res.minor_allele_variants == {
        "s33": ("chrXV_100_T_to_G",), "s10": ("chrXV_100_T_to_G",),
    }
    t10 = (tmp_path / "wd" / "vcf" / "s10.fil.vcf").read_text()
    assert "MAG_AF=0.100" in t10 and "MAG_MINOR" in t10
    t33 = (tmp_path / "wd" / "vcf" / "s33.fil.vcf").read_text()
    assert "MAG_AF=0.325" in t33 and "MAG_MINOR" in t33
    # NEGATIVE control, same table, same code path.
    thom = (tmp_path / "wd" / "vcf" / "sHOM.fil.vcf").read_text()
    assert "MAG_" not in thom, (
        "a clean full call must stay byte-identical; tagging everything would make the "
        "flag uninformative"
    )
    assert "sHOM" not in res.minor_allele_variants


@_needs_bcftools
def test_f11b_on_a_REAL_diploid_bcftools_plate_vcf(tmp_path):
    """The same assertion against bcftools 1.21 output rather than a hand-written file.

    NN 252's F11b was measured on exactly this fixture (three samples at 100/33/10% of
    30x, `call -mv`) and got `rescued=0, minor={}`.
    """
    from magestic.pipelines.wgs.core.vcf_adapter import split_plate_vcf

    ref = None
    bams = []
    for name, n_alt in (("sHOM", 40), ("s33", 13), ("s10", 4)):
        ref, bam = _bam(tmp_path, name, n_total=40, n_alt=n_alt)
        bams.append(bam)
    plate = tmp_path / "plate.vcf"
    _call(ref, bams, "-mv", plate)

    res = split_plate_vcf(plate, tmp_path / "wd")
    assert res.n_flagged_minor_carrying_alt >= 2, (
        f"a real `call -mv` plate VCF carrying 10% and 33% minor alleles must flag "
        f"both; got {res.n_flagged_minor_carrying_alt} "
        f"(minor_allele_variants={dict(res.minor_allele_variants)})"
    )
    assert "s10" in res.minor_allele_variants and "s33" in res.minor_allele_variants
    assert "sHOM" not in res.minor_allele_variants  # negative control, real data
    assert "MAG_MINOR" in (tmp_path / "wd" / "vcf" / "s10.fil.vcf").read_text()


def test_f11b_rescue_needs_a_READ_COUNT_floor_not_only_a_fraction(tmp_path):
    """A fraction threshold with no count floor selects for low coverage.

    At this cohort's 14.4x median depth a 5% fraction floor rounds to ONE read, whose
    per-position false-positive rate from sequencing error alone is ~1-2%. The floor is
    counted, never silently applied.
    """
    from magestic.pipelines.wgs.core.vcf_adapter import split_plate_vcf

    p = tmp_path / "plate.vcf"
    with p.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tone\ttwo\n")
        # `one` has 1 alt read of 14 (7.1%); `two` has 2 of 14 (14.3%). Both clear a
        # 0.05 fraction floor; only `two` clears a 2-read count floor.
        fh.write("chrXV\t100\t.\tT\tG\t50\t.\t.\tGT:AD"
                 "\t0:13,1\t0:12,2\n")
    res = split_plate_vcf(p, tmp_path / "wd", min_allele_fraction=0.05)
    assert res.n_rescued_by_af == 1, "only the 2-read record may be rescued"
    assert res.n_dropped_below_min_alt_reads == 1
    assert res.n_records_out["one"] == 0 and res.n_records_out["two"] == 1
    # POSITIVE control: lowering the floor to 1 read admits it, so the floor is the
    # thing doing the work and not an unrelated filter.
    res1 = split_plate_vcf(p, tmp_path / "wd1", min_allele_fraction=0.05, min_alt_reads=1)
    assert res1.n_rescued_by_af == 2 and res1.n_dropped_below_min_alt_reads == 0


@_needs_bcftools
def test_f11b_caller_ploidy_is_READ_FROM_THE_FILE_not_declared(tmp_path):
    """`##bcftools_callCommand=` answers the question NN 251 answered in a docstring.

    Measured here on real output for both callers, and across a `bcftools merge`, since
    this cohort's plate VCFs are described as 96-sample merges.
    """
    from magestic.pipelines.wgs.core.vcf_adapter import caller_provenance

    ref, bam = _bam(tmp_path, "sA", n_total=40, n_alt=40)
    hap, dip, allsites = (tmp_path / "h.vcf", tmp_path / "d.vcf", tmp_path / "a.vcf")
    _call(ref, [bam], "-mv --ploidy 1", hap)
    _call(ref, [bam], "-mv", dip)
    _call(ref, [bam], "-m --ploidy 1 -A", allsites)

    h = caller_provenance(hap)
    assert (h.ploidy, h.variants_only, h.n_call_commands) == (1, True, 1)
    assert h.recovers_private_minor_alleles is False
    d = caller_provenance(dip)
    assert (d.ploidy, d.variants_only) == (2, True), (
        "no `--ploidy` means bcftools' diploid DEFAULT, which is the premise of defect 1"
    )
    a = caller_provenance(allsites)
    assert (a.ploidy, a.variants_only, a.keep_alts) == (1, False, True)
    assert a.recovers_private_minor_alleles is True, (
        "this is the only one of the three from which a mixture private to one well is "
        "recoverable at all"
    )
    # A header with no call command is UNOBSERVED, never diploid.
    plain = tmp_path / "plain.vcf"
    _synthetic_plate(plain)
    assert caller_provenance(plain).ploidy is None


def test_f11b_step_2_REFUSES_a_haploid_plate_vcf_with_the_rescue_disabled(tmp_path):
    """At `--ploidy 1` the rescue is the ONLY thing that can see a minority allele.

    So `--ploidy 1` + `min_allele_fraction=None` discards every partial edit while
    looking like a clean haploid run. Detected from the VCF's own header, refused.
    """
    common = import_wgs("pipeline._common")
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    cfg = import_wgs("config")

    plate = tmp_path / "plate.vcf"
    plate.write_text(
        "##fileformat=VCFv4.2\n"
        "##bcftools_callCommand=call -mv --ploidy 1 -Ov -o x.vcf; Date=now\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chrXV\t100\t.\tT\tG\t50\t.\t.\tGT:AD\t0:36,4\n"
    )
    config = _config(tmp_path)
    config.donor_precall = cfg.WGSDonorPrecallConfig(
        plate_vcfs=[plate], min_allele_fraction=None
    )
    with pytest.raises(common.StepIncompleteError, match="ploidy 1"):
        step_2._require_recoverable_caller(config, plate)

    # POSITIVE control: with the rescue enabled the same file is accepted, so the
    # refusal is about the COMBINATION and not about the file.
    config.donor_precall = cfg.WGSDonorPrecallConfig(
        plate_vcfs=[plate], min_allele_fraction=0.05
    )
    step_2._require_recoverable_caller(config, plate)


def test_f11b_a_plate_with_two_different_call_commands_is_refused(tmp_path):
    """`bcftools merge` keeps only the first input's command line, so a plate assembled
    from differently-called inputs would hide the disagreement. Two lines is a refusal."""
    common = import_wgs("pipeline._common")
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    config = _config(tmp_path)
    plate = tmp_path / "mixed.vcf"
    plate.write_text(
        "##fileformat=VCFv4.2\n"
        "##bcftools_callCommand=call -mv --ploidy 1 -Ov -o a.vcf; Date=now\n"
        "##bcftools_callCommand=call -mv -Ov -o b.vcf; Date=now\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
    )
    with pytest.raises(common.StepIncompleteError, match="callCommand"):
        step_2._require_recoverable_caller(config, plate)


# --------------------------------------------------------------------------- #
# F5 / F6 / F6b — WITHDRAWN with the GATE-3 alias wiring (planning/268)         #
# --------------------------------------------------------------------------- #
# The four `_require_alias_consistency` tests guarded a ContigAlias encode/decode that
# step_2 no longer performs: the pinned TRACE builds variant ids from structured fields,
# so an underscore-bearing contig needs no alias (its crash-pins in test_known_answer_242
# now assert that). `core/contig_alias.py` and its own tests are unchanged.


def _frame_fixture(tmp_path: Path, *, wt_matches: bool):
    """A one-window table whose `wt_seq` either matches the reference or does not."""
    import pysam

    ref = tmp_path / "ref.fa"
    ref.write_text(">chrXV\n" + "ACGTACGTAC" * 3 + "\n")
    pysam.faidx(str(ref))
    wt = "ACGTACGTAC" if wt_matches else "TTTTTTTTTT"
    table = tmp_path / "windows.tsv"
    table.write_text(
        "oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\n"
        f"o1\tchrXV\t1\t10\t{wt}\tACGTACGTAA\n"
    )
    return ref, table


def test_f10_require_frame_check_False_now_SUPPRESSES_a_frame_MISMATCH(tmp_path, caplog):
    """🔴 BLOCKER. `verify_window_frame` was called UNWRAPPED, so `False` suppressed the
    refusal only for a MISSING/EMPTY `wt_seq`. A `wt_seq` that was present and disagreed
    with the reference still raised `FrameMismatchError` out of step_4 — so §7.2's
    documented escape hatch did not exist at any switch setting, and a cohort whose
    `wt_seq` legitimately diverges could not be run.

    The hatch does NOT mean the call is trusted: the oligo is left out of
    `orientations`, so `designed_positions_by_oligo` omits it and the gate falls back to
    the conservative window minimum rather than gating at a base whose frame is wrong.
    """
    import logging

    import magestic.pipelines.wgs.pipeline.step_4_reconcile as s4
    cfg = import_wgs("config")
    wt = import_wgs("core.window_table")

    ref, table = _frame_fixture(tmp_path, wt_matches=False)
    audit = wt.load_windows_audited(str(table))
    record = type("R", (), {"sample": "s1", "reference": str(ref)})()

    # (a) require_frame_check=True still REFUSES — the positive control, same fixture.
    config = _config(tmp_path)
    config.reconcile = cfg.WGSReconcileConfig(require_frame_check=True)
    with pytest.raises(wt.FrameMismatchError):
        s4._check_reference_frame(config, record, audit)

    # (b) require_frame_check=False accepts, loudly, and discards the designed positions.
    config.reconcile = cfg.WGSReconcileConfig(require_frame_check=False)
    with caplog.at_level(logging.WARNING):
        orientations = s4._check_reference_frame(config, record, audit)
    assert orientations == {}, (
        "an oligo whose frame FAILED must not contribute an orientation; otherwise the "
        "gate would be placed at a base whose frame is known to be wrong"
    )
    assert "NOT VERIFIED" in caplog.text
    assert s4._FRAME_AUDIT["s1"]["n_frame_mismatch_accepted"] == 1
    assert wt.designed_positions_by_oligo(audit, orientations) == {}, (
        "and the omission must propagate, so score_sites falls back to the window minimum"
    )


def test_f10_a_matching_frame_is_still_verified_and_gates_at_the_designed_base(tmp_path):
    """POSITIVE CONTROL for F10: the hatch must not have broken the working path."""
    import magestic.pipelines.wgs.pipeline.step_4_reconcile as s4
    cfg = import_wgs("config")
    wt = import_wgs("core.window_table")

    ref, table = _frame_fixture(tmp_path, wt_matches=True)
    audit = wt.load_windows_audited(str(table))
    record = type("R", (), {"sample": "s2", "reference": str(ref)})()
    config = _config(tmp_path)
    config.reconcile = cfg.WGSReconcileConfig(require_frame_check=True)
    orientations = s4._check_reference_frame(config, record, audit)
    assert orientations == {"o1": "plus"}
    assert s4._FRAME_AUDIT["s2"]["n_frame_mismatch_accepted"] == 0
    assert s4._FRAME_AUDIT["s2"]["n_verified"] == 1
    assert wt.designed_positions_by_oligo(audit, orientations) == {"o1": (10,)}


# --------------------------------------------------------------------------- #
# F7 — the mixture signal never arrived at the cross-stage table               #
# --------------------------------------------------------------------------- #


def test_f7_every_step_2_column_reaches_step_5():
    """🔴 The tripwire the brief asked for: a NEW step-2 column with no step-5 consumer
    must fail this test.

    NN 251 appended six columns and named `minor_allele_variants` as "the only signal
    that it was a mixture"; NN 252 measured that **none of the six appeared anywhere in
    `step_5_join_outcomes.py`**, so the signal died at the join and a reader of the
    cross-stage table saw a bare `designed_variant_installed` for a 10% allele.
    """
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    step_5 = import_wgs("pipeline.step_5_join_outcomes")
    missing = [
        c for c in step_2.OUTCOME_REASON_COLUMNS if c not in step_5.join_columns()
    ]
    assert missing == [], (
        f"step_2 column(s) {missing} have no step_5 consumer. Add them to the join "
        "(step_5.precall_passthrough derives the list FROM step_2, so this should be "
        "impossible unless the derivation was broken)."
    )
    assert step_2.OUTCOME_REASON_COLUMNS, "the guard would be vacuous with no columns"


def test_f7_a_minor_allele_at_the_designed_position_reaches_the_joined_table(tmp_path):
    """The signal arrives, AND it moves `precall_state` EDIT -> mixture_signal with a reason.

    planning/268: master's collapsed `mixed` state was split by the 205 lane into
    `partial_install` and `mixture_signal`; a minority allele fraction is the second."""
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    step_4 = import_wgs("pipeline.step_4_reconcile")
    step_5 = import_wgs("pipeline.step_5_join_outcomes")
    config = _config(tmp_path, samples=("s1", "s2"))
    _stage_4_audited_table(step_4.table_path(config))
    # s1 is `designed_variant_installed` at chrXV_130_G_to_A and carries a MAG_MINOR
    # allele AT THAT POSITION; s2 is unedited and carries none.
    _stage_2_table(step_2.table_path(config), minor=("chrXV_130_G_to_A", ""))

    written = step_5.run_step_5(config)
    rows = [ln.split("\t") for ln in written[0].read_text().splitlines()]
    header, data = rows[0], rows[1:]
    assert "minor_allele_variants" in header and "outcome_reason" in header
    s1 = data[0]
    assert s1[header.index("minor_allele_variants")] == "chrXV_130_G_to_A"
    assert s1[header.index("editing_outcome")] == "designed_variant_installed", (
        "the raw outcome is passed through verbatim; the downgrade is a separate column"
    )
    assert s1[header.index("precall_state")] == "mixture_signal"
    assert s1[header.index("precall_state_reason")] == (
        step_5.PRECALL_STATE_MINOR_ALLELE
    )
    # NEGATIVE control in the same table: the unedited row is untouched.
    s2 = data[1]
    assert s2[header.index("precall_state")] == "no_edit"
    assert s2[header.index("precall_state_reason")] == step_5.PRECALL_STATE_FROM_OUTCOME


def test_f7_a_minority_allele_ELSEWHERE_does_not_downgrade_a_clean_edit(tmp_path):
    """`split_plate_vcf` scans the whole plate VCF, so `minor_allele_variants` is
    genome-wide. A clean 100% designed edit in a cell that also carries a 30%
    SPONTANEOUS variant on another contig must NOT be downgraded — the downgrade is
    restricted to the intersection with that row's own `designed_variants`."""
    step_2 = import_wgs("pipeline.step_2_donor_precall")
    step_4 = import_wgs("pipeline.step_4_reconcile")
    step_5 = import_wgs("pipeline.step_5_join_outcomes")
    config = _config(tmp_path, samples=("s1", "s2"))
    _stage_4_audited_table(step_4.table_path(config))
    _stage_2_table(step_2.table_path(config), minor=("chrIV_9999_C_to_T", ""))
    written = step_5.run_step_5(config)
    rows = [ln.split("\t") for ln in written[0].read_text().splitlines()]
    header, data = rows[0], rows[1:]
    assert data[0][header.index("precall_state")] == "edit"
    assert data[0][header.index("precall_state_reason")] == (
        step_5.PRECALL_STATE_FROM_OUTCOME
    )
    # ...and the id is still CARRIED, so the reader can see it.
    assert data[0][header.index("minor_allele_variants")] == "chrIV_9999_C_to_T"


def test_f7_the_minor_allele_id_matches_the_designed_variant_id_namespace():
    """The intersection only works if the two id strings are built the same way.

    `split_plate_vcf` writes `{CHROM}_{POS}_{REF}_to_{ALT}` from the RAW plate VCF (it
    runs BEFORE `alias_vcf_contigs` in `step_2._run_one_sample`), and TRACE's ids are
    decoded back to raw names by `_decode_outcome_table`, so both land in the same
    namespace. Pinned here because a mismatch would make the downgrade silently never
    fire — a detector returning zero because it was disabled.
    """
    import ast

    step_2 = import_wgs("pipeline.step_2_donor_precall")
    from magestic.pipelines.wgs.core.vcf_adapter import split_plate_vcf

    # (a) Behavioural: the rescued id carries the RAW, underscore-bearing contig name.
    raw = "blank_MAGESTIC3_pF835_guide_donor_plasmid"
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        plate = td / "p.vcf"
        plate.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
            f"{raw}\t130\t.\tG\tA\t50\t.\t.\tGT:AD\t0:36,4\n"
        )
        res = split_plate_vcf(plate, td / "wd", min_allele_fraction=0.05)
    got = res.minor_allele_variants["s1"]
    assert got == (f"{raw}_130_G_to_A",), got

    # (b) step_2 reads those ids back from the split VCF it shipped, unchanged, and
    # (c) no longer aliases, so TRACE's ids carry the same raw contig name. The alias
    # half of the original claim (decode before the downgrade) is WITHDRAWN with the
    # alias wiring (planning/268).
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        plate = td / "p.vcf"
        plate.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
            f"{raw}\t130\t.\tG\tA\t50\t.\t.\tGT:AD\t0:36,4\n"
        )
        res = split_plate_vcf(plate, td / "wd", min_allele_fraction=0.05)
        assert step_2._minor_alleles_from_vcf(Path(res.vcf_paths["s1"])) == got
    tree = ast.parse(Path(step_2.__file__).read_text())
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == "_run_one_sample")
    names = {n.func.id if isinstance(n.func, ast.Name) else n.func.attr
             for n in ast.walk(fn) if isinstance(n, ast.Call)
             and (isinstance(n.func, ast.Name) or hasattr(n.func, "attr"))}
    assert "_minor_alleles_from_vcf" in names
    assert not names & {"alias_vcf_contigs", "alias_fasta", "_decode_outcome_table"}


# --------------------------------------------------------------------------- #
# F2 — a refusal that READ as comprehensive                                    #
# --------------------------------------------------------------------------- #

#: NN 252's exact 12-key denominator. 10 of 12 passed the old guard.
_SV_RATE_KEYS_252 = (
    "large_deletion_rate", "translocation_rate", "translocation_frequency",
    "inversions_rate", "n_inversions_per_cell", "svRate", "duplications_fraction",
    "rearrangement_rate", "cnv_rate", "indel_gt50bp_rate",
    "sv_pct_cells", "fraction_cells_with_inversion",
)

#: The NEGATIVE control, and it is not optional: these are real, measured quantities in
#: this pipeline (`step_3.SHARD_COLUMNS` carries the first four) and a guard that
#: refused them would be useless. A refusal test with no negative control cannot tell a
#: comprehensive guard from one that refuses everything.
_LEGITIMATE_KEYS = (
    "sample", "upstream_count", "downstream_count", "total_count",
    "n_integrated_contigs", "integration_status", "sv_status", "half_b_status",
    "read_count_mid_point", "frags_genome", "frags_plasmid", "frags_chimeric",
    "sv_plasmid_calls", "n_deletions", "deletion_count", "variant_class",
    "designed_variants", "editing_outcome", "agreement", "n_nocov",
)


@pytest.mark.parametrize("key", _SV_RATE_KEYS_252)
def test_f2_every_plausible_sv_rate_key_is_refused(key):
    """All 12 of NN 252's keys, individually. 10 of them used to pass.

    `large_deletion_rate` is the name a well-intentioned author reaches for, and TRACE
    ships three UNCALLED large-deletion helpers, which is what makes the hole realistic.
    """
    sv = import_wgs("core.sv_status")
    with pytest.raises(sv.SvRateRefused):
        sv.require_no_sv_rate({key: 0.07})


@pytest.mark.parametrize("key", _LEGITIMATE_KEYS)
def test_f2_a_legitimate_key_is_NOT_refused(key):
    """POSITIVE CONTROL for the refusal's specificity, in the same test module.

    Includes `sv_status` itself (the refusal's own output column) and `deletion_count` /
    `n_deletions` — a deletion inside the donor window IS measurable
    (`edit_profile.variant_class` records DELETION), so only a SIZE-qualified deletion
    name is an SV.
    """
    sv = import_wgs("core.sv_status")
    sv.require_no_sv_rate({key: 1})  # must not raise
    sv.require_no_sv_rate_in_columns([key])


def test_f2_a_NESTED_sv_rate_is_refused_and_the_path_is_named():
    """`{'summary': {...sv rate...}}` passed: only top-level keys were inspected.

    A provenance dict is nested by construction — `write_provenance` takes `parameters=`
    and `counts=` sub-mappings — so this is the realistic shape, not an adversarial one.
    """
    sv = import_wgs("core.sv_status")
    nested = {"summary": {"sv" + "_rate": 0.07}, "n_cells": 96}
    with pytest.raises(sv.SvRateRefused) as exc:
        sv.require_no_sv_rate(nested)
    assert "summary." in str(exc.value), "the offending PATH must be reported"
    # ...and a nested legitimate dict still passes.
    sv.require_no_sv_rate({"summary": {"upstream_count": 3}, "n_cells": 96})


def test_f2_a_tsv_HEADER_is_refused_too_and_is_wired_into_step_3_and_step_5():
    """The guard was never applied to a written header, so an SV-rate COLUMN reached a
    file unrefused. Now wired at both write sites."""
    sv = import_wgs("core.sv_status")
    step_3 = import_wgs("pipeline.step_3_sv_plasmid")
    step_5 = import_wgs("pipeline.step_5_join_outcomes")
    with pytest.raises(sv.SvRateRefused, match="column"):
        sv.require_no_sv_rate_in_columns(["sample", "inversions_rate"])
    # The real headers must PASS — the negative control for the wiring.
    sv.require_no_sv_rate_in_columns(step_3.SHARD_COLUMNS)
    sv.require_no_sv_rate_in_columns(step_5.join_columns())
    for mod in (step_3, step_5):
        assert "require_no_sv_rate_in_columns(" in Path(mod.__file__).read_text(), (
            f"{mod.__name__} writes a header without the column guard"
        )


def test_f2_the_docstring_names_ONLY_call_sites_that_exist():
    """The old docstring said "Called from the step-3 and step-5 provenance/count
    paths"; `grep` found exactly one call, in step_3. That is the defect-14 class — a
    guard documented as wired where it is not — so the docstring's claim is now ASSERTED
    against the source rather than trusted."""
    sv = import_wgs("core.sv_status")
    doc = sv.require_no_sv_rate.__doc__ or ""
    named = [m for m in ("step_3_sv_plasmid", "step_5_join_outcomes") if m in doc]
    assert named == ["step_3_sv_plasmid", "step_5_join_outcomes"], (
        "the docstring must name its real call sites"
    )
    pipeline_dir = Path(import_wgs("pipeline.step_3_sv_plasmid").__file__).parent
    for module in named:
        src = (pipeline_dir / f"{module}.py").read_text()
        assert "require_no_sv_rate" in src, (
            f"the docstring claims a call site in {module} that does not exist"
        )


# --------------------------------------------------------------------------- #
# F8 / F9 — hardening. NN 252 marked these SUSPECTED, NOT REPRODUCED for a     #
# live route, so these are invariants, not bug fixes.                          #
# --------------------------------------------------------------------------- #


def _profile(start=101, end=110, depth=40, edit_at=None, vaf=1.0):
    """A synthetic `WindowProfile` with flat depth and an optional edit."""
    ep = import_wgs("core.edit_profile")
    positions = []
    for p in range(start, end + 1):
        nonref = int(round(depth * vaf)) if p == edit_at else 0
        positions.append(ep.PositionCount(
            pos=p, ref_base="A", depth=depth, n_substitution=nonref,
        ))
    best = max(positions, key=lambda p: (p.vaf, p.depth))
    return ep.WindowProfile(
        chrom="chrXV", start=start, end=end, positions=tuple(positions),
        max_vaf=best.vaf, cov_at_max=best.depth, pos_at_max=best.pos,
    )


def test_f8_an_out_of_window_position_RAISES_instead_of_reading_as_wt():
    """🔴 `depth_at`/`vaf_at` returned `0`/`0.0` for an unprofiled position, silently.

    Measured consequence (NN 252): a fully-covered 40x window with one out-of-window
    designed position classified `nocov` with
    `reason='designed_position_depth_below_min_cov'` — **the reason named depth, not a
    coordinate-frame error** — so a coordinate bug was indistinguishable from missing
    data. The invariant `designed_positions ⊆ [profile.start, profile.end]` is one assert
    wide and is now enforced.
    """
    ep = import_wgs("core.edit_profile")
    prof = _profile(edit_at=103)

    with pytest.raises(ep.PositionOutsideWindowError, match="depth_at"):
        prof.depth_at(115)
    with pytest.raises(ep.PositionOutsideWindowError, match="vaf_at"):
        prof.vaf_at(115)
    with pytest.raises(ep.PositionOutsideWindowError, match="variant_class_at"):
        prof.variant_class_at(100)
    with pytest.raises(ep.PositionOutsideWindowError, match="outside the profiled"):
        ep.classify_window(prof, designed_positions=(115,))

    # POSITIVE CONTROL, same profile: an IN-window position still works, and both ends
    # of the window are inside it (an off-by-one here would be the next F8).
    assert prof.contains(101) and prof.contains(110)
    assert not prof.contains(100) and not prof.contains(111)
    assert prof.depth_at(101) == 40 and prof.depth_at(110) == 40
    call = ep.classify_window(prof, designed_positions=(103,))
    assert (call.status, call.vaf, call.cov) == ("full", 1.0, 40)


def test_f8_a_wrong_IN_window_position_is_still_a_confident_no_edit_BY_DESIGN():
    """Honest scope statement, kept as a test so it is not mistaken for closed.

    NN 252 measured that gating at a wrong but IN-window base turns a real 100% edit
    into `wt / vaf=0.0`. The containment assert cannot see that — an in-window position
    IS in the window — and NN 252 could build no live route to it
    (`SUSPECTED, NOT REPRODUCED:`). What defends against it is `verify_window_frame`
    (F10) plus the design-orientation check (F9), not this assert, and this test records
    the residual exposure rather than papering over it.
    """
    ep = import_wgs("core.edit_profile")
    prof = _profile(edit_at=103)
    assert ep.classify_window(prof, designed_positions=(103,)).status == "full"
    wrong = ep.classify_window(prof, designed_positions=(108,))
    assert (wrong.status, wrong.vaf) == ("wt", 0.0)
    assert wrong.reason == "designed_position_wt_vaf"


def test_f9_a_mixed_strand_design_pair_RAISES_instead_of_widening_the_gate():
    """🔴 `verify_window_frame` compares only `wt_seq`; `design_donor`'s strand was never
    checked, and the orientation verified from one column was applied to both.

    NN 252: with `wt_seq` stored minus and `design_donor` stored plus,
    `designed_positions` returned **10 positions — the entire window — where truth is
    1**, which degrades the depth gate back to the window minimum and the VAF back to
    TRACE's window MAXIMUM (re-admitting NN 243's A4b: a background polymorphism
    promoted to "the edit"), with no flag distinguishing it from a narrow gate.
    """
    wt_mod = import_wgs("core.window_table")
    trace_do = import_wgs("core._trace").donor_outcome()

    window = "ACGTACGTAC"
    design_plus = window[:2] + "T" + window[3:]          # SNV at offset 2
    wt_minus = trace_do.reverse_complement(window)       # wt stored on the minus strand

    # The pair as it SHOULD be stored (both minus) is accepted...
    same = wt_mod.verify_design_orientation(
        "o1", wt_minus, trace_do.reverse_complement(design_plus)
    )
    assert same.same_strand and same.n_mismatch_as_stored == 1

    # ...and the mixed pair is refused.
    with pytest.raises(wt_mod.DesignOrientationError, match="opposite strands"):
        wt_mod.verify_design_orientation("o1", wt_minus, design_plus)


def test_f9_a_design_that_differs_at_every_base_is_refused():
    """The second, independent refusal: it catches the measured symptom directly."""
    wt_mod = import_wgs("core.window_table")
    with pytest.raises(wt_mod.DesignOrientationError, match="ALL"):
        wt_mod.verify_design_orientation("o1", "AAAAAAAA", "CCCCCCCC")
    # POSITIVE CONTROL: a heavily RECODED donor (many changes, same strand) must pass —
    # synonymous recoding is a real design pattern in this program and refusing it
    # would be a false alarm.
    recoded = wt_mod.verify_design_orientation(
        "o2", "ATGAAACGTTTAGGCATTCCA", "ATGAAGCGCCTGGGTATCCCG"
    )
    assert recoded.same_strand and recoded.n_mismatch_as_stored >= 6


def test_f9_the_check_is_WIRED_into_designed_positions_by_oligo(tmp_path):
    """A guard reached by nothing is defect 14's class. Driven through the real entry
    point `step_4` uses, not by calling the checker directly."""
    wt_mod = import_wgs("core.window_table")
    trace_do = import_wgs("core._trace").donor_outcome()
    window = "ACGTACGTAC"
    design_plus = window[:2] + "T" + window[3:]
    table = tmp_path / "w.tsv"
    table.write_text(
        "oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\n"
        f"o1\tchrXV\t101\t110\t{trace_do.reverse_complement(window)}\t{design_plus}\n"
    )
    audit = wt_mod.load_windows_audited(str(table))
    with pytest.raises(wt_mod.DesignOrientationError):
        wt_mod.designed_positions_by_oligo(audit, {"o1": "minus"})


# --------------------------------------------------------------------------- #
# F3b — `_pg_reference` returned None, and the skip was silent                  #
# --------------------------------------------------------------------------- #

#: NN 252's 11 real `bwa` `@PG CL` forms, with the two that returned `None`.
_PG_LINES = {
    "plain": ("bwa mem -t 8 /ref/comp.fa r1.fq r2.fq", "/ref/comp.fa"),
    "rg_tab": ("bwa mem -R @RG\\tID:s1\\tSM:s1 /ref/comp.fa r1.fq", "/ref/comp.fa"),
    "fasta": ("bwa mem /ref/comp.fasta r1.fq r2.fq", "/ref/comp.fasta"),
    "fsa": ("bwa mem /ref/comp.chrnames.fsa r1.fq", "/ref/comp.chrnames.fsa"),
    "fna": ("bwa mem /ref/comp.fna r1.fq", "/ref/comp.fna"),
    "rg_is_fa": ("bwa mem -R @RG\\tID:comp.fa\\tSM:s1 /ref/real.fa r1.fq", "/ref/real.fa"),
    "mem2": ("bwa-mem2 mem -t 16 /ref/comp.fa r1.fq r2.fq", "/ref/comp.fa"),
    # the two NN 252 measured as returning None:
    "bgzipped": ("bwa mem /ref/comp.fa.gz r1.fq r2.fq", "/ref/comp.fa.gz"),
    "index_prefix": ("bwa mem -t 8 /ref/compidx r1.fq r2.fq", "/ref/compidx"),
}


@pytest.mark.parametrize("case", sorted(_PG_LINES))
def test_f3b_pg_reference_recovers_the_reference_from_every_real_bwa_form(case):
    """0 of 11 returned a WRONG token (NN 252), but 2 returned `None`, which downgrades
    to a warning — a check that quietly did not run. Both are now recovered."""
    step_4 = import_wgs("pipeline.step_4_reconcile")
    cl, expected = _PG_LINES[case]
    bam = type("B", (), {"header": type("H", (), {
        "to_dict": staticmethod(lambda: {"PG": [{"ID": "bwa", "CL": cl}]})})()})()
    assert step_4._pg_reference(bam) == expected


def test_f3b_an_unrecoverable_pg_is_None_and_the_skip_is_LOUD_AND_COUNTED():
    """`None` must stay `None` rather than become a guess — a WRONG reference token
    would compare the manifest against something that is not a reference, which is worse
    than not comparing. What changes is that the skip is now `NOT VERIFIED:` and counted
    into the step's provenance instead of a passing log line."""
    step_4 = import_wgs("pipeline.step_4_reconcile")
    src = Path(step_4.__file__).read_text()
    # NEGATIVE control: a genuinely unobservable command line still yields None.
    for cl in ("bwa mem --some-unknown-flag value /ref/x r1.fq", "samtools sort -o x.bam"):
        bam = type("B", (), {"header": type("H", (), {
            "to_dict": staticmethod(lambda cl=cl: {"PG": [{"ID": "bwa", "CL": cl}]})})()})()
        assert step_4._pg_reference(bam) is None
    assert "NOT VERIFIED" in src and "reference_name_check_skipped" in src
    assert "_REFERENCE_NAME_CHECK" in src


def test_f11b_a_BCF_is_refused_by_name_not_by_unicode_decode_error(tmp_path):
    """The cohort's plates are called with `-Ob` (measured in their own header), so a
    `.bcf` reaching the adapter is a realistic wiring mistake. Opening one as text
    raises a `UnicodeDecodeError` from inside `gzip` that names neither file nor cause."""
    from magestic.pipelines.wgs.core.vcf_adapter import (
        VcfAdapterError,
        caller_provenance,
    )

    bcf = tmp_path / "plate.bcf"
    bcf.write_bytes(b"BCF\x02\x02\x00\x00\x00\x00binary")
    with pytest.raises(VcfAdapterError, match="bcftools view"):
        caller_provenance(bcf)


@pytest.mark.parametrize(
    "command,expected",
    [
        ("call -mv --ploidy 1 -Ov -o x.vcf", 1),
        ("call -mv --ploidy=1 -Ov -o x.vcf", 1),       # getopt `=` form
        ("call -m --ploidy=2 -Ov -o x.vcf", 2),
        ("call -mv -Ov -o x.vcf", 2),                   # bcftools' diploid DEFAULT
        ("call -mv --ploidy GRCh38 -Ov -o x.vcf", None),  # a preset: UNOBSERVED
        ("call -mv --ploidy-file p.txt -Ov -o x.vcf", 2),  # a different flag entirely
    ],
)
def test_f11b_ploidy_is_parsed_from_BOTH_getopt_forms(tmp_path, command, expected):
    """`--ploidy=1` is valid syntax, and reading only `--ploidy 1` would parse it as the
    diploid default — bypassing the ploidy-1 refusal. Same shape of hole as F2's
    `svRate`, so it is closed the same way: by parsing, not by matching a literal."""
    from magestic.pipelines.wgs.core.vcf_adapter import caller_provenance

    vcf = tmp_path / f"{abs(hash(command))}.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        f"##bcftools_callCommand={command}; Date=now\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
    )
    assert caller_provenance(vcf).ploidy == expected


def test_f10_f3b_the_audit_COUNTS_SURVIVE_a_per_sample_fanout(tmp_path, monkeypatch):
    """🔴 The silent zero this nearly shipped, and the reason it is a sidecar on disk.

    The first version of the F10/F3b fix kept the audit in module-global dicts filled by
    `_run_one_sample` and summed in `_combine`. That would have reported **zero** on the
    real path, twice over: `_run_one_sample` returns early when the shard already exists,
    and the Snakemake fan-out runs every sample in its OWN PROCESS, so `_combine` would
    sum a dict nothing in that process had filled. `frame_unverified_no_wt_seq` would
    have read 0 for a cohort in which not one window was verified.

    This test simulates the process boundary by CLEARING the globals between the
    per-sample write and the read, which is what a fresh process does.
    """
    import magestic.pipelines.wgs.pipeline.step_4_reconcile as s4
    cfg = import_wgs("config")
    wt = import_wgs("core.window_table")

    config = _config(tmp_path)
    config.reconcile = cfg.WGSReconcileConfig(require_frame_check=False)
    ref, table = _frame_fixture(tmp_path, wt_matches=False)
    audit = wt.load_windows_audited(str(table))
    record = type("R", (), {"sample": "s1", "reference": str(ref)})()

    s4._FRAME_AUDIT.clear()
    s4._REFERENCE_NAME_CHECK.clear()
    s4._check_reference_frame(config, record, audit)
    s4._REFERENCE_NAME_CHECK["s1"] = "skipped_pg_unobserved"
    s4._write_audit_sidecar(config, "s1")

    # The process boundary.
    s4._FRAME_AUDIT.clear()
    s4._REFERENCE_NAME_CHECK.clear()

    found, missing = s4._read_audit_sidecars(config, ["s1", "s2"])
    assert missing == ["s2"], (
        "a sample with no sidecar must be reported MISSING, not summed as zero"
    )
    assert found["s1"]["frame_check"]["n_frame_mismatch_accepted"] == 1
    assert found["s1"]["reference_name_check"] == "skipped_pg_unobserved"


def test_f10_the_no_wt_seq_branch_RECORDS_before_returning(tmp_path):
    """The real cohort's case: the design table has no `wt_seq` at all.

    The first fix returned from that branch BEFORE recording, so
    `frame_unverified_no_wt_seq` would have shipped 0 for a run in which every window
    was unverified — a detector reading zero because it was never asked.
    """
    import magestic.pipelines.wgs.pipeline.step_4_reconcile as s4
    cfg = import_wgs("config")
    wt = import_wgs("core.window_table")

    table = tmp_path / "nowt.tsv"
    table.write_text(
        "oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\n"
        "o1\tchrXV\t1\t10\no2\tchrXV\t21\t30\n"
    )
    audit = wt.load_windows_audited(str(table))
    assert not audit.has_wt_seq
    config = _config(tmp_path)
    config.reconcile = cfg.WGSReconcileConfig(require_frame_check=False)
    record = type("R", (), {"sample": "s9", "reference": "unused"})()
    s4._FRAME_AUDIT.clear()
    assert s4._check_reference_frame(config, record, audit) == {}
    rec = s4._FRAME_AUDIT["s9"]
    assert rec["n_no_wt_seq"] == 2 and rec["n_verified"] == 0
    assert rec["status"] == "no_wt_seq_column"
    # POSITIVE CONTROL: require_frame_check=True on the same table still REFUSES while
    # any window is scored from read evidence (declares no intended allele)...
    config.reconcile = cfg.WGSReconcileConfig(require_frame_check=True)
    common = import_wgs("pipeline._common")
    with pytest.raises(common.StepIncompleteError, match="no `wt_seq` column"):
        s4._check_reference_frame(config, record, audit, undeclared={"o1"})
    # ...and planning/268: records instead when every window declares REF-validated
    # intended alleles, because TRACE then checks the frame itself.
    s4._FRAME_AUDIT.clear()
    assert s4._check_reference_frame(config, record, audit, undeclared=set()) == {}
    assert s4._FRAME_AUDIT["s9"]["status"] == "no_wt_seq_column"


def test_f11b_the_recommended_command_is_window_restricted_and_normalized():
    """The recommendation has to survive the two things Stage 2 actually enforces.

    (a) `step_2` REFUSES the whole plate on any multi-allelic record, because TRACE
        reads ALT whole — so the recommended command pipes through `bcftools norm -m-`.
    (b) all-sites genome-wide is ~1.2e7 records per plate and `_require_vcf_inputs`
        re-reads each plate ONCE PER SAMPLE, so the recommendation is `-R` restricted.

    Pinned as a test because a recommendation that only lives in a planning doc drifts
    away from the constant the code ships.
    """
    from magestic.pipelines.wgs.core import vcf_adapter as va

    cmd = va.MPILEUP_CALL_WINDOWS_COMMAND
    assert "-R {regions}" in cmd
    assert "--ploidy 1" in cmd and " -A " in cmd
    assert "norm -m-" in cmd
    # no `-v` anywhere in the CALL stage: that is the whole point.
    call_stage = cmd.split("|")[1]
    assert " -mv" not in call_stage and "--variants-only" not in call_stage
    # ...and the caveat names the restricted form rather than the unrestricted one.
    assert "MPILEUP_CALL_WINDOWS_COMMAND" in va.PLOIDY_CAVEAT
    assert "infeasible" in va.PLOIDY_CAVEAT
