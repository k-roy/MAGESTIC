"""Failing-before tests for the screen step 10/11 defaults repair (planning NN 218).

Design and contract IDs (C1a-C9 for step 10, S1-S4 for step 11): REPAIR_DESIGN.md in this folder.
Defects (keys into 216/deliverables/reviewed_candidates_through_batch09.tsv): SCR-D01a, SCR-D01b,
SCR-D01c, TECAN-C01, MAN-C01, FL-H03, FL-C02, CT02, FL-C03, DESEQ-SCOPE-R1.

Run (root session only, under the shared test lock):

    PYTHONPATH=<checkout>/src /Users/kevinroy/.venvs/magestic/bin/python -m pytest -q <this file>

What runs: the ACTUAL numbered scripts ``10_split_comparisons.py`` and ``11_deseq2_single.py``,
loaded with importlib from the 218 checkout (never from the editable install that points at
~/repos/magestic), through their real ``main()`` with a patched ``sys.argv``, or through the real
``build_sample_info`` helper. Inputs are tiny synthetic TSVs in a neutral pytest temp directory
(``tmp_path_factory.mktemp("w")``, so no test name can leak a token into a logged path). No Oak
path, no network, no pipeline run.

Naming:
    test_ctl_*     CONTROL: a valid, fully specified input. Must PASS on HEAD 7393ae6 and after.
    test_fix_*     Must FAIL on HEAD 7393ae6 (the defect) and PASS after the repair.
    test_enable_*  Must FAIL on HEAD and pass after; tests the recommended exclusion option C9.
                   Drop this group if the root decides not to adopt C9.

Refusal checks never rely on a message token alone: each also asserts a behaviour that the
current code does not have (non-zero exit, untouched earlier manifest, model never built). The
refusal text inspected is only the exception text, stderr and ERROR-level log records, never
INFO logs (step 10 logs its own argument values at INFO, e.g. ``10_split_comparisons.py:406``).

Step 11 imports pydeseq2 at module level (``11_deseq2_single.py:39-41``) and pydeseq2 is not
installed in the magestic venv (``216/preflight/model_runtime_availability.json``). The step 11
tests therefore install a recording stub of ``pydeseq2`` in ``sys.modules`` for the duration of
each test. Every refusal they check happens before the model is built, and the stub proves it
(zero ``DeseqDataSet`` constructions). The one control asserts the model WAS built exactly once
with the expected matrix and ``design_factors="condition"``. The stub fits nothing; no numerical
DESeq2 claim is made here.
"""
import importlib.util
import itertools
import logging
import os
import sys
import types
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pytest

# ----------------------------------------------------------------------------- checkout
# Resolve the checkout from this file's location (reviews/<folder>/ -> 218/code/magestic), or
# from MAGESTIC_CHECKOUT if the file is moved. Never via magestic.__file__: the venv's editable
# install points at ~/repos/magestic (canonical), which is exactly the code this effort replaces.
_OVERRIDE = os.environ.get("MAGESTIC_CHECKOUT")
CHECKOUT = (Path(_OVERRIDE) if _OVERRIDE
            else Path(__file__).resolve().parents[3]).resolve()  # <checkout>/tests/test_pipelines/bc1_from_screen/
SNAKEMAKE = CHECKOUT / "src" / "magestic" / "pipelines" / "bc1_from_screen" / "snakemake"
STEP10 = SNAKEMAKE / "10_split_comparisons.py"
STEP11 = SNAKEMAKE / "11_deseq2_single.py"
assert STEP10.is_file() and STEP11.is_file(), f"numbered scripts not found under {SNAKEMAKE}"

_ALIASES = itertools.count()


def load_script(path):
    """Execute one numbered script as a fresh module (the way the 216 probes load them)."""
    alias = f"screen_defaults_under_test_{path.stem}_{next(_ALIASES)}"
    spec = importlib.util.spec_from_file_location(alias, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ----------------------------------------------------------------------------- fixture data
BC = {k: k * 20 for k in "ACGT"}
ORDER = ["A", "C", "G", "T"]
# Construction-derived bc1 -> transformation. trafo_3 and trafo_4 share cV586/V634/SPS in the real
# yL437 trafo key, so the oligo names and the SPS below are deliberately identical across them.
TRAFO = {"A": "trafo_3", "C": "trafo_3", "G": "trafo_4", "T": "trafo_4"}
OLIGO = {"A": "SpCas9_V634_001", "C": "SpCas9_V634_002", "G": "SpCas9_V634_001", "T": "SpCas9_V634_002"}
SPS = "CGATAGACGACTGGACAGCA"
LIB = "yLtest"
DATE = 20250101
STD_SAMPLES = ["t0_r1", "t0_r2", "d6_r1", "d6_r2"]
# Today's manifest columns, in order (10_split_comparisons.py:298-310). The SLURM template reads
# fields by position (run_deseq2_array.sbatch:66-68), so any new field must be APPENDED.
MANIFEST_COLUMNS = ["library", "yeast_sublibrary", "comparison_id", "condition", "generation",
                    "screen_date", "reference_file", "treatment_file", "n_reference_samples",
                    "n_treatment_samples", "n_bc1s"]
SENTINEL = "SENTINEL manifest from an earlier valid run\n"


def count(key, sample):
    """Deterministic, distinct, non-zero synthetic count for one bc1 in one sample."""
    return 10 * (ORDER.index(key) + 1) + sum(ord(c) for c in sample) % 50


def write_rows(path, rows):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def write_counts(path, samples, keys=ORDER):
    """Wide count matrix in step 08's layout: bc1, tier, then one column per sample (08:605)."""
    rows = [{"bc1": BC[k], "tier": "abundance_high", **{s: count(k, s) for s in samples}}
            for k in keys]
    return write_rows(path, rows)


def write_reference(path, trafo=None, extra_rows=()):
    mapping = dict(TRAFO, **(trafo or {}))
    rows = [dict(bc1=BC[k], oligo_name=OLIGO[k], trafo_ID=mapping[k], sps=SPS) for k in ORDER]
    return write_rows(path, rows + list(extra_rows))


def std_key(date=DATE):
    """A complete, correctly labelled key: two t0 aliquots and two generation-6 drug samples."""
    return [
        dict(sample_name="t0_r1", condition="YPD", generations=0, screen_date=date, library_ID=LIB),
        dict(sample_name="t0_r2", condition="YPD", generations=0, screen_date=date, library_ID=LIB),
        dict(sample_name="d6_r1", condition="drug", generations=6, screen_date=date, library_ID=LIB),
        dict(sample_name="d6_r2", condition="drug", generations=6, screen_date=date, library_ID=LIB),
    ]


def tecan_key():
    """Corrected-TECAN shape: `media` and `Tecan_date`, no `condition`/`screen_date` columns."""
    return [
        dict(sample_name="e1_t0", media="YPD", Tecan_date=20240315, generations=0, sample_ID="E1-0"),
        dict(sample_name="e1_g6", media="EGTA", Tecan_date=20240315, generations=6, sample_ID="E1-6"),
        dict(sample_name="e2_t0", media="YPD", Tecan_date=20240320, generations=0, sample_ID="E2-0"),
        dict(sample_name="e2_g6", media="EGTA", Tecan_date=20240320, generations=6, sample_ID="E2-6"),
    ]


def manual_key():
    """Corrected-manual shape: `media` and `outgrowth_used_for_phenotyping`, two outgrowths."""
    return [
        dict(sample_name="o1_t0", media="YPD", outgrowth_used_for_phenotyping="outgrowth_A",
             generations=0, biological_replicate="rep_1"),
        dict(sample_name="o1_g8", media="Caffeine", outgrowth_used_for_phenotyping="outgrowth_A",
             generations=8, biological_replicate="rep_1"),
        dict(sample_name="o2_t0", media="YPD", outgrowth_used_for_phenotyping="outgrowth_B",
             generations=0, biological_replicate="rep_2"),
        dict(sample_name="o2_g8", media="Caffeine", outgrowth_used_for_phenotyping="outgrowth_B",
             generations=8, biological_replicate="rep_2"),
    ]


def june_key():
    """20220607 shape: pooled-culture t0 aliquots plus one individual-library recovery, no dates.

    The pooled rows carry the same library value as --library (LIB), so these fixtures stay valid
    even if an implementation also checks agreement with --library (C2 itself checks uniqueness).
    """
    return [
        dict(sample_name="pool_t0_a", condition="time_zero_combined_culture", generations=0, library=LIB),
        dict(sample_name="pool_t0_b", condition="time_zero_combined_culture", generations=0, library=LIB),
        dict(sample_name="ind_t0_1", condition="YPD_outgrowth_from_plate_wash", generations=0, library="yL901"),
        dict(sample_name="pool_d6_a", condition="drug", generations=6, library=LIB),
        dict(sample_name="pool_d6_b", condition="drug", generations=6, library=LIB),
    ]


def io_args(counts, key, out):
    return ["--count-matrix", counts, "--sample-key", key, "--library", LIB, "--output-dir", out]


def split_args(ref, column="trafo_ID"):
    args = ["--split-by-yeast-sublibrary", "--bc1-reference-table", ref]
    if column is not None:
        args += ["--yeast-sublibrary-column", column]
    return args


# ----------------------------------------------------------------------------- running mains
@dataclass
class Outcome:
    code: int
    refusal: str  # exception text + stderr + ERROR-level log messages (never INFO)


def run_main(module, argv, monkeypatch, caplog, capsys):
    """Call module.main() with argv; normalise return / SystemExit / uncaught exception."""
    monkeypatch.setattr(sys, "argv", [str(module.__file__), *[str(a) for a in argv]])
    caplog.clear()
    capsys.readouterr()
    raised = ""
    with caplog.at_level(logging.INFO):
        try:
            code = module.main()
        except SystemExit as exc:
            if exc.code is None:
                code = 0
            elif isinstance(exc.code, int):
                code = exc.code
            else:
                code, raised = 1, str(exc.code)
        except Exception as exc:  # a refusal may legitimately surface as an uncaught ValueError
            code, raised = 1, f"{type(exc).__name__}: {exc}"
    err = capsys.readouterr().err
    errors = [r.getMessage() for r in caplog.records if r.levelno >= logging.ERROR]
    return Outcome(0 if code is None else int(code), "\n".join([raised, err, *errors]))


def assert_refused(outcome, *tokens, any_of=()):
    assert outcome.code != 0, f"input was accepted (exit {outcome.code}); it must be refused"
    assert "unrecognized arguments" not in outcome.refusal, (
        "refused by argparse because a new flag does not exist, not by validation:\n" + outcome.refusal)
    for token in tokens:
        assert token in outcome.refusal, f"refusal does not name {token!r}:\n{outcome.refusal}"
    if any_of:
        assert any(t in outcome.refusal for t in any_of), (
            f"refusal names none of {any_of}:\n{outcome.refusal}")


def seed_previous(out):
    out.mkdir(parents=True, exist_ok=True)
    (out / "comparison_manifest.tsv").write_text(SENTINEL)


def assert_nothing_published(out):
    assert (out / "comparison_manifest.tsv").read_text() == SENTINEL, (
        "a refused run overwrote the earlier manifest")
    for sub in ("reference", "treatment"):
        folder = out / sub
        written = sorted(p.name for p in folder.iterdir()) if folder.exists() else []
        assert written == [], f"a refused run published {sub} files: {written}"


def manifest_frame(out):
    frame = pd.read_csv(out / "comparison_manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(frame.columns[: len(MANIFEST_COLUMNS)]) == MANIFEST_COLUMNS, list(frame.columns)
    return frame


def manifest_rows(out):
    return manifest_frame(out)[MANIFEST_COLUMNS].to_dict("records")


def expected_row(stratum, date, condition, generation, n_ref, n_treat, n_bc1):
    """One manifest row exactly as 10_split_comparisons.py:298-310 / :379-391 writes it today."""
    suffix = "" if stratum == "all" else f"_{stratum}"
    body = f"{condition}_{generation}" if date == "pooled" else f"{date}_{condition}_{generation}"
    cid = f"{LIB}{suffix}_{body}"
    return dict(library=LIB, yeast_sublibrary=stratum, comparison_id=cid, condition=condition,
                generation=str(generation), screen_date=str(date),
                reference_file=f"{cid}_reference.tsv", treatment_file=f"{cid}_treatment.tsv",
                n_reference_samples=str(n_ref), n_treatment_samples=str(n_treat), n_bc1s=str(n_bc1))


def split_table(out, sub, name):
    return pd.read_csv(out / sub / name, sep="\t", index_col=0)


def expected_table(keys, samples):
    return pd.DataFrame({s: [count(k, s) for k in keys] for s in samples},
                        index=pd.Index([BC[k] for k in keys], name="bc1"))


def standard_split_rows():
    return [expected_row("trafo_3", DATE, "drug", 6, 2, 2, 2),
            expected_row("trafo_4", DATE, "drug", 6, 2, 2, 2)]


@pytest.fixture(scope="module")
def step10():
    return load_script(STEP10)


@pytest.fixture
def work(tmp_path_factory):
    # Neutral directory name: pytest's tmp_path embeds the test name in the path.
    return tmp_path_factory.mktemp("w")


# ============================================================================= step 10 controls
def test_ctl_split_by_declared_metadata_column_is_unchanged(step10, work, monkeypatch, caplog, capsys):
    """CONTROL for SCR-D01a/b, SCR-D01c, FL-C02/CT02 (contract C1b, C1f, C10).

    A complete key and a construction-derived bc1->transformation column (`trafo_ID`), passed
    through the EXISTING explicit split flags. Today (10:192) a non-oligo_name column is used
    verbatim, so trafo_3 and trafo_4 stay separate even though they share oligo names and SPS.
    The repair must reproduce the same 11 manifest columns in the same order, the same rows and
    the same reference/treatment matrices. Every treatment group has a baseline and every count
    column has key metadata, which is the valid counterpart of the FL-C02/CT02 refusals.
    """
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", std_key())
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert manifest_rows(out) == standard_split_rows()
    for stratum, keys in (("trafo_3", ["A", "C"]), ("trafo_4", ["G", "T"])):
        cid = f"{LIB}_{stratum}_{DATE}_drug_6"
        pd.testing.assert_frame_equal(split_table(out, "reference", f"{cid}_reference.tsv"),
                                      expected_table(keys, ["t0_r1", "t0_r2"]))
        pd.testing.assert_frame_equal(split_table(out, "treatment", f"{cid}_treatment.tsv"),
                                      expected_table(keys, ["d6_r1", "d6_r2"]))


def test_ctl_dateless_key_with_explicit_pooled_dates_is_unchanged(step10, work, monkeypatch, caplog, capsys):
    """CONTROL for TECAN-C01/MAN-C01 and FL-C03 (contract C3): the migration path.

    The real 20210724, 20210822, 20211206 and 20240806 cAMP combined keys have no screen_date
    column, so today they silently take the all-gen0 fallback (10:463-466). After the repair the
    same analysis must still be available by asking for it explicitly with `--no-consolidate`
    (the pooled-dates estimand). Here every generation-0 column is an aliquot of one library, so
    pooling them is the requested estimand, not the FL-C03 mix. Output must match today's.
    """
    rows = [{k: v for k, v in r.items() if k != "screen_date"} for r in std_key()]
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    argv = io_args(counts, key, out) + split_args(ref) + ["--no-consolidate"]
    outcome = run_main(step10, argv, monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert manifest_rows(out) == [expected_row("trafo_3", "pooled", "drug", 6, 2, 2, 2),
                                  expected_row("trafo_4", "pooled", "drug", 6, 2, 2, 2)]


def test_ctl_identical_duplicate_key_rows_are_accepted(step10, work, monkeypatch, caplog, capsys):
    """CONTROL for FL-H03 (contract C5).

    Two key rows for d6_r1 that differ only in a column step 10 does not use (PCR plate) are the
    normal shape of technical replicates (step 09 names its output by the grouping value,
    09:221,253). They must be accepted and give exactly the control output. Today last-row-wins
    (10:151-155) happens to give the same answer because the rows agree.
    """
    rows = [dict(r, PCR_plate="P1") for r in std_key()]
    rows.append(dict(rows[2], PCR_plate="P2"))
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert manifest_rows(out) == standard_split_rows()


@pytest.mark.parametrize("shape", ["tecan", "manual"])
def test_ctl_adapter_columns_keep_named_date_specific_comparisons(step10, work, monkeypatch, caplog, capsys, shape):
    """CONTROL for TECAN-C01 and MAN-C01 (contract C3).

    The corrected TECAN and manual keys carry `media` plus `Tecan_date` or
    `outgrowth_used_for_phenotyping`. The existing adapter flags already map them (the register
    records 48 named manual comparisons this way). The repair must keep this path byte-for-byte:
    one comparison per date/outgrowth, so dose-changing robot dates and the two outgrowths stay
    separate contrasts rather than replicates.
    """
    if shape == "tecan":
        rows, date_col, cond = tecan_key(), "Tecan_date", "EGTA"
        groups = [(20240315, 6), (20240320, 6)]
    else:
        rows, date_col, cond = manual_key(), "outgrowth_used_for_phenotyping", "Caffeine"
        groups = [("outgrowth_A", 8), ("outgrowth_B", 8)]
    samples = [r["sample_name"] for r in rows]
    counts = write_counts(work / "counts.tsv", samples)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    argv = io_args(counts, key, out) + split_args(ref) + [
        "--condition-column", "media", "--screen-date-column", date_col]
    outcome = run_main(step10, argv, monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert manifest_rows(out) == [expected_row(stratum, date, cond, gen, 1, 1, 2)
                                  for stratum in ("trafo_3", "trafo_4") for date, gen in groups]


def test_ctl_numeric_decimal_generation_column_parses_as_integer(step10, work, monkeypatch, caplog, capsys):
    """CONTROL for SCR-D01c (contract C4).

    A clean numeric column written as 0.0/6.0 is read by pandas as float64, and int(6.0) == 6
    already works today (10:135). The repair must keep it: generation 6, not 0 and not "6.0".
    (The string form of "6.0", which does fail today, is test_fix_string_decimal_* below.)
    """
    rows = [dict(r, generations=float(r["generations"])) for r in std_key()]
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert manifest_rows(out) == standard_split_rows()


# ============================================================ step 10: transformation scope
def test_fix_undeclared_transformation_scope_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01a (contract C1a, C8).

    Today the default (no split flag, 10:89) analyses every barcode together and exits 0 with
    yeast_sublibrary=all. Pooling must never be the default: with neither
    `--split-by-yeast-sublibrary` nor `--pool-transformations`, step 10 must refuse, name both
    options, leave the earlier manifest untouched and publish no split files.
    """
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", std_key())
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out), monkeypatch, caplog, capsys)
    assert_refused(outcome, "--split-by-yeast-sublibrary", "--pool-transformations")
    assert_nothing_published(out)


def test_fix_columns_from_two_transformations_are_refused_even_when_rows_are_split(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01a, sample axis (contract C1e).

    The audit fixture (modules/screen_downstream/reproduce_downstream.py:28) labels the t0
    columns trafo_3 and the treatment columns trafo_4 in the key's `yeast_sublibrary`. A row
    split cannot separate columns, so every comparison would contrast one transformation's
    culture against another's. Today this exits 0. It must be refused, naming both values.
    """
    rows = [dict(r, yeast_sublibrary="trafo_3" if r["generations"] == 0 else "trafo_4")
            for r in std_key()]
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, "trafo_3", "trafo_4")
    assert_nothing_published(out)


def test_fix_explicit_pooling_is_honoured_and_recorded(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE (new flag): SCR-D01a (contract C1a, C1e, C1f).

    Pooling across transformations stays available, but only when asked for. With
    `--pool-transformations` the same two-transformation input must produce today's pooled
    comparison unchanged (yeast_sublibrary=all, all 4 bc1s), and the manifest must record the
    choice in an APPENDED `transformation_scope` column equal to "pooled". Today the flag does
    not exist (argparse exit 2).
    """
    rows = [dict(r, yeast_sublibrary="trafo_3" if r["generations"] == 0 else "trafo_4")
            for r in std_key()]
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    out = work / "out"
    argv = io_args(counts, key, out) + ["--pool-transformations"]
    outcome = run_main(step10, argv, monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    frame = manifest_frame(out)
    assert frame[MANIFEST_COLUMNS].to_dict("records") == [expected_row("all", DATE, "drug", 6, 2, 2, 4)]
    assert frame["transformation_scope"].tolist() == ["pooled"]


def test_fix_split_manifest_records_declared_transformation_source(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01a/b (contract C1f).

    Downstream steps cannot tell a pooled analysis from a split one by `yeast_sublibrary` alone
    ("all" is today's silent default). A split run must append `transformation_scope` =
    "split:<declared column>" to every manifest row. Today the column does not exist.
    """
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", std_key())
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert manifest_frame(out)["transformation_scope"].tolist() == ["split:trafo_ID"] * 2


@pytest.mark.parametrize("column, token", [(None, "oligo_name"), ("oligo_name", "oligo_name"), ("sps", "sps")],
                         ids=["default_column", "explicit_prefix", "shared_sps"])
def test_fix_design_pool_identifier_never_stands_in_for_transformation(step10, work, monkeypatch, caplog, capsys, column, token):
    """FAILS BEFORE: SCR-D01b (contract C1b).

    Today the split defaults to the oligo_name prefix (10:96-99, 188-190), so the trafo_3 and
    trafo_4 barcodes, which share V634 designs, collapse into one "SpCas9" stratum; an SPS column
    likewise gives one stratum. Oligo prefixes and SPS identify design pools, never
    transformations. The split must refuse such a column (and the old default), naming it.
    """
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", std_key())
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref, column), monkeypatch, caplog, capsys)
    assert_refused(outcome, token)
    assert_nothing_published(out)


def test_fix_absent_split_column_is_refused_not_an_empty_success(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01b family, silent empty success (contract C1b, C8).

    Today a column missing from the reference table only logs an error and returns {} (10:181-183);
    no stratum is processed and a headerless manifest is written with exit 0 (10:435, 513-525).
    A typo ("trafo_id" for "trafo_ID") must be refused, naming the requested column.
    """
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", std_key())
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref, "trafo_id"), monkeypatch, caplog, capsys)
    assert_refused(outcome, "trafo_id")
    assert_nothing_published(out)


def test_fix_unassigned_barcode_is_not_analysed_as_an_unknown_transformation(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01b (contract C1d; exclude-versus-refuse is Decision D3).

    Today a bc1 with no transformation value becomes the stratum "unknown" (10:190,192) and gets
    its own comparisons as if it were a transformation. After the repair it belongs to no
    stratum: the manifest lists trafo_3 (A, C) and trafo_4 (G only), and T appears in no file.
    """
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", std_key())
    ref = write_reference(work / "ref.tsv", trafo={"T": None})
    out = work / "out"
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    rows = manifest_rows(out)
    assert [r["yeast_sublibrary"] for r in rows] == ["trafo_3", "trafo_4"]
    assert [r["n_bc1s"] for r in rows] == ["2", "1"]
    for sub in ("reference", "treatment"):
        for path in (out / sub).iterdir():
            assert BC["T"] not in pd.read_csv(path, sep="\t", index_col=0).index, path.name


def test_fix_barcode_with_two_transformations_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01b (contract C1c).

    Today dict(zip()) keeps the last value for a duplicated bc1 (10:194), so A silently moves to
    trafo_4. A barcode assigned to two transformations must be refused, naming the barcode;
    choosing one assignment is not step 10's decision.
    """
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", std_key())
    ref = write_reference(work / "ref.tsv", extra_rows=[
        dict(bc1=BC["A"], oligo_name=OLIGO["A"], trafo_ID="trafo_4", sps=SPS)])
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, BC["A"])
    assert_nothing_published(out)


# ============================================================ step 10: generation parsing
@pytest.mark.parametrize("raw", ["none", "bad", None, "6.5", 6.5, 0.5],
                         ids=["none_token", "nonnumeric", "missing", "string_fraction", "float_fraction", "float_half"])
def test_fix_unusable_generation_is_refused_by_helper(step10, raw):
    """FAILS BEFORE: SCR-D01c (contract C4), actual build_sample_info.

    Today NaN, "none" and unparsable strings become generation 0, i.e. a baseline (10:130-137);
    float 6.5 is truncated to 6 and 0.5 to 0 by int() (10:135). The code's contract is integer
    generations (docstring 10:115; manifest isdigit() 10:296), so each must raise ValueError
    naming the sample.
    """
    key = pd.DataFrame({"sample_name": ["t0_x", "odd_x"], "condition": ["YPD", "drug"],
                        "generations": [0, raw], "screen_date": ["20250101", "20250101"]})
    with pytest.raises(ValueError, match="odd_x"):
        step10.build_sample_info(key, ["t0_x", "odd_x"], "sample_name", "condition",
                                 "generations", "screen_date")


def test_fix_string_decimal_generation_parses_as_integer_in_helper(step10):
    """FAILS BEFORE: SCR-D01c (contract C4), actual build_sample_info.

    In an object column, "6.0" reaches int("6.0"), raises, and becomes 0 (10:135-137): a
    generation-6 sample silently turns into a baseline. It must parse to the integer 6, and its
    text form must be "6" so the comparison key f"{condition}_{gen}" still round-trips through the
    manifest's isdigit() parser (10:223, 294-296).
    """
    key = pd.DataFrame({"sample_name": ["t0_x", "d_x"], "condition": ["YPD", "drug"],
                        "generations": ["0.0", "6.0"], "screen_date": ["20250101", "20250101"]})
    info = step10.build_sample_info(key, ["t0_x", "d_x"], "sample_name", "condition",
                                    "generations", "screen_date")
    assert info["d_x"]["generation"] == 6
    assert f"{info['d_x']['generation']}" == "6"
    assert info["t0_x"]["generation"] == 0


def test_fix_string_decimal_generation_in_mixed_column_reaches_treatment(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01c (contract C0 scope, C4), actual main.

    A clean numeric column is read as float and works (see the control). Here the key also holds
    a plate-wash row for a sample NOT in the count matrix with generation "none", which keeps the
    whole column as strings, as in the real yL406 key. Today "6.0" then becomes 0, all four
    columns are baselines, and step 10 exits 0 with an empty manifest. After the repair the
    treatment comparison exists, and the unrelated "none" row does not cause a refusal because
    only analysis samples are validated.
    """
    rows = std_key()
    rows[2]["generations"] = rows[3]["generations"] = "6.0"
    rows.append(dict(sample_name="pw_unused", condition="plate wash", generations="none",
                     screen_date=DATE, library_ID=LIB))
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert manifest_rows(out) == standard_split_rows()


def test_fix_none_generation_for_an_analysis_sample_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01c (contract C4).

    The real yL406 key gives 162 plate-wash samples the generation "none". Today "none" becomes 0
    (10:132-133), so a plate wash of one transformation joins the pooled t0 reference and step 10
    exits 0. It must be refused, naming the sample (the explicit escape is C9 or a corrected key;
    whether "none" should instead mean "exclude" is Decision D6).
    """
    rows = std_key() + [dict(sample_name="pw_1", condition="plate wash", generations="none",
                             screen_date=DATE, library_ID=LIB)]
    counts = write_counts(work / "counts.tsv", STD_SAMPLES + ["pw_1"])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, "pw_1")
    assert_nothing_published(out)


def test_fix_missing_generation_column_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: SCR-D01c (contract C3).

    Today a key with neither `generations` nor `generation` gives every row generation 0
    (10:130), so every column is a baseline and step 10 exits 0 with an empty manifest. The
    missing column must be refused by name.
    """
    rows = [{k: v for k, v in r.items() if k != "generations"} for r in std_key()]
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, "generations")
    assert_nothing_published(out)


# ============================================================ step 10: condition and date
def test_fix_tecan_key_without_adapter_flags_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: TECAN-C01 (contract C3).

    On the corrected TECAN key the register records exit 0 with five "unknown" comparisons of
    256 reference and 256 treatment columns: `condition` and `screen_date` are absent, every
    medium and date becomes "unknown" (10:140,147) and the all-gen0 fallback fires (10:463-466).
    Both missing columns must be named in one refusal.
    """
    rows = tecan_key()
    counts = write_counts(work / "counts.tsv", [r["sample_name"] for r in rows])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, "condition", "screen_date")
    assert_nothing_published(out)


def test_fix_manual_key_without_date_column_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: MAN-C01 (contract C3).

    With only the condition adapter, the manual key has no date column, so today every date is
    "unknown" and the fallback pools both outgrowth baselines into one reference (10:463-466),
    exit 0. The missing date column must be refused, and the message must offer
    `--no-consolidate` as the explicit pooled-dates alternative.
    """
    rows = manual_key()
    counts = write_counts(work / "counts.tsv", [r["sample_name"] for r in rows])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    argv = io_args(counts, key, out) + split_args(ref) + ["--condition-column", "media"]
    outcome = run_main(step10, argv, monkeypatch, caplog, capsys)
    assert_refused(outcome, "screen_date", "--no-consolidate")
    assert_nothing_published(out)


def test_fix_treatment_sample_without_condition_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: TECAN-C01, row level (contract C4).

    Today a blank condition becomes "unknown" (10:148-149) and forms its own "unknown_6"
    comparison, exit 0. A treatment sample without a condition must be refused by name.
    """
    rows = std_key()
    rows[3]["condition"] = None
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, "d6_r2")
    assert_nothing_published(out)


def test_fix_sample_without_date_is_refused_in_date_mode(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: TECAN-C01, row level (contract C4, C7).

    Today a blank date becomes "unknown" (10:141-142); that group has no baseline, so d6_r2 is
    skipped with a warning while the rest succeeds (exit 0). In date mode every analysis sample
    needs a date: refuse, naming d6_r2.
    """
    rows = std_key()
    rows[3]["screen_date"] = None
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, "d6_r2")
    assert_nothing_published(out)


def test_fix_missing_condition_column_is_refused_by_helper(step10):
    """FAILS BEFORE: TECAN-C01 (contract C3), actual build_sample_info.

    Today a missing condition column silently labels every sample "unknown" (10:147).
    The helper itself must raise ValueError naming the missing column.
    """
    key = pd.DataFrame(std_key()).drop(columns=["condition"])
    with pytest.raises(ValueError, match="condition"):
        step10.build_sample_info(key, STD_SAMPLES, "sample_name", "condition", "generations", "screen_date")


# ============================================================ step 10: duplicates
def test_fix_conflicting_duplicate_rows_are_refused_by_helper(step10):
    """FAILS BEFORE: FL-H03 (contract C5), actual build_sample_info.

    Today the last row wins (10:151-155); reversing the legacy July key changed 960 of 3,980
    mappings. Two rows for d_x with generations 6 and 12 must raise ValueError naming d_x.
    """
    key = pd.DataFrame({"sample_name": ["t0_x", "d_x", "d_x"], "condition": ["YPD", "drug", "drug"],
                        "generations": [0, 6, 12], "screen_date": ["20250101"] * 3})
    with pytest.raises(ValueError, match="d_x"):
        step10.build_sample_info(key, ["t0_x", "d_x"], "sample_name", "condition",
                                 "generations", "screen_date")


def test_fix_conflicting_duplicate_rows_are_refused_in_either_order(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: FL-H03 (contract C5), actual main.

    The same conflict through the CLI, in both row orders. Today each order exits 0 with a
    different grouping (d6_r1 lands in drug_12 or drug_6). Both must be refused, naming d6_r1.
    """
    rows = std_key() + [dict(sample_name="d6_r1", condition="drug", generations=12,
                             screen_date=DATE, library_ID=LIB)]
    counts = write_counts(work / "counts.tsv", STD_SAMPLES)
    ref = write_reference(work / "ref.tsv")
    for label, ordered in (("forward", rows), ("reversed", rows[::-1])):
        key = write_rows(work / f"key_{label}.tsv", ordered)
        out = work / f"out_{label}"
        seed_previous(out)
        outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
        assert_refused(outcome, "d6_r1")
        assert_nothing_published(out)


# ============================================================ step 10: coverage and eligibility
def test_fix_count_column_without_metadata_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: CT02 (contract C6).

    Step 08's pseudo-replicate columns (t0_pseudorep_i, 08:257) are not in the key. Today they
    default to an "unknown"-date baseline (10:161-167), the dated treatments find no baseline,
    and step 10 exits 0 with a newline manifest and no files. It must refuse, naming the column.
    """
    counts = write_counts(work / "counts.tsv", ["d6_r1", "d6_r2", "t0_pseudorep_1", "t0_pseudorep_2"])
    key = write_rows(work / "key.tsv", std_key())
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, "t0_pseudorep_1")
    assert_nothing_published(out)


def test_fix_count_column_without_metadata_is_refused_by_helper(step10):
    """FAILS BEFORE: CT02 (contract C6), actual build_sample_info.

    The defaulting happens inside the helper (10:161-167), so the helper must raise ValueError
    naming the uncovered column.
    """
    key = pd.DataFrame(std_key())
    with pytest.raises(ValueError, match="t0_pseudorep_1"):
        step10.build_sample_info(key, ["d6_r1", "t0_pseudorep_1"], "sample_name", "condition",
                                 "generations", "screen_date")


@pytest.mark.parametrize("mode", [[], ["--no-consolidate"]], ids=["date_mode", "pooled_dates"])
def test_fix_no_eligible_baseline_is_refused(step10, work, monkeypatch, caplog, capsys, mode):
    """FAILS BEFORE: FL-C02 (contract C7, C8).

    With no generation-0 column, today both modes exit 0 and write a one-byte headerless
    manifest (10:261-263 or 341-343, then 513-525). It must refuse, name the treatment group
    (condition "drug") that has no baseline, and leave the earlier manifest untouched.
    """
    rows = [dict(sample_name=s, condition="drug", generations=g, screen_date=DATE, library_ID=LIB)
            for s, g in (("d6_r1", 6), ("d6_r2", 6), ("d12_r1", 12), ("d12_r2", 12))]
    counts = write_counts(work / "counts.tsv", [r["sample_name"] for r in rows])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref) + mode, monkeypatch, caplog, capsys)
    assert_refused(outcome, "drug")
    assert_nothing_published(out)


def test_fix_treatment_date_without_baseline_is_not_skipped_silently(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: FL-C02, partial case (contract C7).

    e6_r1 is dated 20250202, a date with no generation-0 column. Today it is dropped with a
    warning (10:261-263) while the other date succeeds, exit 0. Every requested treatment group
    must be analysed or refused; refuse, naming the date or the sample.
    """
    rows = std_key() + [dict(sample_name="e6_r1", condition="drug", generations=6,
                             screen_date=20250202, library_ID=LIB)]
    counts = write_counts(work / "counts.tsv", STD_SAMPLES + ["e6_r1"])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome, any_of=("20250202", "e6_r1"))
    assert_nothing_published(out)


def test_fix_input_without_any_treatment_sample_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: FL-C02, zero comparisons (contract C7, C8).

    Only generation-0 columns: today group_samples_by_screen_date finds no treatment group, the
    loop writes nothing, and a one-byte manifest is published with exit 0 (10:513-525). Zero
    eligible comparisons must be a non-zero refusal that says no treatment sample was found
    (checked case-insensitively for "treatment") and leaves the earlier manifest untouched.
    """
    rows = std_key()[:2]
    counts = write_counts(work / "counts.tsv", ["t0_r1", "t0_r2"])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    outcome = run_main(step10, io_args(counts, key, out) + split_args(ref), monkeypatch, caplog, capsys)
    assert_refused(outcome)
    assert "treatment" in outcome.refusal.lower(), f"refusal does not mention treatment:\n{outcome.refusal}"
    assert_nothing_published(out)


def test_fix_baseline_mixing_library_captures_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: FL-C03 (contract C2).

    The real 20220607 key labels individual-library recoveries (library yL304 ... yL315) and the
    pooled-culture aliquots (library yL304_yL315) as generation 0; here they are yL901 and the
    pooled library LIB. Today the pooled reference takes every generation-0 column (10:335),
    mixing them, exit 0 (register: 72 mixed reference columns in all 136 comparisons). A matrix
    whose analysis samples span two values of a library column (step 08's LIBRARY_COL_CANDIDATES,
    08:104) must be refused, naming both values.
    """
    rows = june_key()
    counts = write_counts(work / "counts.tsv", [r["sample_name"] for r in rows])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    out = work / "out"
    seed_previous(out)
    argv = io_args(counts, key, out) + split_args(ref) + ["--no-consolidate"]
    outcome = run_main(step10, argv, monkeypatch, caplog, capsys)
    assert_refused(outcome, "yL901", LIB)
    assert_nothing_published(out)


# ============================================================ step 10: enabling option C9
def test_enable_listed_samples_are_excluded_before_validation(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE (new flag): enabling option C9 for FL-C03 / SCR-D01c keys.

    The refusals above make the June and yL406 keys unrunnable unless some columns can be
    removed explicitly. `--exclude-samples-file` lists them; excluded columns must vanish from
    every output and the remaining pooled-culture analysis must run (2 baselines, not 3).
    """
    rows = june_key()
    counts = write_counts(work / "counts.tsv", [r["sample_name"] for r in rows])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    exclude = work / "exclude.txt"
    exclude.write_text("# individual-library recovery, not a pooled baseline\nind_t0_1\n\n")
    out = work / "out"
    argv = io_args(counts, key, out) + split_args(ref) + ["--no-consolidate", "--exclude-samples-file", exclude]
    outcome = run_main(step10, argv, monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert manifest_rows(out) == [expected_row("trafo_3", "pooled", "drug", 6, 2, 2, 2),
                                  expected_row("trafo_4", "pooled", "drug", 6, 2, 2, 2)]
    for sub in ("reference", "treatment"):
        for path in (out / sub).iterdir():
            assert "ind_t0_1" not in pd.read_csv(path, sep="\t", index_col=0).columns, path.name


def test_enable_unknown_excluded_name_is_refused(step10, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE (new flag): enabling option C9.

    An exclusion list naming a column that is not in the count matrix (a typo) must be refused,
    naming it; silently excluding nothing would re-open the mix the list was meant to prevent.
    """
    rows = june_key()
    counts = write_counts(work / "counts.tsv", [r["sample_name"] for r in rows])
    key = write_rows(work / "key.tsv", rows)
    ref = write_reference(work / "ref.tsv")
    exclude = work / "exclude.txt"
    exclude.write_text("ind_t0_9\n")
    out = work / "out"
    seed_previous(out)
    argv = io_args(counts, key, out) + split_args(ref) + ["--no-consolidate", "--exclude-samples-file", exclude]
    outcome = run_main(step10, argv, monkeypatch, caplog, capsys)
    assert_refused(outcome, "ind_t0_9")
    assert_nothing_published(out)


# ============================================================================= step 11
REF_SAMPLES = ["ref_1", "ref_2"]
TRT_SAMPLES = ["trt_1", "trt_2"]


@pytest.fixture
def step11(monkeypatch):
    """Load the actual step 11 script against a recording pydeseq2 stub (see module docstring)."""
    record = {"dds": [], "deseq2_calls": 0, "contrasts": []}

    class DefaultInference:
        def __init__(self, n_cpus=None, **kwargs):
            self.n_cpus = n_cpus

    class DeseqDataSet:
        def __init__(self, counts=None, metadata=None, design_factors=None, **kwargs):
            record["dds"].append(dict(counts=counts.copy(), metadata=metadata.copy(),
                                      design_factors=design_factors, kwargs=sorted(kwargs)))
            self.counts = counts

        def deseq2(self):
            record["deseq2_calls"] += 1

    class DeseqStats:
        def __init__(self, dds, contrast=None, **kwargs):
            record["contrasts"].append(list(contrast))
            features = list(dds.counts.columns)
            n = len(features)
            self.results_df = pd.DataFrame(
                {"baseMean": [10.0] * n, "log2FoldChange": [0.0] * n, "lfcSE": [1.0] * n,
                 "stat": [0.0] * n, "pvalue": [1.0] * n, "padj": [1.0] * n},
                index=pd.Index(features))

        def summary(self):
            return None

    package = types.ModuleType("pydeseq2")
    package.__path__ = []
    monkeypatch.setitem(sys.modules, "pydeseq2", package)
    for name, attr, obj in (("pydeseq2.dds", "DeseqDataSet", DeseqDataSet),
                            ("pydeseq2.default_inference", "DefaultInference", DefaultInference),
                            ("pydeseq2.ds", "DeseqStats", DeseqStats)):
        submodule = types.ModuleType(name)
        setattr(submodule, attr, obj)
        monkeypatch.setitem(sys.modules, name, submodule)
    # 11_deseq2_single.py:28-32 overwrites these at import; setenv first so teardown restores them.
    for var in ("SLURM_CPUS_PER_TASK", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "LOKY_MAX_CPU_COUNT"):
        monkeypatch.setenv(var, "1")
    return load_script(STEP11), record


def write_split(path, samples):
    """A step 10 reference/treatment table: bc1 index, one column per sample (10:279-283)."""
    frame = pd.DataFrame({s: [count(k, s) for k in ORDER] for s in samples},
                         index=pd.Index([BC[k] for k in ORDER], name="bc1"))
    frame.to_csv(path, sep="\t")
    return path


def step11_args(work, out, *extra):
    return ["--reference", write_split(work / "reference.tsv", REF_SAMPLES),
            "--treatment", write_split(work / "treatment.tsv", TRT_SAMPLES),
            "--output", out, "--comparison", "cmp_6", "--library", LIB, *extra]


def feature_metadata(work, values, column="yeast_sublibrary", keys=ORDER):
    return write_rows(work / "features.tsv", [{"bc1": BC[k], column: values[k]} for k in keys])


def sample_metadata(work, values, column="trafo_ID"):
    return write_rows(work / "samples.tsv", [{"sample_name": s, column: v} for s, v in values.items()])


def test_ctl_step11_without_metadata_fits_once_with_condition_design(step11, work, monkeypatch, caplog, capsys):
    """CONTROL for DESEQ-SCOPE-R1 (contract S3, S4).

    Today's CLI without metadata must keep working unchanged: exactly one model, built on the 4
    samples x 4 bc1s matrix, with design_factors="condition", the treatment-vs-reference
    contrast, and an output table with one row per bc1 plus library/comparison columns.
    """
    module, record = step11
    out = work / "results.tsv"
    outcome = run_main(module, step11_args(work, out), monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert len(record["dds"]) == 1 and record["deseq2_calls"] == 1
    call = record["dds"][0]
    assert call["design_factors"] == "condition"
    assert list(call["counts"].index) == REF_SAMPLES + TRT_SAMPLES
    assert list(call["counts"].columns) == [BC[k] for k in ORDER]
    assert list(call["metadata"]["condition"]) == ["reference"] * 2 + ["treatment"] * 2
    assert record["contrasts"] == [["condition", "treatment", "reference"]]
    result = pd.read_csv(out, sep="\t", index_col=0)
    assert list(result.index) == [BC[k] for k in ORDER]
    assert {"library", "comparison", "padj", "log2FoldChange"} <= set(result.columns)


def test_fix_step11_features_from_two_transformations_are_refused(step11, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: DESEQ-SCOPE-R1, feature axis (contract S1, S2).

    The register's synthetic case: features from trafo_A and trafo_B in one condition-only fit
    made 197/200 and 196/200 rows padj<.05 from a shift confined to one transformation. Given
    bc1 metadata that spans two transformations, step 11 must refuse before building the model,
    name both values and leave the output path untouched. Today the flags do not exist (argparse
    exit 2), which assert_refused rejects as "not a validation refusal".
    """
    module, record = step11
    out = work / "results.tsv"
    out.write_text(SENTINEL)
    features = feature_metadata(work, {"A": "trafo_A", "C": "trafo_A", "G": "trafo_B", "T": "trafo_B"})
    argv = step11_args(work, out, "--feature-metadata", features, "--transformation-column", "yeast_sublibrary")
    outcome = run_main(module, argv, monkeypatch, caplog, capsys)
    assert_refused(outcome, "trafo_A", "trafo_B")
    assert record["dds"] == []
    assert out.read_text() == SENTINEL


def test_fix_step11_columns_from_two_transformations_are_refused(step11, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE: DESEQ-SCOPE-R1, sample axis (contract S1, S2).

    Reference columns from trafo_3 and treatment columns from trafo_4 (per supplied sample
    metadata) would make the contrast a between-transformation comparison. Refuse before the
    model is built, naming both values.
    """
    module, record = step11
    out = work / "results.tsv"
    out.write_text(SENTINEL)
    samples = sample_metadata(work, {"ref_1": "trafo_3", "ref_2": "trafo_3", "trt_1": "trafo_4", "trt_2": "trafo_4"})
    argv = step11_args(work, out, "--sample-metadata", samples, "--transformation-column", "trafo_ID")
    outcome = run_main(module, argv, monkeypatch, caplog, capsys)
    assert_refused(outcome, "trafo_3", "trafo_4")
    assert record["dds"] == []
    assert out.read_text() == SENTINEL


def test_fix_step11_explicit_pooling_fits_all_features(step11, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE (new flag): DESEQ-SCOPE-R1 (contract S2).

    With `--pool-transformations` the same two-transformation input is fitted exactly as today:
    one model on all four bc1s.
    """
    module, record = step11
    out = work / "results.tsv"
    features = feature_metadata(work, {"A": "trafo_A", "C": "trafo_A", "G": "trafo_B", "T": "trafo_B"})
    argv = step11_args(work, out, "--feature-metadata", features, "--transformation-column",
                       "yeast_sublibrary", "--pool-transformations")
    outcome = run_main(module, argv, monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert len(record["dds"]) == 1
    assert list(record["dds"][0]["counts"].columns) == [BC[k] for k in ORDER]


def test_fix_step11_single_transformation_metadata_changes_nothing(step11, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE (new flags): DESEQ-SCOPE-R1 (contract S2, S4).

    Metadata showing one transformation on both axes must be accepted, and the model input and
    the written results must be identical to the no-metadata run: the repair checks scope, it does
    not change the model.
    """
    module, record = step11
    plain = work / "plain.tsv"
    scoped = work / "scoped.tsv"
    assert run_main(module, step11_args(work, plain), monkeypatch, caplog, capsys).code == 0
    features = feature_metadata(work, dict.fromkeys(ORDER, "trafo_A"), column="trafo_ID")
    samples = sample_metadata(work, dict.fromkeys(REF_SAMPLES + TRT_SAMPLES, "trafo_A"))
    argv = step11_args(work, scoped, "--feature-metadata", features, "--sample-metadata", samples,
                       "--transformation-column", "trafo_ID")
    outcome = run_main(module, argv, monkeypatch, caplog, capsys)
    assert outcome.code == 0, outcome.refusal
    assert len(record["dds"]) == 2
    pd.testing.assert_frame_equal(record["dds"][0]["counts"], record["dds"][1]["counts"])
    pd.testing.assert_frame_equal(pd.read_csv(plain, sep="\t", index_col=0),
                                  pd.read_csv(scoped, sep="\t", index_col=0))


def test_fix_step11_feature_without_transformation_is_refused(step11, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE (new flags): DESEQ-SCOPE-R1 (contract S2).

    Scope can only be verified if every fitted bc1 is covered. Metadata that lacks bc1 T (it
    could belong to any transformation) must be refused, naming T, before the model is built.
    """
    module, record = step11
    out = work / "results.tsv"
    out.write_text(SENTINEL)
    features = feature_metadata(work, dict.fromkeys(ORDER, "trafo_A"), keys=["A", "C", "G"])
    argv = step11_args(work, out, "--feature-metadata", features, "--transformation-column", "yeast_sublibrary")
    outcome = run_main(module, argv, monkeypatch, caplog, capsys)
    assert_refused(outcome, BC["T"])
    assert record["dds"] == []
    assert out.read_text() == SENTINEL


def test_fix_step11_metadata_without_declared_column_is_refused(step11, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE (new flags): contract S1, Decision D2.

    Which column is authoritative for transformation is Kevin's decision, so step 11 has no
    default: metadata without `--transformation-column` must be refused, naming the flag.
    """
    module, record = step11
    out = work / "results.tsv"
    features = feature_metadata(work, dict.fromkeys(ORDER, "trafo_A"))
    outcome = run_main(module, step11_args(work, out, "--feature-metadata", features), monkeypatch, caplog, capsys)
    assert_refused(outcome, "--transformation-column")
    assert record["dds"] == []


def test_fix_step11_design_pool_column_is_refused(step11, work, monkeypatch, caplog, capsys):
    """FAILS BEFORE (new flags): SCR-D01b applied at step 11 (contract S1).

    An oligo-name column identifies design pools; it can never serve as the transformation column,
    even when its values look distinct.
    """
    module, record = step11
    out = work / "results.tsv"
    features = feature_metadata(work, OLIGO, column="oligo_name")
    argv = step11_args(work, out, "--feature-metadata", features, "--transformation-column", "oligo_name")
    outcome = run_main(module, argv, monkeypatch, caplog, capsys)
    assert_refused(outcome, "oligo_name")
    assert record["dds"] == []
