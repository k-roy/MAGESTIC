"""
No-silent-success guards for `magestic.pipelines.wgs`.

ARCHITECTURE_POLICY section 7, verbatim: "Orchestrators MUST exit non-zero on a
missing input or failed step, and MUST NOT log 'completed successfully' without
having produced outputs. (Anti-pattern: the bc1 stub returns `True` from every
step.)"

Every test here is TWO-part: the run must fail loudly, **and** it must leave no
output behind. `pytest.raises` on its own misses the failure mode spec 166
section 2.7 names explicitly -- the unguarded `refs[spec.chrom]` KeyError that
fires mid-loop "after some samples have already been written to the open output
file", leaving a partial file on disk that looks like a completed run.

The specific silent failures this module inherits from its algorithm layer, each
measured rather than assumed:

* `bc1_consensus.run` globs `bam_dir + "/*" + bam_suffix` with a default suffix
  of `.sort.dedup.recal2.bam`. Our BAMs are `plate_1B_clone_1.bam`, so the glob
  matches nothing, a header-only TSV is written, and the function returns 0.
* `donor_outcome.precall` drops a sample whose VCF is missing with a bare
  `continue` -- no warning, no stderr line, no count. The output TSV looks
  complete.
* `reconcile_cell` defaults `n_genome_bc=0, n_plasmid_bc=0, bc_confidence="high"`,
  so skipping Stage 1 does not give a neutral default, it gives a wrong one that
  inflates confidence.

Author: Kevin R. Roy
"""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import pytest

from .conftest import WGS_CLI_MAINS, import_wgs, require

pytestmark = pytest.mark.unit


#: `pipeline/step_N_*.py` module suffix -> the `run_step_N` it must expose.
STEP_MODULES = {
    "pipeline.step_1_bc1_call": "run_step_1",
    "pipeline.step_2_donor_precall": "run_step_2",
    "pipeline.step_3_sv_plasmid": "run_step_3",
    "pipeline.step_4_reconcile": "run_step_4",
    "pipeline.step_5_join_outcomes": "run_step_5",
}


def _invoke(fn, argv):
    """Call a CLI main, returning ``("raised", exc)`` or ``("returned", rc)``."""
    try:
        rc = fn(list(argv))
    except BaseException as exc:  # noqa: B036 - SystemExit is a valid loud failure
        if isinstance(exc, ModuleNotFoundError) and "trace_crispr" in str(exc):
            pytest.fail(
                "the magestic[wgs] extra is not installed, so this test proved "
                "nothing about failure behavior. `imports cleanly` is a "
                "precondition of this suite, not one of its criteria. "
                f"Original error: {exc}"
            )
        return "raised", exc
    return "returned", rc


def _assert_loud_failure(outcome, value, cwd: Path, what: str, allowed=()):
    """Assert the run failed AND wrote no output into `cwd`.

    Parameters
    ----------
    allowed : iterable of str
        Relative paths that existed before the run (the test's own input files).
        Everything else found afterwards is output the failed run left behind.
    """
    kind, payload = outcome, value
    if kind == "returned":
        assert payload != 0, (
            f"{what} returned {payload!r} (success) on a missing input. "
            "ARCHITECTURE_POLICY section 7 requires a non-zero exit."
        )
    allowed = set(allowed)
    leftovers = sorted(
        rel
        for rel in (
            p.relative_to(cwd).as_posix() for p in cwd.rglob("*") if p.is_file()
        )
        if rel not in allowed
    )
    assert not leftovers, (
        f"{what} failed but left output on disk: {leftovers}. "
        "A partial file from an aborted run reads downstream as a complete one."
    )


# ---------------------------------------------------------------------------
# CLI-level: a missing config must not succeed
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("main_name", WGS_CLI_MAINS)
def test_cli_fails_on_a_nonexistent_config(main_name, isolated_cwd):
    """Every entry point exits non-zero when `--config` names no file."""
    cli = import_wgs("_cli")
    fn = require(cli, main_name, "166 section 5 (entry points)")
    missing = isolated_cwd / "does_not_exist.yaml"
    assert not missing.exists()
    kind, value = _invoke(fn, ["--config", str(missing)])
    _assert_loud_failure(kind, value, isolated_cwd, f"_cli.{main_name}")


@pytest.mark.parametrize("main_name", WGS_CLI_MAINS)
def test_cli_fails_on_an_empty_config(main_name, isolated_cwd):
    """An empty YAML must raise, not fall back to built-in defaults.

    A config that silently defaults is the same defect class as the empty glob:
    the run completes, the output is wrong, and nothing says so.
    """
    cli = import_wgs("_cli")
    fn = require(cli, main_name, "166 section 5 (entry points)")
    empty = isolated_cwd / "empty.yaml"
    empty.write_text("{}\n")
    kind, value = _invoke(fn, ["--config", str(empty)])
    _assert_loud_failure(
        kind,
        value,
        isolated_cwd,
        f"_cli.{main_name} on an empty config",
        allowed={"empty.yaml"},  # the test's own input, not run output
    )


def test_umbrella_fails_on_the_committed_placeholder_config(
    example_wgs_screen_config, isolated_cwd
):
    """The shipped screen YAML points at `/path/to/...`, which does not exist.

    Running it must fail loudly. This test makes no assumption about the config
    schema -- it uses whatever the module committed -- so it cannot drift away
    from the implementation.
    """
    cli = import_wgs("_cli")
    fn = require(cli, "umbrella_main", "166 section 5 (entry points)")
    kind, value = _invoke(fn, ["--config", str(example_wgs_screen_config)])
    _assert_loud_failure(kind, value, isolated_cwd, "magestic-wgs on placeholder paths")


def test_load_config_rejects_an_empty_document(isolated_cwd):
    """`load_config` raises on a config with no required keys."""
    cfg_mod = import_wgs("config")
    load_config = require(cfg_mod, "load_config", "166 section 2.2")
    empty = isolated_cwd / "empty.yaml"
    empty.write_text("{}\n")
    with pytest.raises((KeyError, ValueError, TypeError, AttributeError)):
        load_config(empty)


def test_config_validate_reports_missing_paths(example_wgs_screen_config):
    """Either `load_config` raises, or `validate()` reports the bad paths.

    What is forbidden is the third outcome: a config built from placeholder
    paths that reports itself as valid.
    """
    cfg_mod = import_wgs("config")
    load_config = require(cfg_mod, "load_config", "166 section 2.2")
    try:
        config = load_config(example_wgs_screen_config)
    except Exception:
        return  # raising on an unusable config is a loud failure; that is fine
    validate = require(config, "validate", "166 section 2.2")
    errors = validate()
    assert isinstance(errors, list), f"validate() must return a list, got {type(errors)}"
    assert errors, (
        "validate() returned no errors for a config whose paths are `/path/to/` "
        "placeholders that do not exist on this machine"
    )


def test_run_rejects_a_step_outside_STEPS():
    """`run()` validates the step range BEFORE it touches the config.

    Asserted against `run()` rather than through a CLI flag on purpose: an
    argparse `SystemExit` proves nothing if the flag name is wrong, because
    argparse rejects an unknown *flag* with the same exception it uses for an
    unknown *value*. That would be a vacuous pass hiding inside a raises-check.

    `TypeError`/`AttributeError` are deliberately NOT accepted: those mean the
    step range was validated only after the config was dereferenced, which on a
    real run means work starts before the arguments are known to be sane.
    """
    rp = import_wgs("pipeline.run_pipeline")
    run = require(rp, "run", "166 section 2.5")
    steps = tuple(require(rp, "STEPS", "166 section 2.5"))
    assert "9" not in steps, "test assumes '9' is not a valid step"
    with pytest.raises((ValueError, KeyError)):
        run(object(), from_step="9")


def test_umbrella_cli_rejects_an_unknown_step():
    """The umbrella's step-range flag rejects a value outside `STEPS`.

    The flag is discovered from `main`'s own parser output rather than assumed,
    so this cannot pass because argparse rejected a flag name we invented.
    """
    rp = import_wgs("pipeline.run_pipeline")
    main = require(rp, "main", "166 section 2.5")
    src = inspect.getsource(rp)
    flags = set(re.findall(r'add_argument\(\s*"(--[\w-]+)"', src))
    step_flags = sorted(f for f in flags if "step" in f or f == "--from" or f == "--to")
    assert step_flags, (
        f"run_pipeline.main defines no step-range flag; flags found: {sorted(flags)}. "
        "Section 2.5 requires --steps and --from-step/--end-step."
    )
    with pytest.raises((SystemExit, ValueError, KeyError)):
        main(["--config", "unused.yaml", step_flags[0], "9"])


# ---------------------------------------------------------------------------
# The "returns True from every step" anti-pattern
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(("suffix", "fn_name"), sorted(STEP_MODULES.items()))
def test_step_function_returns_a_path_not_a_bool(suffix, fn_name):
    """`run_step_N` is annotated to return the path(s) it wrote.

    The named anti-pattern is a step that returns `True` unconditionally: the
    orchestrator then reports success without ever checking that an output
    exists. Returning the written path makes the success claim falsifiable --
    the caller can stat it.
    """
    mod = import_wgs(suffix)
    fn = require(mod, fn_name, "166 section 2.5")
    ann = inspect.signature(fn).return_annotation
    assert ann is not inspect.Signature.empty, (
        f"{fn_name} has no return annotation; section 2.5 fixes it as `-> Path`"
    )
    text = ann if isinstance(ann, str) else getattr(ann, "__name__", str(ann))
    assert "bool" not in text.lower(), (
        f"{fn_name} is annotated `-> {text}`. A boolean return is the "
        "anti-pattern ARCHITECTURE_POLICY section 7 names by name."
    )
    assert "Path" in text or "dict" in text.lower(), (
        f"{fn_name} is annotated `-> {text}`; section 2.5 requires it to return "
        "the path(s) it wrote so the caller can verify the output exists"
    )


@pytest.mark.parametrize(("suffix", "fn_name"), sorted(STEP_MODULES.items()))
def test_step_function_takes_the_config(suffix, fn_name):
    """Each step takes `config` first -- paths come from config, never hardcoded."""
    mod = import_wgs(suffix)
    fn = require(mod, fn_name, "166 section 2.5")
    params = list(inspect.signature(fn).parameters)
    assert params and params[0] == "config", (
        f"{fn_name}{tuple(params)} does not take `config` first (section 2.5)"
    )


# ---------------------------------------------------------------------------
# Structural guards against the two inherited silent failures
# ---------------------------------------------------------------------------
#
# Both checks below assert the ABSENCE of a specific forbidden call. A text
# search answers an absence question correctly; it is not being used here to
# claim that a capability is present.


def test_step_1_does_not_call_the_glob_based_trace_entry_point(wgs_src_dir):
    """Stage 1 must use the library entry, not `bc1_consensus.run()`.

    `run()` discovers BAMs by globbing a suffix and returns 0 over an empty
    match, which writes a header-only TSV. Driving the call from the manifest
    instead makes that failure mode unreachable rather than merely unlikely.
    """
    path = wgs_src_dir / "pipeline" / "step_1_bc1_call.py"
    assert path.is_file(), f"missing {path}"
    text = path.read_text()
    offenders = [
        (i, line.strip())
        for i, line in enumerate(text.splitlines(), 1)
        if re.search(r"bc1_consensus\s*\.\s*run\s*\(", line)
    ]
    assert not offenders, (
        f"step_1 calls bc1_consensus.run(): {offenders}. That entry point globs "
        "`bam_dir + '/*' + bam_suffix` and returns 0 on an empty match "
        "(fail condition F1)."
    )


def test_bam_inventory_is_manifest_driven_not_glob_driven(wgs_src_dir):
    """BAM discovery comes from the manifest (section 2.4).

    A glob cannot carry the per-sample reference that decision D1 makes
    load-bearing, and it silently yields an empty set when the suffix is wrong.
    """
    path = wgs_src_dir / "readers" / "bam_inventory.py"
    assert path.is_file(), f"missing {path}"
    text = path.read_text()
    offenders = [
        (i, line.strip())
        for i, line in enumerate(text.splitlines(), 1)
        if re.search(r"\bglob\s*\.\s*glob\s*\(|\.glob\s*\(|\.rglob\s*\(", line)
    ]
    assert not offenders, (
        f"bam_inventory discovers BAMs by globbing: {offenders}. Section 2.4 "
        "requires discovery from the manifest, which is the only source of truth "
        "for the per-sample reference (D1) and for per-plate n (the rearray "
        "cohort is sparse, not 96-complete)."
    )
