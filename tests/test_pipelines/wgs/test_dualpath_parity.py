"""
Dual-path parity for `magestic.pipelines.wgs` (ARCHITECTURE_POLICY section 4).

Section 4, verbatim: "For the same inputs + config, path (A) and path (B) MUST
produce byte-identical output files." This is the first parity test in the repo;
0 of the 4 existing pipeline modules has one, so there was no template to copy.

The file has two halves, and the split is deliberate:

* **Structural parity** (always runs, no data, no snakemake, no cluster).
  Parity is enforced by *construction* in the REDI pattern: every Snakemake rule
  shells out to the same installed console script the direct path calls, so there
  is exactly one entry point per step and one code path. These tests assert that
  construction holds. They are the half that keeps working in CI.

* **Executable parity** (integration-marked, operator-supplied fixture).
  Runs both paths for real and diffs every output byte. It skips only when its
  external prerequisites are absent, and the skip reason names them.

CARVE-OUT, pre-registered in spec 166 section 4.2 rather than discovered when the
test first fails: ARCHITECTURE_POLICY section 8 requires a PROVENANCE.json per
output, and those sidecars carry a wall-clock timestamp, so byte-identity across
two runs is impossible for them by construction. Section 4 permits exactly this
with an explicit comment naming the source of nondeterminism.

    Excluded from the diff : PROVENANCE.json, *.provenance.json
    Source of nondeterminism: wall-clock time of the run
    Everything else -- every data TSV, every .anno, every BAM -- is diffed with
    no exclusion.

Author: Kevin R. Roy
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from .conftest import (
    TARGET_RULE_NAMES,
    WGS_CONSOLE_SCRIPTS,
    files_under,
    import_wgs,
    read_text,
    require,
    rule_library_path,
)

# --- the pre-registered carve-out, as auditable constants ------------------
PARITY_EXCLUDE_NAMES = ("PROVENANCE.json",)
PARITY_EXCLUDE_GLOBS = ("*.provenance.json", "*.cache.json")
# Completion attestations contain output-specific paths and filesystem identities;
# they are verified by restart tests rather than compared across output roots.

#: Env var naming a prepared fixture directory for the executable half. It must
#: contain `config_a.yaml` and `config_b.yaml`, whose runs write into `smoke_a/`
#: and `smoke_b/` beneath it. Kept out of the repo so no cluster path is
#: committed (ARCHITECTURE_POLICY section 8).
PARITY_FIXTURE_ENV = "MAGESTIC_WGS_PARITY_FIXTURE"

#: The five per-step console scripts. The umbrella is keyed `None`.
STEP_SCRIPTS = {k: v for k, v in WGS_CONSOLE_SCRIPTS.items() if k is not None}


# ===========================================================================
# Snakefile parsing
# ===========================================================================

_RULE_START = re.compile(r"^rule\s+([A-Za-z_]\w*)\s*:")
_STRING_LITERAL = re.compile(
    r'"""(.*?)"""' r"|'''(.*?)'''" r'|"([^"\n]*)"' r"|'([^'\n]*)'", re.S
)


def _rule_blocks(text: str) -> dict:
    """Split a Snakefile into ``{rule_name: block_text}``."""
    blocks, name, buf = {}, None, []
    for line in text.splitlines():
        m = _RULE_START.match(line)
        if m:
            if name is not None:
                blocks[name] = "\n".join(buf)
            name, buf = m.group(1), [line]
        elif name is not None:
            buf.append(line)
    if name is not None:
        blocks[name] = "\n".join(buf)
    return blocks


def _shell_text(block: str) -> str:
    """Concatenate every string literal in a rule's `shell:` directive."""
    idx = block.find("shell:")
    if idx < 0:
        return ""
    parts = []
    for m in _STRING_LITERAL.finditer(block[idx:]):
        parts.append(next(g for g in m.groups() if g is not None))
    return " ".join(parts)


@pytest.fixture(scope="module")
def rules_text(wgs_src_dir: Path) -> str:
    """The packaged Snakemake rule library (either accepted house filename)."""
    return read_text(rule_library_path(wgs_src_dir))


@pytest.fixture(scope="module")
def rules(rules_text: str) -> dict:
    """Executable rules only -- `rule all` is an aggregate target, not a step.

    `rule all` carries `input:` and nothing else by design, so holding it to the
    shell / console-script / log / resources contract would be a false failure.
    """
    blocks = _rule_blocks(rules_text)
    assert blocks, "the packaged rule library declares no rules"
    steps = {n: b for n, b in blocks.items() if n not in TARGET_RULE_NAMES}
    assert steps, (
        f"the rule library declares only target rules {sorted(blocks)}; "
        "no step rule means path (B) runs nothing"
    )
    return steps


@pytest.fixture(scope="module")
def workflow_snakefile(repo_root: Path) -> str:
    """The repo-level entry point that includes the rule library."""
    return read_text(repo_root / "workflow" / "Snakefile.wgs")


@pytest.fixture(scope="module")
def workflow_text(workflow_snakefile: str, rules_text: str) -> str:
    """Entry point plus the library it includes, as one document.

    `rule all` and the `config["config"]` read may legitimately live in either
    file -- the implementation puts `rule all` in the included library with
    `default_target: True`, which is order-independent and arguably better than
    declaring it in the entry point. What must hold is that they exist
    *somewhere* in the workflow that Snakemake actually evaluates.
    """
    return workflow_snakefile + "\n" + rules_text


# ===========================================================================
# Structural parity -- always runs
# ===========================================================================


@pytest.mark.unit
def test_every_rule_has_a_shell_directive(rules):
    """No rule may carry a `run:` block.

    A `run:` block executes inline Python, which is a second code path the direct
    orchestrator does not share. That is precisely the duplication section 3
    forbids ("both orchestrators are thin wiring over identical core functions").
    """
    missing = [name for name, block in rules.items() if "shell:" not in block]
    assert not missing, f"rules without a shell: directive (inline Python path): {missing}"


@pytest.mark.unit
def test_every_rule_invokes_a_declared_console_script(rules, entry_points):
    """Each rule's shell command calls one of the six declared scripts.

    This is what makes parity structural rather than a convention someone has to
    maintain: the rule and the direct path enter the same installed entry point.
    """
    declared = {k for k in entry_points if k.startswith("magestic-wgs")}
    assert declared, "no magestic-wgs console scripts declared in pyproject.toml"
    problems = {}
    for name, block in rules.items():
        shell = _shell_text(block)
        called = set(re.findall(r"magestic-wgs[\w-]*", shell))
        if not called:
            problems[name] = f"invokes no magestic-wgs console script: {shell!r}"
        elif not called <= declared:
            problems[name] = f"invokes undeclared script(s): {sorted(called - declared)}"
    assert not problems, f"rules breaking the single-entry-point rule: {problems}"


@pytest.mark.unit
def test_no_rule_shells_out_to_a_trace_cli(rules):
    """A rule must not call TRACE's per-module `main()` CLIs.

    Section 2.1's corollary: that re-creates the inverted dependency through the
    back door, and it bypasses the console script both paths are supposed to
    share -- so it breaks parity and the architecture at the same time.
    """
    offenders = {
        name: _shell_text(block)
        for name, block in rules.items()
        if re.search(r"-m\s+trace_crispr", _shell_text(block))
    }
    assert not offenders, f"rules shelling out to a TRACE CLI: {offenders}"


@pytest.mark.unit
def test_run_pipeline_declares_the_five_stages():
    """`STEPS` is the five-stage tuple the spec fixes (section 2.5)."""
    rp = import_wgs("pipeline.run_pipeline")
    steps = require(rp, "STEPS", "166 section 2.5")
    assert tuple(steps) == ("1", "2", "3", "4", "5"), f"unexpected STEPS: {steps!r}"


@pytest.mark.unit
def test_every_step_has_a_rule(rules, entry_points):
    """Every step in `STEPS` is reachable from Snakemake, Stage 5 included.

    This is the gap the exemplar has and wgs must not copy: REDI's `rule all`
    stops at step_e while steps f/g/h have console scripts but no rules, so the
    module that defines the standard violates its own sections 2 and 4. A step
    with no rule means path (B) silently produces fewer outputs than path (A),
    and the parity diff below would then compare a subset and pass.
    """
    rp = import_wgs("pipeline.run_pipeline")
    steps = tuple(require(rp, "STEPS", "166 section 2.5"))
    invoked = set()
    for block in rules.values():
        invoked |= set(re.findall(r"magestic-wgs[\w-]*", _shell_text(block)))
    uncovered = [s for s in steps if STEP_SCRIPTS.get(s) not in invoked]
    assert not uncovered, (
        f"steps with no Snakemake rule: {uncovered}. "
        f"Expected scripts {[STEP_SCRIPTS.get(s) for s in uncovered]} to appear in a "
        f"rule's shell directive; rules invoke {sorted(invoked)}"
    )


@pytest.mark.unit
def test_both_paths_read_the_same_config(rules, workflow_text):
    """One config YAML drives both paths.

    The Snakefile may parse it to compute file names and the fan-out list, but
    every rule must hand the same `--config` file to the console script, so no
    analysis decision is made twice. Note `--config config=<yaml>` rather than
    `--configfile`: the latter would make Snakemake a second interpreter of the
    same YAML, with its own defaults and coercions, and would break parity in a
    way that reads as a data difference.
    """
    assert 'config["config"]' in workflow_text or "config['config']" in workflow_text, (
        "the workflow never reads the screen YAML from `config[\"config\"]`; "
        "the two paths could end up reading different configs"
    )
    without_config = {
        name: _shell_text(block)
        for name, block in rules.items()
        if "--config" not in _shell_text(block)
    }
    assert not without_config, (
        f"rules invoking a console script without --config: {without_config}"
    )


@pytest.mark.unit
def test_workflow_entry_point_includes_the_rule_library(
    workflow_snakefile, workflow_text, wgs_src_dir
):
    """`workflow/Snakefile.wgs` includes the packaged rules, resolved relatively.

    Resolving against `workflow.basedir` rather than by importing magestic means
    Snakemake's own environment does not need the package on its import path.
    """
    assert "include:" in workflow_snakefile, "workflow/Snakefile.wgs has no include:"
    assert "workflow.basedir" in workflow_snakefile, (
        "workflow/Snakefile.wgs does not resolve the rule library relative to "
        "`workflow.basedir`; a path that depends on the caller's cwd or on "
        "magestic being importable inside the snakemake env is not reproducible"
    )
    library = rule_library_path(wgs_src_dir).name
    assert library in workflow_snakefile, (
        f"workflow/Snakefile.wgs does not include wgs/snakemake/{library}"
    )
    assert re.search(r"^rule\s+all\s*:", workflow_text, re.M), (
        "no `rule all` in the entry point or the library it includes, so "
        "`snakemake` has no default target and path (B) would produce nothing"
    )


@pytest.mark.unit
def test_every_flag_a_rule_passes_is_defined_by_the_cli(rules, wgs_src_dir):
    """Path (B) may only pass flags path (A)'s argparse actually accepts.

    This is the sharpest failure a structural parity test can catch, because it
    is not a difference in output -- it is path (B) not running at all.
    Snakemake reports the rule as failed, and argparse's message
    ("unrecognized arguments: --combine") appears only in the rule's log file,
    which is exactly where a fan-out over hundreds of samples hides it.

    Flags are read from `_cli.py`'s `add_argument` calls rather than by
    importing and introspecting the parsers, because every `*_main` lazy-imports
    its step module before parsing, so `--help` would require the full
    `magestic[wgs]` environment just to enumerate flag names.
    """
    cli_src = wgs_src_dir / "_cli.py"
    assert cli_src.is_file(), f"missing {cli_src}"
    defined = set(re.findall(r'add_argument\(\s*"(--[\w-]+)"', cli_src.read_text()))
    assert defined, f"no add_argument() flags found in {cli_src}; parse failed"

    unknown = {}
    for name, block in rules.items():
        shell = _shell_text(block)
        if "magestic-wgs" not in shell:
            continue
        used = set(re.findall(r"(?<![\w-])(--[\w-]+)", shell))
        missing = sorted(used - defined)
        if missing:
            unknown[name] = missing
    assert not unknown, (
        f"Snakemake rules pass flags no wgs console script defines: {unknown}. "
        f"Flags defined in _cli.py: {sorted(defined)}. "
        "Path (B) would abort with `unrecognized arguments` while path (A) runs, "
        "so there would be no output to compare at all."
    )


@pytest.mark.unit
def test_rules_log_and_declare_resources(rules):
    """House rule shape: each rule tees to a log and declares resources.

    Without a log, a preempted or failed array task leaves nothing to grep, and
    the completion signal degrades to an exit code -- which this program has
    already been burned by.
    """
    missing_log = [n for n, b in rules.items() if "log:" not in b]
    missing_res = [n for n, b in rules.items() if "resources:" not in b]
    assert not missing_log, f"rules with no log: directive: {missing_log}"
    assert not missing_res, f"rules with no resources: directive: {missing_res}"


# ===========================================================================
# Executable parity -- integration
# ===========================================================================


@pytest.fixture
def parity_fixture_dir() -> Path:
    """Operator-prepared fixture directory, or skip naming what is missing."""
    raw = os.environ.get(PARITY_FIXTURE_ENV)
    if not raw:
        pytest.skip(
            f"{PARITY_FIXTURE_ENV} is unset. Set it to a directory containing "
            "config_a.yaml and config_b.yaml whose runs write to smoke_a/ and "
            "smoke_b/ beneath it (spec 166 section 7.3)."
        )
    d = Path(raw)
    assert d.is_dir(), f"{PARITY_FIXTURE_ENV} does not point at a directory: {d}"
    for name in ("config_a.yaml", "config_b.yaml"):
        assert (d / name).is_file(), f"parity fixture is missing {name}: {d}"
    return d


def _run(cmd, cwd: Path) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    assert proc.returncode == 0, (
        f"command failed ({proc.returncode}): {' '.join(str(c) for c in cmd)}\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
    )
    return proc


@pytest.mark.integration
@pytest.mark.slow
def test_cli_and_snakemake_paths_are_byte_identical(repo_root, parity_fixture_dir):
    """Path (A) and path (B) produce byte-identical outputs (section 4).

    Vacuity guards are first-class here. "Diff every output file" passes
    trivially over zero files, which is exactly the F1 shape -- header-only
    output with exit code 0 from an empty BAM glob. So the test asserts that
    files were produced, that both sides produced the same relative paths, and
    that at least one compared TSV carries a data row, before it asserts
    equality of bytes.
    """
    if shutil.which("magestic-wgs") is None:
        pytest.skip("magestic-wgs console script not on PATH; install magestic[wgs]")
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake not on PATH; install magestic[snakemake]")

    fixture = parity_fixture_dir
    out_a, out_b = fixture / "smoke_a", fixture / "smoke_b"

    _run(
        ["magestic-wgs", "--config", str(fixture / "config_a.yaml"), "--steps", "1"],
        cwd=fixture,
    )
    _run(
        [
            "snakemake",
            "--snakefile",
            str(repo_root / "workflow" / "Snakefile.wgs"),
            "--config",
            f"config={fixture / 'config_b.yaml'}",
            "steps=1",
            "-j",
            "1",
        ],
        cwd=fixture,
    )

    assert out_a.is_dir(), f"path (A) produced no output directory: {out_a}"
    assert out_b.is_dir(), f"path (B) produced no output directory: {out_b}"

    a = files_under(out_a, PARITY_EXCLUDE_NAMES, PARITY_EXCLUDE_GLOBS)
    b = files_under(out_b, PARITY_EXCLUDE_NAMES, PARITY_EXCLUDE_GLOBS)

    # --- vacuity guards, before any equality assertion ---------------------
    assert a, f"path (A) wrote no comparable files under {out_a} (F1 shape)"
    assert b, f"path (B) wrote no comparable files under {out_b} (F1 shape)"
    assert set(a) == set(b), (
        f"the two paths wrote different file sets: "
        f"only in A={sorted(set(a) - set(b))} only in B={sorted(set(b) - set(a))}"
    )
    tsvs = [rel for rel in a if rel.endswith(".tsv")]
    assert tsvs, f"neither path wrote a .tsv; files compared: {sorted(a)}"
    with_rows = [
        rel for rel in tsvs if len(a[rel].read_text().rstrip("\n").splitlines()) > 1
    ]
    assert with_rows, (
        "every compared TSV is header-only. This is fail condition F1 -- most "
        "likely TRACE's bam_suffix default was inherited and the BAM glob "
        "matched nothing. Count rows; the exit code is worthless here."
    )

    # --- the parity assertion ---------------------------------------------
    differing = [rel for rel in sorted(a) if a[rel].read_bytes() != b[rel].read_bytes()]
    assert not differing, (
        f"dual-path parity violated (ARCHITECTURE_POLICY section 4). "
        f"Files differing outside the pre-registered carve-out: {differing}"
    )


@pytest.mark.integration
@pytest.mark.slow
def test_provenance_sidecar_is_written_and_carries_a_sha(parity_fixture_dir):
    """The carve-out must exclude something that actually exists.

    An excluded-but-absent PROVENANCE.json would make the carve-out look benign
    while section 8 goes unsatisfied, and the empty-sha trap (Sherlock's git
    predates `git -C`, exiting 129 with empty stdout) would go unnoticed --
    acceptance criteria P6 and F7.
    """
    import json

    out_a = parity_fixture_dir / "smoke_a"
    if not out_a.is_dir():
        pytest.skip("no smoke_a/ present; run the parity test first")
    sidecars = list(out_a.rglob("PROVENANCE.json")) + list(
        out_a.rglob("*.provenance.json")
    )
    assert sidecars, (
        f"no provenance sidecar under {out_a}; ARCHITECTURE_POLICY section 8 "
        "requires one per output, and the parity carve-out exists only for them"
    )
    empty = []
    for s in sidecars:
        payload = json.loads(s.read_text())
        sha = str(payload.get("code", {}).get("commit") or "").strip()
        if not re.fullmatch(r"[0-9a-f]{40}", sha):
            empty.append(s.name)
        assert re.fullmatch(r"[0-9a-f]{64}",
                           payload.get("code", {}).get("magestic_wgs_source_sha256", ""))
        assert re.fullmatch(r"[0-9a-f]{64}",
                           payload.get("trace", {}).get("trace_wgs_source_sha256", ""))
    assert not empty, (
        f"provenance sidecars record an empty git sha: {empty}. "
        "Sherlock's git predates `git -C` (exit 129, empty stdout), so a "
        "subprocess helper ships a null sha while the run still exits 0."
    )
