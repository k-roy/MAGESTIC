"""
Import and packaging contract for `magestic.pipelines.wgs`.

Covers four things the module cannot be correct without:

1. Every file the layout (spec 166 section 1) names imports cleanly.
2. `__all__` and the house sub-package docstring convention hold.
3. All six console scripts declared in section 5 resolve to real callables, the
   `wgs` extra exists, and the TRACE pin is a full 40-character sha.
4. The dependency stays ONE-DIRECTIONAL: wgs calls TRACE as a library and never
   shells out to TRACE's de-facto module CLIs (ARCHITECTURE_POLICY section 11).

These tests deliberately import at module scope with no `importorskip`. If the
module is absent or the optional `magestic[wgs]` extra is not installed, this
file errors at collection -- which is the correct loud signal. A skip here would
turn "the module does not exist" into a green suite.

Author: Kevin R. Roy
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest

from .conftest import (
    WGS_CLI_MAINS,
    WGS_CONSOLE_SCRIPTS,
    WGS_PKG,
    import_wgs,
    read_text,
    require,
    rule_library_path,
)

pytestmark = pytest.mark.unit


#: Every module the spec's layout names, as a dotted suffix of `WGS_PKG`.
LAYOUT_MODULES = (
    "",
    "_cli",
    "config",
    "core",
    "core.manifest",
    "core.locus_config",
    "core.samplesheet",
    "core.vcf_adapter",
    "core.reconcile_driver",
    "core.emit",
    "readers",
    "readers.bam_inventory",
    "readers.oligo_tables",
    "readers.wgs_calls",
    "pipeline",
    "pipeline.run_pipeline",
    "pipeline.step_1_bc1_call",
    "pipeline.step_2_donor_precall",
    "pipeline.step_3_sv_plasmid",
    "pipeline.step_4_reconcile",
    "pipeline.step_5_join_outcomes",
    "qc",
    "qc.pipeline_qc",
    "qc.provenance",
    "snakemake",
    "slurm",
    "scripts",
)

#: Sub-packages that carry the house one-line docstring and `__all__`.
SUBPACKAGES = ("", "core", "readers", "pipeline", "qc", "snakemake", "slurm", "scripts")



# ---------------------------------------------------------------------------
# 1. Layout imports
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("suffix", LAYOUT_MODULES)
def test_layout_module_imports(suffix):
    """Every module named by the layout imports."""
    mod = import_wgs(suffix)
    assert mod is not None


def test_package_all_is_the_house_list():
    """`__all__` matches the sub-package list section 2.2 fixes."""
    wgs = import_wgs()
    assert getattr(wgs, "__all__", None) == [
        "config",
        "core",
        "readers",
        "pipeline",
        "qc",
    ], f"unexpected __all__: {getattr(wgs, '__all__', None)!r}"


def test_every_all_entry_is_importable():
    """Each `__all__` entry resolves to a real submodule.

    A dangling `__all__` entry does not raise on a plain import -- only on
    `import *` -- so it needs its own check.
    """
    wgs = import_wgs()
    unresolvable = []
    for name in wgs.__all__:
        try:
            importlib.import_module(f"{WGS_PKG}.{name}")
        except ImportError:
            unresolvable.append(name)
    assert not unresolvable, f"__all__ entries that do not import: {unresolvable}"


@pytest.mark.parametrize("suffix", SUBPACKAGES)
def test_subpackage_docstring_convention(suffix):
    """Sub-packages open with `magestic.pipelines.wgs[.<sub>] - <purpose>`."""
    mod = import_wgs(suffix)
    doc = (mod.__doc__ or "").strip()
    assert doc, f"{mod.__name__} has no module docstring (house convention)"
    expected = WGS_PKG if not suffix else f"{WGS_PKG}.{suffix}"
    first = doc.splitlines()[0].strip()
    assert first.startswith(expected), (
        f"{mod.__name__} docstring opens with {first!r}; "
        f"house convention is a first line starting {expected!r}"
    )


def test_package_docstring_records_stage_2_is_gated():
    """The package docstring must say Stage 2 is GATED (section 2.2).

    Stage 2 cannot run today (no gVCF anywhere; the underscore-contig crash
    against the mandated composite reference; no per-oligo step-2 donor). An
    unannotated Stage 2 invites a zero in `synthetic_errors` being read as
    "no synthesis errors found" rather than "axis disabled".
    """
    doc = (import_wgs().__doc__ or "").upper()
    assert "GATED" in doc, (
        "package docstring does not mark Stage 2 as GATED; spec 166 section 2.2 "
        "requires it to name the four stages and Stage 2's gated status"
    )


# ---------------------------------------------------------------------------
# 2. Console scripts
# ---------------------------------------------------------------------------


def test_all_six_console_scripts_are_declared(entry_points):
    """pyproject declares exactly the six scripts section 5 names."""
    declared = {k for k in entry_points if k.startswith("magestic-wgs")}
    expected = set(WGS_CONSOLE_SCRIPTS.values())
    assert declared == expected, (
        f"console script drift: missing={sorted(expected - declared)} "
        f"unexpected={sorted(declared - expected)}"
    )


@pytest.mark.parametrize("script", sorted(set(WGS_CONSOLE_SCRIPTS.values())))
def test_console_script_target_resolves(entry_points, script):
    """Each `magestic-wgs*` target imports and is callable."""
    target = entry_points[script]
    module_path, _, func_name = target.partition(":")
    assert module_path.startswith(WGS_PKG), (
        f"{script} points outside the wgs package: {target}"
    )
    mod = importlib.import_module(module_path)
    fn = require(mod, func_name, "166 section 5 (entry points)")
    assert callable(fn), f"{target} is not callable"


def test_cli_exposes_every_main():
    """`_cli` carries all six `*_main` callables (the REDI pattern)."""
    cli = import_wgs("_cli")
    for name in WGS_CLI_MAINS:
        fn = require(cli, name, "166 section 5 (entry points)")
        assert callable(fn), f"_cli.{name} is not callable"


# ---------------------------------------------------------------------------
# 3. The `wgs` extra and the TRACE pin
# ---------------------------------------------------------------------------


def test_wgs_extra_exists(optional_dependencies):
    """`[project.optional-dependencies]` declares a `wgs` extra."""
    assert "wgs" in optional_dependencies, (
        f"no `wgs` extra declared; extras present: {sorted(optional_dependencies)}"
    )


def test_full_extra_includes_wgs(optional_dependencies):
    """`full` pulls wgs in, so a full install exercises this module."""
    full = " ".join(optional_dependencies.get("full", []))
    assert "wgs" in full, f"`full` extra does not include wgs: {full!r}"


def test_trace_is_pinned_to_one_immutable_release(optional_dependencies):
    """TRACE is pinned to ONE PyPI release (`==X.Y.Z`), not a range, branch or git URL.

    A moving requirement would let the installed algorithm layer change underneath a
    parity test that asserts byte-identical output. A PyPI release cannot be replaced,
    so `==` fixes the code as firmly as the 40-character sha this test used to require;
    a git URL is refused because PyPI will not publish a package that declares one
    (SGTC planning/268 section 5), which kept `magestic[wgs]` off PyPI.
    """
    reqs = optional_dependencies["wgs"]
    trace = [r for r in reqs if "trace" in r.lower()]
    assert trace, f"`wgs` extra does not depend on TRACE: {reqs}"
    assert len(trace) == 1, f"more than one TRACE requirement: {trace}"
    assert "git+" not in trace[0] and "@" not in trace[0], (
        f"TRACE is a direct URL requirement: {trace[0]!r}. PyPI refuses to publish it.")
    pin = re.fullmatch(r"\s*trace-crispr\s*==\s*(\d+)\.(\d+)\.(\d+)\s*", trace[0])
    assert pin, f"TRACE is not pinned to one release with `==X.Y.Z`: {trace[0]!r}"
    assert tuple(map(int, pin.groups())) >= (0, 7, 0), (
        f"TRACE {trace[0]!r} predates 0.7.0, the first release with bc1_consensus, "
        "reconcile and sv_plasmid; older releases give ModuleNotFoundError on 3 of 4 modules.")


def test_rapidfuzz_is_a_hard_dependency_of_the_extra(optional_dependencies):
    """rapidfuzz is required, not optional (spec section 4.2).

    `donor_outcome.levenshtein` is defined twice -- a rapidfuzz wrapper and a
    pure-Python fallback under `except ImportError`. Two environments could take
    different code paths, which would break dual-path parity in a way that looks
    like a data difference.
    """
    reqs = " ".join(optional_dependencies["wgs"])
    assert "rapidfuzz" in reqs, (
        "rapidfuzz is not in the `wgs` extra; donor_outcome would silently fall "
        "back to a pure-Python Levenshtein in some environments (section 4.2)"
    )


def test_extra_documents_the_samtools_binary(pyproject_text):
    """The extras table names `samtools`, which cannot be a Python dependency.

    `sv_plasmid` shells out to it three times and defaults to the bare string
    "samtools". Nothing in a wheel can express that, so a comment is the only
    place a reader learns it, and step 3 preflights it at run time.

    The whole `[project.optional-dependencies]` table is scanned rather than a
    byte window above the `wgs` block, so reordering the extras or moving the
    comment does not fail this for a formatting reason.
    """
    table = re.search(
        r"^\[project\.optional-dependencies\]\s*$(.*?)(?=^\[)",
        pyproject_text,
        re.M | re.S,
    )
    assert table, "pyproject.toml has no [project.optional-dependencies] table"
    assert "samtools" in table.group(1), (
        "the optional-dependencies table carries no comment naming the required "
        "external `samtools` binary (spec section 5)"
    )


# ---------------------------------------------------------------------------
# 4. The one-directional dependency boundary
# ---------------------------------------------------------------------------


def _wgs_source_files(wgs_src_dir: Path):
    """Every committed .py under the module, guaranteed non-empty.

    The `assert files` is load-bearing: an absence check ("no file does X")
    passes trivially over zero files, so a module directory that exists but is
    empty would make three guards below report green.
    """
    assert wgs_src_dir.is_dir(), f"wgs module source directory missing: {wgs_src_dir}"
    files = sorted(
        p for p in wgs_src_dir.rglob("*.py") if "__pycache__" not in p.parts
    )
    assert files, f"no Python source under {wgs_src_dir}; absence checks would be vacuous"
    return files


def test_wgs_consumes_trace_as_a_library(wgs_src_dir):
    """At least one wgs file imports `trace_crispr` -- the orchestration seam."""
    hits = [
        p
        for p in _wgs_source_files(wgs_src_dir)
        if re.search(r"^\s*(from|import)\s+trace_crispr", p.read_text(), re.M)
    ]
    assert hits, (
        "no file under wgs/ imports trace_crispr. This module is the "
        "orchestration layer over TRACE; if it calls nothing, it is not wired up."
    )


def test_wgs_never_shells_out_to_a_trace_cli(wgs_src_dir):
    """No wgs file invokes TRACE's per-module `main()` CLIs.

    Section 2.1's corollary: shelling out to `python -m trace_crispr.wgs.<mod>`
    re-creates the inverted dependency through the back door and bypasses the
    library contract the pin is supposed to freeze.

    This is an ABSENCE check on a forbidden token, which a text search answers
    correctly. It is not being used to assert that a capability is present.
    """
    offenders = {}
    forbidden = re.compile(r"-m\s+trace_crispr|python[\w.]*\s+-m\s+trace_crispr")
    for p in _wgs_source_files(wgs_src_dir):
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if forbidden.search(line):
                offenders.setdefault(p.name, []).append((i, line.strip()))
    assert not offenders, f"wgs shells out to a TRACE CLI: {offenders}"


def test_no_cluster_path_literals_in_module_source(wgs_src_dir, oak_path_pattern):
    """No internal cluster path literal in committed module code (policy 8).

    Committed defaults are placeholders; CI greps the diff for these.
    """
    offenders = {}
    for p in _wgs_source_files(wgs_src_dir):
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if oak_path_pattern.search(line):
                offenders.setdefault(p.name, []).append(i)
    assert not offenders, (
        f"internal cluster path literals in committed wgs source: {offenders}"
    )


def test_snakemake_and_slurm_assets_are_packaged(wgs_src_dir):
    """The rule library and the array template ship inside the package.

    Both paths resolve these relative to the installed package, so a missing
    asset is an installation-time failure rather than a run-time one.
    """
    # Raises AssertionError naming both accepted filenames if neither is present.
    library = rule_library_path(wgs_src_dir)
    assert read_text(library).strip(), f"packaged rule library is empty: {library}"
    template = wgs_src_dir / "slurm" / "template.sbatch"
    assert template.is_file(), f"missing packaged asset: {template}"
    assert read_text(template).strip(), f"slurm template is empty: {template}"
