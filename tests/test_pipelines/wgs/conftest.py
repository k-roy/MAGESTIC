"""
Shared fixtures and helpers for the `magestic.pipelines.wgs` test suite.

Scope note (deliberate deviation from the survey, recorded here so a reviewer does
not "fix" it): the wgs fixtures live in this module-local conftest rather than in
`tests/conftest.py`. pytest collects nested conftests, so they are available to
every test in this directory, and keeping them here avoids editing a file shared
with every other pipeline's suite. `testpaths = ["tests"]` only excludes test
*files* that live under `src/` -- a nested conftest under `tests/` is collected
normally.

Everything here is synthetic. No test in this directory requires a real Oak or
Sherlock path, and no committed file in this directory may contain one
(ARCHITECTURE_POLICY section 8; see `oak_path_pattern`).

Author: Kevin R. Roy
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path
from types import ModuleType

import pytest

# ---------------------------------------------------------------------------
# Package under test
# ---------------------------------------------------------------------------

WGS_PKG = "magestic.pipelines.wgs"

#: The six console scripts declared for this module by spec section 5.
#: Keyed by the pipeline step they drive; `None` = the umbrella entry point.
WGS_CONSOLE_SCRIPTS = {
    None: "magestic-wgs",
    "1": "magestic-wgs-bc1",
    "2": "magestic-wgs-donor-precall",
    "3": "magestic-wgs-sv-plasmid",
    "4": "magestic-wgs-reconcile",
    "5": "magestic-wgs-join",
}

#: `_cli` callables backing those scripts, in the same order as section 5's table.
WGS_CLI_MAINS = (
    "umbrella_main",
    "bc1_main",
    "donor_precall_main",
    "sv_plasmid_main",
    "reconcile_main",
    "join_outcomes_main",
)

#: Accepted basenames for the packaged Snakemake rule library, in preference
#: order. Spec 166 section 1 writes `Snakefile.rules` (the bdb/gdb style) while
#: the REDI exemplar the same spec tells wgs to copy uses a plain `Snakefile`.
#: The contract that matters is that a packaged rule library exists and that
#: every rule shells out to an installed console script -- not which of two
#: house filenames it uses -- so the tests accept either and say so.
RULE_LIBRARY_NAMES = ("Snakefile.rules", "Snakefile")

#: Snakemake's conventional aggregate target. It carries only `input:`, so it is
#: exempt from the per-rule shell/log/resources requirements.
TARGET_RULE_NAMES = frozenset({"all"})


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def require(module: ModuleType, name: str, spec_ref: str):
    """Fetch `name` from `module`, failing with the spec clause that requires it.

    Used for every attribute whose name the specification pins. A plain
    `from ... import name` would surface a name mismatch between two concurrently
    built components as an opaque ImportError at collection time; this reports
    *which clause* is unsatisfied instead.

    Parameters
    ----------
    module : ModuleType
        Already-imported module to read the attribute from.
    name : str
        Attribute name fixed by the specification.
    spec_ref : str
        Human-readable spec citation, e.g. ``"166 section 2.5"``.

    Returns
    -------
    object
        The attribute.
    """
    obj = getattr(module, name, None)
    assert obj is not None, (
        f"{module.__name__}.{name} is missing or None. "
        f"Required by {spec_ref}. "
        f"Public names present: {sorted(n for n in vars(module) if not n.startswith('_'))}"
    )
    return obj


def import_wgs(suffix: str = "") -> ModuleType:
    """Import a wgs submodule by dotted suffix (``""`` = the package root)."""
    return importlib.import_module(f"{WGS_PKG}.{suffix}" if suffix else WGS_PKG)


def read_text(path: Path) -> str:
    """Read a committed repo file, failing with the path when it is absent."""
    assert path.is_file(), f"required file is missing: {path}"
    return path.read_text()


def rule_library_path(wgs_src_dir: Path) -> Path:
    """Locate the packaged Snakemake rule library, accepting either house name."""
    snake_dir = wgs_src_dir / "snakemake"
    for name in RULE_LIBRARY_NAMES:
        candidate = snake_dir / name
        if candidate.is_file():
            return candidate
    raise AssertionError(
        f"no packaged Snakemake rule library under {snake_dir}; "
        f"expected one of {list(RULE_LIBRARY_NAMES)} (spec 166 section 1)"
    )


def files_under(root: Path, exclude_names=(), exclude_globs=()) -> dict:
    """Map relative-path -> absolute Path for every file under `root`.

    Parameters
    ----------
    root : Path
        Directory to walk.
    exclude_names : iterable of str
        Exact basenames to drop.
    exclude_globs : iterable of str
        `Path.match` patterns to drop.

    Returns
    -------
    dict
        ``{relative_posix_path: absolute Path}``.
    """
    out = {}
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        if p.name in exclude_names:
            continue
        if any(p.match(g) for g in exclude_globs):
            continue
        out[p.relative_to(root).as_posix()] = p
    return out


# ---------------------------------------------------------------------------
# Repo-location fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Repository root, derived from this file's location."""
    # tests/test_pipelines/wgs/conftest.py -> repo root is parents[3]
    root = Path(__file__).resolve().parents[3]
    assert (root / "pyproject.toml").is_file(), f"repo root not found from {__file__}"
    return root


@pytest.fixture(scope="session")
def wgs_src_dir(repo_root: Path) -> Path:
    """`src/magestic/pipelines/wgs/` -- the module under test, on disk."""
    return repo_root / "src" / "magestic" / "pipelines" / "wgs"


@pytest.fixture(scope="session")
def pyproject_text(repo_root: Path) -> str:
    """Raw text of pyproject.toml.

    Parsed by regex rather than `tomllib` on purpose: the repo's floor is Python
    3.10 but the local interpreter used to run this suite is 3.9, where
    `tomllib` does not exist. A regex read of a handful of pinned lines is
    strictly more portable than an optional-import dance that could skip.
    """
    return read_text(repo_root / "pyproject.toml")


@pytest.fixture(scope="session")
def entry_points(pyproject_text: str) -> dict:
    """`[project.scripts]` as ``{script_name: "module:function"}``."""
    block = re.search(
        r"^\[project\.scripts\]\s*$(.*?)(?=^\[)", pyproject_text, re.M | re.S
    )
    assert block, "pyproject.toml has no [project.scripts] table"
    return {
        m.group(1): m.group(2)
        for m in re.finditer(r'^([\w-]+)\s*=\s*"([^"]+)"', block.group(1), re.M)
    }


@pytest.fixture(scope="session")
def optional_dependencies(pyproject_text: str) -> dict:
    """`[project.optional-dependencies]` as ``{extra: [requirement, ...]}``."""
    block = re.search(
        r"^\[project\.optional-dependencies\]\s*$(.*?)(?=^\[)",
        pyproject_text,
        re.M | re.S,
    )
    assert block, "pyproject.toml has no [project.optional-dependencies] table"
    extras = {}
    for m in re.finditer(r"^([\w-]+)\s*=\s*\[(.*?)^\]", block.group(1), re.M | re.S):
        reqs = re.findall(r'"([^"]+)"', m.group(2))
        extras[m.group(1)] = reqs
    return extras


@pytest.fixture(scope="session")
def example_wgs_screen_config(repo_root: Path) -> Path:
    """The committed per-cohort screen YAML for wgs (spec section 1).

    Used instead of a hand-written config so the no-silent-success tests make no
    guess about the config schema: whatever the module committed is what gets
    exercised. Its paths are placeholders (`/path/to/...`) that do not exist,
    which is precisely the "missing input" condition those tests need.
    """
    screens = repo_root / "workflow" / "config" / "screens"
    candidates = sorted(screens.glob("*wgs*.yaml"))
    assert candidates, (
        f"no wgs screen config committed under {screens} "
        "(spec 166 section 1 requires workflow/config/screens/<screen>_wgs.yaml)"
    )
    return candidates[0]


# ---------------------------------------------------------------------------
# Policy fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def oak_path_pattern() -> re.Pattern:
    """Regex matching an internal cluster path literal.

    ARCHITECTURE_POLICY section 8 forbids real cluster paths in committed repo
    code and says CI greps the diff for them. The pattern is assembled from
    fragments so that this test file does not itself contain the literal string
    it bans -- otherwise the guard would flag its own guard.
    """
    return re.compile("/" + "oak/stanford/groups/" + "lars" + "ms")


@pytest.fixture(scope="session")
def trace_default_bam_suffix() -> str:
    """TRACE's `bam_suffix` default, which matches none of this program's BAMs.

    `bc1_consensus.run` defaults to this string and discovers BAMs with
    `glob(bam_dir + "/*" + bam_suffix)`. Our BAMs are `plate_1B_clone_1.bam`, so
    an inherited default yields an empty glob, a header-only TSV, and exit code
    0 -- fail condition F1, the single most likely silent failure in the module.
    """
    return ".sort.dedup.recal2.bam"


# ---------------------------------------------------------------------------
# Synthetic data fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def isolated_cwd(tmp_path, monkeypatch) -> Path:
    """Run a test in an empty directory and return it.

    The no-silent-success tests assert that a failing run leaves *no* output
    behind, so they need a directory whose contents are attributable.
    """
    work = tmp_path / "cwd"
    work.mkdir()
    monkeypatch.chdir(work)
    return work


@pytest.fixture
def manifest_rows() -> list:
    """Synthetic per-sample manifest rows (decision D1).

    Reference is a per-SAMPLE column, never a cohort scalar: the measured BAM
    sweep found one directory spanning two references (315 + 342). These rows
    reproduce that shape deliberately -- three samples, two references.
    """
    return [
        {
            "sample": "plate_1B_clone_1",
            "bam": "bam/plate_1B_clone_1.bam",
            "reference": "refA_T2T_full.fasta",
            "cohort": "v11_cohort",
            "plate": "1B",
            "background": "yKR888",
        },
        {
            "sample": "plate_1B_clone_2",
            "bam": "bam/plate_1B_clone_2.bam",
            "reference": "refA_T2T_full.fasta",
            "cohort": "v11_cohort",
            "plate": "1B",
            "background": "yKR888",
        },
        {
            "sample": "rearray_RED_s1",
            "bam": "bam/rearray_RED_s1.bam",
            "reference": "refB_T2T_full.fasta",
            "cohort": "rearray_post",
            "plate": "RED",
            "background": "yKR900",
        },
    ]


@pytest.fixture
def tmp_manifest_tsv(tmp_path, manifest_rows) -> Path:
    """Write `manifest_rows` as a header'd TSV, return the path."""
    cols = ["sample", "bam", "reference", "cohort", "plate", "background"]
    path = tmp_path / "manifest.tsv"
    lines = ["\t".join(cols)]
    lines += ["\t".join(str(r[c]) for c in cols) for r in manifest_rows]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def tmp_empty_bam_dir(tmp_path) -> Path:
    """A directory holding BAM-named files that TRACE's default suffix misses.

    The files are named the way this program names them (`<sample>.bam`), so a
    stage that inherited TRACE's `.sort.dedup.recal2.bam` default globs zero
    files here -- the F1 condition, reproduced without needing a real BAM.
    """
    d = tmp_path / "bam"
    d.mkdir()
    for name in ("plate_1B_clone_1.bam", "plate_1B_clone_2.bam"):
        (d / name).write_bytes(b"")
        (d / (name + ".bai")).write_bytes(b"")
    return d


@pytest.fixture
def donor_samplesheet_columns() -> list:
    """TRACE's seven required samplesheet columns (`donor_outcome._COL`)."""
    return [
        "sample",
        "chrom",
        "donor_start_coord",
        "donor_end_coord",
        "wt_seq",
        "design_donor",
        "synthetic_donor",
    ]


@pytest.fixture
def flat_locus_sidecar_keys() -> list:
    """The six keys `bc1_consensus.tracks_from_config` indexes with `cfg[...]`.

    Every one is a hard `KeyError` if absent, and `load_locus_config` is not a
    YAML parser -- it splits each line on the first colon and ignores nesting
    entirely, so a nested screen YAML is silently mis-read rather than rejected.
    """
    return [
        "barcode_locus",
        "plasmid_bc_locus",
        "barcode_left_sequence",
        "barcode_right_sequence",
        "plasmid_bc_left_sequence",
        "plasmid_bc_right_sequence",
    ]
