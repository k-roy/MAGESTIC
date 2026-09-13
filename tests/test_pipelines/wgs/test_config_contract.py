"""
Config and committed-defaults contract for `magestic.pipelines.wgs`.

Three groups, each anchored on a name or a value the specification fixes:

1. `WGSConfig` / `load_config` expose the section 2.2 API, and the config does
   NOT carry a cohort-level genome. Decision D1 is the dominant structural
   constraint of this module -- reference is a per-SAMPLE attribute carried in a
   manifest, because the measured BAM sweep found one directory spanning two
   references (315 yKR900 + 342 yKR901 in a single `bam/`). The spec states the
   consequence bluntly: "If `config.genome` appears anywhere, the design is
   already wrong." This file is where that gets checked.

2. The committed screen YAML holds placeholders only, and it sets `bam_suffix`
   explicitly -- section 2.6's parameter-ownership rule. TRACE's default is
   `.sort.dedup.recal2.bam`, which matches none of this program's BAMs; an
   inherited default produces a header-only TSV and exit code 0.

3. The SLURM array template carries the measured Sherlock traps and, critically,
   does NOT carry the AVX-512 constraint that was retracted on 2026-08-06.

Author: Kevin R. Roy
"""

from __future__ import annotations

import dataclasses
import inspect
import re
from pathlib import Path

import pytest
import yaml

from .conftest import import_wgs, read_text, require

pytestmark = pytest.mark.unit


@pytest.fixture(scope="module")
def screen_yaml(example_wgs_screen_config: Path) -> dict:
    """The committed per-cohort screen YAML, parsed."""
    data = yaml.safe_load(example_wgs_screen_config.read_text())
    assert isinstance(data, dict), (
        f"{example_wgs_screen_config} does not parse to a mapping"
    )
    return data


@pytest.fixture(scope="module")
def sbatch_text(wgs_src_dir: Path) -> str:
    return read_text(wgs_src_dir / "slurm" / "template.sbatch")


# ---------------------------------------------------------------------------
# 1. WGSConfig API, and the absence of a cohort-level genome
# ---------------------------------------------------------------------------


def test_config_exposes_the_house_api():
    """`WGSConfig` and `load_config` exist with the section 2.2 shape."""
    cfg_mod = import_wgs("config")
    wgs_config = require(cfg_mod, "WGSConfig", "166 section 2.2")
    require(cfg_mod, "load_config", "166 section 2.2")
    assert dataclasses.is_dataclass(wgs_config), "WGSConfig must be a dataclass"
    for method in ("validate", "ensure_directories", "get_working_dir"):
        assert callable(getattr(wgs_config, method, None)), (
            f"WGSConfig.{method}() is missing. Section 2.2 borrows it from "
            "guide_donor_bc0/config.py so `--use-scratch` works."
        )


def test_get_working_dir_takes_a_subdir():
    """`get_working_dir(subdir)` redirects under scratch when asked."""
    wgs_config = require(import_wgs("config"), "WGSConfig", "166 section 2.2")
    params = list(inspect.signature(wgs_config.get_working_dir).parameters)
    params = [p for p in params if p != "self"]
    assert params, "WGSConfig.get_working_dir() takes no subdir argument"


#: Field names that would make reference a cohort-level scalar (decision D1).
#: `reference_root` is explicitly allowed: it is a search root for resolving
#: relative `reference` values out of the manifest, not a reference itself.
BANNED_REFERENCE_FIELDS = frozenset(
    {"genome", "genome_fasta", "reference", "reference_fasta", "genome_path"}
)


def _walk_dataclass_fields(dc, seen=None, prefix=""):
    """Yield ``(dotted_name, field)`` for a dataclass and every nested one."""
    seen = seen if seen is not None else set()
    if dc in seen:
        return
    seen.add(dc)
    for f in dataclasses.fields(dc):
        yield f"{prefix}{f.name}", f
        for candidate in (f.type, getattr(f, "default_factory", None)):
            if dataclasses.is_dataclass(candidate) and isinstance(candidate, type):
                yield from _walk_dataclass_fields(
                    candidate, seen, prefix=f"{prefix}{f.name}."
                )


def test_config_has_no_cohort_level_genome():
    """No scalar genome/reference field anywhere on the config tree (D1).

    A single reference per run cannot describe the rearray cohort, whose one
    `bam/` directory spans two references (315 yKR900 + 342 yKR901). Worse,
    chrVIII differs across the three in-use references (568,236 / 544,282 /
    560,249 -- the CUP1 contraction, with no single offset), so a cohort-level
    genome silently invites coordinate comparisons that are not portable.

    The walk is transitive on purpose: `WGSConfig` is a tree of nested
    dataclasses, and `WGSPaths` is exactly where a cohort-level `genome` would
    look most natural. Checking only the top level would leave the hole open
    precisely where the mistake is likeliest.
    """
    wgs_config = require(import_wgs("config"), "WGSConfig", "166 section 2.2")
    present = sorted(
        name
        for name, _ in _walk_dataclass_fields(wgs_config)
        if name.rsplit(".", 1)[-1] in BANNED_REFERENCE_FIELDS
    )
    assert not present, (
        f"the config tree carries cohort-level reference field(s) {present}. "
        "Decision D1: reference is a per-SAMPLE manifest column. One measured "
        "BAM directory spans two references."
    )


def _keys_to_depth(mapping, depth=2, prefix=""):
    """Yield dotted key names down to `depth` levels of nested mappings."""
    for key, value in mapping.items():
        yield f"{prefix}{key}"
        if depth > 1 and isinstance(value, dict):
            yield from _keys_to_depth(value, depth - 1, prefix=f"{prefix}{key}.")


def test_screen_yaml_declares_a_manifest_not_a_genome(screen_yaml):
    """The screen YAML points at a manifest and declares no cohort genome.

    Both halves check KEY NAMES, not a substring of the stringified document. A
    substring search would match a comment-derived value or a path fragment and
    would answer "is this word present" rather than "is this key declared".
    """
    keys = list(_keys_to_depth(screen_yaml, depth=2))
    leaves = {k.rsplit(".", 1)[-1] for k in keys}

    banned = sorted(leaves & (BANNED_REFERENCE_FIELDS - {"reference_root"}))
    assert not banned, (
        f"screen YAML declares cohort-level reference key(s) {banned} (D1). "
        "Put a `reference` COLUMN in the manifest instead."
    )
    manifest_keys = [k for k in keys if k.rsplit(".", 1)[-1] == "manifest"]
    assert manifest_keys, (
        f"screen YAML declares no `manifest` key; keys found: {sorted(keys)}. "
        "The manifest is the only source of truth for the per-sample reference "
        "(D1) and for per-plate n (the rearray cohort is sparse, not 96-complete)."
    )


# ---------------------------------------------------------------------------
# 2. Committed defaults
# ---------------------------------------------------------------------------


def test_screen_yaml_sets_bam_suffix_explicitly(screen_yaml, trace_default_bam_suffix):
    """`bam_suffix` is owned here, not inherited from TRACE (section 2.6).

    This is the F1 guard. TRACE's default matches none of this program's BAMs,
    an empty glob writes a header-only TSV, and the call returns 0.
    """
    text = yaml.safe_dump(screen_yaml)
    assert "bam_suffix" in text, (
        "screen YAML does not set `bam_suffix`. TRACE defaults it to "
        f"{trace_default_bam_suffix!r}, which matches none of our "
        "`plate_1B_clone_1.bam` files -- the empty glob then writes a "
        "header-only TSV and returns 0 (fail condition F1)."
    )
    found = re.search(r"bam_suffix\s*:\s*['\"]?([^'\"\n]+)", text)
    assert found, "could not read the bam_suffix value from the screen YAML"
    assert found.group(1).strip() != trace_default_bam_suffix, (
        "screen YAML restates TRACE's default bam_suffix verbatim; that value "
        "matches none of this program's BAM filenames"
    )


def test_screen_yaml_holds_placeholders_only(example_wgs_screen_config, oak_path_pattern):
    """No internal cluster path literal in a committed config (policy section 8)."""
    text = example_wgs_screen_config.read_text()
    hits = [
        (i, line.strip())
        for i, line in enumerate(text.splitlines(), 1)
        if oak_path_pattern.search(line)
    ]
    assert not hits, f"committed screen YAML carries cluster path literals: {hits}"
    assert "/path/to" in text, (
        f"{example_wgs_screen_config.name} contains no `/path/to/` placeholder; "
        "committed defaults must be placeholders (policy section 8)"
    )


def test_screen_yaml_declares_cluster_resources(screen_yaml):
    """A `resources:` block sizes each step for the array submission."""
    assert "resources" in screen_yaml, (
        f"screen YAML has no `resources:` block; keys: {sorted(screen_yaml)}"
    )
    resources = screen_yaml["resources"]
    assert isinstance(resources, dict) and resources, "`resources:` is empty"


# ---------------------------------------------------------------------------
# 3. The SLURM array template
# ---------------------------------------------------------------------------


def test_sbatch_has_no_avx512_constraint_or_node_exclude(sbatch_text):
    """No `--constraint`, no `--exclude` (retracted 2026-08-06).

    Measured: the `larsms` partition is exactly two Milan nodes (`sh03-07n10`,
    `sh03-08n22`, both CPU_GEN:MLN). The AVX-512 constraint matches no node in
    the partition, so `-p larsms` hard-fails and `-p larsms,owners` is silently
    evicted to the preemptible owners queue. `--exclude` on those two node names
    excludes 100% of the partition. Both of the repo's existing templates
    predate the retraction; wgs must not copy them.
    """
    offenders = [
        (i, line.strip())
        for i, line in enumerate(sbatch_text.splitlines(), 1)
        if re.search(r"--constraint|--exclude", line)
    ]
    assert not offenders, (
        f"wgs slurm template carries a retracted constraint/exclude: {offenders}"
    )


def test_sbatch_sources_bashrc_before_set_u(sbatch_text):
    """`source ~/.bashrc` comes BEFORE `set -u`.

    `/etc/profile.d/256term.sh` has an unbound COLORTERM; with `set -u` already
    active the job dies in about three seconds with no useful log.
    """
    src = sbatch_text.find("source ~/.bashrc")
    assert src >= 0, "slurm template never sources ~/.bashrc"
    setu = re.search(r"^\s*set\s+-[a-z]*u", sbatch_text, re.M)
    assert setu, "slurm template never sets `set -u`"
    assert src < setu.start(), (
        "`set -u` precedes `source ~/.bashrc`; the unbound COLORTERM in "
        "/etc/profile.d/256term.sh then kills the job in ~3 s with no useful log"
    )


def test_sbatch_skip_checks_its_own_output(sbatch_text):
    """Array tasks are idempotent: each skip-checks its output at task start.

    Policy section 6 requires this because the owners queue preempts and
    requeues routinely. This is a SHAPE check on the template (a `-f` test on an
    output path followed by an early exit) -- it is evidence the idiom is
    present, not proof that the task is idempotent.
    """
    has_file_test = re.search(r"if\s+\[\s+-f\s", sbatch_text) or re.search(
        r"if\s+\[\[\s+-f\s", sbatch_text
    )
    assert has_file_test, (
        "slurm template has no `if [ -f <output> ]` skip-check. Section 6 "
        "requires each array task to skip-check its own output at task start, "
        "because the owners queue requeues and a monolithic task loses all "
        "work in progress."
    )
    assert re.search(r"exit\s+0", sbatch_text), (
        "the skip-check never exits 0 early; a requeued task would redo the work"
    )


def test_sbatch_invokes_a_declared_console_script(sbatch_text, entry_points):
    """The template's final action is one of the module's console scripts."""
    declared = {k for k in entry_points if k.startswith("magestic-wgs")}
    called = set(re.findall(r"magestic-wgs[\w-]*", sbatch_text))
    assert called, "slurm template invokes no magestic-wgs console script"
    assert called <= declared, (
        f"slurm template calls undeclared script(s): {sorted(called - declared)}"
    )
