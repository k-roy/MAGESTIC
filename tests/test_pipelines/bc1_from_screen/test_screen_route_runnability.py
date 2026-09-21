"""Failing-before tests for the screen-route runnability repairs (planning NN 218).

Design and exact code changes: REPAIR_DESIGN.md in this folder (section numbers below refer to it).

Every test asserts the REPAIRED contract. Against HEAD 7393ae6 each "failing-before" case fails, and it
passes once the matching change in REPAIR_DESIGN.md section 3 lands. Positive controls pass on both.
Nothing asserts a current exception, so no test can pass merely because the code still crashes.

What is exercised is the actual production code: numbered scripts run as subprocesses (or, for step 03,
their real ``main()`` in-process with plotting stubbed), the core functions they call, and the command lines
rendered from the shipped ``snakemake/Snakefile.rules``. All inputs are tiny synthetic tables under
``tmp_path``. No Oak path, network, cluster or canonical-repo code is touched.

Run from the root session (the matplotlib / scikit-learn / statsmodels / pydeseq2 cases need the audit's
dependency tree on PYTHONPATH; without it those cases SKIP with a reason, they do not pass):

    CHECKOUT=planning/218_magestic_production_corrections/code/magestic
    DEPS=planning/216_magestic_comprehensive_audit/dependencies
    PYTHONPATH=$CHECKOUT/src:$DEPS PYTHONDONTWRITEBYTECODE=1 \
        /Users/kevinroy/.venvs/magestic/bin/python -m pytest -p no:cacheprovider -q \
        planning/218_magestic_production_corrections/reviews/screen_route_runnability/test_screen_route_runnability.py

The venv's editable install points at the CANONICAL repo (~/repos/magestic); the session fixture
``checkout_provenance`` hard-fails every test if ``magestic`` resolves anywhere but this checkout.
Override discovery with MAGESTIC_CHECKOUT=<checkout root> if the file is moved.
"""
from __future__ import annotations

import ast
import fnmatch
import importlib.util
import inspect
import itertools
import json
import os
import re
import shlex
import string
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("MPLBACKEND", "Agg")
if "MPLCONFIGDIR" not in os.environ:  # keep matplotlib's font cache out of the home directory
    os.environ["MPLCONFIGDIR"] = tempfile.mkdtemp(prefix="screen_route_mpl_")
sys.dont_write_bytecode = True  # never drop __pycache__ into the checkout under test

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------------------------------
# Checkout under test and import provenance
# ---------------------------------------------------------------------------------------------

def _find_checkout() -> Path:
    override = os.environ.get("MAGESTIC_CHECKOUT")
    if override:
        return Path(override).resolve()
    marker = Path("src/magestic/pipelines/bc1_from_screen/snakemake/Snakefile.rules")
    here = Path(__file__).resolve()
    for base in here.parents:
        for candidate in (base, base / "code" / "magestic"):
            if (candidate / marker).is_file():
                return candidate
    raise RuntimeError("Cannot locate the MAGESTIC checkout; set MAGESTIC_CHECKOUT=<checkout root>.")


CHECKOUT = _find_checkout()
CHECKOUT_SRC = CHECKOUT / "src"
if "magestic" not in sys.modules and str(CHECKOUT_SRC) not in sys.path:
    sys.path.insert(0, str(CHECKOUT_SRC))

SCREEN = CHECKOUT_SRC / "magestic" / "pipelines" / "bc1_from_screen"
SNAKEMAKE_DIR = SCREEN / "snakemake"
RULES_FILE = SNAKEMAKE_DIR / "Snakefile.rules"
PYTHON = sys.executable

# Third-party modules that live only in the audit dependency tree. A case that needs one SKIPS (never
# passes) when it is missing.
OPTIONAL_DEPENDENCIES = {"matplotlib", "sklearn", "statsmodels", "pydeseq2", "seaborn", "patsy", "anndata"}


@pytest.fixture(scope="session", autouse=True)
def checkout_provenance():
    """Hard-fail (not skip) unless ``magestic`` is imported from the checkout under test.

    The venv's ``_editable_impl_magestic.pth`` points at ~/repos/magestic/src. Without this guard an
    in-process import could silently bind the canonical core to a checkout script, a mixed build that
    invalidates every result in this file.
    """
    import magestic

    loaded = Path(magestic.__file__).resolve()
    if CHECKOUT_SRC.resolve() not in loaded.parents:
        pytest.fail(
            f"magestic imported from {loaded}, not from the checkout under test {CHECKOUT_SRC}. "
            f"Put {CHECKOUT_SRC} first on PYTHONPATH.",
            pytrace=False,
        )
    return loaded


def _child_env(tmp_path: Path, **extra) -> dict:
    env = dict(os.environ)
    inherited = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join(p for p in (str(CHECKOUT_SRC), inherited) if p)
    mpl = tmp_path / "mplconfig"
    mpl.mkdir(exist_ok=True)
    env.update(
        MPLBACKEND="Agg", MPLCONFIGDIR=str(mpl), PYTHONDONTWRITEBYTECODE="1",
        SLURM_CPUS_PER_TASK="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", LOKY_MAX_CPU_COUNT="1",
        SCREEN_ROUTE_CHECKOUT_SRC=str(CHECKOUT_SRC),
    )
    env.update({key: str(value) for key, value in extra.items()})
    return env


def _run(argv, tmp_path: Path, timeout: int = 300, **env_extra) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(a) for a in argv], capture_output=True, text=True, timeout=timeout,
        cwd=tmp_path, env=_child_env(tmp_path, **env_extra),
    )


def _write_tsv(path, frame: pd.DataFrame, **kwargs) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=kwargs.pop("index", False), **kwargs)
    return path


def _explain(result: subprocess.CompletedProcess, argv=None) -> str:
    head = f"argv: {' '.join(map(str, argv))}\n" if argv else ""
    return f"{head}rc={result.returncode}\n--- stdout ---\n{result.stdout[-3000:]}\n--- stderr ---\n{result.stderr[-3000:]}"


def _skip_if_optional_dependency_missing(result: subprocess.CompletedProcess, what: str) -> None:
    found = re.search(r"MISSING_DEPENDENCY (\S+)", result.stderr) or re.search(
        r"ModuleNotFoundError: No module named '([^']+)'", result.stderr)
    if found and found.group(1).split(".")[0] in OPTIONAL_DEPENDENCIES:
        pytest.skip(f"{what}: optional dependency {found.group(1)!r} not importable; "
                    "add planning/216_magestic_comprehensive_audit/dependencies to PYTHONPATH")


# ---------------------------------------------------------------------------------------------
# Snakemake rule rendering (no Snakemake install needed)
# ---------------------------------------------------------------------------------------------
# Follows tests/test_pipelines/bc1_donor_bc0/test_workflow_guide_link.py (take the shell template from the
# shipped rules file, format it, keep the part before the redirect, shlex-split), adapted to this file's
# `shell:` DIRECTIVE form. input/output/params/log are evaluated from the rule text itself, so a corrected
# rule is rendered from its own declarations. Two Snakemake semantics are reproduced deliberately:
#   * wildcards are applied to input/output/log paths and to `str` params, but a pathlib.Path param is
#     passed through literally (the file's own dry-run note, Snakefile.rules:197-199);
#   * callables (input functions, param lambdas) are called with the job's wildcards and their result is
#     used as returned.

DEFAULT_WILDCARDS = dict(screen="scr", library="yL1", comparison="cmpA", sample="s1", condition="Caffeine")
DEFAULT_COMPARISONS = ("cmpA", "cmpB")
_RULE_HEAD = re.compile(r"^rule\s+(\w+)\s*:\s*$", re.M)
_SECTION_HEAD = re.compile(r"^ {4}(\w+):[ \t]*(.*)$")
_SHELL_OPERATORS = {"|", "||", "&&", ";", ">", ">>", "2>", "&>", "2>&1"}


def _rule_blocks() -> dict:
    text = RULES_FILE.read_text()
    heads = list(_RULE_HEAD.finditer(text))
    return {
        head.group(1): text[head.end(): heads[i + 1].start() if i + 1 < len(heads) else len(text)]
        for i, head in enumerate(heads)
    }


def _rule_sections(block: str) -> dict:
    sections, current = {}, None
    for line in block.splitlines():
        if not line.strip():
            if current is not None:
                sections[current].append("")
            continue
        indent = len(line) - len(line.lstrip(" "))
        if indent == 0:
            break  # a column-0 banner comment ends the rule body
        head = _SECTION_HEAD.match(line) if indent == 4 else None
        if head:
            current = head.group(1)
            sections[current] = [head.group(2)] if head.group(2).strip() else []
        elif current is not None:
            sections[current].append(line)
    return {name: "\n".join(lines) for name, lines in sections.items()}


def _expand(pattern, **wildcards):
    names = list(wildcards)
    choices = [v if isinstance(v, (list, tuple)) else [v] for v in wildcards.values()]
    return [str(pattern).format(**dict(zip(names, combo))) for combo in itertools.product(*choices)]


def _passthrough(value, *args, **kwargs):
    return value


def _flat(values) -> list:
    flat = []
    for value in values:
        flat.extend(value if isinstance(value, (list, tuple)) else [value])
    return [Path(str(v)) for v in flat]


class _ShellFormatter(string.Formatter):
    """Format a shell template as Snakemake does: lists join with spaces; the `:q` spec shell-quotes."""

    def format_field(self, value, spec):
        if isinstance(value, (list, tuple)):
            value = " ".join(str(v) for v in value)
        if spec == "q":
            return shlex.quote(str(value))
        return super().format_field(value, spec)


class RenderedRule:
    """One job of one rule of the shipped Snakefile.rules, evaluated from the rule text."""

    def __init__(self, name: str, root: Path, wildcards: dict | None = None, comparisons=DEFAULT_COMPARISONS,
                 samples=("s1", "s2")):
        blocks = _rule_blocks()
        if name not in blocks:
            raise KeyError(f"rule {name!r} not found in {RULES_FILE}")
        self.name = name
        self.wildcards = SimpleNamespace(**(wildcards or DEFAULT_WILDCARDS))
        self.comparisons = list(comparisons)
        project = Path(root) / "project"
        self.env = dict(
            Path=Path, config={}, COMMON_SCRIPTS_DIR=SNAKEMAKE_DIR, PYTHON=PYTHON,
            SAMPLES=list(samples), LIBRARIES=[self.wildcards.library], expand=_expand,
            directory=_passthrough, temp=_passthrough, protected=_passthrough,
            ancient=_passthrough, touch=_passthrough,
            **{d: project / d.lower() for d in ("BC1_COUNTS_DIR", "DESEQ2_DIR", "HIT_CALLING_DIR",
                                                  "KEYFILE_DIR", "LOGS_DIR", "COLLAPSED_FASTQ_DIR")},
        )
        sections = _rule_sections(blocks[name])
        self.input = self._keyword_section(sections.get("input", ""), "input")
        self.output = self._keyword_section(sections.get("output", ""), "output")
        self.params = self._keyword_section(sections.get("params", ""), "params")
        log = sections.get("log", "")
        self.log = self._evaluate(ast.parse(log.strip(), mode="eval").body, "log") if log.strip() else ""
        shell = sections["shell"]
        self.template = ast.literal_eval(shell[shell.index('"""'):].strip())

    def _keyword_section(self, source: str, kind: str) -> dict:
        if not source.strip():
            return {}
        call = ast.parse("_section(\n" + source + "\n)", mode="eval").body
        values = {f"_{i}": self._evaluate(node, kind) for i, node in enumerate(call.args)}
        values.update({kw.arg: self._evaluate(kw.value, kind) for kw in call.keywords})
        return values

    def _evaluate(self, node, kind: str):
        for name in {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}:
            if name.startswith("get_") and name not in self.env:  # project-Snakefile helper functions
                self.env[name] = lambda *args, **kwargs: list(self.comparisons)
        expression = ast.Expression(body=node)
        ast.fix_missing_locations(expression)
        value = eval(compile(expression, str(RULES_FILE), "eval"), self.env)  # noqa: S307 - rule text only
        if callable(value):
            return self._call(value)
        return self._concretize(value, kind)

    def _call(self, function):
        extra = {"input": SimpleNamespace(**getattr(self, "input", {})),
                 "output": SimpleNamespace(**getattr(self, "output", {})),
                 "threads": 1, "resources": SimpleNamespace(mem_mb=1000)}
        names = list(inspect.signature(function).parameters)
        return function(self.wildcards, *[extra.get(n) for n in names[1:]])

    def _concretize(self, value, kind: str):
        def one(item):
            if kind != "params" or isinstance(item, str):
                return str(item).format(**vars(self.wildcards))
            return item  # Path / number / bool params are not wildcard-expanded by Snakemake
        if isinstance(value, (list, tuple)):
            return [one(item) for item in value]
        return one(value)

    def inputs(self) -> list:
        return _flat(self.input.values())

    def outputs(self) -> list:
        return _flat(self.output.values())

    def argv(self) -> list:
        rendered = _ShellFormatter().format(
            self.template, PYTHON=PYTHON, input=SimpleNamespace(**self.input),
            output=SimpleNamespace(**self.output), params=SimpleNamespace(**self.params),
            wildcards=self.wildcards, log=self.log, threads=1, resources=SimpleNamespace(mem_mb=1000),
        )
        tokens = shlex.split(rendered.replace("\\\n", " "))
        for i, token in enumerate(tokens):
            if token in _SHELL_OPERATORS:
                return tokens[:i]
        return tokens


def _run_rule(rule: RenderedRule, tmp_path: Path, timeout: int = 300, **env_extra):
    """Run one rendered job as Snakemake would: output parent directories exist, then the command runs."""
    for path in rule.outputs():
        path.parent.mkdir(parents=True, exist_ok=True)
    argv = rule.argv()
    return _run(argv, tmp_path, timeout=timeout, **env_extra), argv


def _options(argv) -> dict:
    options, i = {}, 2
    while i < len(argv):
        if argv[i].startswith("--"):
            j, values = i + 1, []
            while j < len(argv) and not argv[j].startswith("--"):
                values.append(argv[j])
                j += 1
            options[argv[i]] = values[0] if len(values) == 1 else values
            i = j
        else:
            i += 1
    return options


# Runs a numbered script as `python script.py ...` would, but stops right after argparse succeeds.
_PARSE_ONLY_DRIVER = r"""
import argparse, json, os, runpy, sys
from pathlib import Path
class _Parsed(BaseException):
    pass
_original = argparse.ArgumentParser.parse_args
def _stop_after_parse(self, args=None, namespace=None):
    raise _Parsed(_original(self, args, namespace))
argparse.ArgumentParser.parse_args = _stop_after_parse
script = sys.argv[1]
sys.argv = sys.argv[1:]
sys.path.insert(0, str(Path(script).resolve().parent))
try:
    runpy.run_path(script, run_name="__main__")
except _Parsed as parsed:
    loaded = sys.modules.get("magestic")
    print(json.dumps({"namespace": {k: str(v) for k, v in vars(parsed.args[0]).items()},
                      "magestic_file": getattr(loaded, "__file__", None)}))
    sys.exit(0)
except ModuleNotFoundError as missing:
    print(f"MISSING_DEPENDENCY {missing.name}", file=sys.stderr)
    sys.exit(4)
print("PARSER_NEVER_CALLED", file=sys.stderr)
sys.exit(3)
"""

RULE_NAMES = sorted(_rule_blocks())


# ---------------------------------------------------------------------------------------------
# Rule harness positive control
# ---------------------------------------------------------------------------------------------

def test_rule_harness_renders_known_good_rule07_like_snakemake(tmp_path):
    """POSITIVE CONTROL for the harness. Rule 07 was already repaired and is known to match its script
    (Snakefile.rules:123-165). Rendering it must reproduce the flags and rule defaults verbatim. A harness
    bug would otherwise show up as a false rule failure elsewhere in this file."""
    rule = RenderedRule("tier_07_classification", tmp_path)
    argv = rule.argv()
    opts = _options(argv)
    assert argv[0] == PYTHON and Path(argv[1]) == SNAKEMAKE_DIR / "07_tier_classification.py"
    assert opts["--bc1-reference"].endswith("keyfile_dir/bc1_reference_table.tsv")
    assert opts["--sample-key"].endswith("keyfile_dir/sample_key.tsv")
    assert opts["--library"] == "yL1"
    assert opts["--output-tier-mapping"].endswith("bc1_counts_dir/scr_bc1_tier_mapping_yL1.tsv")
    assert (opts["--tier-1-min"], opts["--tier-2-min"], opts["--t0-condition"], opts["--chunksize"]) == (
        "100", "20", "YPD", "2000000")
    assert "|" not in argv and "tee" not in argv


# ---------------------------------------------------------------------------------------------
# Shipped rules: argparse acceptance and wildcard expansion (UP-C020, WF-SCREEN12/15/16; NEW rule 10)
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("rule_name", RULE_NAMES)
def test_shipped_rule_argv_is_accepted_by_its_script(rule_name, tmp_path):
    """The command line rendered from each shipped rule must pass the script's own argparse.

    Failing before: qc_03_sample (UP-C020), slope_12_fitness (WF-SCREEN12: `--output` ambiguous, required
    flags absent), aggregate_15_hits (WF-SCREEN15), rarefaction_16 (WF-SCREEN16). Each exits 2 today.
    Positive controls: rules 02, 07, 08, 10, 11, 13, 14. Only the parse step runs (a patched
    ArgumentParser.parse_args stops the real script right after parsing), so this isolates the rule/CLI
    contract from data. Scripts whose module imports need an audit-only dependency (03: matplotlib;
    11: pydeseq2; 12: statsmodels) skip without it."""
    rule = RenderedRule(rule_name, tmp_path)
    argv = rule.argv()
    assert argv[0] == PYTHON, argv
    script = Path(argv[1])
    assert script.parent == SNAKEMAKE_DIR and script.is_file(), f"{rule_name} calls a missing script {script}"
    result = _run([PYTHON, "-c", _PARSE_ONLY_DRIVER, *argv[1:]], tmp_path, timeout=180)
    _skip_if_optional_dependency_missing(result, rule_name)
    assert result.returncode == 0, f"{rule_name}: argparse rejected the shipped command\n" + _explain(result, argv)
    parsed = json.loads(result.stdout.strip().splitlines()[-1])
    if parsed["magestic_file"]:
        assert CHECKOUT_SRC.resolve() in Path(parsed["magestic_file"]).resolve().parents


@pytest.mark.parametrize("rule_name", RULE_NAMES)
def test_shipped_rule_argv_has_no_unexpanded_wildcards(rule_name, tmp_path):
    """No rendered argument may keep a literal `{wildcard}`.

    Snakemake leaves a pathlib.Path param unexpanded (Snakefile.rules:197-199), so a param like
    `BC1_COUNTS_DIR / "counts_for_deseq" / "{library}"` reaches the script as a directory literally named
    `{library}`. Failing before: split_10_comparisons (NEW, not in the register; its declared outputs at
    Snakefile.rules:236-238 would then be missing), slope_12_fitness and aggregate_15_hits (`input_dir`).
    The other rules are positive controls."""
    argv = RenderedRule(rule_name, tmp_path).argv()
    leftovers = [token for token in argv if "{" in token or "}" in token]
    assert not leftovers, f"{rule_name}: unexpanded wildcard(s) in the shell command: {leftovers}"


# ---------------------------------------------------------------------------------------------
# Step 13 (SCR-D03a, SCR-D03b; NEW codon_variant_type contract)
# ---------------------------------------------------------------------------------------------

STATUS = "annotation_status"
ANNOTATED = "annotated"
OLIGO_NOT_IN_REGISTRY = "oligo_not_in_registry"
AMBIGUOUS = "ambiguous_reference"
NO_OLIGO = "no_oligo_in_reference"
NOT_IN_REFERENCE = "bc1_not_in_reference"

REGISTRY_ENTRIES = {
    "O_SYN": dict(codon_variant_type="synonymous", systematic_gene_name_locus="YAL001C",
                  common_gene_name_locus="TFC3", aa_change="A10A", category="codon_changes"),
    "O_MIS": dict(codon_variant_type="missense", systematic_gene_name_locus="YAL002W",
                  common_gene_name_locus="VPS8", aa_change="A20V", category="natural_variants"),
}


def _write_registry(common_dir: Path, entries: dict) -> Path:
    """Harmonized-registry fixture carrying every column create_annotation_lookup() requests.

    Built from the module's own column lists: `usecols=` (utils/oligo_annotations.py:168) raises on any
    missing column, and a hand-written list would later fail in a way that looks like a defect."""
    from magestic.utils import oligo_annotations as oa

    columns = list(oa.COLUMN_MAPPING_BC1_FROM_SCREEN) + list(oa.EXTRA_COLUMNS_BC1_FROM_SCREEN)
    rows = []
    for name, values in entries.items():
        row = {column: np.nan for column in columns}
        row.update(oligo_name=name, **values)
        rows.append(row)
    return _write_tsv(common_dir / "annotation_files" / oa.HARMONIZED_OLIGO_TABLE.name,
                      pd.DataFrame(rows, columns=columns))


def _deseq_results(bc1s, lfcs=None) -> pd.DataFrame:
    """Step-11 output shape: bc1-indexed, index label `bc1` (11_deseq2_single.py:143,168-169,186)."""
    lfcs = np.linspace(-1.0, 1.0, len(bc1s)) if lfcs is None else np.asarray(lfcs, dtype=float)
    return pd.DataFrame(
        {"baseMean": 500.0, "log2FoldChange": lfcs, "lfcSE": 0.2, "stat": lfcs / 0.2,
         "pvalue": 0.01, "padj": 0.05, "library": "yL1", "comparison": "cmpA"},
        index=pd.Index(list(bc1s), name="bc1"),
    )


def _step13_fixture(rule: RenderedRule, tmp_path: Path, tier_rows: pd.DataFrame) -> tuple:
    bc1s = ["bcU1", "bcU2", "bcAMB", "bcNOREF", "bcNOREG", "bcNULL"]
    _write_tsv(rule.input["deseq2_results"], _deseq_results(bc1s), index=True)
    _write_tsv(rule.input["ref_table"], pd.DataFrame({
        "bc1": ["bcU1", "bcU2", "bcU2", "bcAMB", "bcAMB", "bcNOREG", "bcNULL"],
        "oligo_name": ["O_SYN", "O_MIS", "O_MIS", "O_SYN", "O_MIS", "O_ABSENT", np.nan],
    }))
    _write_tsv(rule.input["tier_mapping"], tier_rows)
    _write_tsv(rule.input["oligo_design"], pd.DataFrame({"oligo_name": ["O_SYN", "O_MIS"]}))
    common = tmp_path / "common"
    _write_registry(common, REGISTRY_ENTRIES)
    return bc1s, common


def test_step13_rule_annotates_every_result_row_exactly_once(tmp_path):
    """SCR-D03a + SCR-D03b through the shipped rule `annotate_13_results`.

    Today step 13 raises TypeError at 13_annotate_results.py:112 on every run. Behind that, the core merges on
    `oligo_name` before any bc1 join (core/annotation.py:311-315). After the repair (section 3.1):
    * rc 0 and exactly one output row per DESeq2 row (6 in, 6 out, bc1 unique). The CA02-shaped duplicate
      tier rows for bcU1 and the exact duplicate reference rows for bcU2 multiply nothing;
    * a unique reference link is annotated from the registry, including `codon_variant_type`;
    * every unmatched barcode is KEPT with an explicit `annotation_status`: ambiguous (two oligos, candidates
      kept, no representative chosen), absent from the reference, oligo absent from the registry, or a
      reference row with no oligo."""
    rule = RenderedRule("annotate_13_results", tmp_path)
    tiers = pd.DataFrame({
        "bc1": ["bcU1", "bcU1", "bcU2", "bcAMB", "bcNOREG", "bcNULL"],
        "t0_total": [150, 150, 40, 300, 120, 25],
        "tier": ["abundance_high", "abundance_high", "abundance_medium", "abundance_high",
                 "abundance_high", "abundance_medium"],
        "oligo_name": ["O_SYN", "O_OTHER", "O_MIS", "O_SYN", "O_ABSENT", np.nan],
    })
    bc1s, common = _step13_fixture(rule, tmp_path, tiers)
    result, argv = _run_rule(rule, tmp_path, MAGESTIC_COMMON_DIR=common)
    assert result.returncode == 0, _explain(result, argv)

    out = pd.read_csv(rule.output["annotated"], sep="\t")
    assert len(out) == len(bc1s), "annotation must neither multiply nor drop result rows"
    assert out["bc1"].is_unique and set(out["bc1"]) == set(bc1s)
    assert (out["library"] == "yL1").all()
    by = out.set_index("bc1")

    assert by.loc["bcU1", "oligo_name"] == "O_SYN"
    assert by.loc["bcU1", "gene_systematic"] == "YAL001C"
    assert by.loc["bcU1", "codon_variant_type"] == "synonymous"
    assert by.loc["bcU1", "tier"] == "abundance_high"
    assert by.loc["bcU1", STATUS] == ANNOTATED
    assert by.loc["bcU2", "oligo_name"] == "O_MIS" and by.loc["bcU2", STATUS] == ANNOTATED

    assert pd.isna(by.loc["bcAMB", "oligo_name"]) and by.loc["bcAMB", STATUS] == AMBIGUOUS
    assert pd.isna(by.loc["bcAMB", "gene_systematic"])
    assert set(str(by.loc["bcAMB", "bc1_reference_oligo_candidates"]).split(";")) == {"O_SYN", "O_MIS"}
    assert int(by.loc["bcAMB", "bc1_reference_n_oligos"]) == 2

    assert pd.isna(by.loc["bcNOREF", "oligo_name"]) and by.loc["bcNOREF", STATUS] == NOT_IN_REFERENCE
    assert by.loc["bcNOREF", "tier"] == "unknown"
    assert by.loc["bcNOREG", "oligo_name"] == "O_ABSENT" and by.loc["bcNOREG", STATUS] == OLIGO_NOT_IN_REGISTRY
    assert pd.isna(by.loc["bcNOREG", "gene_systematic"])
    assert pd.isna(by.loc["bcNULL", "oligo_name"]) and by.loc["bcNULL", STATUS] == NO_OLIGO


def test_step13_refuses_conflicting_tier_rows_instead_of_multiplying(tmp_path):
    """Step 13's tier join (13_annotate_results.py:116-125) must be many-to-one.

    A barcode with two DIFFERENT tier rows (upstream corruption: a tier is a function of one t0 total) would
    today duplicate its result row once the arity bug is fixed. After the repair the step refuses, naming the
    barcode, and writes nothing. Failing before: the run dies on the unrelated TypeError, and neither the
    barcode name nor a refusal appears."""
    rule = RenderedRule("annotate_13_results", tmp_path)
    tiers = pd.DataFrame({
        "bc1": ["bcU1", "bcU1", "bcU2"],
        "t0_total": [150, 150, 40],
        "tier": ["abundance_high", "abundance_medium", "abundance_medium"],
    })
    _, common = _step13_fixture(rule, tmp_path, tiers)
    result, argv = _run_rule(rule, tmp_path, MAGESTIC_COMMON_DIR=common)
    assert result.returncode != 0, _explain(result, argv)
    assert "bcU1" in result.stdout + result.stderr, _explain(result, argv)
    assert not Path(rule.output["annotated"]).exists()


def test_core_annotation_joins_bc1_reference_before_design_annotation():
    """SCR-D03b at the core boundary, isolated from D03a's CLI crash.

    `annotate_deseq2_results` with bc1-keyed step-11 results and a reference table currently raises
    KeyError 'oligo_name' (core/annotation.py:311-315). After the repair it joins bc1->oligo through the
    one-row-per-bc1 projection first. The row count is conserved, an ambiguous barcode abstains with a
    status, and `codon_variant_type` is emitted for step 14."""
    from magestic.pipelines.bc1_from_screen.core.annotation import annotate_deseq2_results

    results = _deseq_results(["bcU1", "bcAMB", "bcNOREF"]).reset_index()
    annotation = pd.DataFrame({"oligo_name": ["O_SYN", "O_MIS"], "gene_systematic": ["YAL001C", "YAL002W"],
                               "variant_type": ["synonymous", "missense"]})
    reference = pd.DataFrame({"bc1": ["bcU1", "bcAMB", "bcAMB"], "oligo_name": ["O_SYN", "O_SYN", "O_MIS"]})
    out = annotate_deseq2_results(results, annotation, bc1_ref_df=reference)
    assert len(out) == 3 and out["bc1"].is_unique
    by = out.set_index("bc1")
    assert by.loc["bcU1", "gene_systematic"] == "YAL001C"
    assert by.loc["bcU1", "codon_variant_type"] == "synonymous"
    assert pd.isna(by.loc["bcAMB", "oligo_name"]) and by.loc["bcAMB", STATUS] == AMBIGUOUS
    assert by.loc["bcNOREF", STATUS] == NOT_IN_REFERENCE


def test_core_annotation_of_oligo_labelled_results_is_unchanged():
    """POSITIVE CONTROL. A caller that already supplies `oligo_name` is a valid use of the core contract
    (register votes on SCR-D03b). It works today and must keep order, row count and values after the repair.
    Unmatched oligos stay as rows with null annotation."""
    from magestic.pipelines.bc1_from_screen.core.annotation import annotate_deseq2_results

    results = pd.DataFrame({"bc1": ["b1", "b2", "b3"], "oligo_name": ["O_SYN", "O_MIS", "O_ABSENT"],
                            "log2FoldChange": [0.1, 0.2, 0.3]})
    annotation = pd.DataFrame({"oligo_name": ["O_SYN", "O_MIS"], "gene_systematic": ["YAL001C", "YAL002W"],
                               "variant_type": ["synonymous", "missense"]})
    out = annotate_deseq2_results(results, annotation)
    assert list(out["bc1"]) == ["b1", "b2", "b3"]
    assert list(out["gene_systematic"].fillna("NA")) == ["YAL001C", "YAL002W", "NA"]
    assert list(out["log2FoldChange"]) == [0.1, 0.2, 0.3]


def test_step13_output_gives_step14_its_synonymous_null(tmp_path):
    """NEW (not in the register): the repaired step 13 must hand step 14 the column its null model keys on.

    `create_annotation_lookup` renames `codon_variant_type` to `variant_type` (utils/oligo_annotations.py:84,
    224-229). Step 14 uses the synonymous null only if `codon_variant_type` exists (14_hit_calling.py:130-142);
    otherwise it silently falls back to the middle-50% null (:146-161). The shipped rules chain 13 -> 14 on 40
    uniquely referenced barcodes (20 synonymous with LFC median 0.5, 20 missense near -2.5). Step 14's
    `null_mean` must equal the synonymous median. The fixture is checked to discriminate: the fallback gives
    -1.0. Failing before: step 13 crashes. A repair that fixes only the arity also fails here."""
    r13 = RenderedRule("annotate_13_results", tmp_path)
    r14 = RenderedRule("hit_14_calling", tmp_path)
    assert Path(r14.input["annotated"]) == Path(r13.output["annotated"]), "rules 13 -> 14 are not wired"

    syn = [f"bcS{i:02d}" for i in range(20)]
    mis = [f"bcM{i:02d}" for i in range(20)]
    syn_lfc, mis_lfc = np.linspace(0.0, 1.0, 20), np.linspace(-3.0, -2.0, 20)
    everything = pd.Series(np.concatenate([syn_lfc, mis_lfc]))
    q25, q75 = everything.quantile([0.25, 0.75])
    fallback_mean = everything[(everything >= q25) & (everything <= q75)].median()
    synonymous_mean = float(np.median(syn_lfc))
    assert abs(fallback_mean - synonymous_mean) > 0.5, "fixture does not discriminate the two nulls"

    bc1s = syn + mis
    oligo = {b: f"O_{b}" for b in bc1s}
    _write_tsv(r13.input["deseq2_results"], _deseq_results(bc1s, np.concatenate([syn_lfc, mis_lfc])), index=True)
    _write_tsv(r13.input["ref_table"], pd.DataFrame({"bc1": bc1s, "oligo_name": [oligo[b] for b in bc1s]}))
    _write_tsv(r13.input["tier_mapping"], pd.DataFrame({"bc1": bc1s, "t0_total": 200, "tier": "abundance_high"}))
    _write_tsv(r13.input["oligo_design"], pd.DataFrame({"oligo_name": list(oligo.values())}))
    common = tmp_path / "common"
    _write_registry(common, {
        **{oligo[b]: dict(codon_variant_type="synonymous", systematic_gene_name_locus="YAL001C") for b in syn},
        **{oligo[b]: dict(codon_variant_type="missense", systematic_gene_name_locus="YAL002W") for b in mis},
    })

    first, argv13 = _run_rule(r13, tmp_path, MAGESTIC_COMMON_DIR=common)
    assert first.returncode == 0, "step 13 failed\n" + _explain(first, argv13)
    second, argv14 = _run_rule(r14, tmp_path)
    assert second.returncode == 0, "step 14 failed\n" + _explain(second, argv14)
    hits = pd.read_csv(r14.output["hits"], sep="\t")
    null_means = hits.loc[hits["tier"] == "abundance_high", "null_mean"].unique()
    assert len(null_means) == 1
    assert null_means[0] == pytest.approx(synonymous_mean), (
        f"step 14 used null_mean={null_means[0]} (middle-50% fallback = {fallback_mean}); "
        "step 13 must emit codon_variant_type")


# ---------------------------------------------------------------------------------------------
# Step 15 (SCR-D06, SCR-D09, WF-SCREEN15)
# ---------------------------------------------------------------------------------------------

def _variant_row(variant, lfc, significant=True, tier="abundance_high", comparison="C_6") -> dict:
    """Documented variant-level hit schema, the shape of the register fixture
    (216 modules/screen_downstream/reproduce_downstream.py:52-53)."""
    return {"sublibrary_variant_combo": variant, "comparison": comparison, "mean_log2fc": lfc,
            "concordance_fraction": 1.0, "within_variant_cv": 0.0, "is_significant": significant,
            "z_score": 4.0, "n_bc1": 2, "category": "natural_variants", "codon_variant_type": "missense",
            "tier": tier, "aa_change": "A1V", "common_gene_name_locus": "G"}


def _variant_hits(comparison) -> pd.DataFrame:
    return pd.DataFrame([
        _variant_row("SpG_NGG_v1", 4.0, comparison=comparison),
        _variant_row("SpG_NGG_v1", 3.0, tier="abundance_medium", comparison=comparison),
        _variant_row("SpG_NGG_v2", 0.1, significant=False, comparison=comparison),
    ])


def test_aggregate_by_tier_reports_cross_tier_validation():
    """SCR-D06. The register's two-row fixture (reproduce_downstream.py:54-55) raises KeyError
    'n_abundance_high_significant' today: two unnamed lambdas per column are labelled <lambda_0>/<lambda_1>,
    and the flattening at core/aggregation.py:491-499 never matches them. After the repair: one row, one
    significant row per tier, cross-tier validated."""
    from magestic.pipelines.bc1_from_screen.core.aggregation import aggregate_by_tier

    hits = pd.DataFrame([_variant_row("SpG_NGG_v", 4.0), _variant_row("SpG_NGG_v", -4.0, tier="abundance_medium")])
    out = aggregate_by_tier(hits)
    assert len(out) == 1, "an empty frame here would be a vacuous pass"
    row = out.iloc[0]
    assert (row["n_abundance_high_rows"], row["n_abundance_medium_rows"]) == (1, 1)
    assert (row["n_abundance_high_significant"], row["n_abundance_medium_significant"]) == (1, 1)
    assert bool(row["cross_tier_validated"])


def test_aggregate_by_amino_acid_change_names_its_list_columns():
    """SCR-D09. The register's complete-schema fixture (reproduce_followup.py:22) raises KeyError 'categories'
    today: a bare lambda keeps the label <lambda>, which the flattening at core/aggregation.py:355-362 does not
    list. After the repair the list-valued columns carry their documented names and values."""
    from magestic.pipelines.bc1_from_screen.core.aggregation import aggregate_by_amino_acid_change

    hits = pd.DataFrame({
        "comparison": ["C_6"] * 2, "sublibrary_variant_combo": ["v1", "v2"],
        "shared_aa_missense_change": ["A1V"] * 2, "aa_change": ["A1V"] * 2,
        "category": ["natural_variants", "codon_changes"], "mean_log2fc": [4.0, -4.0],
        "is_significant": [True, True], "z_score": [4.0, -4.0],
        "common_gene_name_locus": ["G", "OTHER"], "n_bc1": [2, 2],
    })
    out = aggregate_by_amino_acid_change(hits)
    assert len(out) == 1, "an empty frame here would be a vacuous pass"
    row = out.iloc[0]
    assert sorted(row["categories"]) == ["codon_changes", "natural_variants"]
    assert sorted(row["variant_designs"]) == ["v1", "v2"]
    assert sorted(row["genes"]) == ["G", "OTHER"]
    assert list(row["aa_changes"]) == ["A1V"]
    assert row["n_variant_designs"] == 2
    assert {"category_concordant", "interesting_codon_effect", "is_significant_aggregated"} <= set(out.columns)


def test_step15_cli_aggregates_documented_variant_level_hits(tmp_path):
    """SCR-D06 + SCR-D09 + the `tier1_significant` barrier (15_aggregate_hits.py:257-263), through the real
    CLI on its documented variant-level input. Today it dies at the AA aggregation (KeyError 'categories').
    After the repair all declared-by-script outputs exist, no row is lost, and the stats count the one
    variant significant in both tiers."""
    hits_dir, out_dir = tmp_path / "hits", tmp_path / "aggregated"
    _write_tsv(hits_dir / "hits_C_6.tsv", _variant_hits("C_6"))
    _write_tsv(hits_dir / "hits_C_12.tsv", _variant_hits("C_12"))
    argv = [PYTHON, SNAKEMAKE_DIR / "15_aggregate_hits.py", "--input-dir", hits_dir, "--output-dir", out_dir,
            "--screen-name", "S"]
    result = _run(argv, tmp_path)
    assert result.returncode == 0, _explain(result, argv)
    assert len(pd.read_csv(out_dir / "S_variant_hits_all.tsv", sep="\t")) == 6
    for name in ("S_variant_hits_significant.tsv", "hits_summary_by_tier.tsv", "hits_cross_tier_validated.tsv",
                 "hits_aggregated_by_nuclease_library.tsv", "hits_aggregated_by_amino_acid.tsv"):
        assert (out_dir / name).stat().st_size > 0, name
    stats = json.loads((out_dir / "aggregation_stats.json").read_text())
    assert stats["n_cross_tier_validated"] == 1
    assert {"n_tier1_only", "n_tier2_only"} <= set(stats)


def test_step15_refuses_bc1_level_hits_with_an_explained_error(tmp_path):
    """WF-SCREEN15 expects unsupported input to be rejected explicitly before execution. Step 14 writes
    BC1-level rows (14_hit_calling.py:211-289), and no step produces the variant-level schema (SCR-D07,
    deferred). Today such input dies in a KeyError traceback at 15_aggregate_hits.py:165. After the repair it
    exits non-zero with no traceback, names every missing variant-level column, and writes no hit table."""
    hits_dir, out_dir = tmp_path / "hits", tmp_path / "aggregated"
    _write_tsv(hits_dir / "hits_C_6.tsv", pd.DataFrame({
        "bc1": ["b1", "b2"], "baseMean": [100.0, 100.0], "log2FoldChange": [1.0, -1.0], "pvalue": [0.01, 0.2],
        "padj": [0.05, 0.5], "comparison": ["C_6"] * 2, "tier": ["abundance_high"] * 2,
        "oligo_name": ["O1", "O2"], "z_score": [3.5, -0.5], "is_significant": [True, False],
    }))
    argv = [PYTHON, SNAKEMAKE_DIR / "15_aggregate_hits.py", "--input-dir", hits_dir, "--output-dir", out_dir]
    result = _run(argv, tmp_path)
    text = result.stdout + result.stderr
    assert result.returncode != 0, _explain(result, argv)
    assert "Traceback" not in text, "refusal must be an explained error, not a crash\n" + _explain(result, argv)
    for column in ("sublibrary_variant_combo", "mean_log2fc", "concordance_fraction", "within_variant_cv", "n_bc1"):
        assert column in text, f"refusal does not name missing column {column}"
    assert not list(out_dir.glob("*_variant_hits_*.tsv"))


def test_rule15_consumes_exactly_the_files_rule14_writes(tmp_path):
    """WF-SCREEN15 input identity. Rule 15's declared inputs must be rule 14's outputs, and each must match
    the glob the script actually reads (taken from 15_aggregate_hits.py, currently `hits_*.tsv` at :120).
    Failing before: rule 14 writes `{comparison}_hits.tsv` (Snakefile.rules:355), which the glob never
    matches, so even a flag-only fix would raise "No hit files found"."""
    rule15 = RenderedRule("aggregate_15_hits", tmp_path)
    rule14_outputs = [RenderedRule("hit_14_calling", tmp_path, wildcards=dict(DEFAULT_WILDCARDS, comparison=c)).output["hits"]
                      for c in DEFAULT_COMPARISONS]
    inputs = rule15.inputs()
    assert sorted(map(str, inputs)) == sorted(map(str, rule14_outputs))
    glob = re.search(r'\.glob\("([^"]+)"\)', (SNAKEMAKE_DIR / "15_aggregate_hits.py").read_text()).group(1)
    unmatched = [p.name for p in inputs if not fnmatch.fnmatch(p.name, glob)]
    assert not unmatched, f"step 15 reads {glob!r}; rule 14 writes {unmatched}"


def test_rule15_runs_on_rule14_outputs_and_reruns_identically(tmp_path):
    """WF-SCREEN15 end to end through the shipped rule, plus the rerun hazard. Variant-level hit tables are
    placed exactly where rule 15 declares its inputs; the rendered command must exit 0 and produce every
    declared output. It is then run again unchanged. If `--output-dir` equalled `--input-dir`, the second run
    would glob step 15's own `hits_summary_by_tier.tsv` / `hits_aggregated_by_*.tsv` back in as inputs.
    Failing before: argparse exit 2."""
    rule = RenderedRule("aggregate_15_hits", tmp_path)
    for comparison, path in zip(DEFAULT_COMPARISONS, rule.input["hits"]):
        _write_tsv(path, _variant_hits(comparison))
    argv = rule.argv()
    options = _options(argv)
    assert "--input-dir" in options and "--output-dir" in options, f"rule 15 argv lacks the script's flags: {argv}"
    assert Path(options["--output-dir"]).resolve() != Path(options["--input-dir"]).resolve()
    all_hits = Path(options["--output-dir"]) / f"{options.get('--screen-name', 'screen')}_variant_hits_all.tsv"
    expected_rows = 3 * len(DEFAULT_COMPARISONS)

    first, _ = _run_rule(rule, tmp_path)
    assert first.returncode == 0, _explain(first, argv)
    for path in rule.outputs():
        assert path.exists(), f"declared output missing: {path}"
    assert len(pd.read_csv(all_hits, sep="\t")) == expected_rows

    second, _ = _run_rule(rule, tmp_path)
    assert second.returncode == 0, _explain(second, argv)
    assert len(pd.read_csv(all_hits, sep="\t")) == expected_rows, "rerun ingested step 15's own outputs"


# ---------------------------------------------------------------------------------------------
# Step 16 (SCR-D10b, WF-SCREEN16)
# ---------------------------------------------------------------------------------------------

def _rarefaction_counts(path: Path) -> Path:
    """Register fixture (reproduce_downstream.py:44-45): four barcodes, 280 reads; with depths 20/50/100 and one
    iteration the saturation fit succeeds, so the prediction branch is actually reached."""
    return _write_tsv(path, pd.DataFrame({"bc1": list("abcd"), "counts": [100, 80, 60, 40]}))


def _step16_argv(counts: Path, out: Path, *extra) -> list:
    return [PYTHON, SNAKEMAKE_DIR / "16_rarefaction.py", "--input", counts, "--output", out,
            "--data-type", "bc1_t0", "--sample-name", "s", "--depths", "20", "50", "100",
            "--n-iterations", "1", *extra]


def test_step16_predict_depths_writes_prediction_rows(tmp_path):
    """SCR-D10b. With `--predict-depths 500` the script appends a prediction row, then applies `~` to an
    `is_prediction` column that is NaN on observed rows (16_rarefaction.py:424-448). TypeError, nothing written.
    After the repair: rc 0, a boolean flag, exactly one prediction row at depth 500 (proving the branch ran,
    not skipped), and a summary that excludes predictions."""
    out = tmp_path / "rare" / "s_rarefaction.tsv"
    argv = _step16_argv(_rarefaction_counts(tmp_path / "counts.tsv"), out, "--predict-depths", "500")
    result = _run(argv, tmp_path)
    assert result.returncode == 0, _explain(result, argv)
    table = pd.read_csv(out, sep="\t")
    assert table["is_prediction"].dtype == bool
    assert list(table.loc[table["is_prediction"], "depth"]) == [500]
    assert not table.loc[~table["is_prediction"], "depth"].eq(500).any()
    summary = pd.read_csv(out.with_suffix(".summary.tsv"), sep="\t")
    assert 500 not in set(summary["depth"]) and {20, 50, 100, 280} <= set(summary["depth"])


def test_step16_without_prediction_writes_both_tables(tmp_path):
    """POSITIVE CONTROL for step 16: the same fixture without `--predict-depths` works today (the fill at
    16_rarefaction.py:443-444 runs when no prediction row exists) and must keep working."""
    out = tmp_path / "rare" / "s_rarefaction.tsv"
    argv = _step16_argv(_rarefaction_counts(tmp_path / "counts.tsv"), out)
    result = _run(argv, tmp_path)
    assert result.returncode == 0, _explain(result, argv)
    table = pd.read_csv(out, sep="\t")
    assert not table["is_prediction"].any()
    assert out.with_suffix(".summary.tsv").exists()


def test_rule16_runs_on_its_declared_input(tmp_path):
    """WF-SCREEN16 end to end through the shipped rule. A single-sample parsed table (the step-02 schema,
    unique per bc1) is written at every declared input; the rendered command must exit 0 and produce every
    declared output. It must also count each barcode ONCE: rarefying the multi-sample long table instead would
    count (bc1, sample) rows as barcodes (16_rarefaction.py:194,229). Failing before: argparse exit 2
    (--bc1-counts/--tier-mapping/--library/--output-plot are not the script's flags)."""
    rule = RenderedRule("rarefaction_16", tmp_path)
    counts = pd.DataFrame({"bc1": [f"bc{i:02d}" for i in range(12)], "sample_name": "s1",
                           "counts": [200 - 15 * i for i in range(12)]})
    for path in rule.inputs():
        _write_tsv(path, counts)
    result, argv = _run_rule(rule, tmp_path)
    assert result.returncode == 0, _explain(result, argv)
    for path in rule.outputs():
        assert path.exists(), f"declared output missing: {path}"
    main_table = next(p for p in rule.outputs() if p.suffix == ".tsv" and not p.name.endswith(".summary.tsv"))
    table = pd.read_csv(main_table, sep="\t")
    actual = table[table["is_actual"]]
    assert int(actual["n_unique_detected"].iloc[0]) == counts["bc1"].nunique()
    assert int(actual["total_reads"].iloc[0]) == int(counts["counts"].sum())


# ---------------------------------------------------------------------------------------------
# Step 03 QC (UP-C017, UP-C018, UP-C019, UP-C020)
# ---------------------------------------------------------------------------------------------

QC_SAMPLES = ["s1", "s2", "s3", "s4"]


def _qc_key() -> pd.DataFrame:
    """The audit's QC fixture (216 modules/screen_upstream/reproduce_qc_orchestration.py:29-31)."""
    return pd.DataFrame({"sample_name": QC_SAMPLES, "library_ID": ["L"] * 4, "condition": ["YPD"] * 4,
                         "generations": [0] * 4, "screen_date": ["d"] * 4,
                         "inner_primers": ["KR1966-KR1967"] * 4})


def _qc_counts(count_column: str) -> pd.DataFrame:
    return pd.DataFrame([
        {"sample_name": s, "bc1": f"bc{i:02d}", count_column: (i + 1) * (j + 1) + ((i + j) % 3)}
        for j, s in enumerate(QC_SAMPLES) for i in range(12)
    ])


def _load_step03():
    pytest.importorskip("matplotlib", reason="03_sample_qc.py imports matplotlib at module top (:52)")
    pytest.importorskip("sklearn", reason="03_sample_qc.py imports scikit-learn at module top (:58-59)")
    spec = importlib.util.spec_from_file_location("screen_route_step03_under_test", SNAKEMAKE_DIR / "03_sample_qc.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _run_step03_main_without_plots(module, monkeypatch, argv) -> None:
    """Run the real main() with the four plot functions stubbed, so metric/merge/flag logic is tested alone."""
    for name in ("plot_pca", "plot_read_count_distribution", "plot_clustering_dendrogram", "plot_correlation_heatmap"):
        monkeypatch.setattr(module, name, lambda *args, **kwargs: None)
    monkeypatch.setattr(sys, "argv", ["03_sample_qc.py", *map(str, argv)])
    module.main()


def _read_flags(path: Path) -> pd.DataFrame:
    try:
        flags = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["sample_name", "flag"])
    return flags if "flag" in flags.columns else pd.DataFrame(columns=["sample_name", "flag"])


def test_step03_main_reads_the_current_counts_schema(tmp_path, monkeypatch):
    """UP-C017. Current producers write `bc1, sample_name, counts` (02_parse_bc1_reads.py:1115;
    04_combine_counts.py:22). Step 03 hard-codes `sample_counts` (03_sample_qc.py:225,324,337,425,437), so
    main() raises KeyError 'Column not found: sample_counts' and writes nothing. After the repair the
    per-sample read totals in the summary equal the fixture's. Plotting is stubbed to isolate this defect from
    UP-C018."""
    qc = _load_step03()
    counts = _write_tsv(tmp_path / "long.tsv", _qc_counts("counts"))
    key = _write_tsv(tmp_path / "key.tsv", _qc_key())
    out = tmp_path / "qc"
    _run_step03_main_without_plots(qc, monkeypatch, ["--input-long-table", counts, "--keyfile", key,
                                                     "--output-dir", out, "--plots-dir", tmp_path / "plots"])
    summary = pd.read_csv(out / "sample_qc_summary.tsv", sep="\t").set_index("sample_name")
    expected = _qc_counts("counts").groupby("sample_name")["counts"].sum()
    assert summary["total_reads"].to_dict() == expected.to_dict()


def test_step03_main_passes_computed_concordance_to_the_flagger(tmp_path, monkeypatch):
    """UP-C019. SampleMetrics carries 0.0 placeholders for the correlation and PCA columns (03_sample_qc.py:121-123).
    main() merges the computed values on top (:1006, :1013), producing `_x`/`_y` names, so detect_outliers reads
    the bare names as missing (:571, :583) and never flags. With the audit fixture and a valid threshold of 0.99
    (all four computed means are 0.967-0.985), four POOR_REPLICATE_CONCORDANCE flags are expected. Today none
    fire. Legacy `sample_counts` input and stubbed plots isolate this defect from UP-C017 and UP-C018."""
    qc = _load_step03()
    counts = _write_tsv(tmp_path / "long.tsv", _qc_counts("sample_counts"))
    key = _write_tsv(tmp_path / "key.tsv", _qc_key())
    out = tmp_path / "qc"
    _run_step03_main_without_plots(qc, monkeypatch, [
        "--input-long-table", counts, "--keyfile", key, "--output-dir", out, "--plots-dir", tmp_path / "plots",
        "--min-reads", "0", "--min-entropy", "0", "--min-correlation", "0.99"])
    summary = pd.read_csv(out / "sample_qc_summary.tsv", sep="\t")
    suffixed = [c for c in summary.columns if c.endswith(("_x", "_y"))]
    assert not suffixed, f"merge left suffixed metric columns: {suffixed}"
    assert summary["mean_replicate_correlation"].between(0.9, 0.99, inclusive="neither").all()
    assert "pca_distance" in summary.columns
    flags = _read_flags(out / "flagged_samples.tsv")
    poor = flags.loc[flags["flag"] == "POOR_REPLICATE_CONCORDANCE", "sample_name"]
    assert sorted(poor) == QC_SAMPLES


def test_step03_cli_writes_reports_and_plots(tmp_path):
    """UP-C018, and step 03 end to end through its own CLI. plot_pca re-merges condition/library_ID/generations
    into a PCA table that already has them (03_sample_qc.py:471-475, 620-624), which raises KeyError
    'condition' at :627. All three reports are written only after the plots (:1041-1064), so nothing is
    written. After the repair (UP-C018 and UP-C019 together; the heatmap at :859-864 also needs UP-C019) the
    three TSVs and five PNGs exist."""
    pytest.importorskip("matplotlib", reason="03_sample_qc.py imports matplotlib at module top (:52)")
    pytest.importorskip("sklearn", reason="03_sample_qc.py imports scikit-learn at module top (:58-59)")
    counts = _write_tsv(tmp_path / "long.tsv", _qc_counts("sample_counts"))
    key = _write_tsv(tmp_path / "key.tsv", _qc_key())
    out, plots = tmp_path / "qc", tmp_path / "plots"
    argv = [PYTHON, SNAKEMAKE_DIR / "03_sample_qc.py", "--input-long-table", counts, "--keyfile", key,
            "--output-dir", out, "--plots-dir", plots]
    result = _run(argv, tmp_path)
    assert result.returncode == 0, _explain(result, argv)
    for name in ("sample_qc_summary.tsv", "flagged_samples.tsv", "pca_coordinates.tsv"):
        assert (out / name).exists(), name
    for name in ("pca_by_condition.png", "pca_by_library.png", "read_count_distribution.png",
                 "clustering_dendrogram.png", "replicate_correlation_heatmap.png"):
        assert (plots / name).stat().st_size > 0, name


def test_rule03_runs_on_its_declared_inputs(tmp_path):
    """UP-C020 end to end through the shipped rule `qc_03_sample`. A current-schema long table and a sample key
    are written at the declared inputs; the rendered command must exit 0 and produce every declared output.
    Failing before: the rule passes --input-dir/--output-summary/--output-plot, which the parser rejects
    (exit 2). This also requires the UP-C017..019 script fixes, because the rule feeds the current schema."""
    pytest.importorskip("matplotlib", reason="03_sample_qc.py imports matplotlib at module top (:52)")
    pytest.importorskip("sklearn", reason="03_sample_qc.py imports scikit-learn at module top (:58-59)")
    rule = RenderedRule("qc_03_sample", tmp_path)
    for name, value in rule.input.items():
        for path in _flat([value]):
            _write_tsv(path, _qc_key() if "key" in name else _qc_counts("counts"))
    result, argv = _run_rule(rule, tmp_path)
    assert result.returncode == 0, _explain(result, argv)
    for path in rule.outputs():
        assert path.exists(), f"declared output missing: {path}"


# ---------------------------------------------------------------------------------------------
# Step 12 (WF-SCREEN12)
# ---------------------------------------------------------------------------------------------

def test_rule12_runs_on_its_declared_inputs(tmp_path):
    """WF-SCREEN12 end to end through the shipped rule `slope_12_fitness`. A wide matrix in the step-08 layout
    (bc1, tier, samples; 08_create_count_matrix.py:553-582), a sample key and a reference are written at the
    declared inputs. Three generations with two replicates each sit under the rule's `condition` wildcard. The
    command must exit 0, produce every declared output and fit at least one barcode (a vacuity guard).
    Failing before: argparse exit 2 (`--output` ambiguous; --count-matrix etc. absent). Note: a failure AFTER
    the rule fix would come from the slope core, not the rule; triage it as a separate defect."""
    pytest.importorskip("statsmodels", reason="core/slope_fitness.py imports statsmodels at module top (:32-33)")
    rule = RenderedRule("slope_12_fitness", tmp_path)
    samples = {f"g{g}_r{r}": g for g in (0, 6, 12) for r in (1, 2)}
    bc1s = [f"bc{i:02d}" for i in range(6)]
    matrix = pd.DataFrame({"bc1": bc1s, "tier": "abundance_high",
                           **{s: [100 + 10 * i + 5 * g for i in range(6)] for s, g in samples.items()}})
    key = pd.DataFrame({"sample_name": list(samples), "condition": DEFAULT_WILDCARDS["condition"],
                        "generations": list(samples.values()), "library": DEFAULT_WILDCARDS["library"]})
    reference = pd.DataFrame({"bc1": bc1s, "oligo_name": [f"O{i // 2}" for i in range(6)],
                              "sublibrary_variant_combo": [f"V{i // 2}" for i in range(6)]})
    for name, value in rule.input.items():
        for path in _flat([value]):
            if "matrix" in name:
                _write_tsv(path, matrix)
            elif "key" in name:
                _write_tsv(path, key)
            else:
                _write_tsv(path, reference)
    result, argv = _run_rule(rule, tmp_path)
    assert result.returncode == 0, _explain(result, argv)
    for path in rule.outputs():
        assert path.exists(), f"declared output missing: {path}"
    bc1_table = next(p for p in rule.outputs() if "bc1" in p.name)
    assert len(pd.read_csv(bc1_table, sep="\t")) > 0


# ---------------------------------------------------------------------------------------------
# Shell launcher (UP-C021) and packaged runner (UP-C012)
# ---------------------------------------------------------------------------------------------

def test_run_full_pipeline_passes_bash_syntax_check(tmp_path):
    """UP-C021. `bash -n` must accept the shipped launcher. Today it fails at run_full_pipeline.sh:140
    (`<=` inside `[[ ]]`); the same construct recurs at :151, :164, :192, :204 and :215, and rc 0 requires all
    six fixed. `-n` parses only; nothing is executed. This does not make the launcher runnable (REPAIR_DESIGN
    section 5, decision 5)."""
    result = subprocess.run(["/bin/bash", "-n", str(SCREEN / "run_full_pipeline.sh")],
                            capture_output=True, text=True, timeout=30, cwd=tmp_path)
    assert result.returncode == 0 and not result.stderr.strip(), _explain(result)


_CLI_SHIM = ("import sys; from magestic.pipelines.bc1_from_screen._cli import main; "
             "sys.argv = ['magestic-screen-analysis'] + sys.argv[1:]; main()")


def test_packaged_screen_runner_does_not_report_success_without_outputs(tmp_path):
    """UP-C012 / docs/ARCHITECTURE_POLICY.md section 7. The `magestic-screen-analysis` entry point
    (pyproject.toml:102 -> _cli.py:10) given a valid but empty project currently logs "PIPELINE COMPLETED
    SUCCESSFULLY" and exits 0 with no outputs: every step returns True (run_pipeline.py:62-166). With nothing
    to analyse, a correct runner must fail. This holds for both repair designs: one that runs the steps
    (missing inputs) and the recommended one that fails loudly (steps not implemented)."""
    project = tmp_path / "project"
    project.mkdir()
    argv = [PYTHON, "-c", _CLI_SHIM, "--project-dir", project]
    result = _run(argv, tmp_path)
    assert result.returncode != 0, _explain(result, argv)
    assert "COMPLETED SUCCESSFULLY" not in result.stdout + result.stderr


def test_packaged_screen_runner_rejects_unknown_steps(tmp_path):
    """UP-C012, second false-success path. An unknown `--steps` name is only warned about and skipped
    (run_pipeline.py:292-294), so a typo still ends in "PIPELINE COMPLETED SUCCESSFULLY", exit 0. After the
    repair it is a usage error that names the step."""
    project = tmp_path / "project"
    project.mkdir()
    argv = [PYTHON, "-c", _CLI_SHIM, "--project-dir", project, "--steps", "no_such_step"]
    result = _run(argv, tmp_path)
    text = result.stdout + result.stderr
    assert result.returncode != 0, _explain(result, argv)
    assert "no_such_step" in text and "COMPLETED SUCCESSFULLY" not in text
