"""REDI's source identity must cover every magestic module REDI imports from outside itself.

`core.io.source_digest()` hashes the REDI package plus `EXTERNAL_SOURCES`. This guard
parses REDI's own code (every module and the Snakefile) and fails when a new import
reaches a magestic module outside `pipelines/redi` that the list does not name, so the
next cross-package dependency cannot silently escape cache invalidation again.
"""
import ast
from pathlib import Path

import magestic
from magestic.pipelines.redi.core import io

PACKAGE = Path(magestic.__file__).resolve().parent
REDI = PACKAGE / 'pipelines' / 'redi'


def _module_file(dotted):
    """The source file a `magestic.*` dotted name resolves to, or None."""
    parts = dotted.split('.')[1:]
    candidate = PACKAGE.joinpath(*parts)
    for path in (candidate.with_suffix('.py'), candidate / '__init__.py'):
        if path.is_file():
            return path
    return None


def _imported_files(tree, source):
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = source.parent
                for _ in range(node.level - 1):
                    base = base.parent
                module = 'magestic.' + '.'.join(base.relative_to(PACKAGE).parts + tuple(filter(None, (node.module or '').split('.'))))
            else:
                module = node.module or ''
            names = [module] + [f'{module}.{alias.name}' for alias in node.names]
        else:
            continue
        for name in names:
            if name == 'magestic' or name.startswith('magestic.'):
                path = _module_file(name)
                if path is not None:
                    found.add(path)
    return found


def test_every_cross_package_import_is_in_the_source_identity():
    imported = set()
    for source in sorted(REDI.rglob('*.py')):
        imported |= _imported_files(ast.parse(source.read_text()), source)
    # A Snakefile is not Python as a whole; its top-level import statements are.
    snakefile = REDI / 'snakemake' / 'Snakefile'
    imports = [line for line in snakefile.read_text().splitlines() if line.startswith(('import ', 'from '))]
    assert imports, 'expected top-level imports in the REDI Snakefile'
    imported |= _imported_files(ast.parse('\n'.join(imports)), snakefile)
    external = {p.relative_to(PACKAGE).as_posix() for p in imported if REDI not in p.parents}
    # Package __init__ files reached only as a parent of a listed module carry no REDI logic.
    external = {p for p in external if not p.endswith('__init__.py')}
    assert external == set(io.EXTERNAL_SOURCES), (
        f'REDI imports {sorted(external)} from outside pipelines/redi, but the source identity lists '
        f'{sorted(io.EXTERNAL_SOURCES)}; update core.io.EXTERNAL_SOURCES')


def test_source_files_and_digest_cover_the_declared_scope():
    files = io.source_files()
    assert set(sorted(REDI.rglob('*.py'))) <= set(files)
    assert [p for p in files if REDI not in p.parents] == [PACKAGE / rel for rel in io.EXTERNAL_SOURCES]
    assert all(p.is_file() for p in files)
