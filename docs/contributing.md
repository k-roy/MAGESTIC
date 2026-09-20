# Contributing

Contributions are welcome. MAGESTIC is a research codebase; the priorities are
**reproducibility** of published libraries and **clarity** for collaborators.

## Development setup

```bash
git clone https://github.com/k-roy/MAGESTIC.git
cd MAGESTIC
pip install -e ".[full,dev]"
```

## Tooling

The repository is configured for:

- **pytest** — `pytest` (markers: `slow`, `integration`, `unit`)
- **ruff** — linting (`ruff check src`)
- **black** — formatting (`black src`)
- **mypy** — type checking (`mypy`)

```bash
pytest -m "not slow"
ruff check src && black --check src
```

## Conventions

- Keep parameters in each pipeline's `config.py` dataclass — no hardcoded paths.
- Pipeline steps must be **idempotent** and skip work already done (clusters
  preempt jobs).
- Preserve `__component_version__` on migrated pipelines; changing analysis
  behavior that affects published libraries needs a version bump and a note in
  the [Changelog](changelog.md).
- Document public functions with Google-style docstrings — they render into the
  [API Reference](api/index.md).

## Reporting issues

Open an issue at
[github.com/k-roy/MAGESTIC/issues](https://github.com/k-roy/MAGESTIC/issues),
or contact Kevin R. Roy ([kevinrjroy@gmail.com](mailto:kevinrjroy@gmail.com)).
