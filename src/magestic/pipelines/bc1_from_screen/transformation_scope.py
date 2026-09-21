"""Transformation scope and generation parsing shared by screen steps 10 and 11.

Kept free of magestic imports so the numbered scripts can load it as a sibling
module under interpreters that lack the installed package, exactly as they load
``config.py`` and ``library_scope.py``.
"""
import math

#: Columns that identify a design pool, not a transformation (normalised names). SPS,
#: V/cV library IDs and oligo-name prefixes name what was synthesized; trafo_3 and
#: trafo_4 share V634, so none of these can stand in for transformation identity.
DESIGN_POOL_COLUMNS = frozenset((
    "oligo_name", "sps", "v_library", "v_harmonized", "cv_library", "corresponding_step_1_v",
    "oligo_pool", "v_library_parent_subpool", "guide_donor_oligo_subpool", "step_1_plasmid_library",
))

#: Per-sample key columns that carry transformation identity.
SAMPLE_TRANSFORMATION_COLUMNS = ("trafo_ID", "yeast_sublibrary")


def normalise_column_name(name) -> str:
    """Lower-case, spaces to underscores, for comparing column names."""
    return str(name).strip().lower().replace(" ", "_")


def is_design_pool_column(name) -> bool:
    """True if the column names a design pool rather than a transformation."""
    return normalise_column_name(name) in DESIGN_POOL_COLUMNS


def parse_generation(raw) -> int:
    """Generation as an integer: 6, 6.0, "6.0" and " 6 " are 6; "-5" is -5.

    Missing values, "none", text, fractional values such as 0.5 and non-finite values
    raise ValueError instead of becoming 0 (a baseline) or being truncated.
    """
    try:
        value = float(str(raw).strip())
    except (TypeError, ValueError):
        raise ValueError(f"not a number: {raw!r}") from None
    if not math.isfinite(value) or not value.is_integer():
        raise ValueError(f"not an integer generation: {raw!r}")
    return int(value)
