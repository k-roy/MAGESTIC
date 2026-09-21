"""Which samples belong to one yeast library: one rule shared by steps 07 and 08.

Kept free of magestic imports so the numbered scripts can load it as a sibling
module under interpreters that lack the installed package, exactly as they load
``config.py``.
"""
import logging
import re
from typing import List, Optional

import pandas as pd

logger = logging.getLogger(__name__)

#: Column names, in priority order, that may carry the library assignment in a
#: sample key. `library_ID` is what the 20250912 combined_key.tsv uses (and what
#: the Feb-2026 production log filtered on).
LIBRARY_COL_CANDIDATES = ("library_ID", "library_id", "library")

#: A yeast-library identifier as written in sample keys (yL406, yL437, ...).
YEAST_LIBRARY_ID = re.compile(r"^yL\d+$")


def resolve_library_samples(
    sample_key: pd.DataFrame,
    library: str,
    sample_name_col: str = "sample_name",
) -> Optional[List[str]]:
    """
    Resolve the set of sample_names belonging to ``library`` (defect D3).

    Tries each of LIBRARY_COL_CANDIDATES in order and uses the first column that
    actually CONTAINS ``library``. A column that exists but never mentions the
    library is not a match -- e.g. yL406's key has a `yL` column holding trafo
    IDs (`trafo_1-1`), which would otherwise resolve to an empty sample set and
    wrongly abort a perfectly good single-library screen.

    If no candidate column contains ``library`` but one of them declares yeast
    libraries (values such as ``yL437``), the key explicitly describes other
    libraries and ``library`` is most likely a typo. That raises instead of
    returning None, because an unfiltered matrix would silently mix libraries.

    Returns:
        List of sample names, or None if no candidate column carries yeast-library
        identifiers at all (caller must warn loudly and skip filtering, NOT filter
        to nothing). This is the single-library-screen case, e.g. yL406.

    Raises:
        ValueError: the key declares yeast libraries, none of which is ``library``.
    """
    if sample_name_col not in sample_key.columns:
        logger.warning(
            f"  Sample key has no '{sample_name_col}' column "
            f"(has {list(sample_key.columns)[:12]}...); cannot apply the library filter."
        )
        return None

    declared = {}
    for col in LIBRARY_COL_CANDIDATES:
        if col not in sample_key.columns:
            continue
        values = sample_key[col].astype(str).str.strip()
        mask = values == library
        if not mask.any():
            logger.info(
                f"  Column '{col}' exists but contains no '{library}' "
                f"(values seen: {sorted(values.unique())[:8]}); "
                f"trying next candidate."
            )
            libraries = sorted({v for v in values if YEAST_LIBRARY_ID.match(v)})
            if libraries:
                declared[col] = libraries
            continue
        samples = sample_key.loc[mask, sample_name_col].dropna().astype(str).unique().tolist()
        logger.info(f"  Samples for library {library}: {len(samples):,} (from column '{col}')")
        return samples

    if declared:
        raise ValueError(
            f"Library '{library}' is not in the sample key, which declares yeast libraries "
            f"{declared}. Refusing to continue without a library filter, because that would "
            f"mix libraries; check the --library value."
        )

    logger.warning(
        f"  NO library column found for '{library}' among {list(LIBRARY_COL_CANDIDATES)} "
        f"in the sample key -- the library filter is NOT being applied and the matrix "
        f"will span every sample present in the counts file. This is correct for a "
        f"single-library screen (e.g. yL406) and WRONG for a multi-library one; if this "
        f"screen has more than one library, add a 'library_ID' column to the sample key."
    )
    return None
