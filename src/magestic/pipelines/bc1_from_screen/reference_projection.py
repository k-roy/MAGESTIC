"""One row per bc1 from the bc1 -> oligo reference table, shared by screen steps 07 and 13.

Kept free of magestic imports so the numbered scripts can load it as a sibling module
under interpreters that lack the installed package, exactly as they load config.py.
"""
import pandas as pd

# Per-row outcome of joining bc1-keyed results to designs. "annotated" means only that
# one legacy reference link named the oligo (and, in step 13, that the registry knows
# it); it is not validated attribution, construct association or genotype.
ANNOTATION_STATUS_COL = "annotation_status"
STATUS_ANNOTATED = "annotated"
STATUS_OLIGO_NOT_IN_REGISTRY = "oligo_not_in_registry"
STATUS_AMBIGUOUS_REFERENCE = "ambiguous_reference"
STATUS_NO_OLIGO_IN_REFERENCE = "no_oligo_in_reference"
STATUS_BC1_NOT_IN_REFERENCE = "bc1_not_in_reference"


def project_bc1_reference(bc1_ref_df: pd.DataFrame, oligo_col: str = "oligo_name") -> pd.DataFrame:
    """Exactly one row per bc1 from a (possibly multi-row) bc1 reference table.

    A bc1 is assigned an oligo only when its rows name exactly one distinct non-null
    oligo. Two or more distinct oligos abstain (oligo NA) and keep the full candidate
    set; no representative row is chosen. Exact duplicate rows do not create ambiguity.
    """
    missing = {"bc1", oligo_col} - set(bc1_ref_df.columns)
    if missing:
        raise ValueError(f"bc1 reference table lacks column(s) {sorted(missing)}")
    ref = bc1_ref_df[["bc1", oligo_col]]
    named = ref.dropna(subset=[oligo_col]).drop_duplicates()
    grouped = named.groupby("bc1")[oligo_col]
    n_oligos = grouped.nunique()
    candidates = grouped.agg(lambda s: ";".join(sorted(s.astype(str).unique())))
    proj = pd.DataFrame({"bc1": pd.unique(ref["bc1"].dropna())})
    proj["bc1_reference_n_oligos"] = proj["bc1"].map(n_oligos).fillna(0).astype(int)
    proj["bc1_reference_oligo_candidates"] = proj["bc1"].map(candidates)
    unique = proj["bc1_reference_n_oligos"] == 1
    proj[oligo_col] = proj["bc1_reference_oligo_candidates"].where(unique)
    status = pd.Series(STATUS_AMBIGUOUS_REFERENCE, index=proj.index)
    status[proj["bc1_reference_n_oligos"] == 0] = STATUS_NO_OLIGO_IN_REFERENCE
    status[unique] = STATUS_ANNOTATED  # provisional; a registry join may demote it
    proj[ANNOTATION_STATUS_COL] = status
    assert proj["bc1"].is_unique
    return proj
