"""
magestic.design.guide_selection — Allowable guide classification from off-target analysis.

This module reads the off-target analysis file produced by
:mod:`magestic.design.pam_search` and partitions every guide into an
*allowed* set (unique in the genome for practical purposes) or a *disallowed*
set (two or more effective target sites).

The allowance criterion mirrors the original ``assign_allowable_guides.py``
script:

    A guide is **allowed** if its total count of effective targets == 1,
    where an effective target is any record whose ``mismatch_position`` is
    either ``"NA"`` (perfect match) or an integer **> 10** (mismatch outside
    the PAM-proximal 10-bp seed window).

    Mismatches at positions 1–10 (seed) are *not* counted because they are
    expected to permit HDR repair with high fidelity, just like the intended
    edit.

Ported from:
  ``common/scripts/guide_donor_design/assign_allowable_guides.py``

Typical usage::

    from magestic.design.guide_selection import assign_allowable_guides, load_allowed_guides

    assign_allowable_guides(
        "guide_off_targets_NGNN.txt",
        "allowed_guides_NGNN.txt",
        "disallowed_guides_NGNN.txt",
    )
    allowed = load_allowed_guides("allowed_guides_NGNN.txt")  # set of guide seqs

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Set, Tuple

logger = logging.getLogger(__name__)

# Output column header shared by both allowed and disallowed files
_OUTPUT_HEADER = "chrom\tPAM_coord\tstrand\tfull_guide\tfull_PAM\n"

# Seed region boundary: mismatches at positions > this threshold are counted
# as effective off-target hits (positions 1–10 are considered within the seed
# and tolerated because they will be corrected by the donor).
_SEED_BOUNDARY = 10


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def assign_allowable_guides(
    off_target_file: str | Path,
    allowed_outfile: str | Path,
    disallowed_outfile: str | Path,
    *,
    progress_interval: int = 100_000,
) -> Tuple[int, int]:
    """
    Read the off-target analysis file and write allowed/disallowed guide lists.

    Parameters
    ----------
    off_target_file : str or Path
        Off-target analysis TSV produced by
        :func:`~pam_search.search_pam_sites_genomewide`.
    allowed_outfile : str or Path
        Destination for guides with exactly one effective target site.
    disallowed_outfile : str or Path
        Destination for guides with two or more effective target sites.
    progress_interval : int
        Log a progress message every N guide records processed. Default 100 000.

    Returns
    -------
    (n_allowed, n_disallowed) : tuple of int
        Counts of guides written to each output file.

    Notes
    -----
    The function makes **two passes** over the data:

    Pass 1 — count effective targets per guide sequence.
        A record contributes to the count if ``mismatch_position`` is
        ``"NA"`` (perfect match) or ``int(mismatch_position) > 10``
        (outside the 10-bp seed).

    Pass 2 — write site records to allowed/disallowed files.
        Every (chrom, PAM_coord, strand, full_guide, full_PAM) record from
        the input is routed to the appropriate output file based on the
        count computed in Pass 1.  Note that the output rows describe the
        *target sites*, not the query guides — one guide may appear on
        multiple rows if it has multiple perfect-match locations (all of
        those rows will go to the disallowed file).
    """
    off_target_file = Path(off_target_file)
    allowed_outfile = Path(allowed_outfile)
    disallowed_outfile = Path(disallowed_outfile)

    allowed_outfile.parent.mkdir(parents=True, exist_ok=True)
    disallowed_outfile.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Pass 1: counting effective targets from %s ...", off_target_file)

    # Pass 1 — count effective targets per guide sequence AND collect all site records
    guides_to_num_targets: dict[str, int] = {}
    # Each element: (chrom, PAM_coord_str, strand, full_target_guide, full_target_PAM)
    site_records: list[tuple[str, str, str, str, str]] = []
    records_read = 0

    with open(off_target_file) as infile:
        _header = infile.readline()  # consume header
        for line in infile:
            records_read += 1
            if records_read % progress_interval == 0:
                logger.info("  %d records read...", records_read)

            parts = line.rstrip("\n").split("\t")
            (
                _full_query_guide,
                _num_mismatches,
                mismatch_position,
                _pam_type,
                _truncated_target,
                chrom,
                pam_coord,
                strand,
                full_target_guide,
                full_target_pam,
            ) = parts

            site_records.append((chrom, pam_coord, strand, full_target_guide, full_target_pam))

            # Initialise counter for this guide if not seen yet
            if full_target_guide not in guides_to_num_targets:
                guides_to_num_targets[full_target_guide] = 0

            # Count only effective off-target hits (perfect match OR mismatch outside seed)
            if mismatch_position == "NA" or int(mismatch_position) > _SEED_BOUNDARY:
                guides_to_num_targets[full_target_guide] += 1

    logger.info(
        "  Pass 1 complete: %d records, %d unique guide sequences.",
        records_read, len(guides_to_num_targets),
    )

    # Pass 2 — partition site records into allowed / disallowed
    logger.info(
        "Pass 2: writing allowed → %s and disallowed → %s ...",
        allowed_outfile, disallowed_outfile,
    )

    n_allowed = 0
    n_disallowed = 0
    records_written = 0

    with open(allowed_outfile, "w") as allowed_fh, open(disallowed_outfile, "w") as disallowed_fh:
        allowed_fh.write(_OUTPUT_HEADER)
        disallowed_fh.write(_OUTPUT_HEADER)

        for chrom, pam_coord, strand, full_guide, full_pam in site_records:
            records_written += 1
            if records_written % progress_interval == 0:
                logger.info(
                    "  %d records written (%d allowed, %d disallowed)...",
                    records_written, n_allowed, n_disallowed,
                )

            row = "\t".join([chrom, pam_coord, strand, full_guide, full_pam]) + "\n"

            if guides_to_num_targets[full_guide] == 1:
                allowed_fh.write(row)
                n_allowed += 1
            else:
                disallowed_fh.write(row)
                n_disallowed += 1

    logger.info(
        "Guide selection complete: %d allowed, %d disallowed site records written.",
        n_allowed, n_disallowed,
    )
    return n_allowed, n_disallowed


# ---------------------------------------------------------------------------
# Loader helpers consumed by oligo_generation
# ---------------------------------------------------------------------------

def load_allowed_guides(allowed_guides_file: str | Path) -> Set[str]:
    """
    Load the set of allowed guide sequences from a file produced by
    :func:`assign_allowable_guides`.

    Parameters
    ----------
    allowed_guides_file : str or Path
        Path to the allowed-guides TSV (columns: chrom, PAM_coord, strand,
        full_guide, full_PAM).

    Returns
    -------
    set of str
        Set of 20-nt guide sequences that passed the uniqueness filter.
    """
    allowed: Set[str] = set()
    with open(allowed_guides_file) as fh:
        _header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                full_guide = parts[3]
                allowed.add(full_guide)
    logger.info("Loaded %d allowed guides from %s", len(allowed), allowed_guides_file)
    return allowed


def load_guide_site_map(
    allowed_guides_file: str | Path,
) -> dict[str, tuple[str, int, str, str]]:
    """
    Load a mapping from guide sequence to its genomic site.

    Parameters
    ----------
    allowed_guides_file : str or Path
        Path to the allowed-guides TSV.

    Returns
    -------
    dict
        ``{full_guide: (chrom, PAM_coord, strand, full_PAM)}``
        mapping for each allowed guide.  If a guide appears on multiple rows
        (which should not happen for an allowed guide by definition), the last
        row wins.
    """
    site_map: dict[str, tuple[str, int, str, str]] = {}
    with open(allowed_guides_file) as fh:
        _header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                chrom, pam_coord_str, strand, full_guide, full_pam = parts[:5]
                site_map[full_guide] = (chrom, int(pam_coord_str), strand, full_pam)
    logger.info(
        "Loaded guide→site map for %d guides from %s", len(site_map), allowed_guides_file,
    )
    return site_map
