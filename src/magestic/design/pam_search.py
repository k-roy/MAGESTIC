"""
magestic.design.pam_search — Genome-wide PAM site enumeration and off-target analysis.

This module produces the pre-computed guide databases consumed by the oligo
design step.  For each PAM site in the genome it:

1. Records the 20-nt guide sequence (full-length spacer).
2. Computes 0-mismatch and 1-mismatch off-target site counts using a truncated
   (19-nt) seed-region search.

The output is written to a plain-text file with the format::

    full_guide  num_mismatches  mismatch_position  PAM_type  truncated_target
    chrom  idx  strand  full_guide  full_PAM

This file is then processed by :mod:`magestic.design.guide_selection` to
split guides into allowed/disallowed sets.

Ported from:
  ``common/scripts/guide_donor_design/process_NGNN_guides_genomewide.py``
  ``common/scripts/guide_donor_design/process_NGG_and_NAG_guides_genomewide.py``

Typical usage::

    from magestic.design.guide_donor_functions import load_genome
    from magestic.design.nuclease import SPG
    from magestic.design.pam_search import search_pam_sites_genomewide

    genome = load_genome("MAGESTIC_background_strain.fasta")
    search_pam_sites_genomewide(genome, SPG, "guide_off_targets_NGNN.txt")

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from magestic.design.guide_donor_functions import hamming_distance, rev_comp
from magestic.design.nuclease import NucleaseConfig

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

# (chrom, coord, strand, full_guide, full_PAM)
_SiteInfo = Tuple[str, int, str, str, str]

# truncated_guide -> list of sites that match it
_TruncatedDict = Dict[str, List[_SiteInfo]]

# Length of the truncated guide used for 1-mismatch off-target lookup
_SEED_LENGTH = 19


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def search_pam_sites_genomewide(
    genome_seq: dict,
    nuclease: NucleaseConfig,
    output_file: str | Path,
    *,
    seed_length: int = _SEED_LENGTH,
    progress_interval: int = 500_000,
) -> Path:
    """
    Enumerate all PAM sites genome-wide and write a 0/1-mismatch off-target
    analysis file.

    This is a long-running computation (~10–60 min depending on PAM breadth).
    Run it once per genome/nuclease combination and store the result as a
    permanent reference file.

    Parameters
    ----------
    genome_seq : dict
        ``{chrom: sequence_string}`` from :func:`~guide_donor_functions.load_genome`.
    nuclease : NucleaseConfig
        Nuclease configuration specifying the PAM pattern and guide orientation.
    output_file : str or Path
        Destination for the off-target analysis TSV.
    seed_length : int
        Length of the truncated guide used for 1-mismatch seed-region search.
        Default 19 (PAM-proximal bases, which matter most for specificity).
    progress_interval : int
        Log a progress message every N guides processed. Default 500 000.

    Returns
    -------
    Path
        The path to the written output file.

    Notes
    -----
    The off-target allowance criterion used downstream (:func:`~guide_selection.
    assign_allowable_guides`) is: a guide is **allowed** if it has exactly
    **one** target site counting both perfect matches and 1-mismatch hits in
    the seed region (positions 1–10 from the PAM).  The 1-mismatch threshold
    is applied to the ``seed_length``-nt truncated guide.
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    pam_pattern = nuclease.pam_pattern
    pam_length = nuclease.pam_length
    guide_length = nuclease.guide_length
    pam_upstream = nuclease.pam_upstream
    all_pams = set(nuclease.all_pams)

    logger.info(
        "Enumerating PAM sites for %s (%s, %d-bp guide, %d chromosomes)...",
        nuclease.name, pam_pattern, guide_length, len(genome_seq),
    )

    # Pass 1 — collect all full-length guides and index by truncated seed
    truncated_dict: _TruncatedDict = {}
    full_length_guides: set[str] = set()

    for chrom, chrom_seq in genome_seq.items():
        logger.debug("  scanning %s (%d bp)...", chrom, len(chrom_seq))
        for idx in range(guide_length + pam_length, len(chrom_seq) - guide_length - pam_length):
            # Forward strand
            fwd_pam = chrom_seq[idx: idx + pam_length]
            if fwd_pam in all_pams:
                if pam_upstream:
                    full_guide = chrom_seq[idx - guide_length: idx]
                else:
                    full_guide = chrom_seq[idx + pam_length: idx + pam_length + guide_length]
                if "N" not in full_guide and len(full_guide) == guide_length:
                    trunc = full_guide[-seed_length:]
                    truncated_dict.setdefault(trunc, []).append(
                        (chrom, idx, "+", full_guide, fwd_pam)
                    )
                    full_length_guides.add(full_guide)

            # Reverse strand (check RC of the window ending at idx)
            rev_pam_rc = rev_comp(chrom_seq[idx - pam_length: idx])
            if rev_pam_rc in all_pams:
                if pam_upstream:
                    full_guide = rev_comp(chrom_seq[idx: idx + guide_length])
                else:
                    full_guide = rev_comp(
                        chrom_seq[idx - pam_length - guide_length: idx - pam_length]
                    )
                if "N" not in full_guide and len(full_guide) == guide_length:
                    trunc = full_guide[-seed_length:]
                    truncated_dict.setdefault(trunc, []).append(
                        (chrom, idx, "-", full_guide, rev_pam_rc)
                    )
                    full_length_guides.add(full_guide)

    logger.info(
        "  Found %d unique full-length guides; writing off-target analysis to %s",
        len(full_length_guides), output_file,
    )

    # Pass 2 — write off-target records
    dash_prefix = "-" * (guide_length - seed_length)
    header_cols = [
        "full_guide", "num_mismatches", "mismatch_position", "PAM_type",
        "truncated_target", "chrom", "idx", "strand", "full_guide", "full_PAM",
    ]

    guides_written = 0
    with open(output_file, "w") as outfile:
        outfile.write("\t".join(header_cols) + "\n")

        for full_guide in full_length_guides:
            guides_written += 1
            if guides_written % progress_interval == 0:
                logger.info("  processed %d guides...", guides_written)

            trunc = full_guide[-seed_length:]

            # 0-mismatch targets
            perfect_targets = truncated_dict.get(trunc, [])
            for chrom, coord, strand, tgt_guide, tgt_pam in perfect_targets:
                row = [
                    full_guide, "0", "NA", tgt_pam,
                    dash_prefix + trunc,
                    chrom, str(coord), strand, tgt_guide, tgt_pam,
                ]
                outfile.write("\t".join(row) + "\n")

            # 1-mismatch targets (mutate each position in the seed)
            trunc_list = list(trunc)
            for pos_idx in range(seed_length):
                orig_base = trunc_list[pos_idx]
                for mismatch_base in "acgt":
                    if mismatch_base.upper() == orig_base:
                        continue
                    mutant = trunc_list[:]
                    mutant[pos_idx] = mismatch_base
                    mutant_seq = "".join(mutant).upper()
                    for chrom, coord, strand, tgt_guide, tgt_pam in truncated_dict.get(
                        mutant_seq, []
                    ):
                        row = [
                            full_guide, "1", str(pos_idx + 1), tgt_pam,
                            dash_prefix + mutant_seq,
                            chrom, str(coord), strand, tgt_guide, tgt_pam,
                        ]
                        outfile.write("\t".join(row) + "\n")

    logger.info("Off-target analysis complete: %d guides written.", guides_written)
    return output_file
