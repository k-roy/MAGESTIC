"""
magestic.design.oligo_generation — VCF-based guide-donor oligo design.

This module is the main orchestrator of the MAGESTIC guide-donor design
pipeline.  Given a processed variant file (output of
:mod:`magestic.design.vcf_processing`), a genome sequence, pre-classified
guide databases (output of :mod:`magestic.design.guide_selection`), and
TIF-seq transcript coverage bedgraphs, it:

1. Builds a homology-directed repair (HDR) donor for each variant.
2. Finds all PAM sites within ``max_dist_to_pam`` bp of the variant.
3. Filters guides for: disruption of re-recognition, off-target uniqueness,
   absence of the internal RE site, poly-T run limits, and exclusion lists.
4. Assigns guides to synthesis subpools based on PAM preference (NGNG > NGNH
   for SpG).
5. Removes residual BspQI/SapI sites at donor junctions.
6. Assembles and writes the final oligo sequences.

Output files
------------
- **oligo_pool** — two-column TSV: ``oligo_name``, ``oligo_sequence``
- **info** — full design metadata (guide, donor, PAM, disruption position, …)
- **untargetable** — variants for which no acceptable guide was found and why

Ported from:
  ``projects/QTL/20240404_SpG_QTL_guide_donor_libraries/final_scripts/
  generate_guide_donor_sequences_from_vcf.py``

  ``projects/QTL/20240404_SpG_QTL_guide_donor_libraries/final_scripts/
  guide_donor_functions_for_vcf.py``

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from magestic.design.guide_donor_functions import rev_comp
from magestic.design.nuclease import NucleaseConfig, expand_pam_motif

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

# (chrom, PAM_coord, strand) → (full_guide, full_PAM)
_CoordsToGuides = Dict[Tuple[str, int, str], Tuple[str, str]]

# (disruption_position, strand, guide, PAM, PAM_coord)
_GuideInfo = Tuple[int, str, str, str, int]

# Default restriction sites for BspQI / SapI Golden Gate assembly
_DEFAULT_RESTRICTION_SITES: Tuple[str, str] = ("GAAGAGC", "GCTCTTC")


# ---------------------------------------------------------------------------
# Guide database loaders
# ---------------------------------------------------------------------------

def load_processed_guides(
    guides_file: str | Path,
    pams_to_use: List[str],
) -> _CoordsToGuides:
    """
    Load a pre-classified guide file (allowed or disallowed) into a coordinate
    lookup dict.

    Parameters
    ----------
    guides_file : str or Path
        Allowed or disallowed guides TSV produced by
        :func:`~guide_selection.assign_allowable_guides`.
        Columns: chrom, PAM_coord, strand, full_guide, full_PAM.
    pams_to_use : list of str
        Concrete PAM sequences to keep (e.g. all expansions of ``"NGNN"``).
        Records whose ``full_PAM`` is not in this set are skipped.

    Returns
    -------
    dict
        ``{(chrom, PAM_coord, strand): (guide, PAM)}``
    """
    pams_set = set(pams_to_use)
    coords_to_guides: _CoordsToGuides = {}
    with open(guides_file) as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pam_coord_str, strand, guide, pam = parts[:5]
            if pam in pams_set:
                coords_to_guides[(chrom, int(pam_coord_str), strand)] = (guide, pam)
    logger.info(
        "Loaded %d guide sites from %s", len(coords_to_guides), guides_file,
    )
    return coords_to_guides


def load_guides_from_fasta(
    fasta_files: List[str | Path],
    guide_length: int,
    pam_length: int,
    pam_upstream: bool,
    pams_to_use: List[str],
) -> Set[str]:
    """
    Collect all guide sequences present in a set of FASTA sequences.

    Used to build an exclusion list of guides that target MAGESTIC backbone
    elements (guide-donor plasmid, Cas9/Cas12a expression cassettes, barcode
    loci).  Any guide that appears in these sequences cannot be used safely
    in the library.

    Parameters
    ----------
    fasta_files : list of str or Path
        FASTA files to scan. Sequences are treated as linear for non-plasmid
        contexts; a 200-bp wrap-around is added to handle circular plasmids.
    guide_length : int
        Length of the protospacer sequence (e.g. 20).
    pam_length : int
        Length of the PAM motif (e.g. 4 for NGNN).
    pam_upstream : bool
        ``True`` if the PAM is 3′ of the guide (Cas9); ``False`` if 5′ (Cas12a).
    pams_to_use : list of str
        Concrete PAM sequences to consider.

    Returns
    -------
    set of str
        Uppercase guide sequences that must be excluded.
    """
    pams_set = set(pams_to_use)
    guides_to_exclude: Set[str] = set()

    for fasta_file in fasta_files:
        seq = ""
        with open(fasta_file) as fh:
            for line in fh:
                if not line.startswith(">"):
                    seq += line.strip().upper()
        # Account for circularity with a 200-bp wrap-around
        extended_seq = seq + seq[:200]
        for idx in range(30, len(extended_seq) - 30):
            fwd_pam = extended_seq[idx: idx + pam_length]
            if fwd_pam in pams_set:
                if pam_upstream:
                    guides_to_exclude.add(extended_seq[idx - guide_length: idx])
                else:
                    guides_to_exclude.add(
                        extended_seq[idx + pam_length: idx + pam_length + guide_length]
                    )
            rev_pam = rev_comp(extended_seq[idx - pam_length: idx])
            if rev_pam in pams_set:
                if pam_upstream:
                    guides_to_exclude.add(rev_comp(extended_seq[idx: idx + guide_length]))
                else:
                    guides_to_exclude.add(
                        rev_comp(
                            extended_seq[idx - pam_length - guide_length: idx - pam_length]
                        )
                    )

    logger.info(
        "Loaded %d guides to exclude from %d FASTA file(s)",
        len(guides_to_exclude), len(fasta_files),
    )
    return guides_to_exclude


# ---------------------------------------------------------------------------
# Donor construction
# ---------------------------------------------------------------------------

def construct_donor(
    chrom: str,
    variant_coord: int,
    ref: str,
    alt: str,
    chrom_length: int,
    individual_variants_lst: List[List[str]],
    guide_length: int,
    donor_length: int,
    minimum_homology: int,
    genome_seq: dict,
    plus_strand_bedgraph: dict,
    minus_strand_bedgraph: dict,
    restriction_sites: Tuple[str, str] = _DEFAULT_RESTRICTION_SITES,
) -> Optional[Tuple]:
    """
    Build a homology-directed repair (HDR) donor for a single variant.

    The donor is centred on the variant with equal homology arms.  The strand
    (sense vs. antisense) is determined by TIF-seq coverage: the donor is
    assembled sense to the majority transcript strand to prevent the donor
    ssDNA from hybridising to the mRNA template.

    If the donor contains a BspQI/SapI recognition site that cannot be
    removed by shifting the homology arms (while maintaining
    ``minimum_homology`` bp on each side), ``None`` is returned.

    Parameters
    ----------
    chrom : str
        Chromosome name (e.g. ``"chrI"``).
    variant_coord : int
        1-based genomic position of the first variant base.
    ref, alt : str
        Reference and alternate allele strings.
    chrom_length : int
        Total length of the chromosome sequence.
    individual_variants_lst : list of list
        Parsed individual variant records: each inner list is
        ``[pos_str, ref, alt, variant_type]``.
    guide_length, donor_length, minimum_homology : int
        Design parameters forwarded from the nuclease/oligo config.
    genome_seq : dict
        ``{chrom: sequence_string}`` from :func:`~guide_donor_functions.load_genome`.
    plus_strand_bedgraph, minus_strand_bedgraph : dict
        Per-position transcript coverage from
        :func:`~guide_donor_functions.load_transcript_coverage`.
    restriction_sites : tuple of str
        RE recognition sequences to avoid in the donor.

    Returns
    -------
    tuple or None
        ``(donor_shift, donor_start_coord, donor_end_coord, donor_strand,
        left_side_donor, right_side_donor, ref, alt, donor, genomic_target)``
        on success, or ``(None, donor_string)`` if no valid donor exists.
    """
    ref_length = len(ref)
    start_of_variant = variant_coord - 1  # convert to 0-based
    end_of_variant = variant_coord + ref_length - 1

    left_arm = (donor_length - len(alt)) // 2
    right_arm = donor_length - len(alt) - left_arm

    # Handle variants near chromosome boundaries
    donor_shift = 0
    at_right_end = end_of_variant + right_arm > chrom_length
    at_left_end = start_of_variant - left_arm < 1
    if at_right_end:
        donor_shift = chrom_length - (end_of_variant + right_arm)  # negative
    if at_left_end:
        donor_shift = left_arm - start_of_variant  # positive

    max_donor_shift = right_arm - minimum_homology

    left_side_donor = genome_seq[chrom][start_of_variant - left_arm + donor_shift: start_of_variant]
    right_side_donor = genome_seq[chrom][end_of_variant: end_of_variant + right_arm + donor_shift]
    donor = left_side_donor + alt + right_side_donor

    middle_donor = donor[minimum_homology: -minimum_homology]
    re_in_left = _restriction_site_in_seq(left_side_donor + alt, restriction_sites)
    re_in_right = _restriction_site_in_seq(alt + right_side_donor, restriction_sites)

    # Determine if donor can be made RE-site free
    if (
        _restriction_site_in_seq(middle_donor, restriction_sites)
        or (re_in_left and re_in_right)
        or (at_right_end and re_in_left)
        or (at_left_end and re_in_right)
    ):
        return None, donor

    elif re_in_left:
        donor_shift += 1
        while (
            _restriction_site_in_seq(left_side_donor + alt + right_side_donor, restriction_sites)
            and donor_shift <= max_donor_shift
            and (end_of_variant + right_arm + donor_shift) <= chrom_length
        ):
            donor_shift += 1
            left_side_donor = genome_seq[chrom][start_of_variant - left_arm + donor_shift: start_of_variant]
            right_side_donor = genome_seq[chrom][end_of_variant: end_of_variant + right_arm + donor_shift]
        if _restriction_site_in_seq(right_side_donor, restriction_sites) or _restriction_site_in_seq(left_side_donor, restriction_sites):
            return None, donor

    elif re_in_right:
        donor_shift -= 1
        while (
            _restriction_site_in_seq(left_side_donor + alt + right_side_donor, restriction_sites)
            and abs(donor_shift) <= max_donor_shift
            and start_of_variant - left_arm >= 1
        ):
            donor_shift -= 1
            left_side_donor = genome_seq[chrom][start_of_variant - left_arm + donor_shift: start_of_variant]
            right_side_donor = genome_seq[chrom][end_of_variant: end_of_variant + right_arm + donor_shift]
        if _restriction_site_in_seq(left_side_donor, restriction_sites) or _restriction_site_in_seq(right_side_donor, restriction_sites):
            return None, donor

    donor = left_side_donor + alt + right_side_donor
    donor_start_coord = start_of_variant - left_arm + donor_shift
    donor_end_coord = end_of_variant + right_arm + donor_shift
    genomic_target = genome_seq[chrom][donor_start_coord: donor_end_coord + ref_length]

    # Determine donor strand from TIF-seq coverage
    total_plus = sum(
        plus_strand_bedgraph[chrom].get(c, 0)
        for c in range(donor_start_coord + guide_length, donor_end_coord - guide_length)
    )
    total_minus = sum(
        minus_strand_bedgraph[chrom].get(c, 0)
        for c in range(donor_start_coord + guide_length, donor_end_coord - guide_length)
    )

    if total_plus >= total_minus:
        donor_strand = "+"
    else:
        donor_strand = "-"
        donor = rev_comp(donor)

    return (
        donor_shift, donor_start_coord, donor_end_coord,
        donor_strand, left_side_donor, right_side_donor,
        ref, alt, donor, genomic_target,
    )


# ---------------------------------------------------------------------------
# PAM and guide filtering helpers
# ---------------------------------------------------------------------------

def find_pams_near_variants(
    chrom: str,
    variant_coord: int,
    ref: str,
    alt: str,
    individual_variants_lst: List[List[str]],
    coords_to_guides: _CoordsToGuides,
    guide_length: int,
    pam_length: int,
    pam_upstream: bool = True,
) -> List[_GuideInfo]:
    """
    Find all guide-PAM sites within ``guide_length * 2`` bp of a variant.

    For each candidate guide, compute the distance from the PAM to the nearest
    disrupted coordinate (i.e., a base that differs between REF and ALT).
    This disruption position tells us how far into the seed region the variant
    falls — closer to the PAM = better disruption.

    Parameters
    ----------
    chrom : str
        Chromosome name.
    variant_coord : int
        1-based coordinate of the leftmost base of the combined variant.
    ref, alt : str
        Combined REF/ALT alleles for the variant cluster.
    individual_variants_lst : list of list
        ``[[pos_str, ref, alt, variant_type], ...]`` for each constituent SNP/INDEL.
    coords_to_guides : dict
        Combined allowed + disallowed guide coordinate dict from
        :func:`load_processed_guides`.
    guide_length, pam_length : int
    pam_upstream : bool
        True for Cas9-family (PAM 3′ of guide); False for Cas12a-family.

    Returns
    -------
    list of (disruption_position, strand, guide, PAM, PAM_coord)
        Sorted by disruption_position (ascending = closer to PAM = better).
    """
    variant_start = variant_coord
    variant_end = variant_coord + len(ref)

    # Identify every genomic position disrupted by the edit
    disrupted_coords: List[int] = []
    for individual_variant in individual_variants_lst:
        vcoord_str, vref, valt, _vtype = individual_variant
        vcoord = int(vcoord_str)
        ref_len = len(vref)
        alt_len = len(valt)
        min_len = min(ref_len, alt_len)

        if alt_len != ref_len:
            disrupted_coords.append(vcoord + ref_len)

        for i in range(min_len):
            if vref[i] != valt[i]:
                disrupted_coords.append(vcoord + i)

        for i in range(1, min_len + 1):
            if vref[-i] != valt[-i]:
                disrupted_coords.append(vcoord + ref_len - i)

    result: set[_GuideInfo] = set()
    search_start = variant_start - guide_length * 2
    search_end = variant_end + guide_length * 2

    for pam_coord in range(search_start, search_end):
        for strand in ("+", "-"):
            if (chrom, pam_coord, strand) not in coords_to_guides:
                continue
            guide, pam = coords_to_guides[(chrom, pam_coord, strand)]

            disruption_position = guide_length + 100  # sentinel = far away

            if pam_upstream:
                if strand == "+":
                    distances = [(pam_coord + 1) - e for e in disrupted_coords]
                else:
                    distances = [e - pam_coord for e in disrupted_coords]
            else:
                if strand == "+":
                    distances = [e - (pam_coord + pam_length) for e in disrupted_coords]
                else:
                    distances = [(pam_coord - pam_length + 1) - e for e in disrupted_coords]

            for dist in distances:
                if -pam_length < dist < disruption_position:
                    disruption_position = dist

            if disruption_position <= guide_length:
                result.add((disruption_position, strand, guide, pam, pam_coord))

    return sorted(result)


def _restriction_site_in_seq(
    seq: str,
    restriction_sites: Tuple[str, str] = _DEFAULT_RESTRICTION_SITES,
) -> bool:
    """Return True if any restriction site appears in *seq* (case-insensitive)."""
    seq_upper = seq.upper()
    return any(site in seq_upper for site in restriction_sites)


def guide_disruption(
    guide: str,
    pam_to_check: str,
    donor: str,
    pam_upstream: bool,
    donor_length: int,
) -> bool:
    """
    Return True if the donor disrupts re-recognition of *guide* by *pam_to_check*.

    The donor is scanned for any window that perfectly matches ``guide + PAM``
    (or ``PAM + guide`` for Cas12a) on either strand.  If a perfect match is
    found, the guide is *not* disrupted (returns False).

    Parameters
    ----------
    guide : str
        20-nt guide sequence.
    pam_to_check : str
        Concrete or IUPAC PAM to test (e.g. ``"NGG"``).
    donor : str
        Donor DNA sequence.
    pam_upstream : bool
        True for Cas9 (guide upstream of PAM); False for Cas12a.
    donor_length : int
        Length of donor window to scan.

    Returns
    -------
    bool
        True if the guide is disrupted (safe to use); False if still targetable.
    """
    from magestic.design.guide_donor_functions import hamming_distance

    if pam_upstream:
        guide_with_pam = guide + pam_to_check
    else:
        guide_with_pam = pam_to_check + guide

    query_len = len(guide_with_pam)
    for idx in range(donor_length - query_len):
        window = donor[idx: idx + query_len]
        if hamming_distance(window, guide_with_pam) == 0:
            return False
        if hamming_distance(rev_comp(window), guide_with_pam) == 0:
            return False
    return True


def alternative_pam_edit(
    guide: str,
    pam_to_check: str,
    donor: str,
    pam_upstream: bool,
    donor_length: int,
) -> bool:
    """
    Return True if the donor creates an alternative PAM site for *guide*.

    Used to check for promiscuous PAMs (e.g. NAG for SpCas9, TTCC for LbCas12a)
    that would allow re-cutting even though the primary PAM is disrupted.  When
    this returns True, the guide count is decremented (the guide still works).

    Parameters
    ----------
    guide : str
    pam_to_check : str
        Alternative/promiscuous PAM to test.
    donor, pam_upstream, donor_length : see :func:`guide_disruption`.

    Returns
    -------
    bool
        True if an alternative PAM site is present in the donor.
    """
    from magestic.design.guide_donor_functions import hamming_distance

    if pam_upstream:
        guide_with_pam = guide + pam_to_check
    else:
        guide_with_pam = pam_to_check + guide

    query_len = len(guide_with_pam)
    for idx in range(donor_length - query_len):
        window = donor[idx: idx + query_len]
        if hamming_distance(window, guide_with_pam) == 0:
            return True
        if hamming_distance(rev_comp(window), guide_with_pam) == 0:
            return True
    return False


def remove_bspqi_sites_from_donor(
    fwd_primer: str,
    guide: str,
    internal_cloning_site: str,
    donor: str,
    rev_primer: str,
    restriction_sites: Tuple[str, str] = _DEFAULT_RESTRICTION_SITES,
) -> Tuple[Optional[str], list]:
    """
    Remove spurious BspQI/SapI sites at donor-junction positions.

    The fully assembled oligo is checked for RE recognition sites.  The
    cloning step requires exactly **one** such site (the intended one in the
    internal cloning linker).  If additional sites are present at the donor
    junctions, single-base substitutions at the first or last donor position
    are tried until the count reaches 1.

    Parameters
    ----------
    fwd_primer, guide, internal_cloning_site, donor, rev_primer : str
        Oligo components.
    restriction_sites : tuple of str
        The two orientations of the RE recognition sequence.

    Returns
    -------
    (new_donor, junction_changes) : tuple
        *new_donor* is ``None`` if no single-base fix was found.
        *junction_changes* is a list of change records (for the info output).
    """
    oligo_upper = (fwd_primer + guide + internal_cloning_site + donor + rev_primer).upper()
    n_sites = sum(oligo_upper.count(site) for site in restriction_sites)
    junction_changes: list = []

    if n_sites == 1:
        return donor, junction_changes

    # Try single-base substitution at each end of the donor
    for end_idx in (0, -1):
        orig_base = donor[end_idx]
        for new_base in "ATCG":
            if end_idx == 0:
                candidate = new_base + donor[1:]
            else:
                candidate = donor[:-1] + new_base
            cand_upper = (
                fwd_primer + guide + internal_cloning_site + candidate + rev_primer
            ).upper()
            n_new = sum(cand_upper.count(site) for site in restriction_sites)
            if n_new < n_sites:
                n_sites = n_new
                donor = candidate
            if n_sites == 1:
                junction_changes.append(
                    ("mutating donor index", end_idx, "from", orig_base, "to", new_base)
                )
                return donor, junction_changes

    logger.warning("Could not remove all extra RE sites from donor junction.")
    return None, None


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def design_oligos_from_vcf(
    processed_vcf_file: str | Path,
    genome_seq: dict,
    nuclease: NucleaseConfig,
    allowed_guides_file: str | Path,
    disallowed_guides_file: str | Path,
    plus_strand_bedgraph: dict,
    minus_strand_bedgraph: dict,
    rev_primer_sequences: List[str],
    output_oligo_pool_file: str | Path,
    output_info_file: str | Path,
    output_untargetable_file: str | Path,
    *,
    fasta_exclusion_files: Optional[List[str | Path]] = None,
    subpool_idx: int = 0,
    pam_preference_order: Optional[List[List[str]]] = None,
    oligo_length: int = 200,
    minimum_homology: int = 30,
    max_guides_per_donor: int = 4,
    max_dist_to_pam: int = 20,
    max_variants_per_position: int = 10,
    global_max_poly_t: int = 8,
    nuclease_label: str = "",
    progress_interval: int = 1_000,
) -> None:
    """
    Design guide-donor oligos for all variants in a processed VCF file.

    Parameters
    ----------
    processed_vcf_file : str or Path
        Tab-separated file produced by
        :func:`~vcf_processing.write_processed_vcf`.
        Columns: CHROM, POS, REF, ALT, variant_description, SAMPLES, ORF, NAME.
    genome_seq : dict
        ``{chrom: sequence}`` from :func:`~guide_donor_functions.load_genome`.
    nuclease : NucleaseConfig
        Nuclease configuration (PAM, guide length, RE sites, etc.).
    allowed_guides_file : str or Path
        Allowed guides TSV from :func:`~guide_selection.assign_allowable_guides`.
    disallowed_guides_file : str or Path
        Disallowed guides TSV from the same source.
    plus_strand_bedgraph, minus_strand_bedgraph : dict
        TIF-seq coverage (per-chromosome, per-position) from
        :func:`~guide_donor_functions.load_transcript_coverage`.
    rev_primer_sequences : list of str
        Reverse subpool priming sequences (one per subpool).
        Must have at least ``subpool_idx + len(pam_preference_order)`` entries.
    output_oligo_pool_file : str or Path
        Two-column TSV output: oligo_name, oligo_sequence.
    output_info_file : str or Path
        Detailed design metadata TSV.
    output_untargetable_file : str or Path
        TSV listing variants for which no acceptable guide was found.
    fasta_exclusion_files : list of str/Path, optional
        FASTA files to screen for reserved guide sequences.  Guides targeting
        these sequences are excluded from the design.
    subpool_idx : int
        Index into *rev_primer_sequences* for the primary subpool. Default 0.
    pam_preference_order : list of list, optional
        Nested list of concrete PAM sequences defining subpool priority.
        Guides with PAMs in ``pam_preference_order[0]`` go to subpool
        ``subpool_idx``; those in ``pam_preference_order[1]`` go to
        ``subpool_idx + 1``; etc.
        If ``None``, all PAMs are placed in a single subpool.
    oligo_length : int
        Total desired oligo length in bases (e.g. 200). Default 200.
    minimum_homology : int
        Minimum homology arm length on each side of the variant. Default 30.
    max_guides_per_donor : int
        Maximum number of guide-donor pairs to generate per variant. Default 4.
    max_dist_to_pam : int
        Maximum distance (bp) from the variant to the PAM site. Default 20.
    max_variants_per_position : int
        Process at most this many overlapping variants at each position. Default 10.
    global_max_poly_t : int
        Absolute maximum consecutive T's allowed in a guide (hard cap). Default 8.
    nuclease_label : str
        Label prepended to oligo names (e.g. ``"SpG_Cas9"``).
    progress_interval : int
        Log a progress message every N variants processed. Default 1 000.
    """
    processed_vcf_file = Path(processed_vcf_file)
    output_oligo_pool_file = Path(output_oligo_pool_file)
    output_info_file = Path(output_info_file)
    output_untargetable_file = Path(output_untargetable_file)

    for p in (output_oligo_pool_file, output_info_file, output_untargetable_file):
        p.parent.mkdir(parents=True, exist_ok=True)

    guide_length = nuclease.guide_length
    pam_length = nuclease.pam_length
    pam_upstream = nuclease.pam_upstream
    fwd_priming_seq = nuclease.fwd_priming_site
    internal_cloning_site = nuclease.internal_cloning_site
    restriction_sites = nuclease.restriction_sites
    pams_to_use = nuclease.all_pams
    max_poly_t = nuclease.max_poly_t

    donor_length = nuclease.donor_length_for_oligo(oligo_length, len(rev_primer_sequences[0]))

    if pam_preference_order is None:
        pam_preference_order = [pams_to_use]

    nuclease_label = nuclease_label or nuclease.name

    # -------------------------------------------------------------------
    # Load guide databases
    # -------------------------------------------------------------------
    logger.info("Loading guide databases...")
    coords_allowed = load_processed_guides(allowed_guides_file, pams_to_use)
    coords_disallowed = load_processed_guides(disallowed_guides_file, pams_to_use)
    coords_to_guides: _CoordsToGuides = {**coords_allowed, **coords_disallowed}

    guides_to_exclude: Set[str] = set()
    if fasta_exclusion_files:
        guides_to_exclude = load_guides_from_fasta(
            fasta_exclusion_files, guide_length, pam_length, pam_upstream, pams_to_use,
        )

    # -------------------------------------------------------------------
    # Load processed variants
    # -------------------------------------------------------------------
    variants: dict = {}  # chrom → {variant_coord: [(ref, alt, desc, desc_lst, samples, orf, name), ...]}
    with open(processed_vcf_file) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_str, ref, alt, desc, samples, orf, name = parts[:8]
            variant_coord = int(pos_str)
            desc_lst = [e.split("_") for e in desc.split(",")]
            variants.setdefault(chrom, {})
            variants[chrom].setdefault(variant_coord, []).append(
                (ref, alt, desc, desc_lst, samples, orf, name)
            )

    total_variants = sum(len(v) for chrom_v in variants.values() for v in chrom_v.values())
    logger.info("Loaded %d variants from %s", total_variants, processed_vcf_file)

    # -------------------------------------------------------------------
    # Write output files
    # -------------------------------------------------------------------
    info_header = [
        "samples", "ORF", "NAME", "chrom", "variant_coord", "ref", "alt",
        "individual_variants", "variant_type", "samples2", "disruption_position",
        "guide", "donor", "PAM_coord", "PAM_strand", "PAM",
        "donor_start_coord", "donor_end_coord", "donor_shift",
        "junction_changes", "donor_strand", "genomic_target",
        "oligo_name", "oligo_sequence",
    ]
    untargetable_header = [
        "samples", "ORF", "NAME", "chrom", "variant_coord", "ref", "alt",
        "individual_variants", "variant_type", "donor", "reason_not_targetable",
    ]

    variants_processed = 0
    guide_donors_created: Set[Tuple[str, str]] = set()
    pam_counts: dict = {}

    with (
        open(output_oligo_pool_file, "w") as pool_fh,
        open(output_info_file, "w") as info_fh,
        open(output_untargetable_file, "w") as untargetable_fh,
    ):
        pool_fh.write("oligo_name\toligo_sequence\n")
        info_fh.write("\t".join(info_header) + "\n")
        untargetable_fh.write("\t".join(untargetable_header) + "\n")

        for chrom in sorted(variants):
            chrom_length = len(genome_seq[chrom])
            chrom_variants = variants[chrom]
            logger.info(
                "Processing %s (%d bp, %d variant positions)...",
                chrom, chrom_length, len(chrom_variants),
            )

            for variant_coord in sorted(chrom_variants):
                for overlapping_variant in chrom_variants[variant_coord][:max_variants_per_position]:
                    variants_processed += 1
                    if variants_processed % progress_interval == 0:
                        logger.info("  %d variants processed...", variants_processed)

                    ref, alt, desc, desc_lst, samples, orf, name = overlapping_variant
                    variant_name = f"{chrom}_{variant_coord}_{ref}_{alt}"

                    if len(desc_lst) > 1:
                        variant_type = "LINKED"
                    elif len(ref) == len(alt) == 1:
                        variant_type = "SNP"
                    elif len(ref) == len(alt):
                        variant_type = "MNP"
                    else:
                        variant_type = "INDEL"

                    # Build donor
                    donor_result = construct_donor(
                        chrom, variant_coord, ref, alt, chrom_length,
                        desc_lst, guide_length, donor_length, minimum_homology,
                        genome_seq, plus_strand_bedgraph, minus_strand_bedgraph,
                        restriction_sites,
                    )

                    if donor_result[0] is None:
                        _, raw_donor = donor_result
                        untargetable_fh.write(
                            "\t".join([
                                samples, orf, name, chrom, str(variant_coord),
                                ref, alt, desc, variant_type, raw_donor,
                                "donor_with_restriction_site_not_allowing_shifting",
                            ]) + "\n"
                        )
                        continue

                    (
                        donor_shift, donor_start_coord, donor_end_coord,
                        donor_strand, left_side_donor, right_side_donor,
                        ref, alt, donor, genomic_target,
                    ) = donor_result

                    # Find nearby PAM sites
                    disruption_guide_info = find_pams_near_variants(
                        chrom, variant_coord, ref, alt, desc_lst,
                        coords_to_guides, guide_length, pam_length, pam_upstream,
                    )

                    num_guides = 0
                    consecutive_t_allowed = max_poly_t
                    edit_dist_allowed = max_dist_to_pam
                    reasons_not_targetable: list = []

                    while (
                        num_guides < max_guides_per_donor
                        and consecutive_t_allowed < global_max_poly_t + 1
                    ):
                        reasons_not_targetable = []
                        passing_guides: List[_GuideInfo] = []

                        for guide_info in disruption_guide_info:
                            disruption_position, pam_strand, guide, pam, pam_coord = guide_info
                            passing = True

                            # 1. Guide must disrupt all compatible PAM patterns
                            for compatible_pam in [nuclease.pam_pattern]:
                                if not guide_disruption(guide, compatible_pam, donor, pam_upstream, donor_length):
                                    passing = False
                                    reasons_not_targetable.append(
                                        (guide, compatible_pam, "variant_does_not_disrupt_PAM")
                                    )

                            # 2. Guide must not target MAGESTIC backbone elements
                            if guide in guides_to_exclude:
                                passing = False
                                reasons_not_targetable.append(guide + "_guide_in_guides_to_exclude")

                            # 3. Guide must not contain internal RE site
                            if _restriction_site_in_seq(guide, restriction_sites):
                                passing = False
                                reasons_not_targetable.append(guide + "_guide_contains_restriction_site")

                            # 4. Guide must be in allowed set (not have off-target sites)
                            if (chrom, pam_coord, pam_strand) in coords_disallowed:
                                passing = False
                                reasons_not_targetable.append(guide + "_guide_has_off_target_site")

                            # 5. Guide must not have excessive poly-T run
                            if "T" * (consecutive_t_allowed + 1) in guide:
                                passing = False
                                reasons_not_targetable.append(guide + "_guide_has_excessive_T_stretch")

                            # 6. Disruption position must be within allowed distance
                            if disruption_position > edit_dist_allowed:
                                passing = False
                                reasons_not_targetable.append(
                                    f"{disruption_position}_bp_disruption_position_beyond_acceptable_range"
                                )

                            if passing:
                                passing_guides.append(guide_info)

                        if not passing_guides:
                            reasons_not_targetable.append("no_acceptable_guide_found")
                            consecutive_t_allowed += 1
                        else:
                            # Assign guides to subpools by PAM preference
                            for subpool_offset, preferred_pams in enumerate(pam_preference_order):
                                current_subpool = subpool_idx + subpool_offset
                                rev_primer = rev_primer_sequences[current_subpool]

                                for guide_info in passing_guides:
                                    disruption_position, pam_strand, guide, pam, pam_coord = guide_info

                                    if pam not in preferred_pams:
                                        continue

                                    oligo_name = "_".join([
                                        nuclease_label, orf, name, variant_name,
                                        pam, "PAM", pam_strand, "PAM_strand",
                                        str(pam_coord), "PAM_coord",
                                        "guide_disruption_position", str(disruption_position),
                                    ])

                                    # Remove any spurious RE sites at donor junctions
                                    new_donor, junction_changes = remove_bspqi_sites_from_donor(
                                        fwd_priming_seq, guide, internal_cloning_site,
                                        donor, rev_primer, restriction_sites,
                                    )

                                    if new_donor is None:
                                        reasons_not_targetable.append(
                                            (oligo_name, "contains_restriction_site_cannot_be_removed")
                                        )
                                        disruption_guide_info = [
                                            g for g in disruption_guide_info if g != guide_info
                                        ]
                                        continue

                                    donor = new_donor
                                    guide_donor_key = (guide, donor)
                                    if guide_donor_key in guide_donors_created:
                                        logger.debug(
                                            "Duplicate guide-donor: %s %s", name, variant_name,
                                        )
                                        reasons_not_targetable.append(
                                            (guide, "guide_donor_already_created")
                                        )
                                        disruption_guide_info = [
                                            g for g in disruption_guide_info if g != guide_info
                                        ]
                                        continue

                                    oligo_sequence = (
                                        fwd_priming_seq + guide + internal_cloning_site
                                        + donor + rev_primer
                                    )
                                    guide_donors_created.add(guide_donor_key)
                                    num_guides += 1
                                    disruption_guide_info = [
                                        g for g in disruption_guide_info if g != guide_info
                                    ]

                                    if num_guides <= max_guides_per_donor:
                                        pam_counts[pam] = pam_counts.get(pam, 0) + 1
                                        info_row = [
                                            samples, orf, name, chrom, str(variant_coord),
                                            ref, alt, desc, variant_type,
                                            samples, str(disruption_position), guide, donor,
                                            str(pam_coord), pam_strand, pam,
                                            str(donor_start_coord), str(donor_end_coord),
                                            str(donor_shift), str(junction_changes),
                                            donor_strand, genomic_target,
                                            oligo_name, oligo_sequence,
                                        ]
                                        info_fh.write("\t".join(info_row) + "\n")
                                        pool_fh.write(f"{oligo_name}\t{oligo_sequence}\n")

                                    if num_guides >= max_guides_per_donor:
                                        break
                                if num_guides >= max_guides_per_donor:
                                    break
                            if num_guides >= max_guides_per_donor:
                                break

                    if num_guides == 0:
                        untargetable_fh.write(
                            "\t".join([
                                samples, orf, name, chrom, str(variant_coord),
                                ref, alt, desc, variant_type, donor,
                                str(reasons_not_targetable),
                            ]) + "\n"
                        )

    logger.info(
        "Design complete: %d variants processed, %d unique guide-donor pairs created.",
        variants_processed, len(guide_donors_created),
    )
    logger.info("PAM usage counts: %s", sorted(pam_counts.items(), key=lambda x: -x[1])[:10])
