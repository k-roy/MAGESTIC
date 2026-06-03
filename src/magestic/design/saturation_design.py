"""
magestic.design.saturation_design — Saturation mutagenesis oligo design.

This module designs guide-donor oligo pairs that introduce every possible
amino acid substitution at every position of a target ORF (saturation
mutagenesis).  Each codon change is accompanied by synonymous codon
substitutions in flanking positions to improve guide-DNA disruption scores —
ensuring that after HDR, the nuclease cannot re-cut the edited locus.

Design logic
------------
For each amino acid position (excluding M1 and stop codon by default):

1. Starting with 0 flanking synonymous codon changes, find the best
   overlapping guide (highest disruption score, computed by
   :func:`~oligo_design.get_disruption_score`).

2. If the best score < ``min_disruption_score``, spread synonymous changes
   one codon at a time — upstream *and* downstream in separate subsets —
   until the threshold is met or the ORF boundary / exon boundary is reached.

3. Build a full donor library for the best guide: all 64 codons at the target
   position, plus ``FOLD_OVERREPRESENTATION_OF_PSEUDO_WT_CONTROLS``
   replicates of the synonymous (no-AA-change) controls.

4. Assign donors to subpools by a sliding chromosomal window (each window
   covers ``subpool_window_size`` bp of the ORF).

Output files
------------
- **oligo_pool** — tab-separated: donor_name, category, oligo_sequence,
  guide_seq, donor_seq, aa_num, oligo_length, disruption_score

Ported from:
  ``common/scripts/guide_donor_design/
  generate_amino_acid_variant_guide_donor_oligos_for_saturation_editing_v5.py``

Author: Kevin R. Roy <kevinrjroy@gmail.com>
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from magestic.design.guide_donor_functions import (
    get_ORF_seq_and_CDS_coords,
    get_aa_num_to_codon_coords_aa,
    hamming_distance,
    rev_comp,
)
from magestic.design.oligo_design import (
    get_disruption_score,
    get_sense_and_antisense_PAMs_for_ORF,
    assemble_mutated_donor_sequence,
)
from magestic.design.nuclease import NucleaseConfig

# ---------------------------------------------------------------------------
# PAM finder for saturation design (uses allowed_guides set from MAGESTIC 3.0
# guide-selection rather than the legacy Azimuth/BLAST dict)
# ---------------------------------------------------------------------------

def _find_orf_pams(
    genome_seq: dict,
    orf_chrom: str,
    orf_exon_coords: list,
    nuclease: NucleaseConfig,
    allowed_guides: set,
    bp_flanking: int = 20,
) -> tuple:
    """
    Find allowed PAM sites in and around an ORF for saturation oligo design.

    Unlike the legacy ``get_sense_and_antisense_PAMs_for_ORF`` (which queries an
    Azimuth/BLAST off-target dict), this function checks candidate guide sequences
    directly against the *allowed_guides* set produced by ``magestic-guide-selection``.

    Parameters
    ----------
    genome_seq : dict
    orf_chrom : str
    orf_exon_coords : list
        Flat list of 0-based genomic coordinates for every nucleotide of the CDS.
    nuclease : NucleaseConfig
    allowed_guides : set of str
        Set of allowed 20-nt guide sequences (uppercase) from ``load_allowed_guides()``.
    bp_flanking : int
        Extra bases to scan beyond the ORF boundaries (default 20).

    Returns
    -------
    (plus_strand_pams, minus_strand_pams) : tuple of dict
        Each dict maps PAM coordinate → (strand, pam_seq, pam_coord, guide_seq, guide_with_pam).
        Coordinates and format are identical to ``get_sense_and_antisense_PAMs_for_ORF``.
    """
    guide_length = nuclease.guide_length
    pam_length = nuclease.pam_length
    pam_upstream = nuclease.pam_upstream
    all_pams = nuclease.all_pams
    max_poly_t = nuclease.max_poly_t
    chrom_seq = genome_seq[orf_chrom]

    start_coord = min(orf_exon_coords)
    end_coord = max(orf_exon_coords)
    scan_start = start_coord - guide_length - pam_length - bp_flanking
    scan_end = end_coord + guide_length + pam_length + bp_flanking

    plus_strand_pams: dict = {}
    minus_strand_pams: dict = {}

    for coord in range(scan_start, scan_end):
        if coord < 0 or coord + pam_length > len(chrom_seq):
            continue

        # --- plus strand ---
        pam_seq = chrom_seq[coord: coord + pam_length].upper()
        if pam_seq in all_pams:
            if pam_upstream:
                guide_seq = chrom_seq[coord + pam_length: coord + pam_length + guide_length].upper()
                guide_with_pam = pam_seq + guide_seq
            else:
                guide_seq = chrom_seq[coord - guide_length: coord].upper()
                guide_with_pam = guide_seq + pam_seq
            if (guide_seq in allowed_guides
                    and max_poly_t * "T" not in guide_seq):
                plus_strand_pams[coord] = ("+", pam_seq, coord, guide_seq, guide_with_pam)

        # --- minus strand ---
        rc_pam_seq = rev_comp(pam_seq)
        if rc_pam_seq in all_pams:
            if pam_upstream:
                rc_guide_seq = rev_comp(chrom_seq[coord - guide_length: coord].upper())
                rc_guide_with_pam = rc_pam_seq + rc_guide_seq
            else:
                rc_guide_seq = rev_comp(chrom_seq[coord + pam_length: coord + pam_length + guide_length].upper())
                rc_guide_with_pam = rc_guide_seq + rc_pam_seq
            if (rc_guide_seq in allowed_guides
                    and max_poly_t * "T" not in rc_guide_seq):
                minus_strand_pams[coord] = ("-", rc_pam_seq, coord, rc_guide_seq, rc_guide_with_pam)

    return plus_strand_pams, minus_strand_pams

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default constants (may be overridden via function arguments)
# ---------------------------------------------------------------------------

_MIN_DISRUPTION_SCORE = 4
_FOLD_PSEUDO_WT_CONTROLS = 5
_SUBPOOL_WINDOW_SIZE = 250

_DEFAULT_RESTRICTION_SITES = ("GAAGAGC", "GCTCTTC")


# ---------------------------------------------------------------------------
# Synonymous codon selection
# ---------------------------------------------------------------------------

def _build_synonymous_codon_preferences(
    codon_to_aa: dict,
    aa_to_codon_fraction: dict,
) -> dict:
    """
    Build a per-codon ranking of synonymous codons (highest → lowest frequency,
    skipping the current codon then wrapping).

    Returns
    -------
    dict
        ``{codon: [preferred_codon_1, preferred_codon_2, ...]}``
        ordered by decreasing frequency, with the current codon last.
    """
    preferences: dict = {}
    for codon, aa in codon_to_aa.items():
        if aa not in aa_to_codon_fraction:
            preferences[codon] = [codon]
            continue
        ranked = [c for _, c in sorted(aa_to_codon_fraction[aa], reverse=True)]
        try:
            cur_idx = ranked.index(codon)
        except ValueError:
            preferences[codon] = ranked
            continue
        # Codon after current first, then those before (wrap around)
        preferences[codon] = ranked[cur_idx + 1:] + ranked[:cur_idx + 1]
    return preferences


def _pick_best_synonymous_codon(
    current_codon: str,
    aa: str,
    suboptimal_removed_aa_to_codon: dict,
    context_left: str,
    context_right: str,
    restriction_sites: Tuple[str, str],
) -> str:
    """
    Select the synonymous codon with the largest Hamming distance from
    *current_codon* that does not introduce a restriction site into the local
    sequence context.

    Parameters
    ----------
    current_codon : str
    aa : str
        Amino acid encoded by *current_codon*.
    suboptimal_removed_aa_to_codon : dict
        High-frequency synonymous codon options (suboptimal codons removed).
    context_left, context_right : str
        Sequence flanking the codon within the restriction-site window (7 bp).
    restriction_sites : tuple of str

    Returns
    -------
    str
        The selected synonymous codon.  Returns *current_codon* unchanged if no
        valid alternative exists.
    """
    best_codon = current_codon
    best_dist = 0
    for alt_codon in suboptimal_removed_aa_to_codon.get(aa, [current_codon]):
        dist = hamming_distance(alt_codon, current_codon)
        if dist > best_dist:
            test_seq = context_left + alt_codon + context_right
            if not any(site in test_seq.upper() for site in restriction_sites):
                best_dist = dist
                best_codon = alt_codon
    return best_codon


# ---------------------------------------------------------------------------
# Donor construction
# ---------------------------------------------------------------------------

def _mutate_synonymous_codons(
    orf_chrom: str,
    orf_strand: str,
    synonymous_codons_to_mutate: int,
    aa_num: int,
    orf_info: dict,
    orf_seq: List[str],
    genome_seq: dict,
    donor_length: int,
    minimum_homology: int,
    suboptimal_removed_aa_to_codon: dict,
    codon_to_aa: dict,
    restriction_sites: Tuple[str, str],
    upstream: bool,
) -> Tuple[str, str, int]:
    """
    Build a synonymous-change-only control donor for a given aa position.

    Returns
    -------
    (control_donor, left_donor, right_donor, ref_allele_start_coord)
    """
    aa_idx = aa_num - 1

    if orf_strand == "+":
        if upstream:
            ref_start = orf_info[aa_num - synonymous_codons_to_mutate][1][0]
        else:
            ref_start = orf_info[aa_num][1][0]
    else:
        if upstream:
            ref_start = orf_info[aa_num][1][2]
        else:
            ref_start = orf_info[aa_num + synonymous_codons_to_mutate][1][2]

    mutated_region_length = 3 * (1 + synonymous_codons_to_mutate)
    homologous_arm = donor_length - mutated_region_length
    left_arm = homologous_arm // 2
    right_arm = homologous_arm - left_arm

    left_donor = genome_seq[orf_chrom][ref_start - left_arm: ref_start]
    right_donor = genome_seq[orf_chrom][ref_start + mutated_region_length: ref_start + mutated_region_length + right_arm]

    temp_orf_seq = orf_seq[:]
    for step in range(1, synonymous_codons_to_mutate + 1):
        if upstream:
            target_aa_num = aa_num - step
        else:
            target_aa_num = aa_num + step
        target_aa_idx = target_aa_num - 1
        codon, aa_coords, aa = orf_info[target_aa_num]

        # Window for restriction site check: 3 bp context on each side
        re_window_size = max(len(r) for r in restriction_sites) - 1
        context_left = "".join(temp_orf_seq[max(0, target_aa_idx - 1): target_aa_idx])
        context_right = "".join(temp_orf_seq[target_aa_idx + 1: target_aa_idx + 2])

        best = _pick_best_synonymous_codon(
            codon, aa, suboptimal_removed_aa_to_codon,
            context_left, context_right, restriction_sites,
        )
        temp_orf_seq[target_aa_idx] = best

    if upstream:
        aa_range = (aa_num - synonymous_codons_to_mutate, aa_num)
    else:
        aa_range = (aa_num, aa_num + synonymous_codons_to_mutate)

    control_donor = assemble_mutated_donor_sequence(
        orf_strand, left_donor, right_donor, temp_orf_seq, aa_range,
        mutated_region_lowercase=True,
    )
    return control_donor, left_donor, right_donor, ref_start, temp_orf_seq


def _generate_donor_library(
    left_donor: str,
    right_donor: str,
    orf: str,
    orf_strand: str,
    synonymous_codons_to_mutate: int,
    aa_num: int,
    codon_to_aa: dict,
    orf_seq: List[str],
    temp_orf_seq: List[str],
    orf_info: dict,
    ref_start: int,
    upstream: bool,
    highest_frequency_codons: List[str],
) -> dict:
    """
    Generate the full donor library for a single amino acid position.

    Each entry in the library corresponds to one of the 64 codons at the
    target position, combined with the synonymous changes in *temp_orf_seq*.

    Returns
    -------
    dict
        ``{donor_name: (codon_pool, full_donor_seq, vcf_fields)}``
    """
    aa_idx = aa_num - 1
    wt_codon, wt_aa_coords, wt_aa = orf_info[aa_num]

    # Collect synonymous change descriptions for the name
    syn_codons: list = []
    syn_aas: list = []
    for step in range(1, synonymous_codons_to_mutate + 1):
        if upstream:
            neighbor_aa_num = aa_num - step
        else:
            neighbor_aa_num = aa_num + step
        neighbor_codon, _, neighbor_aa = orf_info[neighbor_aa_num]
        syn_codon = temp_orf_seq[neighbor_aa_num - 1]
        syn_codons.append(f"{neighbor_codon}{neighbor_aa_num}{syn_codon}")
        syn_aas.append(f"{neighbor_aa}{neighbor_aa_num}{neighbor_aa}")

    syn_codons_name = ",".join(syn_codons)
    syn_aas_name = ",".join(syn_aas)
    direction = "upstream_synonymous_changes:" if upstream else "downstream_synonymous_changes:"

    donor_library: dict = {}

    for mutant_codon in codon_to_aa:
        tmp = temp_orf_seq[:]
        tmp[aa_idx] = mutant_codon
        if upstream:
            aa_range = (aa_num - synonymous_codons_to_mutate, aa_num)
        else:
            aa_range = (aa_num, aa_num + synonymous_codons_to_mutate)

        if upstream:
            wt_region_codons = orf_seq[aa_num - synonymous_codons_to_mutate - 1: aa_idx + 1]
        else:
            wt_region_codons = orf_seq[aa_idx: aa_num + synonymous_codons_to_mutate]

        mutated_region = "".join(tmp[aa_range[0] - 1: aa_range[1]])
        wt_region = "".join(wt_region_codons)
        if orf_strand == "-":
            mutated_region = rev_comp(mutated_region)
            wt_region = rev_comp(wt_region)

        aa_change = f"{wt_aa}{aa_num}{codon_to_aa[mutant_codon]}"
        codon_change = f"{wt_codon}{aa_num}{mutant_codon}"

        if mutant_codon == wt_codon:
            donor_type = "synonymous_codon_control"
        else:
            donor_type = "synonymous_codon_variant"

        donor_name = "_".join([orf, aa_change, donor_type, codon_change, direction, syn_aas_name, syn_codons_name])
        full_donor_seq = assemble_mutated_donor_sequence(
            orf_strand, left_donor, right_donor, tmp, aa_range, mutated_region_lowercase=True,
        )
        vcf_fields = ("", ref_start, wt_region, mutated_region, donor_name)

        if donor_type == "synonymous_codon_variant":
            codon_pool = 0 if mutant_codon in highest_frequency_codons else 1
            donor_library[donor_name] = (codon_pool, full_donor_seq, vcf_fields)
        else:
            # Synonymous control: include in both subpool rows (0 and 1)
            donor_library[donor_name] = (0, full_donor_seq, vcf_fields)
            # Perfect donor (WT control) — no synonymous changes
            perfect_tmp = orf_seq[:]
            perfect_donor_seq = assemble_mutated_donor_sequence(
                orf_strand, left_donor, right_donor, perfect_tmp, aa_range, mutated_region_lowercase=True,
            )
            perfect_name = "_".join([orf, aa_change, "perfect_donor_control", codon_change, "None", "None", "None"])
            donor_library[perfect_name] = (0, perfect_donor_seq, (vcf_fields[0], ref_start, wt_region, wt_region, perfect_name))

    return donor_library


# ---------------------------------------------------------------------------
# Per-ORF guide search
# ---------------------------------------------------------------------------

def _find_best_guide_for_aa(
    aa_num: int,
    orf_info: dict,
    orf_chrom: str,
    orf_strand: str,
    orf_exon_coords: list,
    orf_seq: List[str],
    genome_seq: dict,
    donor_length: int,
    minimum_homology: int,
    suboptimal_removed_aa_to_codon: dict,
    codon_to_aa: dict,
    plus_strand_pams: dict,
    minus_strand_pams: dict,
    restriction_sites: Tuple[str, str],
    upstream: bool,
    max_synonymous_changes: int,
    min_disruption_score: int,
    pam: str,
    aa_num_to_exon: dict,
) -> Optional[Tuple]:
    """
    Find the best guide for a single amino acid position by spreading synonymous
    changes until the disruption score threshold is reached.

    Returns
    -------
    (guide_seq, pam_strand, pam_coord, guide_seq_with_pam, control_donor,
    donor_library, mutated_region_chromosomal_range, best_score)
    or None if no valid guide is found.
    """
    from magestic.design.oligo_design import get_disruption_score

    current_exon = aa_num_to_exon[aa_num][0]
    best_score = -1
    best_output = None
    synonymous_codons_to_mutate = 0

    while synonymous_codons_to_mutate <= max_synonymous_changes:
        # Boundary checks
        if upstream and aa_num - synonymous_codons_to_mutate <= 1:
            break
        if not upstream and aa_num + synonymous_codons_to_mutate > len(orf_seq):
            break

        # Check we haven't crossed an exon boundary
        if synonymous_codons_to_mutate > 0:
            if upstream:
                neighbor = aa_num - synonymous_codons_to_mutate
            else:
                neighbor = aa_num + synonymous_codons_to_mutate
            if len(aa_num_to_exon[neighbor]) > 1 or aa_num_to_exon[neighbor][0] != current_exon:
                break

        control_donor_result = _mutate_synonymous_codons(
            orf_chrom, orf_strand, synonymous_codons_to_mutate,
            aa_num, orf_info, orf_seq, genome_seq, donor_length, minimum_homology,
            suboptimal_removed_aa_to_codon, codon_to_aa, restriction_sites, upstream,
        )
        control_donor, left_donor, right_donor, ref_start, temp_orf_seq = control_donor_result

        if upstream:
            mutated_range = (ref_start, ref_start + 3 * (1 + synonymous_codons_to_mutate))
        else:
            mutated_range = (ref_start, ref_start + 3 * (1 + synonymous_codons_to_mutate))

        # Find overlapping guides
        start_coord, end_coord = mutated_range
        from magestic.design.guide_donor_functions import hamming_distance as _hd
        pam_len = len(pam.replace("N", "A"))  # approximate for search window

        for coord in range(start_coord, end_coord + 20 + 1):
            if coord in plus_strand_pams:
                pam_strand, pam_seq, pam_coord, guide_seq, guide_seq_with_pam = plus_strand_pams[coord]
                score = get_disruption_score(control_donor, guide_seq_with_pam, pam)
                if score > best_score:
                    best_score = score
                    donor_library = _generate_donor_library(
                        left_donor, right_donor, orf_chrom, orf_strand,
                        synonymous_codons_to_mutate, aa_num, codon_to_aa,
                        orf_seq, temp_orf_seq, orf_info, ref_start, upstream,
                        highest_frequency_codons=[],
                    )
                    best_output = (guide_seq, pam_strand, pam_coord, guide_seq_with_pam, control_donor, donor_library, mutated_range, best_score)

        for coord in range(start_coord - 20 - pam_len - 1, end_coord):
            if coord in minus_strand_pams:
                pam_strand, pam_seq, pam_coord, guide_seq, guide_seq_with_pam = minus_strand_pams[coord]
                score = get_disruption_score(control_donor, guide_seq_with_pam, pam)
                if score > best_score:
                    best_score = score
                    donor_library = _generate_donor_library(
                        left_donor, right_donor, orf_chrom, orf_strand,
                        synonymous_codons_to_mutate, aa_num, codon_to_aa,
                        orf_seq, temp_orf_seq, orf_info, ref_start, upstream,
                        highest_frequency_codons=[],
                    )
                    best_output = (guide_seq, pam_strand, pam_coord, guide_seq_with_pam, control_donor, donor_library, mutated_range, best_score)

        if best_score >= min_disruption_score:
            break

        synonymous_codons_to_mutate += 1

    return best_output


# ---------------------------------------------------------------------------
# Subpool assignment
# ---------------------------------------------------------------------------

def _assign_subpools(
    mutated_chromosomal_ranges: List[Tuple[int, int]],
    orf_exon_coords: list,
    subpool_window_size: int,
) -> Tuple[dict, dict]:
    """
    Assign chromosomal ranges to synthesis subpools by sliding window.

    Returns
    -------
    (chromosomal_range_to_subpool, subpool_window_to_subpool_num)
    """
    sorted_ranges = sorted(set(mutated_chromosomal_ranges))
    range_to_subpool: dict = {}
    window_to_num: dict = {}
    current_subpool = 0
    start_pos = min(orf_exon_coords)
    current_window = (start_pos, start_pos + subpool_window_size)
    window_to_num[current_window] = 0

    for chrom_range in sorted_ranges:
        if current_window[0] <= chrom_range[0] <= current_window[1] and current_window[0] <= chrom_range[1] <= current_window[1]:
            range_to_subpool[chrom_range] = current_subpool, current_window
        else:
            current_subpool += 1
            current_window = (chrom_range[0], chrom_range[0] + subpool_window_size)
            window_to_num[current_window] = current_subpool
            range_to_subpool[chrom_range] = current_subpool, current_window

    return range_to_subpool, window_to_num


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def design_saturation_oligos(
    orf_name: str,
    genome_seq: dict,
    orf_dict: dict,
    codon_to_aa: dict,
    aa_to_codon_fraction: dict,
    suboptimal_removed_aa_to_codon: dict,
    nuclease: NucleaseConfig,
    rev_primer_sequences: List[str],
    output_oligo_pool_file: str | Path,
    allowed_guides: set,
    *,
    starting_subpool_idx: int = 0,
    oligo_length: int = 230,
    minimum_homology: int = 20,
    subpool_window_size: int = _SUBPOOL_WINDOW_SIZE,
    min_disruption_score: int = _MIN_DISRUPTION_SCORE,
    fold_pseudo_wt_controls: int = _FOLD_PSEUDO_WT_CONTROLS,
    stop_codons_included: bool = True,
    mutate_initiating_methionine: bool = False,
    mutate_stop_codon: bool = True,
    restriction_sites: Tuple[str, str] = _DEFAULT_RESTRICTION_SITES,
    max_poly_t_in_last_4bp: int = 2,
) -> None:
    """
    Design saturation mutagenesis guide-donor oligos for a target ORF.

    For each amino acid position in *orf_name*, generates guide-donor pairs
    covering all possible amino acid substitutions and synonymous controls.

    Parameters
    ----------
    orf_name : str
        Systematic gene name (e.g. ``"YPR060C"``).
    genome_seq : dict
        ``{chrom: sequence_string}`` from :func:`~guide_donor_functions.load_genome`.
    orf_dict : dict
        ORF coordinate dict from :func:`~guide_donor_functions.load_ORF_coordinates`.
    codon_to_aa : dict
        Full codon→amino-acid mapping (64 entries).
    aa_to_codon_fraction : dict
        ``{aa: [(fraction, codon), ...]}`` sorted by decreasing frequency.
    suboptimal_removed_aa_to_codon : dict
        High-frequency codon options per amino acid (suboptimal codons removed).
    nuclease : NucleaseConfig
        Nuclease configuration (PAM, guide length, RE sites, etc.).
    rev_primer_sequences : list of str
        Reverse subpool priming sequences.
    output_oligo_pool_file : str or Path
        Output TSV path.
    allowed_guides : set of str
        Set of allowed guide sequences (uppercase) from
        :func:`~guide_selection.load_allowed_guides`.  Candidates not in this
        set are silently skipped.
    starting_subpool_idx : int
        Index into *rev_primer_sequences* for the first subpool. Default 0.
    oligo_length : int
        Target total oligo length. Default 230.
    minimum_homology : int
        Minimum homology arm length. Default 20.
    subpool_window_size : int
        Chromosomal window size (bp) for subpool assignment. Default 250.
    min_disruption_score : int
        Minimum guide disruption score to accept. Default 4.
    fold_pseudo_wt_controls : int
        Number of copies of synonymous controls to write. Default 5.
    stop_codons_included : bool
        Include stop codon as a target mutation. Default True.
    mutate_initiating_methionine : bool
        Include M1 in the saturation scan. Default False.
    mutate_stop_codon : bool
        Include mutations of the stop codon. Default True.
    restriction_sites : tuple of str
    max_poly_t_in_last_4bp : int
        Maximum T's allowed in the last 4 bp of the guide seed. Default 2.
    """
    from magestic.design.oligo_design import get_disruption_score
    from magestic.design.guide_donor_functions import (
        get_ORF_seq_and_CDS_coords,
        get_aa_num_to_codon_coords_aa,
    )

    output_oligo_pool_file = Path(output_oligo_pool_file)
    output_oligo_pool_file.parent.mkdir(parents=True, exist_ok=True)

    guide_length = nuclease.guide_length
    pam = nuclease.pam_pattern
    pam_length = nuclease.pam_length
    fwd_priming_seq = nuclease.fwd_priming_site
    internal_cloning_site = nuclease.internal_cloning_site
    max_poly_t = nuclease.max_poly_t

    donor_length = nuclease.donor_length_for_oligo(oligo_length, len(rev_primer_sequences[0]))
    max_synonymous_changes = (donor_length - minimum_homology * 2) // 3

    logger.info("Designing saturation oligos for %s ...", orf_name)

    # -------------------------------------------------------------------
    # Load ORF structure
    # -------------------------------------------------------------------
    orf_chrom, orf_strand, orf_exon_coords, orf_seq_str = get_ORF_seq_and_CDS_coords(
        orf_name, orf_dict, genome_seq,
    )
    orf_seq = [orf_seq_str[i: i + 3] for i in range(0, len(orf_seq_str), 3)]

    orf_info, aa_num_to_exon = get_aa_num_to_codon_coords_aa(
        orf_chrom, orf_strand, orf_exon_coords, orf_seq_str, codon_to_aa,
    )

    # Build highest-frequency codon list for subpool assignment
    highest_frequency_codons = [
        sorted(fracs, reverse=True)[0][1]
        for fracs in aa_to_codon_fraction.values()
    ]

    # -------------------------------------------------------------------
    # Find PAMs in the ORF region using MAGESTIC 3.0 allowed_guides set
    # -------------------------------------------------------------------
    logger.info("Finding PAMs for %s ...", orf_name)
    plus_strand_pams, minus_strand_pams = _find_orf_pams(
        genome_seq, orf_chrom, orf_exon_coords, nuclease, allowed_guides,
    )

    # -------------------------------------------------------------------
    # Iterate over amino acid positions
    # -------------------------------------------------------------------
    combined_guide_donor_library: dict = {}
    mutated_ranges: list = []

    start_aa = 1 if mutate_initiating_methionine else 2
    end_seq = orf_seq if mutate_stop_codon else orf_seq[:-1]

    for aa_idx_0, codon in enumerate(end_seq[start_aa - 1:]):
        aa_num = aa_idx_0 + start_aa
        codon_coords = orf_exon_coords[(aa_num - 1) * 3: aa_num * 3]

        # Skip codons split across exon boundaries
        if abs(codon_coords[2] - codon_coords[0]) != 2:
            logger.debug("Skipping split codon at aa %d", aa_num)
            continue

        if aa_num % 50 == 0:
            logger.info("  amino acid %d / %d ...", aa_num, len(orf_seq))

        for upstream in (True, False):
            result = _find_best_guide_for_aa(
                aa_num=aa_num,
                orf_info=orf_info,
                orf_chrom=orf_chrom,
                orf_strand=orf_strand,
                orf_exon_coords=orf_exon_coords,
                orf_seq=orf_seq,
                genome_seq=genome_seq,
                donor_length=donor_length,
                minimum_homology=minimum_homology,
                suboptimal_removed_aa_to_codon=suboptimal_removed_aa_to_codon,
                codon_to_aa=codon_to_aa,
                plus_strand_pams=plus_strand_pams,
                minus_strand_pams=minus_strand_pams,
                restriction_sites=restriction_sites,
                upstream=upstream,
                max_synonymous_changes=max_synonymous_changes,
                min_disruption_score=min_disruption_score,
                pam=pam,
                aa_num_to_exon=aa_num_to_exon,
            )

            if result is None:
                logger.debug("No guide found for aa %d (upstream=%s)", aa_num, upstream)
                continue

            guide_seq, pam_strand, pam_coord, guide_with_pam, control_donor, donor_library, mutated_range, best_score = result

            for donor_name, (codon_pool, full_donor_seq, vcf_fields) in donor_library.items():
                combined_guide_donor_library[donor_name] = (
                    codon_pool, mutated_range, full_donor_seq, vcf_fields,
                    guide_seq, pam_strand, pam_coord, guide_with_pam[guide_length:],
                    best_score,
                )
            mutated_ranges.append(mutated_range)

    # -------------------------------------------------------------------
    # Assign subpools
    # -------------------------------------------------------------------
    logger.info("Assigning %d donors to subpools...", len(combined_guide_donor_library))
    range_to_subpool, window_to_num = _assign_subpools(
        mutated_ranges, orf_exon_coords, subpool_window_size,
    )
    total_subpool_columns = len(window_to_num)

    # -------------------------------------------------------------------
    # Write output
    # -------------------------------------------------------------------
    header = "\t".join([
        "donor_name", "category", "oligo_sequence", "guide_seq",
        "donor_seq", "aa_num", "oligo_length", "disruption_score",
    ]) + "\n"

    donor_rows = []
    for donor_name, entry in combined_guide_donor_library.items():
        (
            subpool_row_num, mutated_range, full_donor_seq, vcf_fields,
            guide_seq, pam_strand, pam_coord, pam_seq, best_score,
        ) = entry
        disruption_score = get_disruption_score(full_donor_seq, guide_seq + pam_seq, pam)
        if "control" not in donor_name and disruption_score <= 0:
            continue

        subpool_col_num, subpool_window = range_to_subpool.get(mutated_range, (0, (0, 0)))
        subpool_primer_idx = (
            starting_subpool_idx
            + subpool_row_num * total_subpool_columns
            + subpool_col_num
        )

        oligo_sequence = (
            fwd_priming_seq.lower()
            + guide_seq.upper()
            + internal_cloning_site.lower()
            + full_donor_seq
            + rev_primer_sequences[subpool_primer_idx].lower()
        )
        category = (
            f"{subpool_primer_idx}:{orf_strand}_{orf_chrom}_"
            f"mutated_region_{mutated_range}_phenotyping_window_{subpool_window}"
        )

        try:
            aa_num = int(donor_name.split("_")[1][1:-1])
        except (IndexError, ValueError):
            aa_num = 0

        donor_rows.append((
            aa_num, donor_name, category, oligo_sequence,
            guide_seq, full_donor_seq, oligo_length, disruption_score,
        ))

    donor_rows.sort()

    with open(output_oligo_pool_file, "w") as fh:
        fh.write(header)
        for row in donor_rows:
            aa_num, donor_name, category, oligo_sequence, guide_seq, full_donor_seq, oligo_len, disruption_score = row
            n_copies = fold_pseudo_wt_controls if "control" in donor_name else 1
            line = "\t".join([
                donor_name, category, oligo_sequence, guide_seq,
                full_donor_seq, str(aa_num), str(len(oligo_sequence)), str(disruption_score),
            ]) + "\n"
            for _ in range(n_copies):
                fh.write(line)

    logger.info(
        "Saturation design complete: %d oligos written to %s",
        sum(fold_pseudo_wt_controls if "control" in r[1] else 1 for r in donor_rows),
        output_oligo_pool_file,
    )
