# -*- coding: utf-8 -*-
"""
guide_donor_functions.py

Utility functions for CRISPR guide RNA design and donor DNA assembly.
Supports SpCas9 (NGG PAM, guide upstream of PAM) and Cas12a (NNNN PAM, guide downstream).

Author: Kevin R. Roy
"""

import logging
import math

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Reference data loaders
# ---------------------------------------------------------------------------

def load_codon_table(filename: str) -> tuple[dict, dict]:
    """
    Load a two-column codon table (codon<TAB>amino_acid).

    Returns:
        codon_to_aa: dict mapping codon string to single-letter amino acid
        aa_to_codon: dict mapping amino acid to list of synonymous codons
    """
    codon_to_aa: dict = {}
    aa_to_codon: dict = {}
    with open(filename) as f:
        for line in f:
            codon, aa = line.strip().split()
            codon_to_aa[codon] = aa
            aa_to_codon.setdefault(aa, []).append(codon)
    return codon_to_aa, aa_to_codon


def load_processed_guides(
    coords_to_guides_infilename: str, PAMS_TO_USE: tuple
) -> dict:
    """
    Load pre-processed guide database.

    Returns a dict keyed by (chrom, PAM_coord, strand) with values (guide, PAM).
    Only entries whose PAM is in PAMS_TO_USE are retained.
    """
    coords_to_guides: dict = {}
    with open(coords_to_guides_infilename) as infile:
        infile.readline()  # skip header
        for line in infile:
            chrom, PAM_coord, PAM_strand, guide, PAM = line.strip().split()
            if PAM in PAMS_TO_USE:
                coords_to_guides[(chrom, int(PAM_coord), PAM_strand)] = (guide, PAM)
    return coords_to_guides


def _load_sequence_from_file(filepath: str) -> str:
    """
    Load a DNA sequence from a FASTA or GenBank (.gb) file.

    Both formats are detected automatically by file extension and content.
    Returns an uppercase sequence string with no spaces or digits.
    """
    filepath = str(filepath)
    seq_parts: list = []
    if filepath.endswith(".gb") or filepath.endswith(".gbk"):
        in_origin = False
        with open(filepath, "r") as f:
            for line in f:
                if line.startswith("ORIGIN"):
                    in_origin = True
                    continue
                if in_origin:
                    if line.startswith("//"):
                        break
                    # GenBank origin lines: "        1 atgc atgc ..."
                    seq_parts.append("".join(c for c in line if c.isalpha()))
    else:
        with open(filepath, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    seq_parts.append(line.strip())
    return "".join(seq_parts).upper()


def load_guides_from_fasta(
    fasta_infilename_list: list,
    PAM_LENGTH: int,
    ALL_PAMS: tuple,
    GUIDE_LENGTH: int,
    GUIDE_UPSTREAM_PAM: bool,
) -> set:
    """
    Extract all guide sequences present in a list of sequence files (FASTA or GenBank).

    Returns a set of guide sequences to exclude from library design to avoid
    re-cutting any of the delivery vectors or integrated strain components.
    """
    guide_sequences_to_exclude: set = set()
    for fasta_infilename in fasta_infilename_list:
        seq = _load_sequence_from_file(fasta_infilename)
        extended_seq = seq + seq[:200]  # handle circular plasmids
        for idx in range(30, len(extended_seq) - 30):
            PAM_to_check = extended_seq[idx : idx + PAM_LENGTH]
            if PAM_to_check in ALL_PAMS:
                guide = (
                    extended_seq[idx - GUIDE_LENGTH : idx]
                    if GUIDE_UPSTREAM_PAM
                    else extended_seq[idx + PAM_LENGTH : idx + PAM_LENGTH + GUIDE_LENGTH]
                )
                if "N" not in guide:
                    guide_sequences_to_exclude.add(guide)
            PAM_to_check_rc = rev_comp(extended_seq[idx - PAM_LENGTH : idx])
            if PAM_to_check_rc in ALL_PAMS:
                guide = (
                    rev_comp(extended_seq[idx : idx + GUIDE_LENGTH])
                    if GUIDE_UPSTREAM_PAM
                    else rev_comp(
                        extended_seq[idx - PAM_LENGTH - GUIDE_LENGTH : idx - PAM_LENGTH]
                    )
                )
                if "N" not in guide:
                    guide_sequences_to_exclude.add(guide)
    return guide_sequences_to_exclude


def generate_empty_chromosome_dict() -> dict:
    """Return a dict with S. cerevisiae chromosome names as keys and empty dicts as values."""
    chromosome_numerals = "I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI".split(",")
    return {f"chr{num}": {} for num in chromosome_numerals}


def load_genome(genome_file: str, Mito: bool = False) -> dict:
    """
    Load genome sequences from a FASTA file.

    Returns a dict mapping chromosome names to uppercase sequence strings.
    """
    genome_dict: dict = {}
    with open(genome_file, "r") as genome:
        for line in genome:
            if line.startswith(">"):
                current_chromosome = line.strip()[1:]
                genome_dict[current_chromosome] = ""
            else:
                genome_dict[current_chromosome] += line.strip().upper()
    logger.info("Loaded genome from %s (%d sequences)", genome_file, len(genome_dict))
    return genome_dict


def load_ORF_coordinates(gff_filename: str) -> dict:
    """
    Parse CDS entries from a GFF file.

    Returns a dict mapping ORF systematic names to lists of
    (chrom, start, end, strand) exon tuples (1-based GFF coordinates).
    Mitochondrial chromosome (chrmt) entries are skipped.
    """
    ORF_dict: dict = {}
    with open(gff_filename, "r") as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            try:
                chrom, source, region, start, end, period, strand, frame, annotation = (
                    line.split("\t")
                )
                start, end = int(start), int(end)
            except ValueError:
                logger.debug("Skipping malformed GFF line: %s", line.rstrip())
                continue
            if region == "CDS" and chrom != "chrmt":
                gene_name = get_gene_name_from_gff_annotation(annotation)
                if gene_name:
                    gene_name = gene_name[:-4]
                    ORF_dict.setdefault(gene_name, []).append((chrom, start, end, strand))
    return ORF_dict


def load_systematic_to_common_gene_name_dict(gff_filename: str) -> dict:
    """
    Build a mapping from systematic ORF names to common gene names using GFF 'gene' entries.
    """
    systematic_to_common: dict = {}
    with open(gff_filename, "r") as gff_file:
        for line in gff_file:
            info = line.strip().split("\t")
            if len(info) < 9:
                continue
            chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
            if chromosome not in ("chrMito", "chrmt") and sequence_type == "gene":
                systematic_gene_name = get_systematic_gene_name(annotation)
                common_gene_name = get_common_gene_name(annotation)
                systematic_to_common[systematic_gene_name] = common_gene_name
    return systematic_to_common


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def rev_comp(DNA: str) -> str:
    """Return the reverse complement of a DNA string (supports IUPAC ambiguity codes)."""
    comp_bases = {
        "A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
        "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
        "B": "V", "D": "H", "H": "D", "V": "B",
        "a": "t", "t": "a", "g": "c", "c": "g", "n": "n",
        "r": "y", "y": "r", "s": "s", "w": "w", "k": "m", "m": "k",
        "b": "v", "d": "h", "h": "d", "v": "b",
        " ": " ", "": "",
    }
    return "".join(comp_bases.get(base, base) for base in DNA[::-1])


def base_match(base1: str, base2: str, ambiguous_bases_allowed: bool = True) -> bool:
    """Return True if two IUPAC bases are compatible."""
    compatible_bases = {
        "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
        "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
    }
    if base1 == base2:
        return True
    if ambiguous_bases_allowed:
        return (base1 in compatible_bases and base2 in compatible_bases[base1]) or (
            base2 in compatible_bases and base1 in compatible_bases[base2]
        )
    return False


def hamming_distance(
    string1: str, string2: str, ambiguous_bases_allowed: bool = True
) -> int:
    """Return the number of mismatching positions between two equal-length strings."""
    if len(string1) != len(string2):
        raise ValueError(f"Lengths differ: {len(string1)} vs {len(string2)}\n{string1}\n{string2}")
    return sum(not base_match(a, b, ambiguous_bases_allowed) for a, b in zip(string1, string2))


def get_seq_with_largest_hamming_dist(query_seq: str, subject_seq_list: list) -> str:
    """Return the sequence in subject_seq_list with the largest Hamming distance to query_seq."""
    return max(subject_seq_list, key=lambda s: hamming_distance(query_seq, s))


def get_seq_with_smallest_hamming_dist(query_seq: str, subject_seq_list: list) -> str:
    """
    Return the sequence in subject_seq_list with the smallest nonzero Hamming distance to query_seq.
    Falls back to subject_seq_list[0] if all alternatives are identical to the query.
    """
    filtered = [s for s in subject_seq_list if hamming_distance(query_seq, s) != 0]
    return (
        min(filtered, key=lambda s: hamming_distance(query_seq, s))
        if filtered
        else subject_seq_list[0]
    )


def restriction_site_in_seq(seq: str, restriction_sites: tuple) -> bool:
    """Return True if any restriction site is present in seq."""
    return any(site in seq for site in restriction_sites)


# ---------------------------------------------------------------------------
# GFF annotation parsing helpers
# ---------------------------------------------------------------------------

def get_gene_name_from_gff_annotation(string: str) -> str | None:
    """Extract the Name= field from a GFF9 annotation string, or None if absent."""
    if "Name=" not in string:
        return None
    id_start = string.find("Name=") + len("Name=")
    id_end = string.find(";", id_start)
    return string[id_start:id_end] if id_end != -1 else string[id_start:]


def get_systematic_gene_name(string: str) -> str:
    """Extract the ID= field from a GFF9 annotation string."""
    id_start = string.find("ID=") + len("ID=")
    id_end = string.find(";", id_start)
    return string[id_start:id_end] if id_end != -1 else string[id_start:]


def get_common_gene_name(string: str) -> str:
    """Extract common gene name (gene= field) or fall back to Name=."""
    if "gene=" in string:
        id_start = string.find("gene=") + len("gene=")
    else:
        id_start = string.find("Name=") + len("Name=")
    id_end = string.find(";", id_start)
    return string[id_start:id_end] if id_end != -1 else string[id_start:]


# ---------------------------------------------------------------------------
# ORF / codon coordinate mapping
# ---------------------------------------------------------------------------

def get_ORF_seq_and_CDS_coords(
    ORF: str, genome_seq: dict, ORF_dict: dict
) -> tuple[str, str, list, list]:
    """
    Extract the CDS sequence and per-base genomic coordinates for an ORF.

    Returns (chrom, strand, exon_coords, ORF_seq) where:
    - exon_coords: list of 1-based genomic coordinates, one per nucleotide of the CDS
    - ORF_seq: list of codon strings (including stop codon)
    """
    exons = ORF_dict[ORF]
    ORF_seq: list = []
    exon_coords: list = []
    strand = exons[0][3]
    frame_offset = 0
    stop_codon_indices: list = []
    stop_codons = {"TAA", "TAG", "TGA"}

    if strand == "+":
        for exon in exons:
            chrom, start, end, _ = exon
            exon_seq = genome_seq[chrom][start - 1 : end]
            if frame_offset:
                ORF_seq[-1] = ORF_seq[-1] + exon_seq[:frame_offset]
                if ORF_seq[-1] in stop_codons:
                    stop_codon_indices.append(len(ORF_seq))
                exon_coords.extend([start + i for i in range(frame_offset)])
            for idx in range(frame_offset, len(exon_seq), 3):
                frame_offset = 0
                codon = exon_seq[idx : idx + 3]
                ORF_seq.append(codon)
                if codon in stop_codons:
                    stop_codon_indices.append(len(ORF_seq))
                exon_coords.extend([start + i + idx for i in range(len(codon))])
                if len(codon) < 3:
                    frame_offset = 3 - len(codon)
    else:
        for exon in reversed(exons):
            chrom, start, end, _ = exon
            exon_seq = rev_comp(genome_seq[chrom][start - 1 : end])
            if frame_offset:
                ORF_seq[-1] = ORF_seq[-1] + exon_seq[:frame_offset]
                if ORF_seq[-1] in stop_codons:
                    stop_codon_indices.append(len(ORF_seq))
                exon_coords.extend([end - i for i in range(frame_offset)])
            for idx in range(frame_offset, len(exon_seq), 3):
                frame_offset = 0
                codon = exon_seq[idx : idx + 3]
                ORF_seq.append(codon)
                if codon in stop_codons:
                    stop_codon_indices.append(len(ORF_seq))
                exon_coords.extend([end - i - idx for i in range(len(codon))])
                if len(codon) < 3:
                    frame_offset = 3 - len(codon)

    final_codon_length = len(ORF_seq[-1])
    if final_codon_length != 3:
        logger.warning("%s: final codon is not 3 bp (%s), truncating", ORF, ORF_seq[-1])
        ORF_seq = ORF_seq[:-1]
        if final_codon_length in (1, 2):
            exon_coords = exon_coords[:-final_codon_length]

    stop_codons_found = sum(ORF_seq.count(c) for c in stop_codons)
    if stop_codons_found > 1:
        logger.warning(
            "%s: %d stop codons found at positions %s (total aa length %d, exons %s)",
            ORF, stop_codons_found, stop_codon_indices, len(ORF_seq), exons,
        )
    if stop_codons_found == 0:
        logger.warning("%s: no stop codon found in %s", ORF, exons)

    return chrom, strand, exon_coords, ORF_seq


def get_aa_num_to_codon_coords_aa(
    chrom: str,
    ORF_strand: str,
    ORF_exon_coords: list,
    ORF_seq: list,
    codon_to_aa: dict,
) -> tuple[dict, dict]:
    """
    Map 1-based amino acid numbers to (codon, genomic_coords_tuple, amino_acid).

    Returns:
        ORF_info: dict[aa_num] = (codon_str, (coord1, coord2, coord3), aa_single_letter)
        aa_num_to_exon: dict[aa_num] = list of exon numbers spanned by that codon
    """
    aa_num_to_exon: dict = {}
    aa_length = len(ORF_seq)
    previous_aa_end_coord = None
    current_exon_num = 1

    for aa_num in range(1, aa_length + 1):
        current_aa_coords = ORF_exon_coords[(aa_num - 1) * 3 : aa_num * 3]
        if previous_aa_end_coord is not None and abs(current_aa_coords[0] - previous_aa_end_coord) != 1:
            current_exon_num += 1
        aa_num_to_exon[aa_num] = [current_exon_num]
        if abs(current_aa_coords[2] - current_aa_coords[0]) != 2:
            aa_num_to_exon[aa_num] = [current_exon_num, current_exon_num + 1]
            current_exon_num += 1
        previous_aa_end_coord = current_aa_coords[2]

    ORF_info: dict = {}
    for idx, codon in enumerate(ORF_seq):
        aa_num = idx + 1
        aa = codon_to_aa[codon.upper()]
        aa_coords = tuple(ORF_exon_coords[idx * 3 : (idx + 1) * 3])
        ORF_info[aa_num] = (codon, aa_coords, aa)

    return ORF_info, aa_num_to_exon


# ---------------------------------------------------------------------------
# Guide design helpers
# ---------------------------------------------------------------------------

def assign_num_syn_changes_for_Cas9(
    leftmost_codon_coord: int,
    rightmost_codon_coord: int,
    PAM_coord: int,
    PAM_strand: str,
) -> tuple[int, int]:
    """
    Determine how many synonymous codon changes to introduce on each side of the
    target codon(s) to disrupt guide recognition (SpCas9, PAM upstream of guide).

    Returns (syn_changes_to_left, syn_changes_to_right).
    """
    syn_changes_to_left = syn_changes_to_right = 0
    if PAM_strand == "+":
        guide_start = PAM_coord - 20
        guide_midpoint = PAM_coord - 10
        cleavage_site = PAM_coord - 3
        if rightmost_codon_coord <= guide_start:
            syn_changes_to_right = (guide_midpoint - rightmost_codon_coord) // 3 + 1
        elif rightmost_codon_coord <= guide_midpoint:
            syn_changes_to_right = (cleavage_site - rightmost_codon_coord) // 3 + 1
        elif leftmost_codon_coord >= PAM_coord:
            syn_changes_to_left = (leftmost_codon_coord - cleavage_site) // 3 + 1
        else:
            syn_changes_to_left = max((leftmost_codon_coord - guide_midpoint) // 3 + 1, 0)
            syn_changes_to_right = max((PAM_coord - rightmost_codon_coord) // 3 + 1, 0)
    elif PAM_strand == "-":
        guide_start = PAM_coord + 20
        guide_midpoint = PAM_coord + 10
        cleavage_site = PAM_coord + 3
        if leftmost_codon_coord >= guide_start:
            syn_changes_to_left = (leftmost_codon_coord - guide_midpoint) // 3 + 1
        elif leftmost_codon_coord >= guide_midpoint:
            syn_changes_to_left = (leftmost_codon_coord - cleavage_site) // 3 + 1
        elif rightmost_codon_coord <= PAM_coord:
            syn_changes_to_right = (cleavage_site - rightmost_codon_coord) // 3 + 1
        else:
            syn_changes_to_right = max((guide_midpoint - rightmost_codon_coord) // 3 + 1, 0)
            syn_changes_to_left = max((leftmost_codon_coord - PAM_coord) // 3 + 1, 0)
    return syn_changes_to_left, syn_changes_to_right


def assign_num_syn_changes_for_Cas12a(
    leftmost_codon_coord: int,
    rightmost_codon_coord: int,
    PAM_coord: int,
    PAM_strand: str,
) -> tuple[int, int]:
    """
    Determine how many synonymous codon changes to introduce on each side of the
    target codon(s) to disrupt guide recognition (Cas12a, PAM downstream of guide).

    Returns (syn_changes_to_left, syn_changes_to_right).
    """
    syn_changes_to_left = syn_changes_to_right = 0
    if PAM_strand == "+":
        guide_start = PAM_coord + 4
        guide_midpoint = PAM_coord + 14
        guide_end = PAM_coord + 24
        if leftmost_codon_coord >= guide_end:
            syn_changes_to_left = (leftmost_codon_coord - guide_midpoint) // 3 + 1
        elif leftmost_codon_coord >= guide_midpoint:
            syn_changes_to_left = (leftmost_codon_coord - guide_start) // 3 + 1
        elif rightmost_codon_coord <= guide_start:
            syn_changes_to_right = (guide_start - rightmost_codon_coord) // 3 + 2
        else:
            syn_changes_to_right = max((guide_midpoint - rightmost_codon_coord) // 3 + 1, 0)
            syn_changes_to_left = max((leftmost_codon_coord - guide_start) // 3 + 2, 0)
    elif PAM_strand == "-":
        guide_start = PAM_coord - 4
        guide_midpoint = PAM_coord - 14
        guide_end = PAM_coord - 24
        if rightmost_codon_coord <= guide_end:
            syn_changes_to_right = (guide_midpoint - rightmost_codon_coord) // 3 + 1
        elif rightmost_codon_coord <= guide_midpoint:
            syn_changes_to_right = (guide_start - rightmost_codon_coord) // 3 + 1
        elif leftmost_codon_coord >= guide_start:
            syn_changes_to_left = (leftmost_codon_coord - guide_start) // 3 + 2
        else:
            syn_changes_to_left = max((leftmost_codon_coord - guide_midpoint) // 3 + 1, 0)
            syn_changes_to_right = max((guide_start - rightmost_codon_coord) // 3 + 2, 0)
    return syn_changes_to_left, syn_changes_to_right


def guide_disruption(
    guide_seq: str,
    PAM_to_check: str,
    donor: str,
    GUIDE_UPSTREAM_PAM: bool,
    PRINT: bool = False,
) -> bool:
    """
    Check whether the donor sequence disrupts re-recognition by the guide.

    Returns True if the guide+PAM motif is absent from the donor (i.e. disrupted),
    False if it is still present (oligo would be re-cut after HDR).
    """
    donor = donor.upper()
    guide_seq_with_PAM = (
        guide_seq + PAM_to_check if GUIDE_UPSTREAM_PAM else PAM_to_check + guide_seq
    )
    query_length = len(guide_seq_with_PAM)
    for idx in range(len(donor) - query_length):
        query = donor[idx : idx + query_length]
        query_rc = rev_comp(query)
        if (
            hamming_distance(query, guide_seq_with_PAM) == 0
            or hamming_distance(query_rc, guide_seq_with_PAM) == 0
        ):
            if PRINT:
                print("guide not disrupted", query, guide_seq_with_PAM, guide_seq, PAM_to_check, donor)
            return False
    if PRINT:
        print(guide_seq_with_PAM, "is disrupted by the donor:", donor)
    return True


# ---------------------------------------------------------------------------
# Donor assembly
# ---------------------------------------------------------------------------

def remove_RE_site_from_donor_with_syn_change(
    seq: str,
    donor: str,
    codon_to_aa: dict,
    suboptimal_removed_aa_to_codon: dict,
    PRINT: bool = False,
) -> tuple[str, str]:
    """
    Remove a restriction enzyme recognition site from the donor by introducing a
    synonymous codon substitution at the overlapping codon position.

    Returns (new_donor, annotated_change_string).
    """
    RE_site_seq_idx = donor.upper().index(seq)
    aa_idx = math.ceil(RE_site_seq_idx / 3)
    if PRINT:
        print(seq, donor, RE_site_seq_idx, aa_idx)
    codon = donor[aa_idx * 3 : (aa_idx + 1) * 3].upper()
    aa = codon_to_aa[codon]
    other_syn_codons = suboptimal_removed_aa_to_codon[aa]
    replacement = get_seq_with_smallest_hamming_dist(codon, other_syn_codons).lower()
    if PRINT:
        print(aa_idx, codon, aa, other_syn_codons, replacement)
    new_donor = donor[: aa_idx * 3] + replacement + donor[(aa_idx + 1) * 3 :]
    annotated_change = (
        f"restriction site:{seq} at idx in sense donor: {RE_site_seq_idx} "
        f"removed with synonymous codon: {replacement}"
    )
    if PRINT:
        print(new_donor, annotated_change)
    return new_donor, annotated_change


def remove_BspQI_sites_from_donor_in_final_oligo(
    fwd_primer: str,
    guide: str,
    guide_donor_separator: str,
    donor: str,
    subpool_rev_priming_sequence: str,
    donor_info: tuple,
) -> tuple[str | None, list]:
    """
    Resolve extra BspQI (SapI) sites that arise at the junction between the
    donor and the flanking oligo elements by mutating the terminal base(s) of the donor.

    The correctly assembled oligo contains exactly one BspQI site (in the internal
    cloning site). Any additional sites must be removed.

    Returns (donor, junction_changes) where donor is None if sites could not be resolved.
    """
    def count_BspQI(seq: str) -> int:
        upper = seq.upper()
        return upper.count("GAAGAGC") + upper.count("GCTCTTC")

    oligo_sequence = fwd_primer + guide + guide_donor_separator + donor + subpool_rev_priming_sequence
    total_sites = count_BspQI(oligo_sequence)
    junction_changes: list = []

    if total_sites == 1:
        return donor, junction_changes

    logger.warning(
        "Oligo has %d BspQI sites (expected 1); attempting junction repair.\n%s",
        total_sites, oligo_sequence,
    )

    for donor_idx in (0, -1):
        current_donor_base = donor[donor_idx]
        for other_base in "ATCG":
            new_donor = (
                other_base + donor[1:] if donor_idx == 0 else donor[:-1] + other_base
            )
            new_oligo = fwd_primer + guide + guide_donor_separator + new_donor + subpool_rev_priming_sequence
            new_sites = count_BspQI(new_oligo)
            if new_sites < total_sites:
                total_sites = new_sites
                donor = new_donor
            if total_sites == 1:
                junction_change = (
                    "mutating donor index", donor_idx,
                    "from", current_donor_base, "to", other_base,
                    "removed junction restriction site",
                )
                junction_changes.append(junction_change)
                logger.info("Junction BspQI site resolved: %s", junction_change)
                return new_donor, junction_changes

    logger.warning("Could not resolve all BspQI junction sites; donor will be skipped.")
    return None, []


def assemble_donor(
    systematic_ORF_name: str,
    ORF_info: dict,
    chrom: str,
    ORF_strand: str,
    aa_nums: list,
    ref_aas: list,
    alt_aas: list,
    ref_codons: list,
    alt_codons: list,
    syn_changes_to_left: int,
    syn_changes_to_right: int,
    donor_length: int,
    MINIMUM_HOMOLOGY: int,
    codon_to_aa: dict,
    suboptimal_removed_aa_to_codon: dict,
    genome_seq: dict,
    PRINT: bool = False,
    seqs_to_avoid: tuple | None = None,
) -> tuple:
    """
    Build a donor DNA sequence for the requested amino acid substitution(s).

    The donor consists of:
    - Left homology arm (genomic sequence)
    - Upstream synonymous codon changes (lowercase)
    - Altered codon(s) encoding the desired amino acid(s) (lowercase)
    - Downstream synonymous codon changes (lowercase)
    - Right homology arm (genomic sequence)

    Synonymous changes in the seed/PAM-proximal region are introduced to prevent
    Cas9 re-cutting after successful HDR.

    Returns a 13-tuple on success:
        (extended_donor, donor, donor_changed, donor_changes,
         aa_nums, ref_aas, alt_aas, ref_codons, alt_codons,
         num_US_syn_changes, US_syn_codons_str,
         num_DS_syn_changes, DS_syn_codons_str)

    Returns (None, reason_string) on failure.
    """
    ORF_length = len(ORF_info)
    num_US_syn_changes = syn_changes_to_left if ORF_strand == "+" else syn_changes_to_right
    num_DS_syn_changes = syn_changes_to_right if ORF_strand == "+" else syn_changes_to_left

    # Clamp upstream syn changes to not cross the start codon
    first_aa = aa_nums[0]
    if first_aa - num_US_syn_changes < 1:
        num_US_syn_changes = first_aa - 1
        logger.warning(
            "%s %s: clamped num_US_syn_changes to %d (start codon boundary)",
            systematic_ORF_name, aa_nums, num_US_syn_changes,
        )

    # Clamp upstream syn changes to not cross an exon boundary
    for aa in range(first_aa - num_US_syn_changes, first_aa):
        coord_curr = ORF_info[aa][1][0]
        coord_next = ORF_info[aa + 1][1][0]
        if abs(coord_curr - coord_next) != 3:
            num_US_syn_changes = first_aa - (aa + 1)
            logger.warning(
                "%s %s: clamped num_US_syn_changes to %d (exon boundary)",
                systematic_ORF_name, aa_nums, num_US_syn_changes,
            )

    # Clamp downstream syn changes to not cross the stop codon
    last_aa = aa_nums[-1]
    if last_aa + num_DS_syn_changes > ORF_length:
        num_DS_syn_changes = ORF_length - last_aa
        logger.warning(
            "%s %s: clamped num_DS_syn_changes to %d (stop codon boundary)",
            systematic_ORF_name, aa_nums, num_DS_syn_changes,
        )

    # Clamp downstream syn changes to not cross an exon boundary
    for aa in range(last_aa + num_DS_syn_changes, last_aa, -1):
        coord_curr = ORF_info[aa][1][2]
        coord_prev = ORF_info[aa - 1][1][2]
        if abs(coord_prev - coord_curr) != 3:
            num_DS_syn_changes = (aa - 1) - last_aa
            logger.warning(
                "%s %s: clamped num_DS_syn_changes to %d (exon boundary)",
                systematic_ORF_name, aa_nums, num_DS_syn_changes,
            )

    total_homology_length = donor_length - (num_US_syn_changes + num_DS_syn_changes + len(aa_nums)) * 3
    if total_homology_length < MINIMUM_HOMOLOGY * 2:
        return None, "donor does not meet minimum homology requirement"

    # Build upstream synonymous codons
    US_syn_start_coord = ORF_info[first_aa - num_US_syn_changes][1][0]
    US_syn_codons = []
    for aa in range(first_aa - num_US_syn_changes, first_aa):
        codon, _, current_aa = ORF_info[aa]
        US_syn_codons.append(
            get_seq_with_largest_hamming_dist(codon, suboptimal_removed_aa_to_codon[current_aa])
        )

    # Build downstream synonymous codons
    DS_syn_end_coord = ORF_info[last_aa + num_DS_syn_changes][1][2]
    DS_syn_codons = []
    for aa in range(last_aa + 1, last_aa + 1 + num_DS_syn_changes):
        codon, _, current_aa = ORF_info[aa]
        DS_syn_codons.append(
            get_seq_with_largest_hamming_dist(codon, suboptimal_removed_aa_to_codon[current_aa])
        )

    # Extract homology arms from genome
    extended_arm = 500
    if ORF_strand == "+":
        US_hom_len = round(total_homology_length / 6) * 3
        DS_hom_len = total_homology_length - US_hom_len
        left_start = US_syn_start_coord - US_hom_len - 1
        left_stop = US_syn_start_coord - 1
        left_homology = genome_seq[chrom][left_start:left_stop]
        ext_left = genome_seq[chrom][left_start - extended_arm : left_start]
        right_start = DS_syn_end_coord
        right_stop = DS_syn_end_coord + DS_hom_len
        right_homology = genome_seq[chrom][right_start:right_stop]
        ext_right = genome_seq[chrom][right_stop : right_stop + extended_arm]
    else:
        DS_hom_len = round(total_homology_length / 6) * 3
        US_hom_len = total_homology_length - DS_hom_len
        left_start = US_syn_start_coord
        left_stop = US_syn_start_coord + DS_hom_len
        left_homology = rev_comp(genome_seq[chrom][left_start:left_stop])
        ext_left = rev_comp(genome_seq[chrom][left_stop : left_stop + extended_arm])
        right_start = DS_syn_end_coord - US_hom_len - 1
        right_stop = DS_syn_end_coord - 1
        right_homology = rev_comp(genome_seq[chrom][right_start:right_stop])
        ext_right = rev_comp(genome_seq[chrom][right_start - extended_arm : right_start])

    if PRINT:
        print(
            "DS_hom_len", DS_hom_len, "US_hom_len", US_hom_len,
            "left_homology", left_homology, "US_syn_codons", US_syn_codons,
            "alt_codons", alt_codons, "DS_syn_codons", DS_syn_codons,
            "right_homology", right_homology,
        )

    donor = (
        left_homology
        + "".join(US_syn_codons).lower()
        + "".join(alt_codons).lower()
        + "".join(DS_syn_codons).lower()
        + right_homology
    )

    if len(donor) != donor_length:
        logger.warning("Donor length mismatch: expected %d, got %d", donor_length, len(donor))

    donor_changed = False
    donor_changes: list = []
    if seqs_to_avoid:
        for _ in range(5):
            for seq in seqs_to_avoid:
                if seq in donor.upper():
                    donor_changed = True
                    donor, change = remove_RE_site_from_donor_with_syn_change(
                        seq, donor, codon_to_aa, suboptimal_removed_aa_to_codon, PRINT
                    )
                    donor_changes.append(change)

    extended_donor = ext_left + donor + ext_right

    if seqs_to_avoid:
        for seq in seqs_to_avoid:
            if seq in donor.upper():
                return None, "donor could not be constructed due to unresolvable restriction sites"

    return (
        extended_donor, donor, donor_changed, donor_changes,
        aa_nums, ref_aas, alt_aas, ref_codons, alt_codons,
        num_US_syn_changes, ",".join(US_syn_codons),
        num_DS_syn_changes, ",".join(DS_syn_codons),
    )
