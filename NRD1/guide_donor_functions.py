import math


def load_processed_guides(coords_to_guides_infilename, PAMS_TO_USE):
    """
    Loads pre-processed guides from a file.
    Returns a dictionary with (chrom, PAM_coord, strand) as keys and (guide, PAM) as values.
    """
    coords_to_guides = {}
    with open(coords_to_guides_infilename) as infile:
        infile.readline()  # Skip header
        for line in infile:
            chrom, PAM_coord, PAM_strand, guide, PAM = line.strip().split()
            if PAM in PAMS_TO_USE:
                coords_to_guides[(chrom, int(PAM_coord), PAM_strand)] = (guide, PAM)
    return coords_to_guides


def load_guides_from_fasta(
    fasta_infilename_list, PAM_LENGTH, ALL_PAMS, GUIDE_LENGTH, GUIDE_UPSTREAM_PAM
):
    """
    Loads guide sequences from fasta files, checking for all possible guides with given PAMs.
    Returns a set of guide sequences to exclude from design.
    """
    guide_sequences_to_exclude = set()
    for fasta_infilename in fasta_infilename_list:
        seq = ""
        with open(fasta_infilename, "r") as fasta:
            for line in fasta:
                if not line.startswith(">"):
                    seq += line.strip().upper()
        extended_seq = seq + seq[:200]  # Account for circular plasmids
        for idx in range(30, len(extended_seq) - 30):
            PAM_to_check = extended_seq[idx : idx + PAM_LENGTH]
            if PAM_to_check in ALL_PAMS:
                guide = (
                    extended_seq[idx - GUIDE_LENGTH : idx]
                    if GUIDE_UPSTREAM_PAM
                    else extended_seq[
                        idx + PAM_LENGTH : idx + PAM_LENGTH + GUIDE_LENGTH
                    ]
                )
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
                guide_sequences_to_exclude.add(guide)
    return guide_sequences_to_exclude


def generate_empty_chromosome_dict():
    """
    Returns a dictionary with chromosome names as keys and empty dicts as values.
    """
    chromosome_numerals = "I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI".split(
        ","
    )
    return {f"chr{num}": {} for num in chromosome_numerals}


def load_genome(genome_file, Mito=False):
    """
    Loads genome sequence from fasta file.
    Returns a dictionary with chromosome names as keys and sequences as values.
    """
    genome_dict = {}
    with open(genome_file, "r") as genome:
        for line in genome:
            if line.startswith(">"):
                current_chromosome = line.strip()[1:]
                genome_dict[current_chromosome] = ""
            else:
                genome_dict[current_chromosome] += line.strip().upper()
    print(f"{genome_file} genome file loaded")
    return genome_dict


def rev_comp(DNA):
    """
    Returns the reverse complement of a DNA string.
    """
    comp_bases = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "n": "n",
        "r": "y",
        "y": "r",
        "s": "s",
        "w": "w",
        "k": "m",
        "m": "k",
        "b": "v",
        "d": "h",
        "h": "d",
        "v": "b",
        " ": " ",
        "": "",
    }
    return "".join(comp_bases.get(base, base) for base in DNA[::-1])


def base_match(base1, base2, ambiguous_bases_allowed=True):
    """
    Returns True if two IUPAC bases are compatible, False otherwise.
    """
    compatible_bases = {
        "R": "AG",
        "Y": "CT",
        "S": "GC",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "ACGT",
    }
    if base1 == base2:
        return True
    if ambiguous_bases_allowed:
        return (base1 in compatible_bases and base2 in compatible_bases[base1]) or (
            base2 in compatible_bases and base1 in compatible_bases[base2]
        )
    return False


def hamming_distance(string1, string2, ambiguous_bases_allowed=True):
    """
    Returns the number of mismatches between two strings of equal length.
    """
    if len(string1) != len(string2):
        raise ValueError(f"Lengths of strings not equal! \n{string1}\n{string2}")
    return sum(
        not base_match(a, b, ambiguous_bases_allowed) for a, b in zip(string1, string2)
    )


def get_gene_name_from_gff_annotation(string):
    """
    Extracts gene name from GFF annotation string.
    """
    if "Name=" not in string:
        return None
    id_start = string.find("Name=") + len("Name=")
    id_end = string.find(";", id_start)
    return string[id_start:id_end] if id_end != -1 else string[id_start:]


def load_ORF_coordinates(gff_filename):
    """
    Loads ORF coordinates from GFF file.
    Returns a dictionary with ORF names as keys and lists of exon tuples as values.
    """
    ORF_dict = {}
    with open(gff_filename, "r") as gff:
        for line in gff:
            if not line.startswith("#"):
                try:
                    (
                        chrom,
                        source,
                        region,
                        start,
                        end,
                        period,
                        strand,
                        frame,
                        annotation,
                    ) = line.split("\t")
                    start, end = int(start), int(end)
                except Exception:
                    print(line)
                    continue
                if region == "CDS" and chrom != "chrmt":
                    gene_name = get_gene_name_from_gff_annotation(annotation)
                    if gene_name:
                        gene_name = gene_name[:-4]
                        ORF_dict.setdefault(gene_name, []).append(
                            (chrom, start, end, strand)
                        )
    return ORF_dict


def get_ORF_seq_and_CDS_coords(ORF, genome_seq, ORF_dict):
    """
    Returns chromosome, strand, exon coordinates, and ORF sequence for a given ORF.
    """
    exons = ORF_dict[ORF]
    ORF_seq, exon_coords = [], []
    strand = exons[0][3]
    frame_offset, stop_codon_indices = 0, []
    stop_codons = {"TAA", "TAG", "TGA"}
    if strand == "+":
        for exon in exons:
            chrom, start, end, _ = exon
            exon_seq = genome_seq[chrom][start - 1 : end]
            if frame_offset:
                previous_partial_codon = ORF_seq[-1]
                rest_of_codon = exon_seq[:frame_offset]
                ORF_seq[-1] = previous_partial_codon + rest_of_codon
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
                previous_partial_codon = ORF_seq[-1]
                rest_of_codon = exon_seq[:frame_offset]
                ORF_seq[-1] = previous_partial_codon + rest_of_codon
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
        print(f"{ORF} has ORF not multiple of 3, with final codon {ORF_seq[-1]}")
        ORF_seq = ORF_seq[:-1]
        if final_codon_length in (1, 2):
            exon_coords = exon_coords[:-final_codon_length]
    stop_codons_found = sum(ORF_seq.count(codon) for codon in stop_codons)
    if stop_codons_found > 1:
        print(
            f"{ORF} has {stop_codons_found} stop codons at positions {stop_codon_indices} in the ORF with total aa length: {len(ORF_seq)} and exons {exons}"
        )
    if stop_codons_found == 0:
        print(f"{ORF} has no stop codon {exons}")
    return chrom, strand, exon_coords, ORF_seq


def get_aa_num_to_codon_coords_aa(
    chrom, ORF_strand, ORF_exon_coords, ORF_seq, codon_to_aa
):
    """
    Maps amino acid numbers to codon, coordinates, and amino acid.
    Returns ORF_info and aa_num_to_exon dictionaries.
    """
    current_exon_num = 1
    aa_num_to_exon = {}
    aa_length = len(ORF_seq)
    previous_aa_end_coord = None
    for aa_num in range(1, aa_length + 1):
        current_aa_coords = ORF_exon_coords[(aa_num - 1) * 3 : aa_num * 3]
        if (
            previous_aa_end_coord is not None
            and abs(current_aa_coords[0] - previous_aa_end_coord) != 1
        ):
            current_exon_num += 1
        aa_num_to_exon[aa_num] = [current_exon_num]
        if abs(current_aa_coords[2] - current_aa_coords[0]) != 2:
            aa_num_to_exon[aa_num] = [current_exon_num, current_exon_num + 1]
            current_exon_num += 1
        previous_aa_end_coord = current_aa_coords[2]
    ORF_info = {}
    for idx, codon in enumerate(ORF_seq):
        aa_num = idx + 1
        aa = codon_to_aa[codon.upper()]
        aa_coords = tuple(ORF_exon_coords[idx * 3 : (idx + 1) * 3])
        ORF_info[aa_num] = (codon, aa_coords, aa)
    return ORF_info, aa_num_to_exon


def get_systematic_gene_name(string):
    """
    Extracts systematic gene name from GFF annotation string.
    """
    id_start = string.find("ID=") + len("ID=")
    id_end = string.find(";", id_start)
    return string[id_start:id_end] if id_end != -1 else string[id_start:]


def get_common_gene_name(string):
    """
    Extracts common gene name from GFF annotation string, or systematic if not present.
    """
    if "gene=" not in string:
        id_start = string.find("Name=") + len("Name=")
        id_end = string.find(";", id_start)
        return string[id_start:id_end] if id_end != -1 else string[id_start:]
    id_start = string.find("gene=") + len("gene=")
    id_end = string.find(";", id_start)
    return string[id_start:id_end] if id_end != -1 else string[id_start:]


def load_systematic_to_common_gene_name_dict(gff_filename):
    """
    Loads mapping from systematic to common gene names from GFF file.
    """
    systematic_to_common_gene_name_dict = {}
    with open(gff_filename, "r") as gff_file:
        for line in gff_file:
            info = line.strip().split("\t")
            if len(info) < 9:
                continue
            (
                chromosome,
                region,
                sequence_type,
                start,
                end,
                period,
                strand,
                frame,
                annotation,
            ) = info
            if chromosome not in ("chrMito", "chrmt") and sequence_type == "gene":
                systematic_gene_name = get_systematic_gene_name(annotation)
                common_gene_name = get_common_gene_name(annotation)
                systematic_to_common_gene_name_dict[systematic_gene_name] = (
                    common_gene_name
                )
    return systematic_to_common_gene_name_dict


def assign_num_syn_changes_for_Cas9(
    leftmost_codon_coord, rightmost_codon_coord, PAM_coord, PAM_strand
):
    """
    Returns the number of synonymous changes to left and right of the target range for Cas9.
    """
    syn_changes_to_right = syn_changes_to_left = 0
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
            syn_changes_to_left = max(
                (leftmost_codon_coord - guide_midpoint) // 3 + 1, 0
            )
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
            syn_changes_to_right = max(
                (guide_midpoint - rightmost_codon_coord) // 3 + 1, 0
            )
            syn_changes_to_left = max((leftmost_codon_coord - PAM_coord) // 3 + 1, 0)
    return syn_changes_to_left, syn_changes_to_right


def assign_num_syn_changes_for_Cas12a(
    leftmost_codon_coord, rightmost_codon_coord, PAM_coord, PAM_strand
):
    """
    Returns the number of synonymous changes to left and right of the target range for Cas12a.
    """
    syn_changes_to_right = syn_changes_to_left = 0
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
            syn_changes_to_right = max(
                (guide_midpoint - rightmost_codon_coord) // 3 + 1, 0
            )
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
            syn_changes_to_left = max(
                (leftmost_codon_coord - guide_midpoint) // 3 + 1, 0
            )
            syn_changes_to_right = max(
                (guide_start - rightmost_codon_coord) // 3 + 2, 0
            )
    return syn_changes_to_left, syn_changes_to_right


def get_seq_with_largest_hamming_dist(query_seq, subject_seq_list):
    """
    Returns the sequence from subject_seq_list with the largest Hamming distance to query_seq.
    """
    return max(
        subject_seq_list, key=lambda subject: hamming_distance(query_seq, subject)
    )


def get_seq_with_smallest_hamming_dist(query_seq, subject_seq_list):
    """
    Returns the sequence from subject_seq_list with the smallest nonzero Hamming distance to query_seq.
    """
    filtered = [
        subject
        for subject in subject_seq_list
        if hamming_distance(query_seq, subject) != 0
    ]
    return (
        min(filtered, key=lambda subject: hamming_distance(query_seq, subject))
        if filtered
        else subject_seq_list[0]
    )


def remove_RE_site_from_donor_with_syn_change(
    seq, donor, codon_to_aa, suboptimal_removed_aa_to_codon, PRINT=False
):
    """
    Removes a restriction enzyme site from donor by introducing a synonymous codon change.
    """
    RE_site_seq_idx = donor.upper().index(seq)
    aa_idx = math.ceil(RE_site_seq_idx / 3)
    syn_codons_to_add = ""
    if PRINT:
        print(seq, donor, RE_site_seq_idx, aa_idx)
    codon = donor[aa_idx * 3 : (aa_idx + 1) * 3].upper()
    aa = codon_to_aa[codon]
    other_syn_codons = suboptimal_removed_aa_to_codon[aa]
    seq_with_smallest_hamming_dist = get_seq_with_smallest_hamming_dist(
        codon, other_syn_codons
    ).lower()
    syn_codons_to_add += seq_with_smallest_hamming_dist
    if PRINT:
        print(
            aa_idx,
            codon,
            aa,
            other_syn_codons,
            seq_with_smallest_hamming_dist,
            syn_codons_to_add,
        )
    new_donor = donor[: aa_idx * 3] + syn_codons_to_add + donor[(aa_idx + 1) * 3 :]
    annotated_change = f"restriction site:{seq} at idx in sense donor: {RE_site_seq_idx} removed with synonymous codons: {syn_codons_to_add}"
    if PRINT:
        print(new_donor, annotated_change)
    return new_donor, annotated_change


def remove_BspQI_sites_from_donor_in_final_oligo(
    fwd_primer,
    guide,
    guide_donor_separator,
    donor,
    subpool_rev_priming_sequence,
    donor_info,
):
    """
    Removes BspQI restriction sites from donor sequence in final oligo.
    """
    oligo_sequence = (
        fwd_primer
        + guide
        + guide_donor_separator
        + donor
        + subpool_rev_priming_sequence
    )
    oligo_sequence_upper = oligo_sequence.upper()
    total_BspQI_sites = oligo_sequence_upper.count(
        "GAAGAGC"
    ) + oligo_sequence_upper.count("GCTCTTC")
    junction_changes = []
    if total_BspQI_sites == 1:
        return donor, junction_changes
    print(
        "oligo has BspQI site at junction between donor and cloning site or donor and subpool priming site"
    )
    print(donor_info)
    print(oligo_sequence)
    print(oligo_sequence_upper)
    for donor_idx in (0, -1):
        current_donor_base = donor[donor_idx]
        for other_base in "ATCG":
            new_donor = (
                other_base + donor[1:] if donor_idx == 0 else donor[:-1] + other_base
            )
            oligo_sequence = (
                fwd_primer
                + guide
                + guide_donor_separator
                + new_donor
                + subpool_rev_priming_sequence
            )
            oligo_sequence_upper = oligo_sequence.upper()
            modified_donor_total_BspQI_sites = oligo_sequence_upper.count(
                "GAAGAGC"
            ) + oligo_sequence_upper.count("GCTCTTC")
            if modified_donor_total_BspQI_sites < total_BspQI_sites:
                total_BspQI_sites = modified_donor_total_BspQI_sites
                donor = new_donor
            if total_BspQI_sites == 1:
                junction_change = (
                    "mutating donor index",
                    donor_idx,
                    "from",
                    current_donor_base,
                    "to",
                    other_base,
                    "removed junction restriction site",
                )
                junction_changes.append(junction_change)
                print(junction_change)
                return new_donor, junction_changes
    return None, None


def assemble_donor(
    systematic_ORF_name,
    ORF_info,
    chrom,
    ORF_strand,
    aa_nums,
    ref_aas,
    alt_aas,
    ref_codons,
    alt_codons,
    syn_changes_to_left,
    syn_changes_to_right,
    donor_length,
    MINIMUM_HOMOLOGY,
    codon_to_aa,
    suboptimal_removed_aa_to_codon,
    genome_seq,
    PRINT=False,
    seqs_to_avoid=None,
):
    """
    Assembles donor DNA sequence for mutagenesis, introducing synonymous changes and checking for restriction sites.
    Returns donor info tuple.
    """
    ORF_length = len(ORF_info)
    num_US_syn_changes = (
        syn_changes_to_left if ORF_strand == "+" else syn_changes_to_right
    )
    num_DS_syn_changes = (
        syn_changes_to_right if ORF_strand == "+" else syn_changes_to_left
    )
    first_aa_num_to_mutate = aa_nums[0]
    if first_aa_num_to_mutate - num_US_syn_changes < 1:
        num_US_syn_changes = first_aa_num_to_mutate - 1
        print(
            systematic_ORF_name,
            aa_nums,
            "US synonymous changes require traversing the start codon, restricted num_US_syn_changes to:",
            num_US_syn_changes,
        )
    for current_aa_num in range(
        first_aa_num_to_mutate - num_US_syn_changes, first_aa_num_to_mutate
    ):
        current_codon, current_aa_coords, _ = ORF_info[current_aa_num]
        next_aa_num = current_aa_num + 1
        next_codon, next_aa_coords, _ = ORF_info[next_aa_num]
        if abs(current_aa_coords[0] - next_aa_coords[0]) != 3:
            num_US_syn_changes = first_aa_num_to_mutate - next_aa_num
            print(
                systematic_ORF_name,
                aa_nums,
                "US synonymous changes require traversing an exon, restricted num_US_syn_changes to:",
                num_US_syn_changes,
            )
    last_aa_num_to_mutate = aa_nums[-1]
    if last_aa_num_to_mutate + num_DS_syn_changes > ORF_length:
        num_DS_syn_changes = ORF_length - last_aa_num_to_mutate
        print(
            systematic_ORF_name,
            aa_nums,
            "DS synonymous changes require traversing the stop codon, restricted num_DS_syn_changes to:",
            num_DS_syn_changes,
        )
    for current_aa_num in range(
        last_aa_num_to_mutate + num_DS_syn_changes, last_aa_num_to_mutate, -1
    ):
        current_codon, current_aa_coords, _ = ORF_info[current_aa_num]
        previous_aa_num = current_aa_num - 1
        previous_codon, previous_aa_coords, _ = ORF_info[previous_aa_num]
        if abs(previous_aa_coords[2] - current_aa_coords[2]) != 3:
            num_DS_syn_changes = previous_aa_num - last_aa_num_to_mutate
            print(
                systematic_ORF_name,
                aa_nums,
                "DS synonymous changes require traversing an exon, restricted num_DS_syn_changes to:",
                num_DS_syn_changes,
            )
    total_homology_length = (
        donor_length - (num_US_syn_changes + num_DS_syn_changes + len(aa_nums)) * 3
    )
    if total_homology_length < MINIMUM_HOMOLOGY * 2:
        return None, "donor does not meet minimum homology requirement"
    US_syn_codons = []
    US_syn_start_coord = ORF_info[first_aa_num_to_mutate - num_US_syn_changes][1][0]
    for current_aa_num in range(
        first_aa_num_to_mutate - num_US_syn_changes, first_aa_num_to_mutate
    ):
        current_codon, _, current_aa = ORF_info[current_aa_num]
        synonymous_codons = suboptimal_removed_aa_to_codon[current_aa]
        US_syn_codons.append(
            get_seq_with_largest_hamming_dist(current_codon, synonymous_codons)
        )
    DS_syn_codons = []
    DS_syn_end_coord = ORF_info[last_aa_num_to_mutate + num_DS_syn_changes][1][2]
    for current_aa_num in range(
        last_aa_num_to_mutate + 1, last_aa_num_to_mutate + 1 + num_DS_syn_changes
    ):
        current_codon, _, current_aa = ORF_info[current_aa_num]
        synonymous_codons = suboptimal_removed_aa_to_codon[current_aa]
        DS_syn_codons.append(
            get_seq_with_largest_hamming_dist(current_codon, synonymous_codons)
        )
    extended_donor_arm_length = 500
    if ORF_strand == "+":
        US_homology_length = round(total_homology_length / 6) * 3
        DS_homology_length = total_homology_length - US_homology_length
        left_donor_start_coord = US_syn_start_coord - US_homology_length - 1
        left_donor_stop_coord = US_syn_start_coord - 1
        left_donor_homology = genome_seq[chrom][
            left_donor_start_coord:left_donor_stop_coord
        ]
        extended_left_donor_homology = genome_seq[chrom][
            left_donor_start_coord - extended_donor_arm_length : left_donor_start_coord
        ]
        right_donor_start_coord = DS_syn_end_coord
        right_donor_stop_coord = DS_syn_end_coord + DS_homology_length
        right_donor_homology = genome_seq[chrom][
            right_donor_start_coord:right_donor_stop_coord
        ]
        extended_right_donor_homology = genome_seq[chrom][
            right_donor_stop_coord : right_donor_stop_coord + extended_donor_arm_length
        ]
    else:
        DS_homology_length = round(total_homology_length / 6) * 3
        US_homology_length = total_homology_length - DS_homology_length
        left_donor_start_coord = US_syn_start_coord
        left_donor_stop_coord = US_syn_start_coord + DS_homology_length
        left_donor_homology = rev_comp(
            genome_seq[chrom][left_donor_start_coord:left_donor_stop_coord]
        )
        extended_left_donor_homology = rev_comp(
            genome_seq[chrom][
                left_donor_stop_coord : left_donor_stop_coord
                + extended_donor_arm_length
            ]
        )
        right_donor_start_coord = DS_syn_end_coord - US_homology_length - 1
        right_donor_stop_coord = DS_syn_end_coord - 1
        right_donor_homology = rev_comp(
            genome_seq[chrom][right_donor_start_coord:right_donor_stop_coord]
        )
        extended_right_donor_homology = rev_comp(
            genome_seq[chrom][
                right_donor_start_coord
                - extended_donor_arm_length : right_donor_start_coord
            ]
        )
    if PRINT:
        print(
            "DS_homology_length",
            DS_homology_length,
            "US_homology_length",
            US_homology_length,
            "left_donor_homology",
            left_donor_homology,
            "US_syn_codons",
            US_syn_codons,
            "alt_codons",
            alt_codons,
            "DS_syn_codons",
            DS_syn_codons,
            "right_donor_homology",
            right_donor_homology,
        )
    donor = (
        left_donor_homology
        + "".join(US_syn_codons).lower()
        + "".join(alt_codons).lower()
        + "".join(DS_syn_codons).lower()
        + right_donor_homology
    )
    if len(donor) != donor_length:
        print("donor is the wrong length", donor, len(donor), donor_length)
    donor_changed = False
    donor_changes = []
    if seqs_to_avoid:
        for _ in range(5):
            for seq in seqs_to_avoid:
                if seq in donor.upper():
                    donor_changed = True
                    donor, annotated_change = remove_RE_site_from_donor_with_syn_change(
                        seq, donor, codon_to_aa, suboptimal_removed_aa_to_codon, PRINT
                    )
                    donor_changes.append(annotated_change)
    extended_donor = (
        extended_left_donor_homology + donor + extended_right_donor_homology
    )
    if seqs_to_avoid:
        for seq in seqs_to_avoid:
            if seq in donor.upper():
                return (
                    None,
                    "donor could not be constructed due to unresolvable restriction sites",
                )
    US_syn_codons_str = ",".join(US_syn_codons)
    DS_syn_codons_str = ",".join(DS_syn_codons)
    donor_info = (
        extended_donor,
        donor,
        donor_changed,
        donor_changes,
        aa_nums,
        ref_aas,
        alt_aas,
        ref_codons,
        alt_codons,
        num_US_syn_changes,
        US_syn_codons_str,
        num_DS_syn_changes,
        DS_syn_codons_str,
    )
    return donor_info


def guide_disruption(guide_seq, PAM_to_check, donor, GUIDE_UPSTREAM_PAM, PRINT=False):
    """
    Checks if guide recognition is disrupted by the donor DNA.
    Returns True if disrupted, False otherwise.
    """
    donor = donor.upper()
    guide_seq_with_PAM = (
        guide_seq + PAM_to_check if GUIDE_UPSTREAM_PAM else PAM_to_check + guide_seq
    )
    query_length = len(guide_seq_with_PAM)
    for idx in range(0, len(donor) - query_length):
        query = donor[idx : idx + query_length]
        query_rc = rev_comp(query)
        if (
            hamming_distance(query, guide_seq_with_PAM) == 0
            or hamming_distance(query_rc, guide_seq_with_PAM) == 0
        ):
            if PRINT:
                print(
                    "guide not disrupted",
                    query,
                    guide_seq_with_PAM,
                    guide_seq,
                    PAM_to_check,
                    donor,
                )
            return False
    if PRINT:
        print(guide_seq_with_PAM, "is disrupted by the donor:", donor)
    return True


def restriction_site_in_seq(seq, RESTRICTION_SITES):
    """
    Checks if any restriction site is present in the sequence.
    """
    return any(site in seq for site in RESTRICTION_SITES)
