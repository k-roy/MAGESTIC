# -*- coding: utf-8 -*-
"""
generate_guide_donors_for_aa_mutations_Cas9_and_Cas12a.py

Designs guides and donors for desired ORF amino acid substitutions.
Input file: tab-separated, no header. Supports 3, 5, 7, 9, or 11 fields (see docstring for formats).
"""

import os
from datetime import date
from guide_donor_functions import *

DATE = date.today().strftime("%Y%m%d")

# Directories and filenames
BASE_DIR = "/Users/kevinroy/Library/CloudStorage/GoogleDrive-kevinroy@stanford.edu/My Drive/SGTC_genome_editing_group/scripts_and_reference_files/guide_donor_design/"
PROJECT_DIR = BASE_DIR + "20240209_NNS_complex_mutagenesis/"
os.chdir(PROJECT_DIR)

genome_filename = BASE_DIR + "MAGESTIC_background_strain.fasta"
gff_filename = BASE_DIR + "MAGESTIC_background_strain_annotations.gff"
codon_table_file = BASE_DIR + "DNA_codon_to_AA.txt"
codon_table_subopt_file = BASE_DIR + "DNA_codon_to_AA_suboptimal_codons_removed.txt"
rev_primer_file = BASE_DIR + "rev_subpool_priming_sequences.txt"
allowed_guides_file = (
    BASE_DIR + "MAGESTIC_base_strain_allowed_Cas9_20mer_guides_NGG.txt"
)
disallowed_guides_file = (
    BASE_DIR + "MAGESTIC_base_strain_disallowed_Cas9_20mer_guides_NGG.txt"
)

fasta_files = [
    BASE_DIR + "MARVEL_round_1_REDI_locus.fasta",
    BASE_DIR + "pS1380.fasta",
    BASE_DIR + "pS1381.fasta",
    BASE_DIR + "yT158.fasta",
    BASE_DIR + "yT159.fasta",
    BASE_DIR + "pF433.fasta",
]

# Cas9 parameters
GUIDE_UPSTREAM_PAM = True
NUCLEASE = "SpCas9"
DEGENERATE_PAM = "NGG"
PAMS_TO_USE = ("GGG", "AGG", "CGG", "TGG")
ALL_PAMS = PAMS_TO_USE
PAM_LENGTH = len(DEGENERATE_PAM)
GUIDE_LENGTH = 20
OLIGO_LENGTH = 200
MIN_HOMOLOGY = 24
MAX_GUIDE_DONORS_PER_EDIT = 2
NUM_CONSEC_T_ALLOWED = 4
RESTRICTION_SITES = ("GAAGAGC", "GCTCTTC")
FWD_PRIMING_SITE = "ctgcgattggcaggcgcgcc"
INTERNAL_CLONING_SITE = "gtttgaagagc"


# Load reference data
def load_codon_table(filename):
    codon_to_aa, aa_to_codon = {}, {}
    with open(filename) as f:
        for line in f:
            codon, aa = line.strip().split()
            codon_to_aa[codon] = aa
            aa_to_codon.setdefault(aa, []).append(codon)
    return codon_to_aa, aa_to_codon


codon_to_aa, aa_to_codon = load_codon_table(codon_table_file)
subopt_codon_to_aa, subopt_aa_to_codon = load_codon_table(codon_table_subopt_file)

with open(rev_primer_file) as f:
    REV_PRIMERS = [line.strip().lower() for line in f]
REV_PRIMING_SEQUENCES = [rev_comp(e) for e in REV_PRIMERS]

genome_seq = load_genome(genome_filename)
coords_to_guides_allowed = load_processed_guides(allowed_guides_file, PAMS_TO_USE)
coords_to_guides_disallowed = load_processed_guides(disallowed_guides_file, PAMS_TO_USE)
coords_to_guides = {**coords_to_guides_allowed, **coords_to_guides_disallowed}
guide_sequences_to_exclude = load_guides_from_fasta(
    fasta_files, PAM_LENGTH, ALL_PAMS, GUIDE_LENGTH, GUIDE_UPSTREAM_PAM
)
systematic_to_common_gene_name = load_systematic_to_common_gene_name_dict(gff_filename)
ORF_dict = load_ORF_coordinates(gff_filename)

project_names = ["NRD1_mut", "NAB3_mut", "SEN1_mut", "SSU72_mut"]
SUBPOOL_IDX = -1

for project_name in project_names:
    SUBPOOL_IDX += 1
    infile = PROJECT_DIR + project_name + ".tsv"
    file_prefix = f"{DATE}_{project_name}"

    outfilename = f"{PROJECT_DIR}{file_prefix}_{NUCLEASE}_{DEGENERATE_PAM}_codon_change_guide_donor_{OLIGO_LENGTH}_bp_oligos_info.txt"
    oligo_pool_outfilename = f"{PROJECT_DIR}{file_prefix}_{NUCLEASE}_{DEGENERATE_PAM}_codon_change_guide_donor_{OLIGO_LENGTH}_bp_oligo_pool.txt"
    untargetable_outfilename = f"{PROJECT_DIR}{file_prefix}_{NUCLEASE}_{DEGENERATE_PAM}_codon_change_guide_donor_{OLIGO_LENGTH}_bp_untargetable.txt"

    systematic_name_to_ORF_info = {}
    ORFs_to_codons_to_mutate = {}

    with open(infile) as fin:
        for line in fin:
            info = line.strip().split("\t")
            ref_codons, alt_codons = "", ""
            if len(info) == 11:
                (
                    systematic_ORF_name,
                    gene_name,
                    aa_nums,
                    ref_aas,
                    alt_aas,
                    ref_codons,
                    alt_codons,
                    *_,
                ) = info
            elif len(info) == 7:
                (
                    systematic_ORF_name,
                    gene_name,
                    aa_nums,
                    ref_aas,
                    alt_aas,
                    ref_codons,
                    alt_codons,
                ) = info
            elif len(info) == 5:
                systematic_ORF_name, gene_name, aa_nums, ref_aas, alt_aas = info
            elif len(info) == 3:
                systematic_ORF_name, gene_name, aa_substitution = info
                ref_aas = [e[0] for e in aa_substitution.split(",")]
                aa_nums = [e[1:-1] for e in aa_substitution.split(",")]
                alt_aas = [e[-1] for e in aa_substitution.split(",")]
            else:
                raise ValueError("Input file format not recognized")

            common_ORF_name = systematic_to_common_gene_name.get(
                systematic_ORF_name, systematic_ORF_name
            )
            if systematic_ORF_name not in systematic_name_to_ORF_info:
                chrom, strand, exons, seq = get_ORF_seq_and_CDS_coords(
                    systematic_ORF_name, genome_seq, ORF_dict
                )
                ORF_info, aa_num_to_exon = get_aa_num_to_codon_coords_aa(
                    chrom, strand, exons, seq, codon_to_aa
                )
                systematic_name_to_ORF_info[systematic_ORF_name] = (
                    chrom,
                    strand,
                    exons,
                    seq,
                    ORF_info,
                    aa_num_to_exon,
                )

            chrom, strand, exons, seq, ORF_info, aa_num_to_exon = (
                systematic_name_to_ORF_info[systematic_ORF_name]
            )

            if len(info) in [5, 7, 11]:
                ref_aas = ref_aas.split(",")
                aa_nums = aa_nums.split(",")
                alt_aas = alt_aas.split(",")
                ref_codons = ref_codons.split(",") if ref_codons else ""
                alt_codons = alt_codons.split(",") if alt_codons else ""
            aa_nums = [int(e) for e in aa_nums]
            aa_num_to_ref_aa = dict(zip(aa_nums, ref_aas))
            aa_num_to_alt_aa = dict(zip(aa_nums, alt_aas))
            derived_ref_codons = [ORF_info[e][0] for e in aa_nums]
            ref_codons = derived_ref_codons

            aa_num_to_ref_codon = dict(zip(aa_nums, ref_codons))
            aa_num_to_alt_codon = dict(zip(aa_nums, alt_codons)) if alt_codons else {}

            ref_aas_temp, alt_aas_temp, aa_nums_temp, alt_codons_temp = [], [], [], []
            aa_num_range = range(min(aa_nums), max(aa_nums) + 1)
            for aa_num in aa_num_range:
                ORF_derived_ref_aa = ORF_info[aa_num][2]
                ORF_derived_ref_codon = ORF_info[aa_num][0]
                aa_nums_temp.append(aa_num)
                ref_aas_temp.append(ORF_derived_ref_aa)
                alt_aas_temp.append(aa_num_to_alt_aa.get(aa_num, ORF_derived_ref_aa))
                if alt_codons and aa_num in aa_nums:
                    alt_codons_temp.append(aa_num_to_alt_codon[aa_num])
                else:
                    possible_codons = subopt_aa_to_codon[alt_aas_temp[-1]]
                    alt_codon = get_seq_with_largest_hamming_dist(
                        ORF_derived_ref_codon, possible_codons
                    )
                    alt_codons_temp.append(alt_codon)
            ref_aas, alt_aas, aa_nums, alt_codons = (
                ref_aas_temp,
                alt_aas_temp,
                aa_nums_temp,
                alt_codons_temp,
            )

            ORFs_to_codons_to_mutate.setdefault(systematic_ORF_name, []).append(
                (aa_nums, ref_aas, alt_aas, ref_codons, alt_codons)
            )

    # Output files
    with open(untargetable_outfilename, "w") as untargetable_out, open(
        oligo_pool_outfilename, "w"
    ) as oligo_pool_out, open(outfilename, "w") as out:

        untargetable_out.write(
            "\t".join(
                [
                    "systematic_ORF_name",
                    "common_ORF_name",
                    "aa_nums",
                    "ref_aas",
                    "alt_aas",
                    "ref_codons",
                    "alt_codons",
                    "ORF_strand",
                    "chrom",
                    "reasons_not_targetable",
                ]
            )
            + "\n"
        )
        oligo_pool_out.write("\t".join(["oligo_name", "oligo_seq"]) + "\n")
        out.write(
            "\t".join(
                [
                    "oligo_name",
                    "oligo_sequence",
                    "subpool_ID",
                    "subpool_number",
                    "common_ORF_name",
                    "systematic_ORF_name",
                    "aa_nums_str",
                    "ref_aas_str",
                    "alt_aas_str",
                    "ref_codons_str",
                    "alt_codons_str",
                    "num_US_syn_changes",
                    "US_syn_codons",
                    "num_DS_syn_changes",
                    "DS_syn_codons",
                    "ORF_strand",
                    "chrom",
                    "PAM_strand",
                    "PAM_coord",
                    "PAM",
                    "PAM_proximal_seed_region_to_codon_changes_midpoint",
                    "leftmost_codon_coord",
                    "rightmost_codon_coord",
                    "FWD_PRIMING_SITE",
                    "guide",
                    "INTERNAL_CLONING_SITE",
                    "donor",
                    "extended_donor",
                    "REV_PRIMING_SITE",
                    "donor_changed",
                    "donor_changes",
                    "junction_changes",
                ]
            )
            + "\n"
        )

        guide_donors_created = set()
        for systematic_ORF_name, mutations in ORFs_to_codons_to_mutate.items():
            common_ORF_name = systematic_to_common_gene_name.get(
                systematic_ORF_name, systematic_ORF_name
            )
            chrom, strand, exons, seq = get_ORF_seq_and_CDS_coords(
                systematic_ORF_name, genome_seq, ORF_dict
            )
            ORF_info, aa_num_to_exon = get_aa_num_to_codon_coords_aa(
                chrom, strand, exons, seq, codon_to_aa
            )
            for aa_nums, ref_aas, alt_aas, ref_codons, alt_codons in mutations:
                variant_targetable = False
                reasons_not_targetable = []
                first_aa, last_aa = aa_nums[0], aa_nums[-1]
                first_codon_coord = ORF_info[first_aa][1][0]
                last_codon_coord = ORF_info[last_aa][1][2]
                if (abs(first_codon_coord - last_codon_coord) + 1) // 3 != len(aa_nums):
                    reasons_not_targetable.append(
                        "aa codons to mutate span an intron, donor cannot be made"
                    )
                else:
                    leftmost, rightmost = (
                        (first_codon_coord, last_codon_coord)
                        if strand == "+"
                        else (last_codon_coord, first_codon_coord)
                    )
                    codon_mid = (first_codon_coord + last_codon_coord) // 2
                    potential_PAMs = []
                    for PAM_strand in "+-":
                        for PAM_coord in range(
                            leftmost - GUIDE_LENGTH * 6, rightmost + GUIDE_LENGTH * 6
                        ):
                            PAM_seed = (
                                PAM_coord - 3 if PAM_strand == "+" else PAM_coord + 3
                            )
                            dist_to_mid = codon_mid - PAM_seed
                            if strand == "-":
                                dist_to_mid = -dist_to_mid
                            if (chrom, PAM_coord, PAM_strand) in coords_to_guides:
                                guide, PAM = coords_to_guides[
                                    (chrom, PAM_coord, PAM_strand)
                                ]
                                if (
                                    guide in guide_sequences_to_exclude
                                    or (chrom, PAM_coord, PAM_strand)
                                    in coords_to_guides_disallowed
                                    or restriction_site_in_seq(guide, RESTRICTION_SITES)
                                    or "T" * (NUM_CONSEC_T_ALLOWED + 1) in guide
                                ):
                                    continue
                                potential_PAMs.append(
                                    (
                                        abs(dist_to_mid),
                                        dist_to_mid,
                                        chrom,
                                        PAM_coord,
                                        PAM_strand,
                                        guide,
                                        PAM,
                                    )
                                )
                    if not potential_PAMs:
                        reasons_not_targetable.append("no acceptable guides found")
                    else:
                        for (
                            _,
                            dist_to_mid,
                            chrom,
                            PAM_coord,
                            PAM_strand,
                            guide,
                            PAM,
                        ) in sorted(potential_PAMs):
                            if GUIDE_UPSTREAM_PAM:
                                syn_left, syn_right = assign_num_syn_changes_for_Cas9(
                                    leftmost, rightmost, PAM_coord, PAM_strand
                                )
                            else:
                                syn_left, syn_right = assign_num_syn_changes_for_Cas12a(
                                    leftmost, rightmost, PAM_coord, PAM_strand
                                )
                            donor_info = assemble_donor(
                                systematic_ORF_name,
                                ORF_info,
                                chrom,
                                strand,
                                aa_nums,
                                ref_aas,
                                alt_aas,
                                ref_codons,
                                alt_codons,
                                syn_left,
                                syn_right,
                                OLIGO_LENGTH
                                - len(FWD_PRIMING_SITE)
                                - GUIDE_LENGTH
                                - len(INTERNAL_CLONING_SITE)
                                - len(REV_PRIMING_SEQUENCES[SUBPOOL_IDX]),
                                MIN_HOMOLOGY,
                                codon_to_aa,
                                subopt_aa_to_codon,
                                genome_seq,
                                PRINT=False,
                                seqs_to_avoid=RESTRICTION_SITES,
                            )
                            if donor_info[0] is None:
                                reasons_not_targetable.append(donor_info)
                            else:
                                (
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
                                    US_syn_codons,
                                    num_DS_syn_changes,
                                    DS_syn_codons,
                                ) = donor_info
                                REV_PRIMING_SITE = REV_PRIMING_SEQUENCES[SUBPOOL_IDX]
                                donor, junction_changes = (
                                    remove_BspQI_sites_from_donor_in_final_oligo(
                                        FWD_PRIMING_SITE,
                                        guide,
                                        INTERNAL_CLONING_SITE,
                                        donor,
                                        REV_PRIMING_SITE,
                                        donor_info,
                                    )
                                )
                                oligo_sequence = (
                                    FWD_PRIMING_SITE
                                    + guide
                                    + INTERNAL_CLONING_SITE
                                    + donor
                                    + REV_PRIMING_SITE
                                )
                                if donor is None:
                                    reasons_not_targetable.append(
                                        "donor could not be assembled into oligo due to unresolvable restriction sites"
                                    )
                                elif not guide_disruption(
                                    guide,
                                    DEGENERATE_PAM,
                                    extended_donor,
                                    GUIDE_UPSTREAM_PAM,
                                    PRINT=False,
                                ):
                                    reasons_not_targetable.append(
                                        "guide not disrupted by donor"
                                    )
                                else:
                                    guide_donor = (guide, donor)
                                    if guide_donor in guide_donors_created:
                                        continue
                                    guide_donors_created.add(guide_donor)
                                    variant_targetable = True
                                    aa_nums_str = ",".join(map(str, aa_nums))
                                    ref_aas_str = ",".join(ref_aas)
                                    alt_aas_str = ",".join(alt_aas)
                                    ref_codons_str = ",".join(ref_codons)
                                    alt_codons_str = ",".join(alt_codons)
                                    oligo_name = "_".join(
                                        map(
                                            str,
                                            [
                                                NUCLEASE,
                                                systematic_ORF_name,
                                                common_ORF_name,
                                                aa_nums_str,
                                                ref_aas_str,
                                                alt_aas_str,
                                                ref_codons_str,
                                                alt_codons_str,
                                                num_US_syn_changes,
                                                "US_syn_changes",
                                                num_DS_syn_changes,
                                                "DS_syn_changes",
                                                strand,
                                                chrom,
                                                PAM_strand,
                                                PAM_coord,
                                                PAM,
                                                dist_to_mid,
                                                leftmost,
                                                rightmost,
                                            ],
                                        )
                                    )
                                    subpool_ID = NUCLEASE + "_codon_changes"
                                    subpool_number = f"subpool_{SUBPOOL_IDX + 1}"
                                    out.write(
                                        "\t".join(
                                            map(
                                                str,
                                                [
                                                    oligo_name,
                                                    oligo_sequence,
                                                    subpool_ID,
                                                    subpool_number,
                                                    common_ORF_name,
                                                    systematic_ORF_name,
                                                    aa_nums_str,
                                                    ref_aas_str,
                                                    alt_aas_str,
                                                    ref_codons_str,
                                                    alt_codons_str,
                                                    num_US_syn_changes,
                                                    US_syn_codons,
                                                    num_DS_syn_changes,
                                                    DS_syn_codons,
                                                    strand,
                                                    chrom,
                                                    PAM_strand,
                                                    PAM_coord,
                                                    PAM,
                                                    dist_to_mid,
                                                    leftmost,
                                                    rightmost,
                                                    FWD_PRIMING_SITE,
                                                    guide,
                                                    INTERNAL_CLONING_SITE,
                                                    donor,
                                                    extended_donor,
                                                    REV_PRIMING_SITE,
                                                    donor_changed,
                                                    donor_changes,
                                                    junction_changes,
                                                ],
                                            )
                                        )
                                        + "\n"
                                    )
                                    oligo_pool_out.write(
                                        f"{oligo_name}\t{oligo_sequence}\n"
                                    )
                                    if (
                                        len(guide_donors_created)
                                        >= MAX_GUIDE_DONORS_PER_EDIT
                                    ):
                                        break
                if not variant_targetable:
                    untargetable_out.write(
                        "\t".join(
                            map(
                                str,
                                [
                                    systematic_ORF_name,
                                    common_ORF_name,
                                    aa_nums,
                                    ref_aas,
                                    alt_aas,
                                    ref_codons,
                                    alt_codons,
                                    strand,
                                    chrom,
                                    reasons_not_targetable,
                                ],
                            )
                        )
                        + "\n"
                    )
