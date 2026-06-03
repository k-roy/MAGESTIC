# -*- coding: utf-8 -*-
"""
generate_guide_donors_for_aa_substitutions_Cas9.py

Design SpCas9 guide RNA / donor DNA oligonucleotides for introducing specific
amino acid substitutions into a yeast ORF via CRISPR-mediated HDR.

This script was used to design the NRD1 alanine-scan / PTC / synonymous-control
oligonucleotide library described in:

    Roy KR et al. (2024) [manuscript in preparation]

Usage
-----
    python -m magestic.projects.nrd1.generate_guide_donors_for_aa_substitutions_Cas9 \\
        --ref-dir /path/to/large/reference/files/ \\
        --rev-primer-file /path/to/NRD1/rev_subpool_priming_sequences.txt \\
        --sequence-map-dir /path/to/NRD1/NRD1_MAGESTIC_3.0_sequence_maps/

All small reference files (codon tables, NRD1 mutation list) are bundled in the
magestic.projects.nrd1.data package and accessed via get_nrd1_data_path().

The three large reference files must be downloaded from Zenodo (see README.md):
    - MAGESTIC_background_strain.fasta
    - MAGESTIC_background_strain_annotations.gff
    - MAGESTIC_base_strain_allowed_Cas9_20mer_guides_NGG.txt
    - MAGESTIC_base_strain_disallowed_Cas9_20mer_guides_NGG.txt

Outputs (written to --output-dir)
---------
    NRD1_SpCas9_NGG_guide_donor_200bp_oligo_pool.tsv  — oligo names + sequences
    NRD1_SpCas9_NGG_guide_donor_200bp_info.tsv        — full per-oligo design details
    NRD1_SpCas9_NGG_guide_donor_200bp_untargetable.tsv — mutations that could not be targeted

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

from magestic.design.guide_donor_functions import (
    assign_num_syn_changes_for_Cas9,
    assign_num_syn_changes_for_Cas12a,
    assemble_donor,
    get_aa_num_to_codon_coords_aa,
    get_ORF_seq_and_CDS_coords,
    get_seq_with_largest_hamming_dist,
    guide_disruption,
    load_codon_table,
    load_genome,
    load_guides_from_fasta,
    load_ORF_coordinates,
    load_processed_guides,
    load_systematic_to_common_gene_name_dict,
    remove_BspQI_sites_from_donor_in_final_oligo,
    restriction_site_in_seq,
    rev_comp,
)
from magestic.projects.nrd1 import get_nrd1_data_path

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s %(name)s: %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# SpCas9 / NGG design parameters
# ---------------------------------------------------------------------------
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
RESTRICTION_SITES = ("GAAGAGC", "GCTCTTC")  # BspQI / SapI recognition sequences
FWD_PRIMING_SITE = "ctgcgattggcaggcgcgcc"
INTERNAL_CLONING_SITE = "gtttgaagagc"

# NRD1 occupies subpool index 0 in the multiplexed library
SUBPOOL_IDX = 0

# Output filename prefix
OUTPUT_PREFIX = "NRD1_SpCas9_NGG_guide_donor_200bp"

# Bundled reference filenames (accessed via get_nrd1_data_path)
CODON_TABLE_FILE = "DNA_codon_to_AA.txt"
CODON_TABLE_SUBOPT_FILE = "DNA_codon_to_AA_suboptimal_codons_removed.txt"
INPUT_FILE = "NRD1_alanine_scan_with_PTC_and_syn_controls.tsv"

# GenBank sequence files for the MAGESTIC 3.0 system components and recoded NRD1 strain.
# Guides that would target any of these integrated sequences are excluded from the library.
SEQUENCE_FILES = [
    "yT138.gb",          # pRNR2-tetR-TUP1 / pTDH3-tetO-tetR @ YORW∆17
    "yT170.gb",          # pGalL-SaCas9 / pRPR1-tetO-HHR-SaCas9_X1 sgRNA @ YNRC∆9
    "yT196.gb",          # MAGESTIC 3.0 Tir1 / Nrd1recoded-AID barcode landing pad @ LEU2
    "YLS30.gb",          # pADH1-OsTIR1-FLAG-AID-NRD1rec @ YPRC∆15
    "pS1439.gb",         # Step-2 cloning vector (2µ, SaCas9 cut sites, guide-donor scaffold)
    "pF834_pF835.gb",    # MAGESTIC 3.0 barcoded insert (pF835)
]

# Large reference file names (in --ref-dir)
GENOME_FILE = "MAGESTIC_background_strain.fasta"
GFF_FILE = "MAGESTIC_background_strain_annotations.gff"
ALLOWED_GUIDES_FILE = "MAGESTIC_base_strain_allowed_Cas9_20mer_guides_NGG.txt"
DISALLOWED_GUIDES_FILE = "MAGESTIC_base_strain_disallowed_Cas9_20mer_guides_NGG.txt"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--ref-dir",
        type=Path,
        required=True,
        help=(
            "Directory containing the four large reference files "
            "(genome FASTA, GFF, allowed/disallowed guide databases)."
        ),
    )
    parser.add_argument(
        "--rev-primer-file",
        type=Path,
        required=True,
        help=(
            "Path to rev_subpool_priming_sequences.txt "
            "(one reverse priming sequence per line, ordered by subpool index)."
        ),
    )
    parser.add_argument(
        "--sequence-map-dir",
        type=Path,
        required=True,
        help=(
            "Directory containing the MAGESTIC 3.0 system GenBank sequence maps "
            "(yT138.gb, yT170.gb, yT196.gb, YLS30.gb, pS1439.gb, pF834_pF835.gb). "
            "Guides targeting these sequences are excluded from the library."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Directory for output TSV files (default: current directory).",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Load codon tables (bundled data) ---
    codon_to_aa, aa_to_codon = load_codon_table(get_nrd1_data_path(CODON_TABLE_FILE))
    subopt_codon_to_aa, subopt_aa_to_codon = load_codon_table(
        get_nrd1_data_path(CODON_TABLE_SUBOPT_FILE)
    )

    # --- Load reverse priming sequences ---
    rev_primers = args.rev_primer_file.read_text().splitlines()
    REV_PRIMING_SEQUENCES = [rev_comp(seq.lower()) for seq in rev_primers if seq.strip()]
    REV_PRIMING_SITE = REV_PRIMING_SEQUENCES[SUBPOOL_IDX]

    # --- Load genome and annotations ---
    genome_seq = load_genome(args.ref_dir / GENOME_FILE)
    coords_to_guides_allowed = load_processed_guides(
        args.ref_dir / ALLOWED_GUIDES_FILE, PAMS_TO_USE
    )
    coords_to_guides_disallowed = load_processed_guides(
        args.ref_dir / DISALLOWED_GUIDES_FILE, PAMS_TO_USE
    )
    coords_to_guides = {**coords_to_guides_allowed, **coords_to_guides_disallowed}
    guide_sequences_to_exclude = load_guides_from_fasta(
        [args.sequence_map_dir / f for f in SEQUENCE_FILES],
        PAM_LENGTH, ALL_PAMS, GUIDE_LENGTH, GUIDE_UPSTREAM_PAM,
    )
    systematic_to_common_gene_name = load_systematic_to_common_gene_name_dict(
        args.ref_dir / GFF_FILE
    )
    ORF_dict = load_ORF_coordinates(args.ref_dir / GFF_FILE)

    # --- Parse input mutations (bundled data) ---
    infile = get_nrd1_data_path(INPUT_FILE)
    systematic_name_to_ORF_info: dict = {}
    ORFs_to_codons_to_mutate: dict = {}

    with open(infile) as fin:
        for line in fin:
            info = line.strip().split("\t")
            if not info or not info[0]:
                continue
            if info[0] == "systematic_ORF_name":  # skip header row
                continue
            ref_codons: list = []
            alt_codons: list = []

            if len(info) == 11:
                systematic_ORF_name, gene_name, aa_nums, ref_aas, alt_aas, ref_codons, alt_codons, *_ = info
            elif len(info) == 7:
                systematic_ORF_name, gene_name, aa_nums, ref_aas, alt_aas, ref_codons, alt_codons = info
            elif len(info) == 5:
                systematic_ORF_name, gene_name, aa_nums, ref_aas, alt_aas = info
            elif len(info) == 3:
                systematic_ORF_name, gene_name, aa_substitution = info
                ref_aas = [e[0] for e in aa_substitution.split(",")]
                aa_nums = [e[1:-1] for e in aa_substitution.split(",")]
                alt_aas = [e[-1] for e in aa_substitution.split(",")]
            else:
                raise ValueError(f"Unrecognized input format ({len(info)} fields): {line.rstrip()}")

            if systematic_ORF_name not in systematic_name_to_ORF_info:
                chrom, strand, exons, seq = get_ORF_seq_and_CDS_coords(
                    systematic_ORF_name, genome_seq, ORF_dict
                )
                ORF_info, aa_num_to_exon = get_aa_num_to_codon_coords_aa(
                    chrom, strand, exons, seq, codon_to_aa
                )
                systematic_name_to_ORF_info[systematic_ORF_name] = (
                    chrom, strand, exons, seq, ORF_info, aa_num_to_exon
                )

            chrom, strand, exons, seq, ORF_info, aa_num_to_exon = systematic_name_to_ORF_info[systematic_ORF_name]

            if len(info) in (5, 7, 11):
                ref_aas = ref_aas.split(",")
                aa_nums = aa_nums.split(",")
                alt_aas = alt_aas.split(",")
                ref_codons = ref_codons.split(",") if ref_codons else []
                alt_codons = alt_codons.split(",") if alt_codons else []

            aa_nums = [int(e) for e in aa_nums]
            aa_num_to_alt_aa = dict(zip(aa_nums, alt_aas))
            aa_num_to_alt_codon = dict(zip(aa_nums, alt_codons)) if alt_codons else {}

            # Fill in the full contiguous range of aa_nums and derive alt codons
            ref_aas_out, alt_aas_out, aa_nums_out, alt_codons_out = [], [], [], []
            for aa_num in range(min(aa_nums), max(aa_nums) + 1):
                ref_codon = ORF_info[aa_num][0]
                ref_aa = ORF_info[aa_num][2]
                alt_aa = aa_num_to_alt_aa.get(aa_num, ref_aa)
                aa_nums_out.append(aa_num)
                ref_aas_out.append(ref_aa)
                alt_aas_out.append(alt_aa)
                if alt_codons and aa_num in aa_num_to_alt_codon:
                    alt_codons_out.append(aa_num_to_alt_codon[aa_num])
                else:
                    possible = subopt_aa_to_codon[alt_aa]
                    alt_codons_out.append(get_seq_with_largest_hamming_dist(ref_codon, possible))
            ref_codons_out = [ORF_info[n][0] for n in aa_nums_out]

            ORFs_to_codons_to_mutate.setdefault(systematic_ORF_name, []).append(
                (aa_nums_out, ref_aas_out, alt_aas_out, ref_codons_out, alt_codons_out)
            )

    # --- Write outputs ---
    info_path = out_dir / f"{OUTPUT_PREFIX}_info.tsv"
    oligo_pool_path = out_dir / f"{OUTPUT_PREFIX}_oligo_pool.tsv"
    untargetable_path = out_dir / f"{OUTPUT_PREFIX}_untargetable.tsv"

    with (
        open(untargetable_path, "w") as untargetable_out,
        open(oligo_pool_path, "w") as oligo_pool_out,
        open(info_path, "w") as info_out,
    ):
        untargetable_out.write(
            "\t".join([
                "systematic_ORF_name", "common_ORF_name",
                "aa_nums", "ref_aas", "alt_aas", "ref_codons", "alt_codons",
                "ORF_strand", "chrom", "reasons_not_targetable",
            ]) + "\n"
        )
        oligo_pool_out.write("oligo_name\toligo_seq\n")
        info_out.write(
            "\t".join([
                "oligo_name", "oligo_sequence",
                "subpool_ID", "subpool_number",
                "common_ORF_name", "systematic_ORF_name",
                "aa_nums_str", "ref_aas_str", "alt_aas_str",
                "ref_codons_str", "alt_codons_str",
                "num_US_syn_changes", "US_syn_codons",
                "num_DS_syn_changes", "DS_syn_codons",
                "ORF_strand", "chrom",
                "PAM_strand", "PAM_coord", "PAM",
                "PAM_proximal_seed_region_to_codon_changes_midpoint",
                "leftmost_codon_coord", "rightmost_codon_coord",
                "FWD_PRIMING_SITE", "guide", "INTERNAL_CLONING_SITE",
                "donor", "extended_donor", "REV_PRIMING_SITE",
                "donor_changed", "donor_changes", "junction_changes",
            ]) + "\n"
        )

        guide_donors_created: set = set()

        for systematic_ORF_name, mutations in ORFs_to_codons_to_mutate.items():
            common_ORF_name = systematic_to_common_gene_name.get(
                systematic_ORF_name, systematic_ORF_name
            )
            chrom, strand, exons, seq = get_ORF_seq_and_CDS_coords(
                systematic_ORF_name, genome_seq, ORF_dict
            )
            ORF_info, _ = get_aa_num_to_codon_coords_aa(chrom, strand, exons, seq, codon_to_aa)

            for aa_nums, ref_aas, alt_aas, ref_codons, alt_codons in mutations:
                variant_targetable = False
                reasons_not_targetable: list = []
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
                            PAM_seed = PAM_coord - 3 if PAM_strand == "+" else PAM_coord + 3
                            dist_to_mid = codon_mid - PAM_seed
                            if strand == "-":
                                dist_to_mid = -dist_to_mid
                            if (chrom, PAM_coord, PAM_strand) not in coords_to_guides:
                                continue
                            guide, PAM = coords_to_guides[(chrom, PAM_coord, PAM_strand)]
                            if (
                                guide in guide_sequences_to_exclude
                                or (chrom, PAM_coord, PAM_strand) in coords_to_guides_disallowed
                                or restriction_site_in_seq(guide, RESTRICTION_SITES)
                                or "T" * (NUM_CONSEC_T_ALLOWED + 1) in guide
                            ):
                                continue
                            potential_PAMs.append(
                                (abs(dist_to_mid), dist_to_mid, chrom, PAM_coord, PAM_strand, guide, PAM)
                            )

                    if not potential_PAMs:
                        reasons_not_targetable.append("no acceptable guides found")
                    else:
                        for _, dist_to_mid, chrom, PAM_coord, PAM_strand, guide, PAM in sorted(potential_PAMs):
                            syn_left, syn_right = assign_num_syn_changes_for_Cas9(
                                leftmost, rightmost, PAM_coord, PAM_strand
                            )
                            donor_result = assemble_donor(
                                systematic_ORF_name, ORF_info, chrom, strand,
                                aa_nums, ref_aas, alt_aas, ref_codons, alt_codons,
                                syn_left, syn_right,
                                OLIGO_LENGTH
                                    - len(FWD_PRIMING_SITE)
                                    - GUIDE_LENGTH
                                    - len(INTERNAL_CLONING_SITE)
                                    - len(REV_PRIMING_SITE),
                                MIN_HOMOLOGY, codon_to_aa, subopt_aa_to_codon,
                                genome_seq, PRINT=False, seqs_to_avoid=RESTRICTION_SITES,
                            )

                            if donor_result[0] is None:
                                reasons_not_targetable.append(donor_result)
                                continue

                            (
                                extended_donor, donor, donor_changed, donor_changes,
                                aa_nums, ref_aas, alt_aas, ref_codons, alt_codons,
                                num_US_syn_changes, US_syn_codons,
                                num_DS_syn_changes, DS_syn_codons,
                            ) = donor_result

                            donor, junction_changes = remove_BspQI_sites_from_donor_in_final_oligo(
                                FWD_PRIMING_SITE, guide, INTERNAL_CLONING_SITE,
                                donor, REV_PRIMING_SITE, donor_result,
                            )

                            if donor is None:
                                reasons_not_targetable.append(
                                    "donor could not be assembled into oligo due to unresolvable restriction sites"
                                )
                                continue

                            oligo_sequence = (
                                FWD_PRIMING_SITE + guide + INTERNAL_CLONING_SITE
                                + donor + REV_PRIMING_SITE
                            )

                            if not guide_disruption(
                                guide, DEGENERATE_PAM, extended_donor, GUIDE_UPSTREAM_PAM, PRINT=False
                            ):
                                reasons_not_targetable.append("guide not disrupted by donor")
                                continue

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

                            oligo_name = "_".join(map(str, [
                                NUCLEASE, systematic_ORF_name, common_ORF_name,
                                aa_nums_str, ref_aas_str, alt_aas_str,
                                ref_codons_str, alt_codons_str,
                                num_US_syn_changes, "US_syn_changes",
                                num_DS_syn_changes, "DS_syn_changes",
                                strand, chrom, PAM_strand, PAM_coord, PAM, dist_to_mid,
                                leftmost, rightmost,
                            ]))

                            info_out.write(
                                "\t".join(map(str, [
                                    oligo_name, oligo_sequence,
                                    NUCLEASE + "_codon_changes",
                                    f"subpool_{SUBPOOL_IDX + 1}",
                                    common_ORF_name, systematic_ORF_name,
                                    aa_nums_str, ref_aas_str, alt_aas_str,
                                    ref_codons_str, alt_codons_str,
                                    num_US_syn_changes, US_syn_codons,
                                    num_DS_syn_changes, DS_syn_codons,
                                    strand, chrom,
                                    PAM_strand, PAM_coord, PAM, dist_to_mid,
                                    leftmost, rightmost,
                                    FWD_PRIMING_SITE, guide, INTERNAL_CLONING_SITE,
                                    donor, extended_donor, REV_PRIMING_SITE,
                                    donor_changed, donor_changes, junction_changes,
                                ])) + "\n"
                            )
                            oligo_pool_out.write(f"{oligo_name}\t{oligo_sequence}\n")

                            if len(guide_donors_created) >= MAX_GUIDE_DONORS_PER_EDIT:
                                break

                if not variant_targetable:
                    untargetable_out.write(
                        "\t".join(map(str, [
                            systematic_ORF_name, common_ORF_name,
                            aa_nums, ref_aas, alt_aas, ref_codons, alt_codons,
                            strand, chrom, reasons_not_targetable,
                        ])) + "\n"
                    )

    logger.info("Done. Outputs written to %s", out_dir)


if __name__ == "__main__":
    main()
