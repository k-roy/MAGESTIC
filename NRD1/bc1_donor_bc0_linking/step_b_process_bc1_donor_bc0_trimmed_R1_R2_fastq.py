"""
This script takes collapsed bc1-donor-bc0 reads from 2 x 150 bp AVITI reads, and processes the bc1, donor fragments, and bc0 sequences.

@AV240704:298-11_run1:2449525515:1:10102:0883:0682 1:N:0:TAGCGGAG+CCCAAGCG
TCGACAGGTCATGCTC CCGGTCACCGCAGACTGTAC CCCTAGGATACGTGCATAACTTCGTATAATGTATGCTATACGAACGGTATCCCCTATTTCATAGGGCACTATCCATCCTTACCCACTCCTCCGT
+
NNNMNNNNNNNNMNMNNNNJNMNNNMNJNNNNNNNKNCNNNMNMNNNMNKNNMKNNINGNNNNBMMNLFMMLMIMKMLMMLKMMMM=;MLLMMMMMMLLMMKMJFLJLLLLKLLLLLHKKLLAFKLGKGE

expected PCR products after merging the 2 x 300 bp reads:
KR1967-KR1590 from yeast_plasmid_assembly
from read 1 after trimming 20 bp:
TCGACAGGTCATGCTCNNNNNNNNNNNNNNNNNNNNCCCTAggatAcgtgcATAACTTCGTATAATGTATGCTATACGAAcggtaTCCCCTATTTCATAGGGCACTATCCATcctTACCCACTCCTCCGTGGGTTTCTTCCACTCAGGTGCGCACCCTTACGTCAGAAGAAACGGGTTTCCTgtttgaagagcNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGATAGACGACTGGACAGCANNNNNNNNNNAGGAAAACAGAC
TCGACAGGTCATGCTCNNNNNNNNNNNNNNNNNNNNCCCTAggatAcgtgcATAACTTCGTATAATGTATGCTATACGAAcggtaTCCCCTATTTCATAGGGCACTATCCATcctTACCCACTCCTCCGT--------------------------------------------------------------------------------------------------------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGATAGACGACTGGACAGCANNNNNNNNNNAGGAAAACAGAC

For the 2 x 150 bp reads, we will not get the hyphenated portion of the insert, missing 41 bp of the donor.
read 1:
TCGACAGGTCATGCTCNNNNNNNNNNNNNNNNNNNNCCCTAggatAcgtgcATAACTTCGTATAATGTATGCTATACGAAcggtaTCCCCTATTTCATAGGGCACTATCCATcctTACCCACTCCTCCGT

read 2:
GTCTGTTTTCCTNNNNNNNNNNTGCTGTCCAGTCGTCTATCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

output:
counts, PS1, bc1, MAGESTIC_3_homology, RE_site, donor, SPS, bc0, msd
"""

import pysam
import timeit
from sys import argv
from tqdm import tqdm
import os

R1_infilename = argv[1]
R2_infilename = argv[2]
outfilename = argv[3]
overwrite = bool(argv[4])  # True
max_reads_per_sample = 10e10
SaCas9_SceI_fragment_seq = "CCCTAggatAcgtgc".upper()  # This is just downstream of bc1.
SaCas9_SceI_fragment_length = len(SaCas9_SceI_fragment_seq)

bc1_length = 20
fwd_bc1_prefix = "CAGGTCATGCTC"
fwd_bc1_suffix = "CCCTAGGATACG"
fwd_bc1_prefix_length = len(fwd_bc1_prefix)

if not overwrite and os.path.exists(outfilename):
    print("Output file already exists. Exiting the program.")
    exit()


def rev_comp(seq):
    return seq.translate(
        str.maketrans("ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn ")
    )[::-1]


def process_R1_seq(seq):
    if fwd_bc1_prefix in seq:
        # in the bc1-donor-bc0 data, bc1 are defined as the 20 bp after TCGACAGGTCATGCTC
        bc1_start = seq.index(fwd_bc1_prefix) + fwd_bc1_prefix_length
        bc1_end = bc1_start + bc1_length
        bc1 = seq[bc1_start:bc1_end]
        PS1 = seq[:bc1_start]
    else:
        bc1 = seq[16 : 16 + bc1_length]
        PS1 = seq[:16]

    SaCas9_SceI_fragment_found = True
    SaCas9_SceI_fragment_idx = seq.find(SaCas9_SceI_fragment_seq)
    if SaCas9_SceI_fragment_idx == -1:
        SaCas9_SceI_fragment_idx = 36
        SaCas9_SceI_fragment_found = False
    SaCas9_SceI_fragment = seq[
        SaCas9_SceI_fragment_idx : SaCas9_SceI_fragment_idx
        + SaCas9_SceI_fragment_length
    ]

    bc1_proximal_MAGESTIC_3_homology = seq[
        SaCas9_SceI_fragment_idx:
    ]  # The whole MAGESTIC homology cannot be captured with the 2 x 150 bp reads.
    # The RE site will be captured by the KR1883 read 2 primer.
    return [
        PS1,
        bc1,
        SaCas9_SceI_fragment,
        SaCas9_SceI_fragment_found,
        bc1_proximal_MAGESTIC_3_homology,
    ]


def process_R2_seq(seq):
    donor_fragment = seq[42:]
    SPS = seq[22:42]
    bc0 = seq[12:22]
    msd = seq[:12]
    donor_bc0_fragment = seq[12:72]
    return [rev_comp(e) for e in [donor_fragment, SPS, bc0, msd, donor_bc0_fragment]]


R1_R2_seqs_to_counts = {}
start_time = timeit.default_timer()
reads_processed_for_sample = 0
with (
    pysam.FastxFile(R1_infilename) as R1,
    pysam.FastxFile(R2_infilename) as R2,
):
    for R1_read, R2_read in tqdm(zip(R1, R2), desc="Processing reads"):
        reads_processed_for_sample += 1
        R1_R2_seqs = process_R1_seq(R1_read.sequence) + process_R2_seq(R2_read.sequence)
        R1_R2_seqs_to_counts[tuple(R1_R2_seqs)] = (
            R1_R2_seqs_to_counts.get(tuple(R1_R2_seqs), 0) + 1
        )

elapsed = timeit.default_timer() - start_time
print(R1_infilename + " took " + str(elapsed) + " seconds to process\n")
print(
    R1_infilename
    + " had "
    + str(reads_processed_for_sample)
    + " total_seqs_processed\n"
)

with open(outfilename, "w") as outfile:
    header = (
        "\t".join(
            "counts, PS1, bc1, SaCas9_SceI_fragment, SaCas9_SceI_fragment_found, bc1_proximal_MAGESTIC_3_homology, donor_fragment, SPS, bc0, msd, donor_bc0_fragment".split(
                ", "
            )
        )
        + "\n"
    )
    outfile.write(header)
    # Sort by count, descending
    for seqs, count in sorted(
        R1_R2_seqs_to_counts.items(), key=lambda x: x[1], reverse=True
    ):
        str_output = f"{count}\t" + "\t".join([str(e) for e in seqs]) + "\n"
        _ = outfile.write(str_output)
