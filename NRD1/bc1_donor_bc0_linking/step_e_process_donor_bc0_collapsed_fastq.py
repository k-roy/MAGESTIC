"""
This script takes collapsed bc1-donor-bc0 reads from 2 x 150 bp AVITI reads, and processes the bc1, donor fragments, and bc0 sequences.

collapsed fastq input:
@443
TGTTTTCCT GGTTTTTTTA TGCTGTCCAGTCGTCTATCG TTGGTTAAGCATATTCATCAAAGAATTCAATTGAGCAGTGGGGTCAAATTGTTGCTGAGGGACTGGTGGAGCCGCTGCGGCAGGTCCCTGAGAAGGTAATGGTTGGTTGGGGGCATAACCATACGGTTG GCTCTTCAAACAGGA
+
MNNMNMNLLJNNLNMNMMMMGMMIMKMMKMLIM5LMNNJNNF8NNMNNNNNNNSSSSSSSSSSSSSSSSSSJRSSSPSSSSSSSSSSSSKSSSSRSSSSSSSSSSSSS@SSQSSSSSSSSPSSSSSSSSSNMMLNHFMFMNMLLICDILBJFDMMMMMKJFM>MCJGIMLD>DHIEEML>+K;
@350
TGTTTTCCT TTCACTTCTC TGCTGTCCAGTCGTCTATCG ACTCTTTTGAAACAAGCCGGACCTGTCCCAAATGTCTAAAAGCATACGAATTTTTTCTTTATGATCGGCGTTACTTTTTGCAATAGCATCGGATAGCAGTTCTTGAATTACTTCGCCCAACGTGTTTAT GCTCTTCAAACAGGA
+
KLLJKLLM=KNKMMINLCKMKMMGMKMKKLEIK5MLMNMMNMMGFNNNMMNNKRSLJSSSSQSSSSSSMSSSSSFSSSSSSPSSPRSSSSDOKSISSSSSSSSSOOPQSSSRSKSSSPSRSRSSIR:SQRG@CJGHJBD7L>;J>=MDCMBE=J=@3GC>K6EJE:NDMLLLMHJLLLKLHK?


expected PCR products:
KR1882-KR1883
donor_bc0
rev
SpCas9_QTL_step_2_or_yeast_plasmid_assembly
view merged read from read 1 direction:
HNHGAATCTGAGTTACTGTCTGTTTTCCTNNNNNNNNNNCGTCTCTCGTCCACAAGTTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNgctcttcaaacAGGAAACCCGTTTCTTCTGACDND
after trimming:
TGTTTTCCTNNNNNNNNNNCGTCTCTCGTCCACAAGTTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNgctcttcaaacAGG

view merged read from read 2 direction:
HNHGTCAGAAGAAACGGGTTTCCTgtttgaagagcNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAACTTGTGGACGAGAGACGNNNNNNNNNNAGGAAAACAGACAGTAACTCAGATTCDND
after trimming:
TCCTgtttgaagagcNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAACTTGTGGACGAGAGACGNNNNNNNNNNAGGAAAACA

# Reads will be processed as viewed from the read 2 direction.

donor and bc0 should be oriented 5' to 3', guide to donor
bc0, counts, seq_without_bc0, yL

"""

import pysam
import timeit
from sys import argv
from tqdm import tqdm
import os

infilename = argv[1]
outfilename = argv[2]
overwrite = True  # bool(argv[3])
max_reads_per_sample = 10e10

RE_site_seq = "gtttgaagagc".upper()
RE_site_length = len(RE_site_seq)

if not overwrite and os.path.exists(outfilename):
    print("Output file already exists. Exiting the program.")
    exit()


def rev_comp(seq):
    return seq.translate(
        str.maketrans("ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn ")
    )[::-1]


def process_merged_seq(seq):
    seq = rev_comp(seq)
    RE_site_found = True
    RE_site_idx = seq.find(RE_site_seq)
    if RE_site_idx == -1:
        RE_site_idx = 3
        RE_site_found = False

    RE_site = seq[RE_site_idx : RE_site_idx + RE_site_length]
    donor = seq[RE_site_idx + RE_site_length : -39]
    SPS = seq[-39:-19]
    bc0 = seq[-19:-9]
    donor_bc0_fragment = seq[-69:-9]
    return [RE_site_found, RE_site, donor, SPS, bc0, donor_bc0_fragment]


start_time = timeit.default_timer()
reads_processed_for_sample = 0

with open(outfilename, "w") as outfile:
    header = (
        "\t".join(
            "counts, RE_site_found, RE_site, donor, SPS, bc0, donor_bc0_fragment".split(
                ", "
            )
        )
        + "\n"
    )
    outfile.write(header)

    with pysam.FastxFile(infilename) as merged_fastq:
        for merged_read in tqdm(merged_fastq, desc="Processing reads"):
            reads_processed_for_sample += 1
            counts = int(merged_read.name)
            seqs = process_merged_seq(merged_read.sequence)
            str_output = f"{counts}\t" + "\t".join([str(e) for e in seqs]) + "\n"
            _ = outfile.write(str_output)

elapsed = timeit.default_timer() - start_time
print(infilename + " took " + str(elapsed) + " seconds to process\n")
print(
    infilename + " had " + str(reads_processed_for_sample) + " total_seqs_processed\n"
)
