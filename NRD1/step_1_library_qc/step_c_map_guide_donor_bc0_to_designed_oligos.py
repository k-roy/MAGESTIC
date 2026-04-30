"""
Created on Fri Jun 19 02:23:35 2020

@author: kevinroy

This script will process each step 1 library separately based on perfect matching.

"""

import pandas as pd
import timeit
import os
from tqdm import tqdm

PROJECT_PATH = "by_project/NNS/202250628_repeat_step_1_library_cloning/"

PROCESSED_DATA_DIR = (
    "/path/to/processed_data/" + PROJECT_PATH
)
os.chdir(PROCESSED_DATA_DIR)
COLLAPSED_FASTQ_DIR = PROCESSED_DATA_DIR + "collapsed_fastq/"
COUNTS_DIR = PROCESSED_DATA_DIR + "guide_donor_bc0/"
for SUB_DIR in "guide_donor_bc0/", "log/":
    if not os.path.exists(PROCESSED_DATA_DIR + SUB_DIR):
        os.mkdir(PROCESSED_DATA_DIR + SUB_DIR)

KEYFILE_DIR = (
    "/path/to/scripts_and_keyfiles/"
    + PROJECT_PATH
    + "keyfiles/"
)
sample_keyfile = KEYFILE_DIR + "step_1_sample_key.tsv"

oligo_pool_filename = KEYFILE_DIR + "20240411_Twist_200mer_oligo_array_order.tsv"


def rev_comp(seq):
    return seq.translate(
        str.maketrans("ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn ")
    )[::-1]


print("loading guide-donors...")

SpCas9_leader_length = 20
SUBPOOL_PRIMING_SEQ_LENGTH = 20
SPS_length = 20
GUIDE_LENGTH = 20
SpCas9_donor_length = 129
SpCas9_cloning_site = "GTTTGAAGAGC"
lines_to_test = 100000000000000000000
debug_print = False  # True #

oligo_pool_df = pd.read_csv(oligo_pool_filename, sep="\t")

oligo_pool_df.columns
oligo_pool_df["guide"] = oligo_pool_df["oligo_seq"].str.slice(20, 40)
oligo_pool_df["donor"] = oligo_pool_df["oligo_seq"].str.slice(51, 180).str.upper()
oligo_pool_df["middle_donor"] = (
    oligo_pool_df["oligo_seq"].str.slice(51 + 30, 51 + 90).str.upper()
)

middle_donor_start_idx = 30
middle_donor_end_idx = 90
middle_donor_length = middle_donor_end_idx - middle_donor_start_idx

# Create dictionaries that will handle cases where guides might have more than one donor, and vice versa.
guide_to_oligo_names = {}
perfect_donor_to_oligo_names = {}
partial_donor_to_oligo_names = {}
guide_donor_to_oligo_names = {}
for index, row in oligo_pool_df.iterrows():
    guide = row["guide"]
    donor = row["donor"]
    oligo_name = row["oligo_name"]
    guide_donor = guide + SpCas9_cloning_site + donor
    guide_donor_to_oligo_names[guide_donor] = oligo_name
    if guide not in guide_to_oligo_names:
        guide_to_oligo_names[guide] = []
    if oligo_name not in guide_to_oligo_names[guide]:
        guide_to_oligo_names[guide] += [oligo_name]
    if donor not in perfect_donor_to_oligo_names:
        perfect_donor_to_oligo_names[donor] = []
    if oligo_name not in perfect_donor_to_oligo_names[donor]:
        perfect_donor_to_oligo_names[donor] += [oligo_name]
    for idx in range(0, len(donor) - middle_donor_length):
        partial_donor = donor[idx : idx + middle_donor_length]
        if partial_donor not in partial_donor_to_oligo_names:
            partial_donor_to_oligo_names[partial_donor] = []
        if oligo_name not in partial_donor_to_oligo_names[partial_donor]:
            partial_donor_to_oligo_names[partial_donor] += [oligo_name]


def match_guide_donor_bc0(bc0_seq_infilename, outfilename):
    start_time = timeit.default_timer()
    lines_processed = 0
    header = "library_ID, scaffold, bc0, counts, perfect_leader, leader, perfect_guide, guide, perfect_donor, perfect_middle_donor_region, donor, donor_bc0_fragment, unambiguous_guide, unambiguous_donor, oligo_name".split(
        ", "
    )
    with open(bc0_seq_infilename, "r") as infile, open(outfilename, "w") as outfile:
        outfile.write("\t".join(header) + "\n")
        infile_header = infile.readline()
        for line in infile:
            lines_processed += 1
            if lines_processed > lines_to_test:
                break
            if lines_processed % 10000 == 0:
                elapsed = timeit.default_timer() - start_time
                print(
                    "Perfect match analysis took: "
                    + str(elapsed / 60)
                    + " minutes for "
                    + str(lines_processed)
                    + " sequences\n"
                )
            (
                library_ID,
                total_counts,
                total_fwd_counts,
                total_rev_counts,
                fwd_0707,
                rev_0707,
                leader,
                guide,
                donor,
                SPS,
                bc0,
                msd,
                donor_bc0_fragment,
                seq,
            ) = line.strip("\n").split("\t")
            counts = total_counts
            cloning_site_found = False
            perfect_leader = False
            perfect_guide = False
            perfect_donor = False
            perfect_middle_donor_region = False
            unambiguous_guide = False
            unambiguous_donor = False
            scaffold, oligo_name = "NA", "NA"
            guide_oligo_names = None
            donor_oligo_names = None
            if guide != "NA":
                cloning_site_found = True
                if leader == "GCAGGCGCGCC":  # GCAGGCGCGCC
                    perfect_leader = True
                scaffold = "WT_SpCas9"
                middle_donor = donor[middle_donor_start_idx:middle_donor_end_idx]
                guide_donor = guide + SpCas9_cloning_site + donor
            if debug_print:
                print(library_ID, bc0, counts, seq)
            if cloning_site_found:
                if guide_donor in guide_donor_to_oligo_names:
                    if debug_print:
                        print(
                            guide_donor, "\nguide_donor in guide_donor_to_oligo_names"
                        )
                    oligo_name = guide_donor_to_oligo_names[guide_donor]
                    perfect_guide = True
                    perfect_donor = True
                    perfect_middle_donor_region = True
                    unambiguous_guide = True
                    unambiguous_donor = True
                else:
                    if debug_print:
                        print(
                            "guide_donor not in guide_donor_to_oligo_names", guide_donor
                        )
                    if guide in guide_to_oligo_names:
                        perfect_guide = True
                        guide_oligo_names = guide_to_oligo_names[guide]
                        if debug_print:
                            print(
                                "guide in guide_to_oligo_names:",
                                guide in guide_to_oligo_names,
                            )
                        if debug_print:
                            print("guide_oligo_names:", guide_to_oligo_names[guide])
                        if len(guide_oligo_names) == 1:
                            unambiguous_guide = True
                    if donor in perfect_donor_to_oligo_names:
                        perfect_donor = True
                        perfect_middle_donor_region = True
                        donor_oligo_names = perfect_donor_to_oligo_names[donor]
                        if debug_print:
                            print(
                                "donor in donor_oligo_names:",
                                donor in perfect_donor_to_oligo_names,
                            )
                        if debug_print:
                            print(
                                "donor_oligo_names:",
                                perfect_donor_to_oligo_names[donor],
                            )
                        if len(donor_oligo_names) == 1:
                            unambiguous_donor = True
                    elif middle_donor in partial_donor_to_oligo_names:
                        perfect_middle_donor_region = True
                        donor_oligo_names = partial_donor_to_oligo_names[middle_donor]
                        if len(donor_oligo_names) == 1:
                            unambiguous_donor = True
                    if guide_oligo_names != None and donor_oligo_names != None:
                        intersection = [
                            oligo_name
                            for oligo_name in guide_oligo_names
                            if oligo_name in donor_oligo_names
                        ]
                        if len(intersection) == 1:
                            oligo_name = intersection[0]
                            unambiguous_guide = True
                            unambiguous_donor = True
                    if guide_oligo_names == None:
                        if unambiguous_donor:
                            oligo_name = donor_oligo_names[0]
                    if donor_oligo_names == None:
                        if unambiguous_guide:
                            oligo_name = guide_oligo_names[0]

                output = [
                    library_ID,
                    scaffold,
                    bc0,
                    counts,
                    perfect_leader,
                    leader,
                    perfect_guide,
                    guide,
                    perfect_donor,
                    perfect_middle_donor_region,
                    donor,
                    donor_bc0_fragment,
                    unambiguous_guide,
                    unambiguous_donor,
                    oligo_name,
                ]
                if debug_print:
                    print(output)
                outfile.write("\t".join([str(e) for e in output]) + "\n")
    outfile.close()

    elapsed = timeit.default_timer() - start_time
    print(
        "Perfect match analysis took: "
        + str(elapsed / 60)
        + " minutes for "
        + str(lines_processed)
        + " sequences\n"
    )


# process files for each library_ID according to the sample key
library_ID_to_sample_info = {}
with open(sample_keyfile) as infile:
    header = infile.readline()
    for line in infile:
        sample_number, library_ID, inner_primers, outer_primers, dates_of_sequencing = (
            line.strip().split("\t")
        )
        if library_ID not in library_ID_to_sample_info:
            library_ID_to_sample_info[library_ID] = [
                (sample_number, inner_primers, outer_primers, dates_of_sequencing)
            ]
        else:
            library_ID_to_sample_info[library_ID] += [
                (sample_number, inner_primers, outer_primers, dates_of_sequencing)
            ]

for library_ID in tqdm(library_ID_to_sample_info, desc="Processing libraries"):
    bc0_seq_infilename = COUNTS_DIR + library_ID + "_guide_donor_bc0_counts.tsv"
    outfilename = COUNTS_DIR + library_ID + "_matched_guide_donor_bc0_counts.tsv"
    match_guide_donor_bc0(bc0_seq_infilename, outfilename)
