"""
Created on Fri Jun 19 02:23:35 2020

@author: kevinroy

This script will process each step 1 library separately based on perfect matching.

"""

import pandas as pd
import timeit
import os

BASE_DIR = '/path/to'
PROJECT_DIR = BASE_DIR + '/projects/NNS/20250628_repeat_step_1_library_cloning'
PROCESSED_DATA_DIR = PROJECT_DIR + '/processed_data/'
COUNTS_DIR = PROCESSED_DATA_DIR + 'guide_donor_bc0/'
KEYFILE_DIR = PROJECT_DIR + '/keyfiles/'
sample_keyfile = KEYFILE_DIR + 'step_1_sample_key.tsv'
oligo_pool_filename = KEYFILE_DIR + '20240411_Twist_200mer_oligo_array_order.tsv'

# need a keyfile matching each V library_ID to the expected SPS and subpool index.
library_ID_subpool_summary_keyfile = KEYFILE_DIR + 'library_ID_subpool_summary_key.tsv'

def rev_comp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdNn ', 'TGCAtgcaYRKMyrkmBVDHbvdhNn '))[::-1]

print('loading guide-donors...')

SpCas9_leader_length = 20
SUBPOOL_PRIMING_SEQ_LENGTH = 20
SPS_length = 20
GUIDE_LENGTH = 20
SpCas9_donor_length = 129
SpCas9_cloning_site = 'GTTTGAAGAGC'
lines_to_test = 100000000000000000000
debug_print = False # True # 

oligo_info = pd.read_csv(oligo_pool_filename, sep = '\t')
oligo_info.columns
oligo_info['guide'] = oligo_info['oligo_seq'].str.slice(20, 40)
oligo_info['donor'] = oligo_info['oligo_seq'].str.slice(51, 180).str.upper()
oligo_info['middle_donor'] = oligo_info['oligo_seq'].str.slice(51+30, 51+90).str.upper()
oligo_info['SPS'] = oligo_info['oligo_seq'].str.slice(180, 201).str.upper()

middle_donor_start_idx = 30
middle_donor_end_idx = 90
middle_donor_length = middle_donor_end_idx - middle_donor_start_idx

library_ID_subpool_summary_key = pd.read_csv(library_ID_subpool_summary_keyfile, sep='\t')
library_ID_subpool_summary_key.columns
library_ID_subpool_summary_key['SPS'] = library_ID_subpool_summary_key['rev_priming_seq'].str.upper()
library_ID_to_expected_SPS = {}
for index, row in library_ID_subpool_summary_key.iterrows():
    print(index, row)
    library_ID = row['library_ID']
    SPS = row['SPS']
    library_ID_to_expected_SPS[library_ID] = SPS


# process files for each library_ID according to the sample key
library_ID_to_sample_info = {}
with open(sample_keyfile) as infile:
    header = infile.readline()
    for line in infile:
        sample_number, library_ID, inner_primers, outer_primers, dates_of_sequencing = line.strip().split('\t')
        if library_ID not in library_ID_to_sample_info:
            library_ID_to_sample_info[library_ID] = [(sample_number, inner_primers, outer_primers, dates_of_sequencing)]
        library_ID_to_sample_info[library_ID] += [(sample_number, inner_primers, outer_primers, dates_of_sequencing)]

library_ID = 'V536'

NRD1_oligo_info_file = BASE_DIR + '/software/MAGESTIC/NRD1/NRD1_SpCas9_NGG_guide_donor_200bp_info.tsv'
for library_ID in 'V536',: # library_ID_to_sample_info: #  
    expected_SPS = library_ID_to_expected_SPS[library_ID]
    bc0_seq_matched_oligo_infilename = COUNTS_DIR + library_ID + '_matched_guide_donor_bc0_counts.tsv'
    bc0_seq_matched_oligo_df = pd.read_csv(bc0_seq_matched_oligo_infilename, sep='\t')
    # AGACGTGAACGTGAACGTGAAAGAGAAAGA CGATAGACGACTGGACAGCA GGCAATGCAA
    bc0_seq_matched_oligo_df['SPS'] = bc0_seq_matched_oligo_df['donor_bc0_fragment'].str.slice(30, 50)
    # left join based on SPS
    merged_df = pd.merge(bc0_seq_matched_oligo_df, library_ID_subpool_summary_key, on='SPS', how='left')
    merged_df.columns
    merged_df['library_ID_y'].value_counts()
    designed_oligos = oligo_info[oligo_info['SPS'] == expected_SPS]
    oligo_info_df = pd.read_csv(NRD1_oligo_info_file, sep='\t')
    df_for_plotting = merged_df.query('oligo_name.str.contains("NRD1", na=False)')
    merged_df[merged_df['oligo_name'].str.contains('NRD1')]['oligo_name'].value_counts()
    merged_df.query('oligo_name.str.contains("NRD1", na=False)')['oligo_name'].value_counts()
    merged_df.columns
    merged_df.groupby
