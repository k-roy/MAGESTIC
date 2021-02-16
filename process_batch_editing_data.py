#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:50:53 2017

@author: kevinroy
"""
import operator
import subprocess
import pandas as pd
from Bio import SeqIO

## note demultiplexing from Novogene means no S numbers in filenames
GOOGLE_DRIVE_DIR= '/Volumes/SPxDrive/' # '/Users/kevinroy/Google_Drive/'

DIR = '/Volumes/Samsung_T5/sequencing_data/191022_Church_method/'
COLLAPSED_FASTQ_DIR = DIR + 'collapsed_fastq/' # argv[1]  #
FASTA_REF_DIR = DIR +  'fasta/'
ALIGNMENTS_DIR = DIR +  'alignments/'
PROCESSED_COUNTS_DIR = DIR + 'processed_counts/'

## load and process keyfile
plate_ID_map_keyfile = DIR + 'plate_ID_map.xlsx'
guide_donor_plasmid_keyfile = DIR + 'plasmid_keyfile.xlsx'
trafo_ID_keyfile = DIR + 'TRAFO_plate_layout.xlsx'
FASTA_REF_DIR = GOOGLE_DRIVE_DIR + 'SGTC_genome_editing_group/Projects/single_guide_donor_editing/natural_variant_guide_donors_fasta/'

plate_ID_map_key = pd.read_excel(plate_ID_map_keyfile)
sample_numbers = pd.Series(['sample_' + str(i) for i in range(1,97) ])
data_dict = plate_ID_map_key.to_dict('list')
new_df = {}
new_df['plate_ID'] = data_dict['plate_ID']
new_df['sample_number'] = ['sample_' + str(i) for i in range(1,97) ]
import itertools
def expand_grid(data_dict):
    rows = itertools.product(*data_dict.values())
    return pd.DataFrame.from_records(rows, columns=data_dict.keys())

plate_ID_map = expand_grid(new_df)
plate_ID_map['plate_ID'] = plate_ID_map['plate_ID'].astype('category')
plate_ID_map['sample_number'] = plate_ID_map['sample_number'].astype('category')
plate_ID_map.plate_ID = pd.Categorical(plate_ID_map.plate_ID, categories = list(plate_ID_map_key.plate_ID), ordered=True)
plate_ID_map.sample_number = pd.Categorical(plate_ID_map.sample_number, categories = list(sample_numbers), ordered=True)

keyfile = plate_ID_map.sort_values(by = ['plate_ID', 'sample_number']).reset_index(drop=True).merge(plate_ID_map_key)
keyfile.to_excel(DIR + 'revised_keyfile.xlsx', index = False )

guide_donor_plasmid_key = pd.read_excel(guide_donor_plasmid_keyfile)
trafo_ID_key = pd.read_excel(trafo_ID_keyfile)

keyfile.columns
# 191112_plate_1D_sample_1_collapsed.fastq
combined_keyfile = keyfile.merge(trafo_ID_key, on = 'sample_number').merge(guide_donor_plasmid_key, on = 'guide_number')
combined_keyfile['full_sample_name'] = '191112_' +  combined_keyfile['plate_ID'] + '_' + combined_keyfile['sample_number']
combined_keyfile.to_excel(DIR + 'combined_keyfile.xlsx', index = False )
# sample_info_fasta_keyfile = DIR +  'keyfile.xlsx' # argv[2]
# combined_keyfile.to_
combined_keyfile = combined_keyfile.sort_values(by=['sample_number'])
for column in combined_keyfile:
    print(column)
    
combined_keyfile.sample_number

for column in guide_donor_plasmid_key:
    print(column)
 
    
def parse_cigar(cigar):
    '''
    input: CIGAR string from SAM alignment
    output: list of 2 item tuples: mapped segment length, CIGAR operation
    '''
    total_cigar_operation_bp = {'N':0, 'X':0, 'I':0, 'D':0, 'S':0, '=':0, 'M':0 }
    segment_length = ''
    parsed_cigar = []
    for char in cigar:
        if char in '0123456789':
            segment_length += char
        else:
            operation = char
           # print(segment_length, cigar)
            try:
                int(segment_length)
            except:
                print(segment_length, cigar)
            mapped_segment_length = int(segment_length)
            segment_length = ''
            parsed_cigar.append( (mapped_segment_length, operation) )
            total_cigar_operation_bp[operation] += mapped_segment_length
    return parsed_cigar, total_cigar_operation_bp

def rev_comp(DNA):
    '''
    input: DNA string
    output: reverse complement string
    '''
    rev_comp_DNA = ''
    comp_bases = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', 'R':'Y', 'Y':'R', ' ':' ', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'n':'n', 'r':'y', 'y':'r', '':''}
    for base in DNA[::-1]:
        rev_comp_DNA += comp_bases[base]
    return rev_comp_DNA

## map reads and process cigar alignments to characterize donor edits and NHEJ indels
import math
for index, row in combined_keyfile.iterrows():
    print(index, row)
    guide_donor_plasmid = row["plasmid_ID"]
    sample_num = int(row["sample_number"].split('_')[1])
    if sample_num not in (1, 2) and sample_num not in list(range(9,22) ): # True: # row['plasmid_ID'] == 'pKR797': #  guide_donor_plasmid != 'pT518_pool_of_32' and row['plasmid_ID'] == 'pKR997': # and int(row["sample_number"][1:]) > 1156 : 
        sample_prefix = row["full_sample_name"]
        sample_filename = COLLAPSED_FASTQ_DIR + sample_prefix +  "_collapsed.fastq"
        print('sample_filename:', sample_filename)
        fasta_reference_filename = FASTA_REF_DIR + row["fasta_reference_filename"]
        sam = ALIGNMENTS_DIR + sample_prefix +  '_bbmap.sam'
        bbmap_cmd = 'bbmap.sh local=t in=' + sample_filename + ' out=' + sam + ' ref=' + fasta_reference_filename
        print(bbmap_cmd)
        print ('bbmap running for:', sample_prefix)
        try:
            bbmap_output = subprocess.check_output(bbmap_cmd, stderr=subprocess.STDOUT, shell=True).decode("utf-8")
               # print(bbmap_output)
        except subprocess.CalledProcessError as e: 
            print ("error>",e.output,'<')
        guide_seq = row["guide_sequence"]
        designed_donor = row["designed_donor"]
        fasta_sequences = SeqIO.parse(open(fasta_reference_filename),'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if 'WT' in name:
                WT_seq = sequence.upper()
        
        print ('processing alignments for:', sample_prefix)
        seq_to_counts = {}
        seq_to_ref_coord_cigar = {}
        
        total_counts = 0
        with open(sam, 'r') as infile:
            for line in infile:
                if line[0] != '@':
                    total_counts += 1
                    info = line.split('\t')
                    read_ID, flag, chromosome, coordinate, mapping_quality, cigar, equals_sign, paired_end_location, observed_template_length, seq, quality_scores = info[:11]
                    coordinate = int(coordinate)
                    counts = int(read_ID)
                    if cigar != '*':
                        total_counts += counts
                        seq_to_counts[seq] = counts
                        seq_to_ref_coord_cigar[seq] = chromosome, coordinate, cigar
            infile.close()
        
        sorted_seq = sorted(seq_to_counts.items(), key=operator.itemgetter(1), reverse=True)    
                
                
            #guide_seq = argv[2]
            #WT_target_seq = argv[3]    
                
            
        if guide_seq in WT_seq:
            guide_start_idx = WT_seq.index(guide_seq)
            guide_end_idx = guide_start_idx + 23
            PAM_idx = guide_start_idx + 20
            PAM_strand = '+'
        else:
            rc_guide_seq = rev_comp(guide_seq)
            guide_start_idx = WT_seq.index(rc_guide_seq) - 3
            guide_end_idx = guide_start_idx + 20
            PAM_idx = guide_start_idx + 3
            PAM_strand = '-'
        
        FLANKING_DONOR_SEQ_TO_EXAMINE = 2
        
        donor_left_portion = designed_donor[:20]
        if donor_left_portion in WT_seq:
            donor_start_idx = max( WT_seq.index(donor_left_portion) - FLANKING_DONOR_SEQ_TO_EXAMINE, 0 )
        elif rev_comp(donor_left_portion ) in WT_seq:
            donor_start_idx = max( WT_seq.index(rev_comp(donor_left_portion )) - FLANKING_DONOR_SEQ_TO_EXAMINE, 0 )
        else:
            donor_start_idx = 0
            
        
        donor_right_portion = designed_donor[-20:]
        if donor_right_portion in WT_seq:
            donor_end_idx = min( WT_seq.index(donor_right_portion) + FLANKING_DONOR_SEQ_TO_EXAMINE, len(WT_seq))
        elif rev_comp(donor_right_portion ) in WT_seq:
            donor_start_idx = min( WT_seq.index(rev_comp(donor_right_portion )) + FLANKING_DONOR_SEQ_TO_EXAMINE, len(WT_seq) )
        else:
            donor_start_idx = len(WT_seq)
              
        outfile = open(PROCESSED_COUNTS_DIR + sample_prefix + '_counts_cigar_seq.txt' , 'w')
            
        outfile.write('sample_name\tcounts\treference_sequence\tperfect_match\tcigar\tindel_in_guide_site\tindel_PAM_dist\tmismatch_in_guide_site\tclosest_mismatch_PAM_dist\tmismatch_in_donor_region\tfraction\tsoft_clipped_bp\tinserted_bp\tnum_inserted_residues\tinsertion_sites\tdeletion_sites\tseq\n')    
            
        for seq_counts in sorted_seq:
            seq, counts = seq_counts
            chromosome, coordinate, cigar = seq_to_ref_coord_cigar[seq]
            parsed_cigar, total_cigar_operation_bp = parse_cigar(cigar)
        #    print(cigar, parsed_cigar)
            seq_idx = 0
            ref_idx = coordinate
            inserted_bp = []
            insertion_sites = []
            deletion_sites = []
            indel_in_guide_site = False
            mismatch_in_guide_site = False
            perfect_match = True
            closest_indel_PAM_dist = 1000
            closest_mismatch_PAM_dist = 1000
            mismatch_in_donor_region = False
            total_soft_clipped_bp = 0
            for cigar_idx in range(len(parsed_cigar)):
        #        print(cigar_idx, parsed_cigar[cigar_idx]  )
                bp, operation = parsed_cigar[cigar_idx]
                if operation in 'ID':
                    perfect_match = False
                    for indel_idx in range(bp):
                        indel_PAM_dist = abs( PAM_idx - ref_idx + indel_idx )
                        if indel_PAM_dist < closest_indel_PAM_dist:
                            closest_indel_PAM_dist = indel_PAM_dist
                        if indel_PAM_dist < 20:
                            indel_in_guide_site = True
                if operation in 'IDXM':
                    perfect_match = False
                    for indel_mismatch_idx in range(bp):
                        indel_mismatch_coord_in_ref = ref_idx + indel_mismatch_idx
                        if donor_end_idx > indel_mismatch_coord_in_ref > donor_start_idx:
                            mismatch_in_donor_region = True
                if operation in 'S':
                    seq_idx += bp
                    total_soft_clipped_bp += bp
                elif operation in 'XM':
                    seq_idx += bp
                    ref_idx += bp
                    perfect_match = False
                    for mismatch_idx in range(bp):
                        mismatch_PAM_dist = abs( PAM_idx - ref_idx + mismatch_idx )
                        if mismatch_PAM_dist < closest_mismatch_PAM_dist:
                            closest_mismatch_PAM_dist = mismatch_PAM_dist
                        if closest_mismatch_PAM_dist < 20:
                            mismatch_in_guide_site = True
    #                        print(seq, coordinate, PAM_idx, ref_idx, cigar)
    #                        break
    #                        break
    #                        break
    #                        break
    #                        break
    #                        break 
                elif operation in '=':
                    seq_idx += bp
                    ref_idx += bp
                elif operation == 'I':
                    inserted_bp +=  [ seq[seq_idx:seq_idx+bp] ]
                    insertion_sites += [ref_idx]
                    seq_idx += bp          
                elif operation == 'D':
                    deletion_sites += [ref_idx]
                    ref_idx += bp
            fraction = float(counts/total_counts)
            
            output = [sample_prefix, counts, chromosome, perfect_match, cigar, indel_in_guide_site , closest_indel_PAM_dist , mismatch_in_guide_site, closest_mismatch_PAM_dist, mismatch_in_donor_region,   fraction, total_soft_clipped_bp,  inserted_bp, inserted_bp, insertion_sites, deletion_sites, seq] 
            output = [str(e) for e in output]    
            outfile.write( '\t'.join(output) + '\n')
        outfile.close()