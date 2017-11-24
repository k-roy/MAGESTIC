# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 16:20:58 2016

@author: kevinroy

Takes in a gene name and window in bp for each subpool, and mutates each codon to either all 20 top frequency codons per amino acid, or all 64 possible codons, with the option of excluding stop codons

-every amino acid change is compared to a control in which only the synonymous changes are made ***
-every change is conducted with synonymous changes downstream using a downstream NGG or CCN and with upstream changes using an upstream NGG or CCN  (guides with 6 or more T's are disallowed these are potent Pol III terminators)
-synonymous codons are spread towards the PAM until the alignment between the guide and donor results in a disruption score of 6 or greater

writes an outfile with the following fields: donor_name, category, oligo_sequence, guide_seq, full_donor_seq,  str(aa_num), str(oligo_length), str(disruption_score)
"""

import GTF_GFF_manipulation, bedgraph_computation, guide_donor_oligo_functions

REF_DIR = '/Users/kevinroy/Dropbox/'
DIR = '/Users/kevinroy/Dropbox/steinmetz_lab/CRISPR_saturation_genome_editing/'
STOP_CODONS_INCLUDED = True
MUTATE_INITIATING_METHIONINE = False
MUTATE_STOP_CODON = True
FOLD_OVERREPRESENTATION_OF_PSEUDO_WT_CONTROLS = 10

OLIGO_LENGTH = 210
GUIDE_LENGTH = 20
SUBPOOL_WINDOW_SIZE = 110

NUM_CONSECUTIVE_T_DISALLOWED = 6
MIN_GUIDE_DONOR_DISRUPTION_SCORE = 6
PAM = 'NGG'
rc_PAM = bedgraph_computation.rev_comp(PAM)
PAM_length = len(PAM)

END_RESTRICTION_SITES =  'GGCGCGCC', 'GCGGCCGC'
INTERNAL_RESTRICTION_SITE ='GCTCTTC'  ## , 'GGTCTC' BsaI site in AmpR ORF would need to be removedr
rc_INTERNAL_RESTRICTION_SITE = bedgraph_computation.rev_comp(INTERNAL_RESTRICTION_SITE)
INTERNAL_RESTRICTION_SITE_LENGTH = len(INTERNAL_RESTRICTION_SITE)
MINIMUM_HOMOLOGY = 20


## fwd priming sequence: GGACTTTggcgcgcc
FWD_PRIMING_SEQUENCE = 'GGACTTTggcgcgcc'.upper()
## internal cloning site will be BspQI_cloning_site = 'GTTTGAAGAGC', a backup would be 'GGTCTC'  GTTTAgagacc for regions that contain the BspQI site
INTERNAL_CLONING_SITE = 'gtttgaagagc'.upper() ## 'GTTTGAAGAGC'  ##  'gtttGAAGAGCGCTCTTCacga' 
## a future alternative would be the BsaI site: GGTCTC, which would give the following internal cloning site: gtttaGAGACC  (BsaI cuts 1 and 5 nt away from its recognition site)
INTERNAL_CLONING_SITE_INDEX = len(FWD_PRIMING_SEQUENCE) + GUIDE_LENGTH + 4
FWD_RESTRICTION_SITE_INDEX = 8

subpool_priming_seq_filename = 'subpool_priming_sequences.txt'
## reverse priming sequences previously validated (all should be the same length)
RIGHT_SUBPOOL_AMPLIFICATION_SEQUENCES  = []
with open(DIR + subpool_priming_seq_filename, 'r') as infile:
    for line in infile:
        RIGHT_SUBPOOL_AMPLIFICATION_SEQUENCES.append(line.strip())
infile.close()
REVERSE_PRIMING_SEQUENCES = [bedgraph_computation.rev_comp(e) for e in RIGHT_SUBPOOL_AMPLIFICATION_SEQUENCES]
STARTING_SUBPOOL_PRIMER_IDX = 0

DONOR_LENGTH = OLIGO_LENGTH - len(FWD_PRIMING_SEQUENCE) - GUIDE_LENGTH - len(INTERNAL_CLONING_SITE) - len(REVERSE_PRIMING_SEQUENCES[0])
MAX_VARIANT_LENGTH = DONOR_LENGTH - MINIMUM_HOMOLOGY*2

RESTRICTION_SITES_SCREENED_OUT = {}
RESTRICTION_SITES_SCREENED_OUT[INTERNAL_RESTRICTION_SITE] = 0
RESTRICTION_SITES_SCREENED_OUT[rc_INTERNAL_RESTRICTION_SITE] = 0
#for site in END_RESTRICTION_SITES:
#    RESTRICTION_SITES_SCREENED_OUT[site] = 0  
## for now don't screen out end restriction sites as Gibson may well be used for the first step


## Load codon table
infilename=REF_DIR + 'DNA_codon_to_AA.txt'
codon_to_aa = {}
aa_to_codon = {}
with open(infilename, 'r') as infile:
    for line in infile:
        codon, aa = line.strip().split()
        codon_to_aa[codon] = aa
        if aa in aa_to_codon:
            aa_to_codon[aa].append(codon)
        else:
            aa_to_codon[aa] = [codon]
infile.close()

## Load codon table with suboptimal codons removed
infilename=REF_DIR + 'DNA_codon_to_AA_suboptimal_codons_removed.txt'
## we removed all codons with equal to or less than 0.1 frequency as defined by 
## http://downloads.yeastgenome.org/unpublished_data/codon/ysc.gene.cod
suboptimal_removed_codon_to_aa = {}
suboptimal_removed_aa_to_codon = {}
with open(infilename, 'r') as infile:
    for line in infile:
        codon, aa = line.strip().split()
        suboptimal_removed_codon_to_aa[codon] = aa
        if aa in suboptimal_removed_aa_to_codon:
            suboptimal_removed_aa_to_codon[aa].append(codon)
        else:
            suboptimal_removed_aa_to_codon[aa] = [codon]
infile.close()


infilename=REF_DIR+ 'codon_frequencies_SGD_aa_abbreviated.txt'
codon_fraction_to_aa = {}
aa_to_codon_fraction = {}
codon_to_fraction = {}


with open(infilename, 'r') as infile:
    for line in infile:
        aa, codon, total_occurrences, frequency, fraction  = line.strip().split()
        fraction = float(fraction)
        codon_fraction_to_aa[(codon, fraction)] = aa
        codon_to_fraction[codon] = fraction
        if aa in aa_to_codon_fraction:
            aa_to_codon_fraction[aa].append( (fraction, codon) )
        else:
            aa_to_codon_fraction[aa] = [ (fraction, codon) ]
    for aa in aa_to_codon_fraction:
        aa_to_codon_fraction[aa] = sorted( aa_to_codon_fraction[aa], reverse = True   )
infile.close()
print(aa_to_codon_fraction)

highest_frequency_codons = []
for aa in aa_to_codon_fraction:
    info = aa_to_codon_fraction[aa]
    highest_frequency, highest_frequency_codon = info[0]
    highest_frequency_codons.append(highest_frequency_codon)
    
synonymous_codon_preferences = {}

for codon in codon_to_aa:
    order_of_preference = []
    aa = codon_to_aa[codon]
    fraction = codon_to_fraction[codon]
    synonymous_info = aa_to_codon_fraction[aa]
    for idx in range(len(synonymous_info)):
        if synonymous_info[idx][1] == codon:
            codon_idx = idx
    for idx in range(codon_idx + 1,len(synonymous_info) ):
        order_of_preference.append( synonymous_info[idx][1]  )
    for idx in range(codon_idx + 1):
        order_of_preference.append( synonymous_info[idx][1]  )
    synonymous_codon_preferences[codon] = order_of_preference
print(synonymous_codon_preferences)


## gff file format: chrI	SGD	CDS	335	649	.	+	0	Parent=YAL069W_mRNA;Name=YAL069W_CDS;orf_classification=Dubious
gff_filename = REF_DIR + 'saccharomyces_cerevisiae_annotations.gff'

genome_filename = REF_DIR + 'saccharomyces_cerevisiae_sequence.fasta'
## Load reference sequence
try:
    genome_seq
except:
    genome_seq = bedgraph_computation.load_genome(genome_filename)
    
## load ORF coordinates    
ORF_dict= GTF_GFF_manipulation.load_ORF_coordinates(gff_filename)


def assign_aa_num_to_exons(ORF_strand, ORF_exon_coords, ORF_seq):
    '''
    takes the strand and exon coords
    returns a tuple of tuples [(1,2), (2,3,4,5,6)], where each tuple represents an exon, and the integers inside represent aa_nums in that exon
    '''
    current_exon_num = 1
    aa_num_to_exon = {}
    aa_length = len(ORF_seq)
    previous_aa_end_coord = None
    for aa_num in range(1, aa_length + 1):
        current_aa_coords = ORF_exon_coords[(aa_num - 1)*3 : aa_num*3 ]
        if previous_aa_end_coord != None and abs( current_aa_coords[0] - previous_aa_end_coord ) != 1: ##exon found in between amino acids
            current_exon_num += 1
            aa_num_to_exon[aa_num] = [current_exon_num]
        else:
            aa_num_to_exon[aa_num] = [current_exon_num]
        if abs( current_aa_coords[2] - current_aa_coords[0] ) != 2:  ##aa is split across exons
            aa_num_to_exon[aa_num] = [current_exon_num, current_exon_num + 1]
            current_exon_num += 1
        previous_aa_end_coord = current_aa_coords[2]
    return aa_num_to_exon

def assemble_mutated_donor_sequence(ORF_strand, left_donor, right_donor, temp_ORF_seq, aa_num, synonymous_codons_to_mutate, upstream):
    aa_idx = aa_num - 1
    if upstream:
        mutated_region = ''.join( temp_ORF_seq[aa_idx - synonymous_codons_to_mutate: aa_idx + 1] )
    else:
        mutated_region = ''.join( temp_ORF_seq[aa_idx: aa_idx + synonymous_codons_to_mutate + 1] )
    if ORF_strand == '-':
        mutated_region = bedgraph_computation.rev_comp(mutated_region)
    if len(mutated_region) % 3 != 0:
        print(temp_ORF_seq[aa_idx - synonymous_codons_to_mutate: aa_idx + 1])
    query_donor = left_donor + mutated_region.lower() + right_donor
    return query_donor
    
def generate_donors(left_donor, right_donor, ORF, ORF_strand, synonymous_codons_to_mutate, aa_num, codon_to_aa, ORF_seq, temp_ORF_seq, ORF_info, ORF_chrom, ref_allele_start_coord, upstream):
    '''
    generates a vcf-format  donor dictionary with the form: 
    donor_library [ YDR025C_Y4F_variant_ACG4GAC_C5C,G6G,H7H,I8I_CCG5CCA,... ] =  codon_pool, chrVII,	18036,   ACG..., 	GAC...	
    donor_library [ YDR025C_Y4Y_control_ACG4ACG_C5C,G6G,H7H,I8I_CCG5CCA,... ] =  codon_pool, chrVII,	18036,   ACG...,	ACG...	
    '''
    aa_idx = aa_num - 1
    donor_library = {}
    synonymous_aa = []
    synonymous_codons = []
    WT_target_codon, WT_target_aa_coords, WT_target_aa = ORF_info[aa_num] 
    if upstream:
        for upstream_aa_num in range(aa_num - synonymous_codons_to_mutate, aa_num):
             upstream_aa_idx = upstream_aa_num - 1
             codon, aa_coords, aa = ORF_info[upstream_aa_num] 
             synonymous_codon = temp_ORF_seq[upstream_aa_idx]
             synonymous_codons.append(codon + str(upstream_aa_num) + synonymous_codon)
             synonymous_aa.append(aa + str(upstream_aa_num) + aa)
    else:
         for downstream_aa_num in range(aa_num + 1, aa_num  + synonymous_codons_to_mutate + 1):
             downstream_aa_idx = downstream_aa_num - 1
             codon, aa_coords, aa = ORF_info[downstream_aa_num] 
             synonymous_codon = temp_ORF_seq[downstream_aa_idx]
             synonymous_codons.append(codon + str(downstream_aa_num) + synonymous_codon)
             synonymous_aa.append(aa + str(downstream_aa_num) + aa)
    synonymous_codons_name = ','.join(synonymous_codons) 
    synonymous_aa_name = ','.join(synonymous_aa)
    ## generate donor and donor names
    for mutant_target_codon in codon_to_aa: 
        temp_ORF_seq_for_target_codon = temp_ORF_seq[:]
        temp_ORF_seq_for_target_codon[aa_idx] = mutant_target_codon
        if upstream:
            mutated_region = ''.join( temp_ORF_seq_for_target_codon[aa_idx - synonymous_codons_to_mutate: aa_idx + 1] )
            WT_region = ''.join( ORF_seq[aa_idx - synonymous_codons_to_mutate: aa_idx + 1] )
        else:
            mutated_region = ''.join( temp_ORF_seq_for_target_codon[aa_idx: aa_idx + synonymous_codons_to_mutate + 1] )
            WT_region = ''.join( ORF_seq[aa_idx: aa_idx + synonymous_codons_to_mutate + 1])
        if ORF_strand == '-':
            mutated_region = bedgraph_computation.rev_comp(mutated_region)
            WT_region = bedgraph_computation.rev_comp(WT_region)
        if mutant_target_codon ==  WT_target_codon:
            donor_type = 'control'
        else:
            donor_type = 'variant'
        aa_change = WT_target_aa + str(aa_num) + codon_to_aa[mutant_target_codon]
        codon_change = WT_target_codon + str(aa_num) + mutant_target_codon
        if upstream:
            spreading_direction = 'upstream_synonymous_changes:'
        else:
            spreading_direction = 'downstream_synonymous_changes:'
        donor_name = '_'.join([ORF, aa_change, donor_type, codon_change,spreading_direction, synonymous_aa_name, synonymous_codons_name] )
        full_donor_seq = assemble_mutated_donor_sequence(ORF_strand, left_donor, right_donor, temp_ORF_seq_for_target_codon, aa_num, synonymous_codons_to_mutate, upstream)
        vcf_fields = (ORF_chrom, ref_allele_start_coord, WT_region, mutated_region, donor_name)
        if donor_type == 'variant':
            if mutant_target_codon in highest_frequency_codons:
                codon_pool = 0
            else:
                codon_pool = 1
            donor_library[donor_name] = codon_pool, full_donor_seq, vcf_fields
        else:
            donor_library[donor_name] = 0, full_donor_seq, vcf_fields
            donor_library[donor_name] = 1, full_donor_seq, vcf_fields
    return donor_library

def combine_donor_libraries(donor_library_1, donor_library_2):
    combined_library = {}
    for donor_name in donor_library_1:
        combined_library[donor_name] = donor_library_1[donor_name]
    for donor_name in donor_library_2:
        combined_library[donor_name] = donor_library_2[donor_name]
    return combined_library
            
def mutate_synonymous_codons(ORF_chrom, ORF_strand, synonymous_codons_to_mutate, aa_num, aa_to_codon, suboptimal_removed_aa_to_codon, ORF_info, ORF_seq, upstream):
    '''
    takes an aa_num and the number of codons to mutate upstream or downstream (designated by upstream = False)
    returns plus_strand_coord start pos of change, ref, alt alleles, and the name of ORF and the synonymous changes
    does not mutate the actual target amino acid

    ## also checks synonymous donor for inadvertent introduction of BspQI sites, and changes to other synonymous codons if possible
    
    ORF_info[aa_num] = ( codon, aa_coords, aa )
    returns a list of donors and donor_infos
    
    chrVII	18036	   ACG...	GAC...	YDR025C_(C5C,G6G,H7H,I8I)_[CCG5CCA,...]_{Y,ACG,4,F,GAC})_variant_subpool_1
    '''
    aa_idx = aa_num - 1
    if ORF_strand == '+':
        if upstream:
            ref_allele_start_coord = ORF_info[aa_num - synonymous_codons_to_mutate][1][0]
        else:
            ref_allele_start_coord = ORF_info[aa_num][1][0]
    else:
        if upstream:
            ref_allele_start_coord = ORF_info[aa_num][1][2]
        else:
            ref_allele_start_coord = ORF_info[aa_num + synonymous_codons_to_mutate][1][2]
    mutated_region_length = 3 + synonymous_codons_to_mutate*3  ## includes the amino acid to be mutated
    homologous_arm_length = DONOR_LENGTH - mutated_region_length
    left_donor_length = homologous_arm_length // 2
    right_donor_length = homologous_arm_length - left_donor_length
    left_donor = genome_seq[ORF_chrom][ref_allele_start_coord - left_donor_length: ref_allele_start_coord]
    right_donor = genome_seq[ORF_chrom][ref_allele_start_coord + mutated_region_length : ref_allele_start_coord + mutated_region_length + right_donor_length]
    extended_left_donor = genome_seq[ORF_chrom][ref_allele_start_coord - left_donor_length - 200: ref_allele_start_coord]
    extended_right_donor = genome_seq[ORF_chrom][ref_allele_start_coord + mutated_region_length : ref_allele_start_coord + mutated_region_length + right_donor_length + 200]
    temp_ORF_seq = ORF_seq[:]
    for idx in range(1, synonymous_codons_to_mutate + 1):
        if upstream:
            current_aa_num = aa_num - idx 
            current_aa_idx = current_aa_num - 1
        else:
            current_aa_num = aa_num + idx 
            current_aa_idx = current_aa_num - 1
        codon, aa_coords, aa = ORF_info[current_aa_num]
        synonymous_codon_with_largest_hamming_dist = codon
        largest_hamming_dist_for_this_aa = 0
        for other_codon in suboptimal_removed_aa_to_codon[aa]:
            temp_ORF_seq[current_aa_idx] = other_codon
            hamming_distance_for_this_codon = bedgraph_computation.hamming_distance(other_codon, codon)
            query_donor = assemble_mutated_donor_sequence(ORF_strand, left_donor, right_donor, temp_ORF_seq, aa_num, synonymous_codons_to_mutate, upstream)
            ## checking for internal restriction sites with synonymoous changes
            if hamming_distance_for_this_codon > largest_hamming_dist_for_this_aa and INTERNAL_RESTRICTION_SITE not in query_donor and rc_INTERNAL_RESTRICTION_SITE not in query_donor:
                temp_ORF_seq_check_all_target_codons = temp_ORF_seq[:]
                all_target_codons_clear_of_BspQI_introduction = True
                for target_codon in codon_to_aa:
                    temp_ORF_seq_check_all_target_codons[aa_idx] = target_codon
                    check_all_target_codons_query_donor = assemble_mutated_donor_sequence(ORF_strand, left_donor[-(INTERNAL_RESTRICTION_SITE_LENGTH- 4):], right_donor[: INTERNAL_RESTRICTION_SITE_LENGTH - 4], temp_ORF_seq_check_all_target_codons, aa_num, synonymous_codons_to_mutate, upstream)
                    ## if possible, remove synonymous codon changes that would generate a restriction site with the target codon and proximal upstream/downstream sequence
                    if INTERNAL_RESTRICTION_SITE in check_all_target_codons_query_donor or rc_INTERNAL_RESTRICTION_SITE in check_all_target_codons_query_donor:
                        all_target_codons_clear_of_BspQI_introduction = False
                if not all_target_codons_clear_of_BspQI_introduction:
                    print('restriction site prevented')
                if all_target_codons_clear_of_BspQI_introduction:
                    largest_hamming_dist_for_this_aa = hamming_distance_for_this_codon
                    synonymous_codon_with_largest_hamming_dist = other_codon
        temp_ORF_seq[current_aa_idx] = synonymous_codon_with_largest_hamming_dist
    control_donor = assemble_mutated_donor_sequence(ORF_strand, extended_left_donor, extended_right_donor, temp_ORF_seq, aa_num, synonymous_codons_to_mutate, upstream)
    if upstream:
        mutated_region = ''.join( temp_ORF_seq[aa_idx - synonymous_codons_to_mutate: aa_idx + 1] )
    else:
        mutated_region = ''.join( temp_ORF_seq[aa_idx : aa_idx + synonymous_codons_to_mutate + 1] )
    if ORF_strand == '-':
        mutated_region = bedgraph_computation.rev_comp(mutated_region)
    mutated_region_chromosomal_range =   ref_allele_start_coord,  ref_allele_start_coord + mutated_region_length
    donor_library = generate_donors(left_donor, right_donor,  ORF, ORF_strand, synonymous_codons_to_mutate, aa_num, codon_to_aa, ORF_seq, temp_ORF_seq, ORF_info, ORF_chrom, ref_allele_start_coord, upstream)
    return control_donor, donor_library, mutated_region_chromosomal_range

def get_proximal_sense_and_antisense_PAMs(genome_seq, ORF_chrom, ORF_strand, ORF_exon_coords, coord, upstream):
    '''
    takes a chromosome, strand, and coordinate
    
    for minus strand, coord is last coord of codon, for plus strand it is the first coord of codon    
    finds 5 PAMs for each strand
    returns the PAM_coordinate, guide, and guide with full PAM sequence (outlaws guides with excessive T stretch as defined by NUM_CONSECUTIVE_T_DISALLOWED )
    '''
    if (ORF_strand == '+' and upstream) or (ORF_strand == '-' and not upstream):
        plus_strand_coord = coord + GUIDE_LENGTH
        minus_strand_coord = coord + GUIDE_LENGTH
    else:
        plus_strand_coord = coord - GUIDE_LENGTH
        minus_strand_coord = coord - GUIDE_LENGTH
        
    plus_strand_PAMs = []
    shift = - GUIDE_LENGTH
    while len(plus_strand_PAMs) < 10 and shift < DONOR_LENGTH - MINIMUM_HOMOLOGY*2 :
        query_seq = genome_seq[ORF_chrom][plus_strand_coord: plus_strand_coord + PAM_length]
        plus_strand_PAM_guide_seq = genome_seq[ORF_chrom][plus_strand_coord - GUIDE_LENGTH: plus_strand_coord ]
        plus_strand_PAM_guide_seq_with_PAM = genome_seq[ORF_chrom][plus_strand_coord - GUIDE_LENGTH: plus_strand_coord + PAM_length ]
        if bedgraph_computation.hamming_distance( query_seq, PAM ) == 0 and NUM_CONSECUTIVE_T_DISALLOWED*'T' not in  plus_strand_PAM_guide_seq:
            plus_strand_PAMs.append( ('+', query_seq, plus_strand_coord, plus_strand_PAM_guide_seq, plus_strand_PAM_guide_seq_with_PAM)   )
        if ORF_strand == '+':
            if upstream:
                plus_strand_coord -= 1
            else:
                plus_strand_coord += 1
        else:
            if upstream:
                plus_strand_coord += 1
            else:
                plus_strand_coord -= 1
        shift += 1
    shift = - GUIDE_LENGTH
    minus_strand_PAMs = []
    
    while len(minus_strand_PAMs) < 10 and shift < DONOR_LENGTH - MINIMUM_HOMOLOGY*2 :
        query_seq = bedgraph_computation.rev_comp( genome_seq[ORF_chrom][minus_strand_coord: minus_strand_coord + PAM_length] )
        minus_strand_PAM_guide_seq = bedgraph_computation.rev_comp( genome_seq[ORF_chrom][minus_strand_coord + PAM_length : minus_strand_coord + PAM_length + GUIDE_LENGTH ] )
        minus_strand_PAM_guide_seq_with_PAM = bedgraph_computation.rev_comp( genome_seq[ORF_chrom][minus_strand_coord : minus_strand_coord + PAM_length + GUIDE_LENGTH ] )
        if bedgraph_computation.hamming_distance( query_seq, PAM ) == 0 and NUM_CONSECUTIVE_T_DISALLOWED*'T' not in minus_strand_PAM_guide_seq:
            minus_strand_PAMs.append ( ('-', query_seq, minus_strand_coord, minus_strand_PAM_guide_seq, minus_strand_PAM_guide_seq_with_PAM) )
        if ORF_strand == '+': 
            if upstream:
                minus_strand_coord -= 1
            else:
                minus_strand_coord += 1
        else:
            if upstream:
                minus_strand_coord += 1
            else:
                minus_strand_coord -= 1
        shift += 1
    return plus_strand_PAMs + minus_strand_PAMs

def generate_guide_donor_oligos_with_synonymous_variants(aa_num_to_exon, ORF_info, aa_num, codon_coords, ORF_chrom, ORF_strand, ORF_exon_coords, ORF_seq, upstream):
    '''
    takes an ORF, its chromosome and strand, and exon coords ( a list of coords of each nt in the ORF, start to stop), the ORF seq, the aa to mutate and the codon idx for that aa
    searches for the nearest upstream CCN and NGG
    mutates upstream codons with synonymous changes (picking codon involving largest hamming distance change) until a disruption score of 6 is achieved
    disruption scores are calculated by considering that mismatches in the 10-20 nt region count as 1, in the seed region count as 2, and in the PAM 3
    
    determine appropriate subpool based on window and codon list
    '''
    current_exon = aa_num_to_exon[aa_num][0]
    if ORF_strand == '+':
        coord = codon_coords[0]
    else:
        coord = codon_coords[2]
    query_PAMs = get_proximal_sense_and_antisense_PAMs(genome_seq, ORF_chrom, ORF_strand, ORF_exon_coords, coord, upstream)
    ## print(query_PAMs)
    synonymous_codons_to_mutate = 0
   
    control_donor, donor_library, mutated_region_chromosomal_range = mutate_synonymous_codons(ORF_chrom, ORF_strand, synonymous_codons_to_mutate, aa_num, aa_to_codon, suboptimal_removed_aa_to_codon, ORF_info, ORF_seq, upstream)
    
    best_disruption_score = -1
    output = None
    for PAM_info in query_PAMs: 
        PAM_strand, PAM_seq, PAM_coord, guide_seq, guide_seq_with_PAM = PAM_info
        disruption_score = guide_donor_oligo_functions.get_disruption_score(control_donor, guide_seq_with_PAM, PAM)
        if disruption_score > best_disruption_score:
            best_disruption_score = disruption_score
            output = guide_seq, PAM_strand, PAM_coord, guide_seq_with_PAM, control_donor, donor_library, mutated_region_chromosomal_range, disruption_score
            # print(best_score, num_exons_in_proximal_codon, proximal_codon_exon, current_exon )
            
    synonymous_codons_to_mutate += 1
    if (upstream and aa_num - synonymous_codons_to_mutate <= 1) or (not upstream and aa_num + synonymous_codons_to_mutate > len(ORF_seq) ):
        return output
    if upstream:
        proximal_codon_exon = aa_num_to_exon[aa_num - synonymous_codons_to_mutate][0]
        num_exons_in_proximal_codon = len(aa_num_to_exon[aa_num - synonymous_codons_to_mutate])
    else:
        proximal_codon_exon = aa_num_to_exon[aa_num + synonymous_codons_to_mutate][0]  
        num_exons_in_proximal_codon = len(aa_num_to_exon[aa_num + synonymous_codons_to_mutate])
 
    while best_disruption_score < MIN_GUIDE_DONOR_DISRUPTION_SCORE and num_exons_in_proximal_codon == 1 and proximal_codon_exon == current_exon:

        ## spread 1 synonymous change (upstream or downstream)
        control_donor, donor_library, mutated_region_chromosomal_range = mutate_synonymous_codons(ORF_chrom, ORF_strand, synonymous_codons_to_mutate, aa_num, aa_to_codon, suboptimal_removed_aa_to_codon, ORF_info, ORF_seq, upstream)
        for PAM_info in query_PAMs: 
            PAM_strand, PAM_seq, PAM_coord, guide_seq, guide_seq_with_PAM = PAM_info
            disruption_score = guide_donor_oligo_functions.get_disruption_score(control_donor, guide_seq_with_PAM, PAM)
            if disruption_score > best_disruption_score:
                best_disruption_score = disruption_score
                output = guide_seq, PAM_strand, PAM_coord, guide_seq_with_PAM, control_donor, donor_library, mutated_region_chromosomal_range, disruption_score
        synonymous_codons_to_mutate += 1
        ## cannot spread synonymous changes to the M1 or beyond stop codon
        if (upstream and aa_num - synonymous_codons_to_mutate <= 1) or (not upstream and aa_num + synonymous_codons_to_mutate > len(ORF_seq) ):
            print('spread beyond ORF')
            break
        else:
            if upstream:
                proximal_codon_exon = aa_num_to_exon[aa_num - synonymous_codons_to_mutate][0]
                num_exons_in_proximal_codon = len(aa_num_to_exon[aa_num - synonymous_codons_to_mutate])
            else:
                proximal_codon_exon = aa_num_to_exon[aa_num + synonymous_codons_to_mutate][0]  
                num_exons_in_proximal_codon = len(aa_num_to_exon[aa_num + synonymous_codons_to_mutate])

    if output == None:
        print('no guides found for aa num:', aa_num, query_PAMs)
    return output
    
def assign_mutated_chromosomal_ranges_to_subpools(mutated_region_chromosomal_ranges, ORF_exon_coords, ORF_strand):
    '''
    takes the chromosomal ranges mutated by the donors, and the global variable SUBPOOL_WINDOW_SIZE
    returns a dictionary of chromosomal ranges to subpool nums
    '''
    sorted_chromosomal_ranges  = sorted( mutated_region_chromosomal_ranges[:] )
    chromosomal_range_to_subpool_num_and_window = {}
    subpool_window_to_subpool_num = {}
    current_subpool_num = 0
    current_subpool_window = ORF_exon_coords[0], ORF_exon_coords[0] + SUBPOOL_WINDOW_SIZE ## move from lowest coord to highest coord
    subpool_window_to_subpool_num[current_subpool_window] = current_subpool_num
    print(sorted_chromosomal_ranges )
    for chromosomal_range in sorted_chromosomal_ranges:
        if current_subpool_window[1] >= chromosomal_range[1] >= current_subpool_window[0] and current_subpool_window[1] >= chromosomal_range[0] >= current_subpool_window[0]:  ## check if range is completely contained in existing window
            chromosomal_range_to_subpool_num_and_window[ chromosomal_range ] = current_subpool_num, current_subpool_window
        else:
            current_subpool_num += 1
            current_subpool_window = chromosomal_range[0], chromosomal_range[0] + SUBPOOL_WINDOW_SIZE
            subpool_window_to_subpool_num [current_subpool_window] = current_subpool_num
            if current_subpool_window[1] >= chromosomal_range[1] >= current_subpool_window[0] and current_subpool_window[1] >= chromosomal_range[0] >= current_subpool_window[0]:
                chromosomal_range_to_subpool_num_and_window[ chromosomal_range ] = current_subpool_num, current_subpool_window
            else:
                print(current_subpool_window, chromosomal_range)
                print('chromosomal range exceeds subpool window size')
    return chromosomal_range_to_subpool_num_and_window, subpool_window_to_subpool_num

def combine_guide_info_with_donor_library(donor_library, guide_seq, PAM_strand, PAM_coord, guide_seq_with_PAM, mutated_region_chromosomal_range, best_score):
    guide_donor_library = {}
    for donor_name in donor_library:
        codon_pool, full_donor_seq, vcf_fields = donor_library[donor_name]
        guide_donor_library[donor_name] = codon_pool, mutated_region_chromosomal_range, full_donor_seq, vcf_fields, guide_seq, PAM_strand, PAM_coord, guide_seq_with_PAM[GUIDE_LENGTH:], best_score
    return guide_donor_library

def generate_guide_donor_library_for_aa_saturation_editing(ORF):
    '''
    takes an ORF name, sequence, window size for subpools, and a codon to aa dictionary
    writes variants to a vcf format file
    '''
    ## load ORF info for a given ORF,
    mutated_region_chromosomal_ranges = []
    ORF_chrom, ORF_strand, ORF_exon_coords, ORF_seq = GTF_GFF_manipulation.get_ORF_seq_and_CDS_coords(ORF, ORF_dict, genome_seq)
    ORF_info =  GTF_GFF_manipulation.get_aa_num_to_codon_coords_aa(ORF_chrom, ORF_strand, ORF_exon_coords, ORF_seq, codon_to_aa)
    aa_num_to_exon =  assign_aa_num_to_exons(ORF_strand, ORF_exon_coords, ORF_seq)
    ## ORF_info[aa_num] = ( codon, aa_coords, aa )
    combined_guide_donor_library = {}
    aa_num = 2
    codon_idx = 3
    if MUTATE_INITIATING_METHIONINE:
        aa_num = 1
        codon_idx = 0
    if MUTATE_STOP_CODON:
        codons_to_mutate = ORF_seq[aa_num-1:]
    else:
        codons_to_mutate = ORF_seq[aa_num-1:-1]
    for codon in codons_to_mutate: ## [:10]
        if aa_num % 10 == 0:      
            print('mutating amino acid number: ',aa_num )
        codon_coords = ORF_exon_coords[codon_idx: codon_idx + 3]  
        ## only mutate codons that are not split across introns
        if abs(codon_coords[2] - codon_coords[0]) == 2: ## don't mutate codons that are split across exons
            for upstream in True, False:
                output = generate_guide_donor_oligos_with_synonymous_variants(aa_num_to_exon, ORF_info, aa_num, codon_coords, ORF_chrom, ORF_strand, ORF_exon_coords, ORF_seq, upstream)
                if output != None:
                    guide_seq, PAM_strand, PAM_coord, guide_seq_with_PAM, control_donor, donor_library, mutated_region_chromosomal_range, best_score = output
                    guide_donor_library = combine_guide_info_with_donor_library(donor_library, guide_seq, PAM_strand, PAM_coord, guide_seq_with_PAM, mutated_region_chromosomal_range, best_score)             
                    combined_guide_donor_library = combine_donor_libraries(combined_guide_donor_library, guide_donor_library)
                    mutated_region_chromosomal_ranges.append(mutated_region_chromosomal_range)
        else:
            print('split codon')
        aa_num += 1
        codon_idx += 3
    chromosomal_range_to_subpool_num_and_window, subpool_window_to_subpool_num = assign_mutated_chromosomal_ranges_to_subpools(mutated_region_chromosomal_ranges, ORF_exon_coords, ORF_strand)
    return combined_guide_donor_library, chromosomal_range_to_subpool_num_and_window, subpool_window_to_subpool_num, ORF_chrom, ORF_strand 

                
def generate_guide_donor_oligos_for_aa_saturation_editing(ORF):
    '''
    takes the ORF to be mutated, and uses global variables FWD_PRIMING_SEQUENCE,  INTERNAL_CLONING_SITE, REVERSE_PRIMING_SEQUENCES and STARTING_SUBPOOL_PRIMER_IDX
    first examines the subpool row as directed by the codon frequency, and then the subpool column as directed by  chromosomal_range_to_subpool_num_and_window
    returns a dictionary of guide_donor_oligos [ donor_name ] = category, oligo_sequence, vcf_fields, subpool_primer_idx
    '''
    combined_guide_donor_library, chromosomal_range_to_subpool_num_and_window, subpool_window_to_subpool_num, ORF_chrom, ORF_strand  = generate_guide_donor_library_for_aa_saturation_editing(ORF)
    total_subpool_columns = len(subpool_window_to_subpool_num)
    guide_donor_oligos = {}
    for donor_name in combined_guide_donor_library:
        subpool_row_num, mutated_region_chromosomal_range, full_donor_seq, vcf_fields, guide_seq, PAM_strand, PAM_coord, PAM_seq, best_score = combined_guide_donor_library[donor_name]
        disruption_score = guide_donor_oligo_functions.get_disruption_score(full_donor_seq, guide_seq + PAM_seq, PAM)
        donor_guide_name = '_'.join( [donor_name, 'PAM', PAM_strand, ORF_chrom, str(PAM_coord), PAM_seq, 'disruption_score', str(disruption_score)] )
        subpool_column_num, subpool_window = chromosomal_range_to_subpool_num_and_window[mutated_region_chromosomal_range]
        subpool_primer_idx = subpool_row_num * total_subpool_columns + STARTING_SUBPOOL_PRIMER_IDX + subpool_column_num
        #print('subpool_row_num', subpool_row_num,'subpool_column_num',  subpool_column_num, 'total_subpool_columns', total_subpool_columns, )
        oligo_sequence = FWD_PRIMING_SEQUENCE.lower() + guide_seq.upper() + INTERNAL_CLONING_SITE.lower() + full_donor_seq + REVERSE_PRIMING_SEQUENCES[ subpool_primer_idx ]
        category = str(subpool_primer_idx) + ':' + '_'.join([ORF_strand, ORF_chrom, 'mutated_region', str(mutated_region_chromosomal_range), 'phenotyping_window', str(subpool_window) ] )
        guide_donor_oligos[ donor_guide_name ] = category, oligo_sequence, vcf_fields, subpool_primer_idx, guide_seq, PAM_seq, full_donor_seq, best_score, disruption_score
    return guide_donor_oligos

ORF= 'YMR079W' ## 
guide_donor_oligos = generate_guide_donor_oligos_for_aa_saturation_editing(ORF)
oligo_library_outfilename = DIR + 'guide_donor_oligos_for_' + ORF + '_codon_variants_with_synonymous_codon_spreading_' + str(OLIGO_LENGTH) + '_bp_oligos_with_' + str(SUBPOOL_WINDOW_SIZE) + '_bp_subpooled_windows.tab'

#def write_vcf(vcf_outfilename, )

def write_guide_donor_oligos(oligo_library_outfilename, guide_donor_oligos):
    with open(oligo_library_outfilename, 'w') as outfile:
        header = 'donor name (ORF, mutated aa, variant or control, mutated codon, synonymous changes upstream or downstream, PAM strand, coordinates, and sequence, disruption score of donor against guide and PAM)	subpool_priming_sequence_number, chromosomal coordinates of mutated region,chromosomal coordinates of subpool phenotyping window	full oligo sequence (fwd priming sequence (lower case), guide sequence (upper case), internal cloning site (lower case), left homology of donor (upper case), mutated region (lower case), right homology of donor (upper case), rev subpool priming sequence (lower case))	guide sequence	donor sequence	amino acid number\toligo length\tdisruption score' + '\n'
        outfile.write(header) 
        guide_donors = []
        for donor_name in guide_donor_oligos:
            info = donor_name.split('_')
            aa_change = info[1]
            WT_aa, aa_num, mut_aa = aa_change[0], int(aa_change[1:-1]), aa_change[-1]
            category, oligo_sequence, vcf_fields, subpool_primer_idx, guide_seq, PAM_seq, full_donor_seq, best_score, disruption_score = guide_donor_oligos[donor_name]
            guide_donors.append( (aa_num, donor_name, category, oligo_sequence, vcf_fields, subpool_primer_idx, guide_seq, PAM_seq, full_donor_seq, best_score, disruption_score) )
        sorted_guide_donors = sorted(guide_donors)
        for oligo_info in sorted_guide_donors:    
            aa_num, donor_name, category, oligo_sequence, vcf_fields, subpool_primer_idx, guide_seq, PAM_seq, full_donor_seq, best_score, disruption_score = oligo_info
            oligo_length = len(oligo_sequence)  
            if oligo_length != OLIGO_LENGTH:
                print('WARNING: oligo_length != OLIGO_LENGTH: ', oligo_info )
            if 'control' in donor_name:
                for idx in range(FOLD_OVERREPRESENTATION_OF_PSEUDO_WT_CONTROLS):
                    output = '\t'.join([donor_name, category, oligo_sequence, guide_seq, full_donor_seq, str(aa_num), str(oligo_length), str(disruption_score) ]) + '\n'
                    outfile.write(output)
            else:
                output = '\t'.join([donor_name, category, oligo_sequence, guide_seq, full_donor_seq,  str(aa_num), str(oligo_length), str(disruption_score) ]) + '\n'
                outfile.write(output)
                
    outfile.close()

write_guide_donor_oligos(oligo_library_outfilename, guide_donor_oligos)     

