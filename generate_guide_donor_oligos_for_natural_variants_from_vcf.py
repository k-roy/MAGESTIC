# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 19:39:28 2016

@author: kevinroy
"""
import  numpy

DIR = '/Users/kevinroy/Dropbox/steinmetz_lab/Scer_natural_variants/'
genome_filename = '/Users/kevinroy/Dropbox/saccharomyces_cerevisiae_sequence.fasta'

OLIGO_LENGTH = 170
MAX_DIST_TO_PAM = 20
MAX_GUIDES_PER_DONOR = 300
GUIDE_LENGTH = 20

strain = 'rm11' # 'sk1' #   'RM11-1A' # 'SK1' #  'SK1' # 'RM11-1A'

MAX_BP_BETWEEN_VARIANTS_TO_COMBINE = 5
linked_variant_outfilename = DIR + strain + '_' + str(MAX_BP_BETWEEN_VARIANTS_TO_COMBINE) + 'bp_linked_SICtools.vcf'
SNP_variant_outfilename = DIR + strain + '_SNPs_SICtools.vcf'
indel_variant_outfilename = DIR + strain + '_indels_SICtools.vcf'

### all guide_donors can be generated from a common file format (.vcf)
indel_vcf_infilename =  DIR + strain + '_indels_WAVE.vcf' ## '/Users/kevinroy/Dropbox/steinmetz_lab/variants_Wave_LinGEN/rm_indel_test.txt' ## 
SNP_vcf_infilename = DIR + strain + '_snps_WAVE.vcf' ## '/Users/kevinroy/Dropbox/steinmetz_lab/variants_Wave_LinGEN/rm11_snps_test.txt'   ## 


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
    
def generate_empty_nonstranded_bedgraph_dict():
    '''
    takes no arguments
    returns a dictionary of strands of chromosomes of empty dictionaries
    '''
    dictionary = {}
    chromosome_numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI'.split(',')
    for chromosome_numeral in chromosome_numerals:
        dictionary['chr' + chromosome_numeral] = {}
    return dictionary

NUM_CONSECUTIVE_T_DISALLOWED = 6
PAM = 'NGG'
END_RESTRICTION_SITES =  'GGCGCGCC', 'GCGGCCGC'
INTERNAL_RESTRICTION_SITE ='GCTCTTC'
rc_INTERNAL_RESTRICTION_SITE = rev_comp(INTERNAL_RESTRICTION_SITE)
MINIMUM_HOMOLOGY = 20


## fwd priming sequence: GGACTTTggcgcgcc
FWD_PRIMING_SEQUENCE = 'GGACTTTggcgcgcc'
## internal cloning site will be BspQI_cloning_site = 'GTTTGAAGAGC'
INTERNAL_CLONING_SITE = 'gtttgaagagc' ## 'GTTTGAAGAGC'  ##  'gtttGAAGAGCGCTCTTCacga' 
## a future alternative would be the BsaI site: GGTCTC, which would give the following internal cloning site: gtttaGAGACC
## BsaI cuts 1 and 5 nt away from its recognition site
INTERNAL_CLONING_SITE_INDEX = len(FWD_PRIMING_SEQUENCE) + GUIDE_LENGTH + 4
FWD_RESTRICTION_SITE_INDEX = 8
## 7 reverse priming sequences previously validated (all should be the same length)
RIGHT_SUBPOOL_AMPLIFICATION_SEQUENCES  = 'gtccttggactttgc', 'cctatccgaactacg', 'gtaagatcacctgcg', 'ggtacgaggcgactg', 'cggagcgttgaggac', 'agaccggtgtgcagg', 'aggtctgcggagagc'
REVERSE_PRIMING_SEQUENCES = [rev_comp(e) for e in RIGHT_SUBPOOL_AMPLIFICATION_SEQUENCES]

DISTANCE_FROM_PAM_PARTITIONS = 0, 5, 10, 15, 20

DONOR_LENGTH = OLIGO_LENGTH - len(FWD_PRIMING_SEQUENCE) - GUIDE_LENGTH - len(INTERNAL_CLONING_SITE) - len(REVERSE_PRIMING_SEQUENCES[0])
MAX_VARIANT_LENGTH = DONOR_LENGTH - MINIMUM_HOMOLOGY*2

RESTRICTION_SITES_SCREENED_OUT = {}
RESTRICTION_SITES_SCREENED_OUT[INTERNAL_RESTRICTION_SITE] = 0
RESTRICTION_SITES_SCREENED_OUT[rc_INTERNAL_RESTRICTION_SITE] = 0
#for site in END_RESTRICTION_SITES:
#    RESTRICTION_SITES_SCREENED_OUT[site] = 0

## modified design for donor-guide orientation
GIBSON_FWD_PRIMING_SEQ = 'GCGAATGGGACTTTgg'.lower()
GIBSON_REV_PRIMING_SEQ = 'Gtttaagagctatgctggaaa'.lower()
primer_length_difference = len(FWD_PRIMING_SEQUENCE) + len(REVERSE_PRIMING_SEQUENCES[0]) - len(GIBSON_FWD_PRIMING_SEQ) - len(GIBSON_REV_PRIMING_SEQ)
total_length_difference = len(INTERNAL_CLONING_SITE) + primer_length_difference 
LEFT_EXTENSION_FOR_GIBSON = total_length_difference // 2
RIGHT_EXTENSION_FOR_GIBSON = total_length_difference - LEFT_EXTENSION_FOR_GIBSON

header = '\t'.join('oligo_name (number at the end corresponds to the position of the disruption that is closest to the PAM),category (each category will have different subpool priming sequences at the 3-prime end of the oligos),chromosome	position_of_PAM,strand_of_PAM,guide_only,donor,ref_allele,variant_allele,chromosomal_pos_of_change,guide_with_PAM,highest_scoring_aligned_donor_segment,alignment ("_"=match,"*"=mismatch,"-"=indel),bp_from_PAM_for_PAM_disruption(s) (+1 = first nucleotide upstream of PAM; 0 = N of the NGG; -1 = first G of NGG; -2 = second G of NGG)'.split(',') ) + '\n'
oligo_header = '\t'.join('oligo_name,category,oligo_sequence'.split(',')  ) + '\n'

convert_chromosome_number_to_numeral = {}
nums = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17'.split(',')
numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI,M'.split(',')
for idx in range(len(nums)):
    convert_chromosome_number_to_numeral['chr' + nums[idx]] = 'chr' + numerals[idx]

def load_vcf(vcf_infilename):
    '''
    takes a vcf file and returns a dictionary of chromosomes of coordinates of ref seq and its variant replacement
    '''
    variant_set = set([])
    variants = generate_empty_nonstranded_bedgraph_dict()
    with open(vcf_infilename) as infile:
        for line in infile:
            if line[0] != '#':
                info = line.strip().split()
                CHROM, POS, REF, ALT, P_VALUE, D_VALUE, TYPE, STRAIN = info  ## CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SK1 = info
                CHROM = convert_chromosome_number_to_numeral[CHROM]
                try:
                    int(POS)
                except:
                    print('error', info)                
                POS = int(POS)
                ## remove constant sequence at the beginning of the REF sequence, as these are not the start positions of the variants
                ## it is assumed that this constant sequence assists in resolving repetitive sequences and ambiguous coordinates (e.g. 0- or 1- based)
                while len(REF) > 1 and len(ALT) > 1 and REF[-1] == ALT[-1]:
                    POS += 1
                    REF = REF[1:]
                    ALT = ALT[1:]
                ## remove constant sequence at the end of the REF and ALT sequence
                while len(REF) > 1 and len(ALT) > 1 and REF[-1] == ALT[-1]:
                    REF = REF[:-1]
                    ALT = ALT[:-1]
                if REF == '' or ALT == '':
                    print(info)
                
                if 'M' not in CHROM and ',' not in REF and ',' not in ALT:  ##   if 'filter' not in FILTER and 'M' not in CHROM and ',' not in REF and ',' not in ALT:
                    variant_name = '_'.join( [CHROM, str(POS), REF, ALT] )
                    variant_info = POS, REF, ALT
                    if variant_name not in variant_set:
                        variant_set.add(variant_name)
                        if POS not in variants[CHROM]:
                            variants[CHROM][POS] = [ (REF, ALT, variant_info), ]
                        else:
                            print('POS clash at CHROM, POS, REF, ALT:', CHROM, POS, REF, ALT, 'with', variants[CHROM][POS])
                            variants[CHROM][POS].append( (REF, ALT, variant_info) )
                    else:
                        print('variant duplicate for', variant_name)
    return variants
    
def load_processed_variants(processed_vcf_infilename, variant_set):
    '''
    takes a vcf file and returns a dictionary of chromosomes of coordinates of ref seq and its variant replacement
    '''
    variants = generate_empty_nonstranded_bedgraph_dict()
    with open(processed_vcf_infilename) as infile:
        for line in infile:
            if line[0] != '#':
                info = line.strip().split('\t')
                CHROM, POS, REF, ALT, individual_variants = info
                variant_name = '_'.join( [CHROM, str(POS), REF, ALT] )
                POS = int(POS)
                if variant_name not in variant_set:
                    variant_set.add(variant_name)
                    if POS not in variants[CHROM]:
                        variants[CHROM][POS] = [ (REF, ALT, individual_variants), ]
                    else:
                        print('POS clash at CHROM, POS, REF, ALT:', CHROM, POS, REF, ALT, 'with', variants[CHROM][POS])
                        variants[CHROM][POS].append( (REF, ALT, individual_variants) )
                else:
                    print('variant duplicate for', variant_name)
    return variants, variant_set
    

def combine_variants(variant_list):
    '''
    takes a variant dict list and returns a dictionary of chromosomes of coordinates of overlapping variants of ref seq and its variant replacement
    '''
    variant_set = set([])
    combined_variants = generate_empty_nonstranded_bedgraph_dict()
    conflicted_variants = generate_empty_nonstranded_bedgraph_dict()
    for variant_dict in variant_list:
        for CHROM in variant_dict:
            for POS in variant_dict[CHROM]:
                for variant in variant_dict[CHROM][POS]:
                    REF, ALT, variant_info = variant
                    variant_name = '_'.join([CHROM, str(POS), REF, ALT])
                    if variant_name not in variant_set:
                        variant_set.add(variant_name)
                        if POS in combined_variants[CHROM]:  ## first variant dict in list gets prioritized
                            print('POS clash at CHROM, POS, REF, ALT:', CHROM, POS, REF, ALT, 'with', combined_variants[CHROM][POS])
                            conflicted_variants[CHROM][POS] = variant_dict[CHROM][POS], combined_variants[CHROM][POS]
                            combined_variants[CHROM][POS].append( (REF, ALT, variant_info) )
                        else:
                            combined_variants[CHROM][POS] = [ (REF, ALT, variant_info), ]
                    else:
                        print('variant duplicate for', variant_name)
    return combined_variants, conflicted_variants
    
def link_close_variants(genome_seq, combined_variants, variant_set):
    '''
    takes the genome sequence, variant dictionaries, and a variant set, and combines the variants in the dictionaries
    in a such a way that closely spaced variants in either dictionary are included as a single variant,
    variants are linked to the variants to the right as long as the linked variants are within MAX_BP_BETWEEN_VARIANTS_TO_COMBINE 
    and the total variant length does not exceed MAX_VARIANT_LENGTH
    this is important for coding regions for ex., so that two close indels or SNPs result in the correct amino acids
    returns a linked variant dictionary with only those variants that get linked together (all single variants will be constructed as well separately)
    '''
    linked_variants = generate_empty_nonstranded_bedgraph_dict()
    for CHROM in combined_variants:
        current_cluster_start_pos = 0
        current_cluster_end_pos = 0
        current_cluster = []
        variant_positions = sorted(combined_variants[CHROM])
        num_variant_positions = len(variant_positions)
        for POS_idx in range( num_variant_positions ):
            POS = variant_positions[POS_idx]
            queue = []
            for overlapping_variant in combined_variants[CHROM][POS]:  ## allow for discordant variants to be combined with other variants to the right
                REF, ALT, variant_info = overlapping_variant
                current_cluster_start_pos = POS
                current_cluster_end_pos = POS +  len(REF)
                ## initialize each potential cluster
                current_cluster = [ variant_info, ]
                queue.append( [current_cluster,  current_cluster_start_pos, current_cluster_end_pos, POS_idx]  )
                while queue != []: 
                    current_cluster_info = queue.pop()
                    current_cluster,  current_cluster_start_pos, current_cluster_end_pos, POS_idx = current_cluster_info
                    ## add cluster if it contains 2 or more variants
                    if len(current_cluster) > 1:
                        ref_seq = genome_seq[CHROM][current_cluster_start_pos - 1: current_cluster_end_pos - 1]            
                        variant_seq = list(ref_seq)
                        offset = current_cluster_start_pos
                        for individual_variant in current_cluster:
                            individual_variant_POS, individual_variant_REF, individual_variant_ALT = individual_variant
                            individual_variant_offset = len(individual_variant_REF) - len(individual_variant_ALT)
                            variant_seq[individual_variant_POS - offset: individual_variant_POS - offset + len(individual_variant_REF) ] = individual_variant_ALT
                            offset = offset + individual_variant_offset
                        variant_seq = ''.join(variant_seq)
                        variant_name = '_'.join([CHROM, str(current_cluster_start_pos), ref_seq, variant_seq])
                        if variant_name not in variant_set:
                            variant_set.add(variant_name)
                            if current_cluster_start_pos in linked_variants[CHROM]:  ## first variant dict in list gets prioritized
                                linked_variants[CHROM][current_cluster_start_pos].append( (ref_seq, variant_seq, current_cluster) )
                            else:
                                linked_variants[CHROM][current_cluster_start_pos] = [ (ref_seq, variant_seq, current_cluster), ]
                        else:
                            print('variant duplicate for', variant_name)
                    next_POS_idx = POS_idx
                    while next_POS_idx < num_variant_positions - 1:
                        next_POS_idx += 1
                        next_POS = variant_positions[next_POS_idx]
                        if next_POS - current_cluster_end_pos < MAX_BP_BETWEEN_VARIANTS_TO_COMBINE: ## expand the cluster to the right if variants are close together
                            for next_overlapping_variant in combined_variants[CHROM][next_POS]:
                                REF, ALT, variant_info = next_overlapping_variant
                                next_current_cluster = current_cluster[:]
                                next_current_cluster.append( (next_POS, REF, ALT)  )
                                next_current_cluster_end_pos = next_POS +  len(REF)
                                queue.append( [next_current_cluster,  current_cluster_start_pos, next_current_cluster_end_pos, next_POS_idx]  )
    return linked_variants, variant_set
    
def write_linked_variants(variants, variant_outfilename):
    '''
    takes a variants dictionary and writes an outfile to store the variants
    returns None
    '''
    total_variants  = 0
    total_individual_variants  = 0
    with open(variant_outfilename, 'w') as outfile:
        for chrom in variants:
            variant_pos = set(list(variants[chrom]))
            unique_variant_pos = list(variant_pos)
            for pos in sorted(unique_variant_pos):
                for variant in variants[chrom][pos]:
                    total_variants += 1
                    variant_info = []
                    lst = variant
                    ref, alt, individual_variants = lst[0], lst[1], lst[2:]
                    for lst in individual_variants:
                        total_individual_variants += 1
                        info = ','.join( [str(e) for e in lst] )
                        variant_info.append(info)
                    output = '\t'.join([chrom, str(pos), ref, alt, ';'.join(variant_info)]) + '\n'
                    outfile.write(output)
    outfile.close()
    print('total_variants', total_variants)
    print('total_individual_variants', total_individual_variants)
    return None
    
def write_variants(variants, variant_outfilename):
    '''
    takes a variants dictionary and writes an outfile to store the variants
    returns None
    '''
    total_variants  = 0
    total_individual_variants  = 0
    with open(variant_outfilename, 'w') as outfile:
        for chrom in variants:
            variant_pos = set(list(variants[chrom]))
            unique_variant_pos = list(variant_pos)
            for pos in sorted(unique_variant_pos):
                for variant in variants[chrom][pos]:
                    total_variants += 1
                    ref, alt, individual_variant = variant[0], variant[1], variant[2:]
                    info = ','.join( [str(e) for e in individual_variant] )
                    output = '\t'.join([chrom, str(pos), ref, alt, info ]) + '\n'
                    outfile.write(output)
    outfile.close()
    print('total_variants', total_variants)
    print('total_individual_variants', total_individual_variants)
    return None

def T_stretch_present(guide_seq_with_PAM):
    '''
    takes the guide sequence and checks for a stretch of T's as specified by the global variable NUM_CONSECUTIVE_T_DISALLOWED
    returns True if T stretch equal to or greater than NUM_CONSECUTIVE_T_DISALLOWED is found
    '''
    if 'T'*NUM_CONSECUTIVE_T_DISALLOWED in guide_seq_with_PAM:
        return True
    else:
        return False
                    
def check_guide_disruption(guide_seq_with_PAM, donor_DNA, left_side_donor_length, right_side_donor_length, individual_variants, target):
    '''
    takes a selected guide RNA sequence and its PAM, and the donor DNA intended to disrupt the recognition, and length of WT sequence on left side of variant
    returns False if guide recognition is not disrupted, and True if there is a disruption
    '''
    info = individual_variants[1:-1].split(',')
    ref = info[1]
    parsed_ref = ''
    for char in ref:
        if char in 'TACG':
            parsed_ref += char
    alt = info[2]
    parsed_alt = ''
    for char in alt:
        if char in 'TACG':
            parsed_alt += char
    left_idx = max ( left_side_donor_length - GUIDE_LENGTH - len(PAM) - len(parsed_ref) - len(parsed_alt), 0 )
    right_idx = min( DONOR_LENGTH - right_side_donor_length + GUIDE_LENGTH + len(PAM) + len(parsed_ref) + len(parsed_alt), DONOR_LENGTH )
    query_length = len(guide_seq_with_PAM)
    for idx in range(left_idx, right_idx): ## len(donor_DNA)):   ## no need to check the entire donor for whether or not disruption occurred
        query = donor_DNA[idx: idx + query_length]
        query_rc = rev_comp(query)
        if hamming_distance(query, guide_seq_with_PAM) == 0:
            return False
        elif hamming_distance(query_rc, guide_seq_with_PAM) == 0:
            return False
    return True

def build_backtrack_for_fitting_alignment(v, w):
    '''
      Input: Two nucleotide strings v and w, where v has length at most 1000 and w has length at most 100.
     Output: A highest-scoring fitting alignment between v and w. Use the simple scoring method in which
     matches count +1 and both the mismatch and indel penalties are 1.
    '''
    indel_penalty = 3
    max_index_on_rightmost_column = (0,0)
    max_score_on_rightmost_column = 0
    num_rows = len(v)+1
    num_columns = len(w)+1
    ## need to use all bases of w in the alignment, and any number of bases in v to maximize the score
    path_array = numpy.zeros((num_rows,num_columns))
    backtracking_array = numpy.zeros((num_rows,num_columns))
    for i in range(1,num_rows):
        ## Build the path array where first column is all zeros, allowing any 5' truncation of w to not count towards the alignment score
        backtracking_array[i][0] = 1
    for j in range(1,num_columns):
        path_array[0][j] = path_array[0][j-1] - indel_penalty
        backtracking_array[0][j] = 2
        ## Proceed with building path array as per the global alignment
    for i in range(1,num_rows):
        for j in range(1,num_columns):
            if [v[i-1]] == [w[j-1]]:
                diagonal_value = 1
            else:
                diagonal_value = -1
            path_array[i][j] = max(path_array[i-1][j] - indel_penalty, path_array[i][j-1] - indel_penalty, path_array[i-1][j-1] + diagonal_value)
            if path_array[i][j] == path_array[i-1][j] - indel_penalty:
                backtracking_array[i][j] = 1
            elif path_array[i][j] == path_array[i][j-1] - indel_penalty:
                backtracking_array[i][j] = 2
            else:
                backtracking_array[i][j] = 3
    ## find the maximum score on the far right column, which represents all of w and a possible 3' truncation of v
    j = num_columns - 1
    for i in range(num_rows):
        if path_array[i][j] > max_score_on_rightmost_column:
            max_score_on_rightmost_column = path_array[i][j]
            max_index_on_rightmost_column = (i, j)

    i = max_index_on_rightmost_column[0]
    j = max_index_on_rightmost_column[1]
    alignment_score = int(path_array[i][j])
    return i, j, backtracking_array, alignment_score

def fitting_alignment(i, j, backtracking_array, v, w):
    '''
     Fitting Alignment Problem: Construct a highest-scoring fitting alignment between two strings.
     Input: Strings v (donor) and w (guide_with_PAM) as well as a matrix backtracking_array.
     Output: A highest-scoring fitting alignment of v and w as defined by the scoring matrix backtracking_array.
     All of w is needed, with possible end truncations of v
     '''
    query = ''
    subject = ''
    alignment = ''
    initial_j =  j
    snp_positions = []
    indel_positions = []
    while j > 0:
        if backtracking_array[i][j] == 3:
            query = v[i-1] + query
            subject = w[j-1] + subject
            if v[i-1] == w[j-1]:
                alignment = '_' + alignment
            else:
                alignment = '*' + alignment
                snp_positions.append( initial_j - j )
            i -= 1
            j -= 1

        elif backtracking_array[i][j] == 1:
            query = v[i-1] + query
            subject = '-' + subject
            alignment = '-' + alignment
            indel_positions.append( initial_j - j )
            i -= 1
            
        else:
            query = '-' + query
            subject = w[j-1] + subject
            alignment = '-' + alignment
            indel_positions.append( initial_j - j )
            j -= 1
    if snp_positions == [] and indel_positions == []:
        print('query, subject, alignment', query, subject, alignment)
    return query, subject, alignment, snp_positions, indel_positions

def process_disruption_positions_for_snps_indels(snp_positions, indel_positions):
    '''
    takes a list of snp positions, and indel positions
    returns the position closest to the PAM that disrupts guide RNA recognition 
    (i.e., if there is a snp in the N of the NGG, there is no disruption, whereas an indel does cause a disruption)
    variants disrupting the PAM give
    '''
    disruptions = []
    for position in snp_positions:
        disruption_position = position - len(PAM) + 1  ## convert to 1-based, so that the first nt upstream the PAM is denoted as position 1
        if disruption_position != 0:  
            disruptions.append(disruption_position)
    for position in indel_positions:
        disruption_position = position - len(PAM) + 1 
        disruptions.append(disruption_position)
    return disruptions
    
def get_disruption_position_in_guide_donor_alignment(guide_with_full_PAM, donor):
    '''
    takes the guide RNA sequence plus full PAM, the PAM, and the donor DNA
    aligs the guide and full PAM to the donor DNA sequence
    returns the best fitting alignment between the guide_with_full_PAM and donor, 
    and the closest position of disruption with 0 denoting PAM disruption, and positive integers denoting distance of SNP/indel from the PAM
    only an indel can have a position of 0, as a snp at position 0 does not disrupt the NGG
    '''
    sense_i, sense_j, sense_backtracking_array, sense_alignment_score = build_backtrack_for_fitting_alignment(donor, guide_with_full_PAM)
    rc_donor = rev_comp(donor)
    antisense_i, antisense_j, antisense_backtracking_array, antisense_alignment_score = build_backtrack_for_fitting_alignment(rc_donor, guide_with_full_PAM)
    if sense_alignment_score > antisense_alignment_score:
        query, subject, alignment, snp_positions, indel_positions = fitting_alignment(sense_i, sense_j, sense_backtracking_array, donor, guide_with_full_PAM)
    else:
        query, subject, alignment, snp_positions, indel_positions = fitting_alignment(antisense_i, antisense_j, antisense_backtracking_array, rc_donor, guide_with_full_PAM)
    disruptions = process_disruption_positions_for_snps_indels(snp_positions, indel_positions)
    return disruptions, query, subject, alignment
  
def generate_guide_donor_pairs(genome_seq, variants, PAM, outfile, oligo_outfile, variant_category, variant_statistics_outfile):
    '''
    takes a linked variant dictionary and searches for guide RNAs within a distance MAX_DIST_TO_PAM from the indicated PAM
    SNPs in the N of the NGG cannot be used!
    returns a dictionary of chromosomes of coordinates for the N in NGG of the guide RNAs of guide, donor tuples
    '''
    guides_per_category = {}
    guides_per_variant = {}
    guides_with_restriction_sites = 0
    donors_with_restriction_sites_not_allowing_shifting = 0  
    oligos_with_restriction_sites_made_with_junction_sequences = []
    ## a restriction site near the edge of the donor DNA should not prevent the variant from being included!
    ## shift the homologous arms to be asymmetric lengths, as long as the shorter arm retains a minimum homology length
    ## keep track of variants without guide
    variants_with_no_PAM = generate_empty_nonstranded_bedgraph_dict()
    PAM_rc = rev_comp(PAM)
    num_variants_with_no_PAM = 0
    num_variants_with_multiple_PAMS = 0
    num_variants_with_single_PAM = 0
    total_variants = 0
    T_stretch_guides = 0
    T_stretch_donors = 0
    for CHROM in variants: ## ('chrIV',): ## 
        print('CHROM', CHROM)
        chrom_length = len(genome_seq[CHROM])       
        variant_pos = set(list(variants[CHROM]))
        unique_variant_pos = list(variant_pos)
        for POS in sorted( unique_variant_pos ):
            for overlapping_variant in variants[CHROM][POS]:
                ref_seq, variant_seq, individual_variants = overlapping_variant
                total_variants += 1
                num_guides = 0
                pos_guide_donor = []
                start_of_variant = POS - 1
                end_of_variant = POS + len(ref_seq) - 1
                variant_length = len(variant_seq)
                left_side_donor_length = (DONOR_LENGTH - variant_length) // 2
                right_side_donor_length = DONOR_LENGTH - variant_length - left_side_donor_length
                ## check if variant is near end of chromosome such that donor DNA will be clipped                
                end_of_chromosome_shift = 0
                at_end_of_chromosome = False
                at_beginning_of_chromosome = False                
                if end_of_variant + right_side_donor_length > chrom_length:
                    end_of_chromosome_shift = end_of_variant + right_side_donor_length - chrom_length
                    right_side_donor_length -= end_of_chromosome_shift
                    left_side_donor_length += end_of_chromosome_shift
                    at_end_of_chromosome = True
                if start_of_variant - left_side_donor_length < 1:
                    end_of_chromosome_shift = left_side_donor_length - start_of_variant
                    left_side_donor_length -= end_of_chromosome_shift
                    right_side_donor_length += end_of_chromosome_shift
                    at_beginning_of_chromosome = True
                donor = genome_seq[CHROM][start_of_variant - left_side_donor_length: start_of_variant] + variant_seq + genome_seq[CHROM][end_of_variant: end_of_variant + right_side_donor_length]
                
                ## check if donor DNA contains restriction sites, and shift homologous arms to asymmetric lengths if necessary            
                donor_shift = 0
                donor_contains_restriction_site_and_shift_will_not_allow_keeping_minimum_homology = False
                left_side_restriction_sites_present = []
                right_side_restriction_sites_present = []
                for restriction_site in RESTRICTION_SITES_SCREENED_OUT:
                    if restriction_site in donor:
                        print('donor_harbors_restriction_site: ', restriction_site)
                        restriction_site_idx = donor.index(restriction_site)
                        left_boundary = left_side_donor_length - MINIMUM_HOMOLOGY
                        right_boundary = left_side_donor_length + len(variant_seq) + MINIMUM_HOMOLOGY
                        if left_side_donor_length - MINIMUM_HOMOLOGY > 0:
                                ## could add donor-guide orientation designs here ##
                            if restriction_site_idx >= right_boundary and not at_beginning_of_chromosome:
                                donor_shift = DONOR_LENGTH - restriction_site_idx  ## shift donor to the left requires positive value
                                right_side_restriction_sites_present.append( (restriction_site, donor_shift) )
                            elif restriction_site_idx <= right_boundary and not at_end_of_chromosome:
                                donor_shift = restriction_site_idx - left_boundary  ## shift donor to the right requires negative value
                                left_side_restriction_sites_present.append( (restriction_site, donor_shift) )
                            else:
                                print('donor_shift_will_not_allow_keeping_minimum_homology: ', restriction_site, restriction_site_idx, donor)
                                ## can include a subset of donors in front of the gRNA if the internal cloning site is in the middle of the donor
                                left_side_donor_length -= LEFT_EXTENSION_FOR_GIBSON
                                right_side_donor_length += RIGHT_EXTENSION_FOR_GIBSON
                if donor_shift > 0 and not donor_contains_restriction_site_and_shift_will_not_allow_keeping_minimum_homology:
                    if right_side_restriction_sites_present == [] or left_side_restriction_sites_present == []:
                        max_donor_shift = 0
                        if right_side_restriction_sites_present != []:
                            for restriction_site_info in right_side_restriction_sites_present:
                                restriction_site, donor_shift = restriction_site_info
                                if donor_shift > max_donor_shift:
                                    max_donor_shift = donor_shift
                        if left_side_restriction_sites_present != []:
                            for restriction_site_info in left_side_restriction_sites_present:
                                restriction_site, donor_shift = restriction_site_info
                                if donor_shift < max_donor_shift:
                                    max_donor_shift = donor_shift
                        left_side_donor_length = left_side_donor_length + donor_shift
                        right_side_donor_length = right_side_donor_length - donor_shift
                    else:
                        print('two restriction sites present on either side of donor', left_side_restriction_sites_present, right_side_restriction_sites_present)
                        donors_with_restriction_sites_not_allowing_shifting += 1
                        donor_contains_restriction_site_and_shift_will_not_allow_keeping_minimum_homology = True
                        left_side_donor_length -= LEFT_EXTENSION_FOR_GIBSON
                        right_side_donor_length += RIGHT_EXTENSION_FOR_GIBSON
                donor = genome_seq[CHROM][start_of_variant - left_side_donor_length: start_of_variant] + variant_seq + genome_seq[CHROM][end_of_variant: end_of_variant + right_side_donor_length]
                target =  genome_seq[CHROM][start_of_variant - left_side_donor_length: start_of_variant] + ' ' +  ref_seq + ' ' +  genome_seq[CHROM][end_of_variant: end_of_variant + right_side_donor_length ]
                ## determine whether variant that disrupts a PAM    
                ## search for a PAM upstream of the variant end position, and downstream of variant start position
                PAMs_found = 0
                upstream_PAMs_coords = set([])
                for idx in range( end_of_variant + len(PAM) + 1, start_of_variant - MAX_DIST_TO_PAM - 1, - 1):
                    query = genome_seq[CHROM][idx - len(PAM): idx ]
                    if hamming_distance(query, PAM_rc) == 0:
                        PAMs_found += 1
                        upstream_PAMs_coords.add( idx )
    
                downstream_PAMs_coords = set([])
                for idx in range( start_of_variant - len(PAM) - 1, end_of_variant + MAX_DIST_TO_PAM + 1):
                    query = genome_seq[CHROM][idx: idx + len(PAM) ]
                    if hamming_distance(query, PAM) == 0:
                        PAMs_found += 1
                        downstream_PAMs_coords.add( idx )
                ## select the guide RNAs with the donor DNA variants closest to the PAMs until MAX_GUIDES_PER_DONOR is filled
                if PAMs_found > 0 :
                    upstream_PAM_coords_selected = set([])
                    downstream_PAM_coords_selected = set([])
                    for coord in range(start_of_variant, end_of_variant + 1):
                        if num_guides <= MAX_GUIDES_PER_DONOR:
                            if coord in upstream_PAMs_coords:
                                guide_seq = rev_comp( genome_seq[CHROM][coord : coord + GUIDE_LENGTH ] )
                                guide_seq_with_PAM = guide_seq + PAM
                                guide_seq_with_full_PAM = rev_comp( genome_seq[CHROM][coord - 3: coord + GUIDE_LENGTH ] )
                                ## analyze selected guide sequences and PAM for mismatches to the donor
                                if coord not in upstream_PAM_coords_selected and check_guide_disruption(guide_seq_with_PAM, donor, left_side_donor_length, right_side_donor_length, individual_variants, target):
                                    upstream_PAM_coords_selected.add(coord)                                    
                                    num_guides += 1
                                    pos_guide_donor.append( (CHROM, coord, '-', guide_seq_with_full_PAM, guide_seq, donor)  )
                        if num_guides <= MAX_GUIDES_PER_DONOR:
                            if coord in downstream_PAMs_coords:
                                guide_seq = genome_seq[CHROM][coord - GUIDE_LENGTH: coord]
                                guide_seq_with_PAM = guide_seq + PAM
                                guide_seq_with_full_PAM = genome_seq[CHROM][coord - GUIDE_LENGTH: coord + 3]
                                if coord not in downstream_PAM_coords_selected and check_guide_disruption(guide_seq_with_PAM, donor, left_side_donor_length, right_side_donor_length, individual_variants, target):
                                    downstream_PAM_coords_selected.add(coord)                                    
                                    num_guides += 1
                                    pos_guide_donor.append( (CHROM, coord, '+', guide_seq_with_full_PAM, guide_seq, donor)  )
                    for idx in range(MAX_DIST_TO_PAM):
                        if num_guides <= MAX_GUIDES_PER_DONOR:
                            coord = start_of_variant - idx - 1
                            if coord in upstream_PAMs_coords:    
                                guide_seq = rev_comp( genome_seq[CHROM][coord : coord + GUIDE_LENGTH ] )
                                guide_seq_with_PAM = guide_seq + PAM
                                guide_seq_with_full_PAM = rev_comp( genome_seq[CHROM][coord - 3: coord + GUIDE_LENGTH ] )
                                if coord not in upstream_PAM_coords_selected and check_guide_disruption(guide_seq_with_PAM, donor, left_side_donor_length, right_side_donor_length, individual_variants, target):
                                    upstream_PAM_coords_selected.add(coord)                                                                        
                                    num_guides += 1
                                    pos_guide_donor.append( (CHROM, coord, '-', guide_seq_with_full_PAM, guide_seq, donor)  )
                        if num_guides <= MAX_GUIDES_PER_DONOR:
                            coord = end_of_variant + idx + 1
                            if coord in downstream_PAMs_coords:
                                guide_seq = genome_seq[CHROM][coord - GUIDE_LENGTH: coord]
                                guide_seq_with_PAM = guide_seq + PAM
                                guide_seq_with_full_PAM = genome_seq[CHROM][coord - GUIDE_LENGTH: coord + 3]
                                if coord not in downstream_PAM_coords_selected and check_guide_disruption(guide_seq_with_PAM, donor, left_side_donor_length, right_side_donor_length, individual_variants, target):
                                    downstream_PAM_coords_selected.add(coord)                                    
                                    num_guides += 1
                                    pos_guide_donor.append( (CHROM, coord, '+', guide_seq_with_full_PAM, guide_seq, donor)  )
                if num_guides == 0:
                    variants_with_no_PAM[CHROM][POS] = ref_seq, variant_seq, individual_variants
                    num_variants_with_no_PAM += 1
                elif num_guides == 1:
                    num_variants_with_single_PAM += 1
                elif num_guides > 1:
                    num_variants_with_multiple_PAMS += 1
                if num_guides not in guides_per_variant:
                    guides_per_variant[num_guides] = 1
                else:
                    guides_per_variant[num_guides] += 1
                for info in pos_guide_donor:
                    CHROM, coord, strand, guide_seq_with_full_PAM, guide_seq, donor = info
                    if strand == '-':
                        donor = rev_comp(donor) ## optional: change the donor to be sense with the guide to avoid donor and guide RNA interactions for donors transcribed with guide RNA
                    disruptions, query, subject, alignment = get_disruption_position_in_guide_donor_alignment(guide_seq_with_full_PAM, donor)
                    if disruptions == []:
                        disruption_position = GUIDE_LENGTH + 1  ## donor has another perfect match to the guide-PAM outside of the middle region of the donor
                                                        ## these oligos are included as controls as they should result in cleavage of donor and target
                    else:
                        disruption_position = min(disruptions)
                    oligo_name = '>' + '_'.join([CHROM, str(POS),ref_seq, 'to', variant_seq, 'PAM', str(coord), strand, variant_category, str( disruption_position )])
                    if disruption_position <= DISTANCE_FROM_PAM_PARTITIONS[0]:
                        reverse_priming_sequence = REVERSE_PRIMING_SEQUENCES[0]
                        category = '1:PAM_disruption_or_1_to_' + str(DISTANCE_FROM_PAM_PARTITIONS[1]) + '_bp_from_PAM_disruption'
                    elif disruption_position <= DISTANCE_FROM_PAM_PARTITIONS[1]:
                        reverse_priming_sequence = REVERSE_PRIMING_SEQUENCES[0]
                        category = '1:PAM_disruption_or_1_to_' + str(DISTANCE_FROM_PAM_PARTITIONS[1]) + '_bp_from_PAM_disruption'
                    elif disruption_position <= DISTANCE_FROM_PAM_PARTITIONS[2]:
                        reverse_priming_sequence = REVERSE_PRIMING_SEQUENCES[1]
                        category = '2:' + str(DISTANCE_FROM_PAM_PARTITIONS[1]) + '_to_' + str(DISTANCE_FROM_PAM_PARTITIONS[2]) + '_bp_from_PAM_disruption'
                    elif disruption_position <= DISTANCE_FROM_PAM_PARTITIONS[3]:
                        reverse_priming_sequence = REVERSE_PRIMING_SEQUENCES[2]
                        category = '3:' + str(DISTANCE_FROM_PAM_PARTITIONS[2]) +  '_to_' + str(DISTANCE_FROM_PAM_PARTITIONS[3]) + '_bp_from_PAM_disruption'
                    elif disruption_position <= DISTANCE_FROM_PAM_PARTITIONS[4]:
                        reverse_priming_sequence = REVERSE_PRIMING_SEQUENCES[3]
                        category = '4:' + str(DISTANCE_FROM_PAM_PARTITIONS[3]) + '_to_' + str(DISTANCE_FROM_PAM_PARTITIONS[4]) + '_bp_from_PAM_disruption'
                    else:  ## included in the distal PAM variants as controls
                        reverse_priming_sequence = REVERSE_PRIMING_SEQUENCES[3]
                        category = '4:' + str(DISTANCE_FROM_PAM_PARTITIONS[3]) + '_to_' + str(DISTANCE_FROM_PAM_PARTITIONS[4]) + '_bp_from_PAM_disruption'
                    restriction_site_in_guide = False
                    restriction_site_in_donor = False
                    restriction_site_in_junction = False
                    oligo_sequence = FWD_PRIMING_SEQUENCE + guide_seq + INTERNAL_CLONING_SITE + donor + reverse_priming_sequence
                    oligo_seq_upper = oligo_sequence.upper()
                    for restriction_site in RESTRICTION_SITES_SCREENED_OUT:
                        if restriction_site in guide_seq:
                            RESTRICTION_SITES_SCREENED_OUT[restriction_site] += 1
                            restriction_site_in_guide = True
                        if restriction_site in donor:
                            RESTRICTION_SITES_SCREENED_OUT[restriction_site] += 1
                            restriction_site_in_donor = True
                        elif restriction_site in oligo_seq_upper:
                            site_index = oligo_seq_upper.index(restriction_site) 
                            if site_index > 9  and site_index != INTERNAL_CLONING_SITE_INDEX:
                                restriction_site_in_junction = True
                                # print('restriction site made with junction between sequences!')
                                # print(oligo_sequence)
                                oligos_with_restriction_sites_made_with_junction_sequences.append(oligo_name)
                    if not restriction_site_in_donor and not restriction_site_in_guide and not restriction_site_in_junction:
                        if T_stretch_present(guide_seq):
                            reverse_priming_sequence = REVERSE_PRIMING_SEQUENCES[4]
                            category = '5:guide_with_T_stretch'
                            oligo_sequence = FWD_PRIMING_SEQUENCE + guide_seq + INTERNAL_CLONING_SITE + donor + reverse_priming_sequence
                            output = '\t'.join( [oligo_name, category] + [str(e) for e in info] + [ref_seq, variant_seq, individual_variants, query, subject, alignment, ','.join([str(e) for e in disruptions]) ] ) +  '\n'
                            oligo_outfile.write(oligo_name + '\t' + category + '\t' + oligo_sequence + '\n')
                            outfile.write(output)
                            T_stretch_guides += 1
                        else:
                            output = '\t'.join( [oligo_name, category] + [str(e) for e in info] + [ref_seq, variant_seq, individual_variants, query, subject, alignment, ','.join([str(e) for e in disruptions]) ] ) +  '\n'
                            oligo_outfile.write(oligo_name + '\t' + category + '\t' + oligo_sequence + '\n')
                            outfile.write(output)
                    if donor_contains_restriction_site_and_shift_will_not_allow_keeping_minimum_homology:
                        category = '6:donor_guide'
                        if not T_stretch_present(donor + guide_seq):
                            oligo_sequence = GIBSON_FWD_PRIMING_SEQ + donor + guide_seq + GIBSON_REV_PRIMING_SEQ
                            info = CHROM, coord, strand, guide_seq_with_full_PAM, guide_seq, donor                            
                            output = '\t'.join( [oligo_name, category] + [str(e) for e in info] + [ref_seq, variant_seq, individual_variants, query, subject, alignment, ','.join([str(e) for e in disruptions]) ] ) +  '\n'
                            oligo_outfile.write(oligo_name + '\t' + category + '\t' + oligo_sequence + '\n')
                            outfile.write(output)
                        else:
                            print('donor_guide has T-stretch', donor)
                            T_stretch_donors += 1
                    if restriction_site_in_guide:
                        guides_with_restriction_sites += 1
                    if restriction_site_in_donor:
                        donors_with_restriction_sites_not_allowing_shifting += 1
                    if category not in guides_per_category:
                        guides_per_category[category] = 0
                    else:
                        guides_per_category[category] += 1
    variant_statistics_outfile.write('\t'.join( [str(e) for e in [variant_category,total_variants,num_variants_with_no_PAM,num_variants_with_single_PAM,num_variants_with_multiple_PAMS,T_stretch_guides,T_stretch_donors,guides_per_variant,RESTRICTION_SITES_SCREENED_OUT,guides_with_restriction_sites,donors_with_restriction_sites_not_allowing_shifting,guides_per_category]]) + '\n')
    return variants_with_no_PAM
    
## Load reference sequence
try:
    genome_seq
except:
    genome_seq = load_genome(genome_filename)

process_oligos = True 
process_variants = True

if process_variants:
    SNP_variants = load_vcf(SNP_vcf_infilename)
    write_variants(SNP_variants, SNP_variant_outfilename)
    
    indel_variants = load_vcf(indel_vcf_infilename)
    write_variants(indel_variants, indel_variant_outfilename)

    
    combined_variants, conflicted_variants = combine_variants( (indel_variants, SNP_variants) )
    variant_set = set([])
    linked_variants, variant_set = link_close_variants(genome_seq, combined_variants, variant_set)
    write_linked_variants(linked_variants, linked_variant_outfilename)


if process_oligos:
    variant_set = set([])
    print('LOADING LINKED VARIANTS...')
    linked_processed_variants, variant_set = load_processed_variants(linked_variant_outfilename, variant_set)
    print('LOADING SNP VARIANTS...')
    SNP_processed_variants, variant_set =  load_processed_variants(SNP_variant_outfilename, variant_set)
    print('LOADING INDEL VARIANTS...')
    indel_processed_variants, variant_set = load_processed_variants(indel_variant_outfilename, variant_set)
    
    
    variant_guide_donor_outfilename = DIR + strain + '_SNPs_indels_with_' + str(MAX_BP_BETWEEN_VARIANTS_TO_COMBINE) + '_bp_linked_' + str(MAX_DIST_TO_PAM) + '_bp_from_' + PAM + '_PAM_guide_donor_' + str(OLIGO_LENGTH) + '_bp_oligos_info.txt'
    variant_guide_donor_oligos_outfilename = DIR + strain + '_SNPs_indels_with_' + str(MAX_BP_BETWEEN_VARIANTS_TO_COMBINE) + '_bp_linked_' + str(MAX_DIST_TO_PAM) + '_bp_from_' + PAM + '_PAM_guide_donor_' + str(OLIGO_LENGTH) + '_bp_oligos.txt'
    variant_statistics_outfilename = DIR + strain + '_SNPs_indels_with_' + str(MAX_BP_BETWEEN_VARIANTS_TO_COMBINE) + '_bp_linked_' + str(MAX_DIST_TO_PAM) + '_bp_from_' + PAM + '_PAM_guide_donor_' + str(OLIGO_LENGTH) + '_bp_oligos_statistics.txt'
    
    outfile = open(variant_guide_donor_outfilename, 'w')
    outfile.write(header)
    oligo_outfile = open(variant_guide_donor_oligos_outfilename, 'w')
    oligo_outfile.write(oligo_header)
    statistics_header_info = ['variant_category','total variants','variants with no PAM found','variants with single PAM','variants with multiple PAMS','guides with a ' + str(NUM_CONSECUTIVE_T_DISALLOWED) + ' or more T stretch', 'donors with restriction sites with a T-stretch','guides per variant'] + [e for e in RESTRICTION_SITES_SCREENED_OUT.keys()] + ['guides_with_restriction_sites' ,'donors_with_restriction_sites_not_allowing_shifting','guides_per_category']
    variant_statistics_header =   '\t'.join(statistics_header_info ) + '\n'
    variant_statistics_outfile = open(variant_statistics_outfilename, 'w')
    variant_statistics_outfile.write(variant_statistics_header )
    
    linked_variants_with_no_PAM = generate_guide_donor_pairs(genome_seq, linked_processed_variants, PAM, outfile, oligo_outfile, 'linked', variant_statistics_outfile)
    indel_variants_with_no_PAM = generate_guide_donor_pairs(genome_seq, indel_processed_variants, PAM, outfile, oligo_outfile, 'indel', variant_statistics_outfile)
    SNP_variants_with_no_PAM = generate_guide_donor_pairs(genome_seq, SNP_processed_variants, PAM, outfile, oligo_outfile, 'SNP', variant_statistics_outfile)
    
    
    
    outfile.close()
    oligo_outfile.close()
    variant_statistics_outfile.close()