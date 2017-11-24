# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 10:36:29 2016

@author: kevinroy
"""
import bedgraph_computation, numpy, GTF_GFF_manipulation

def restriction_sites_in_donors(restriction_sites, donors):
    for site in restriction_sites:
        for donor in donors:
            if site in donor:
                return True
    return False


def check_guide_disruption(log_outfile, individual_variants, variant_type, disruption_position, m0, PAM_chrom, PAM_strand, PAM_coord, guide_seq, donor_shift, donor_start, donor_end, donor_strand, left_side_donor, right_side_donor, ref_seq, variant_seq,  donor, genomic_target):
    guide_disrupted = True
    rc_donor = bedgraph_computation.rev_comp(donor)
    for nt in 'ACTG':
        for seq in donor, rc_donor:    
            if (guide_seq + nt + 'GG') in seq:
                guide_disrupted = False
    if not guide_disrupted: 
        info = 'guide not disrupted',variant_type, individual_variants, '\n', str(disruption_position), '\n', guide_seq, '\n', bedgraph_computation.rev_comp(guide_seq), '\n', left_side_donor + ' ' + ref_seq + ' ' + right_side_donor, '\n', left_side_donor + ' ' + variant_seq + ' ' + right_side_donor, '\n', donor, '\n',genomic_target,'\n', PAM_chrom, PAM_strand, str(PAM_coord), 'donor_shift_' + str(donor_shift), str(donor_start), str(donor_end),  donor_strand
        log_outfile.write('\t'.join(info) + '\n')
        #raise ValueError        
        return False
        
    guide_in_target = False
    rc_target = bedgraph_computation.rev_comp(genomic_target)
    for nt in 'ACTG':
        for seq in genomic_target, rc_target:
            if (guide_seq + nt + 'GG') in seq:
                guide_in_target = True

    if not guide_in_target:
        info = 'guide not in target',variant_type, individual_variants, '\n', str(disruption_position), '\n', guide_seq, '\n', bedgraph_computation.rev_comp(guide_seq), '\n', left_side_donor + ' ' + ref_seq + ' ' + right_side_donor, '\n', left_side_donor + ' ' + variant_seq + ' ' + right_side_donor, '\n', donor, '\n',genomic_target,'\n', PAM_chrom, PAM_strand, str(PAM_coord), 'donor_shift_' + str(donor_shift), str(donor_start), str(donor_end),  donor_strand
        log_outfile.write('\t'.join(info) + '\n')
        return False
    
    return True

def assemble_mutated_donor_sequence(ORF_strand, left_donor, right_donor, ORF_seq, aa_range, mutated_region_lowercase = False):
    '''
    takes the homologous arms and mutated middle sequence,
    returns a complete donor with mutated sequence in lower case
    left donor and right donor are GENOMIC sequence, and thus may contain UTR and intron regions
    temp_ORF_seq is a tuple of codons, from initiating methionine to stop codon, and must be joined prior to assembly
    '''    
    N_terminal_aa_num, C_terminal_aa_num = aa_range
    mutated_region = ''.join( ORF_seq[N_terminal_aa_num - 1: C_terminal_aa_num] )
    if ORF_strand == '-':
        mutated_region = bedgraph_computation.rev_comp(mutated_region)
    if len(mutated_region) % 3 != 0:
        print(aa_range, ORF_seq[N_terminal_aa_num - 1: C_terminal_aa_num],  'is not a multiple of three')
    if mutated_region_lowercase:
        mutated_region = mutated_region.lower()
    query_donor = left_donor + mutated_region + right_donor
    return query_donor

def load_Azimuth_guide_scores_and_BLAST_mismatch_counts(infilename):
    print('loading genome-wide gRNA efficacy and off-target scores...')
    PAM_coord_to_seq_score_mismatch_info = {}
    azimuth_target_seq_to_BLAST_mismatch_counts = {}
    with open(infilename, 'r') as infile:
        infile.readline()
        for line in infile:
            pcid, oseqs, aseqs, azsc, m0, m1, m2, m3, m4, m5 = line.strip().split()
            chromosome, strand, position = pcid.split('_')
            if 'NA' not in line:
                position, m0, m1, m2, m3, m4, m5 = [int(e) for e in (position, m0, m1, m2, m3, m4, m5) ]
            else:    
                position, m0, m1, m2, m3, m4, m5 = [int(e) for e in (position, 20, 20, 20, 20, 20, 20) ]
            chromosome = GTF_GFF_manipulation.convert_chromosome_number_to_numeral(chromosome)
            position -= 1
            
            PAM_coord_to_seq_score_mismatch_info[pcid] = oseqs, aseqs, azsc, m0, m1, m2, m3, m4, m5
            azimuth_target_seq_to_BLAST_mismatch_counts[aseqs] = m0, m1, m2, m3, m4, m5
        infile.close()
    print('finished loading genome-wide gRNA efficacy and off-target scores')
    return PAM_coord_to_seq_score_mismatch_info, azimuth_target_seq_to_BLAST_mismatch_counts

def get_sense_and_antisense_PAMs_for_ORF(genome_seq, ORF_chrom, ORF_exon_coords, GUIDE_LENGTH, PAM_length, NUM_CONSECUTIVE_T_DISALLOWED, MOST_T_IN_LAST_X_BP, azimuth_target_seq_to_BLAST_mismatch_counts, max_m0 =0, max_m1 = 0, max_m2 = 0, BP_FLANKING_ORF = 20, RESTRICTION_SITES_DISALLOWED = ['GGCGCGCC', 'GCGGCCGC', 'GCTCTTC'], FWD_RESTRICTION_SITE_FRAGMENT = 'CGCC' ):
    '''
    takes a chromosome, strand, and coordinate
    
    for minus strand, coord is last coord of codon, for plus strand it is the first coord of codon    
    finds 5 PAMs for each strand
    returns the PAM_coordinate, guide, and guide with full PAM sequence (outlaws guides with excessive T stretch as defined by NUM_CONSECUTIVE_T_DISALLOWED )
    '''
    start_coord = ORF_exon_coords[0]
    end_coord = ORF_exon_coords[-1]
    if end_coord < start_coord:
        temp_end_coord = end_coord
        end_coord = start_coord
        start_coord = temp_end_coord      
    plus_strand_PAMs = {}
    minus_strand_PAMs = {}
    ## PAM_strand, PAM_seq, PAM_coord, guide_seq, guide_seq_with_PAM
    for coord in range(start_coord - GUIDE_LENGTH  - PAM_length - BP_FLANKING_ORF, end_coord + GUIDE_LENGTH + PAM_length + BP_FLANKING_ORF):
        query_seq = genome_seq[ORF_chrom][coord: coord + PAM_length]
        plus_strand_PAM_guide_seq = genome_seq[ORF_chrom][coord - GUIDE_LENGTH: coord ]
        plus_strand_PAM_guide_seq_with_PAM = genome_seq[ORF_chrom][coord - GUIDE_LENGTH: coord + PAM_length ]
        restriction_site_present = False
        for site in RESTRICTION_SITES_DISALLOWED:
            if site in FWD_RESTRICTION_SITE_FRAGMENT + plus_strand_PAM_guide_seq + 'GTTT' or bedgraph_computation.rev_comp(site) in FWD_RESTRICTION_SITE_FRAGMENT + plus_strand_PAM_guide_seq + 'GTTT':
                restriction_site_present = True
        if not restriction_site_present and bedgraph_computation.hamming_distance( query_seq, PAM ) == 0 and NUM_CONSECUTIVE_T_DISALLOWED*'T' not in plus_strand_PAM_guide_seq and plus_strand_PAM_guide_seq[-MOST_T_IN_LAST_X_BP[1]:].count('T') <= MOST_T_IN_LAST_X_BP[0]:
            plus_strand_azimuth_seq = genome_seq[ORF_chrom][coord - GUIDE_LENGTH - 4: coord + PAM_length + 3]
            m0, m1, m2, m3, m4, m5 = azimuth_target_seq_to_BLAST_mismatch_counts[plus_strand_azimuth_seq]   
            if m0 <= max_m0 and m1 <= max_m1 and m2 <= max_m2:
                plus_strand_PAMs[ coord  ] = '+', query_seq, coord, plus_strand_PAM_guide_seq, plus_strand_PAM_guide_seq_with_PAM
        
        rc_query_seq = bedgraph_computation.rev_comp(query_seq)
        minus_strand_PAM_guide_seq = bedgraph_computation.rev_comp(genome_seq[ORF_chrom][coord + PAM_length: coord + PAM_length + GUIDE_LENGTH])
        minus_strand_PAM_guide_seq_with_PAM = bedgraph_computation.rev_comp(genome_seq[ORF_chrom][coord : coord + PAM_length+ GUIDE_LENGTH ])
        restriction_site_present = False
        for site in RESTRICTION_SITES_DISALLOWED:
            if site in FWD_RESTRICTION_SITE_FRAGMENT + minus_strand_PAM_guide_seq + 'GTTT' or bedgraph_computation.rev_comp(site) in FWD_RESTRICTION_SITE_FRAGMENT + minus_strand_PAM_guide_seq + 'GTTT':
                restriction_site_present = True
        if not restriction_site_present and bedgraph_computation.hamming_distance( rc_query_seq, PAM ) == 0 and NUM_CONSECUTIVE_T_DISALLOWED*'T' not in  minus_strand_PAM_guide_seq and minus_strand_PAM_guide_seq[-MOST_T_IN_LAST_X_BP[1]:].count('T') <= MOST_T_IN_LAST_X_BP[0]:
            minus_strand_azimuth_seq = bedgraph_computation.rev_comp(genome_seq[ORF_chrom][coord - 3 : coord + PAM_length+ GUIDE_LENGTH + 4 ])
            m0, m1, m2, m3, m4, m5 = azimuth_target_seq_to_BLAST_mismatch_counts[minus_strand_azimuth_seq]   
            if m0 <= max_m0 and m1 <= max_m1 and m2 <= max_m2:
                minus_strand_PAMs[ coord  ] = '-', rc_query_seq, coord, minus_strand_PAM_guide_seq, minus_strand_PAM_guide_seq_with_PAM
    
    return plus_strand_PAMs, minus_strand_PAMs

def T_stretch_present(guide_seq_with_PAM, NUM_CONSECUTIVE_T_DISALLOWED):
    '''
    takes the guide sequence and checks for a stretch of T's as specified by the global variable NUM_CONSECUTIVE_T_DISALLOWED
    returns True if T stretch equal to or greater than NUM_CONSECUTIVE_T_DISALLOWED is found
    '''
    if 'T'*NUM_CONSECUTIVE_T_DISALLOWED in guide_seq_with_PAM:
        return True
    else:
        return False

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
#    if snp_positions == [] and indel_positions == []:
#        print('query, subject, alignment', query, subject, alignment)
    return query, subject, alignment, snp_positions, indel_positions

def process_disruption_positions_for_snps_indels(PAM, snp_positions, indel_positions):
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
    
def get_disruption_position_in_guide_donor_alignment(PAM, guide_with_full_PAM, donor):
    '''
    takes the guide RNA sequence plus full PAM, the PAM, and the donor DNA
    aligs the guide and full PAM to the donor DNA sequence
    returns the best fitting alignment between the guide_with_full_PAM and donor, 
    and the closest position of disruption with 0 denoting PAM disruption, and positive integers denoting distance of SNP/indel from the PAM
    only an indel can have a position of 0, as a snp at position 0 does not disrupt the NGG
    '''
    sense_i, sense_j, sense_backtracking_array, sense_alignment_score = build_backtrack_for_fitting_alignment(donor, guide_with_full_PAM)
    rc_donor = bedgraph_computation.rev_comp(donor)
    antisense_i, antisense_j, antisense_backtracking_array, antisense_alignment_score = build_backtrack_for_fitting_alignment(rc_donor, guide_with_full_PAM)
    if sense_alignment_score > antisense_alignment_score:
        query, subject, alignment, snp_positions, indel_positions = fitting_alignment(sense_i, sense_j, sense_backtracking_array, donor, guide_with_full_PAM)
    else:
        query, subject, alignment, snp_positions, indel_positions = fitting_alignment(antisense_i, antisense_j, antisense_backtracking_array, rc_donor, guide_with_full_PAM)
    disruptions = process_disruption_positions_for_snps_indels(snp_positions, indel_positions)
    return disruptions, query, subject, alignment
    
def global_alignment(v, w, PAM):
    '''
    takes two sequences v, and w, and best fit alignment of all of w against any part of v
    for the best fit alignment, returns the positions of the disruptions, the two sequences, and the alignment
    '''
    v= v.upper()
    w= w.upper()
    sense_i, sense_j, sense_backtracking_array, sense_alignment_score = build_backtrack_for_fitting_alignment(v, w)
    rc_v = bedgraph_computation.rev_comp(v)
    antisense_i, antisense_j, antisense_backtracking_array, antisense_alignment_score = build_backtrack_for_fitting_alignment(rc_v, w)
    if sense_alignment_score > antisense_alignment_score:
        query, subject, alignment, snp_positions, indel_positions = fitting_alignment(sense_i, sense_j, sense_backtracking_array, v, w)
    else:
        query, subject, alignment, snp_positions, indel_positions = fitting_alignment(antisense_i, antisense_j, antisense_backtracking_array, rc_v, w)
    disruptions = process_disruption_positions_for_snps_indels(PAM, snp_positions, indel_positions)
    return disruptions, query, subject, alignment
    
def get_disruption_score(v, w, PAM):
    '''
    disruption scores are calculated by considering that mismatches in the 10-20 nt region count as 1, in the seed region count as 2, and in the PAM 3
    '''
    disruptions, query, subject, alignment = global_alignment(v, w, PAM)
    score = 0
    for position in disruptions:
        if position in range(11,21):
            score += 1
        elif position in range(1,11):
            score += 2
        elif position in (-1, -2):
            score += 3
    return score


        
v ='TATCGCTCAGGTGCCAGCCGGTTTGAAGAGCAGCGGCCGTTTCTGTCACTGCTAGTTCCGATGAAGATATCGCTCAGGTGCCAGCTGAGGCCATTATTGGATACTTGGATTTCTGAGGTGATCATGACATAGCTTTTTTAGCAAAGTCCAAGGACCCGCCGTGGGTTTCTTCCACTCAAACATGTCATGCATCACGTGCTAGCTGACATGACGTCGTTCAGTACCCAGTATACAATAGACAATGGCAGATC'
y ='GCCATTGTCTATTGTATACTGGGTACTGAACGACGTCATGTCAGCTAGCACGTGATGCATGACATGTTTGAGTGGAAGAAACCCACGGCGGGTCCTTGGACTTTGCTAAAAAAGCTATGTCATGATCACCTCAGAAATCCAAGTATCCAATAATGGCCTCAGCTGGCACCTGAGCGATATCTTCATCGGAACTAGCAGTGACAGAAACGGCCGCTGCTCTTCAAACCGGCTGGCACCTGAGCGATCGGCGC'  #y='TGTTTTACCGAAAACATCGGATGGG'
PAM='NGG'
w=y.upper()
disruptions, query, subject, alignment = global_alignment(v, w, PAM)
score = get_disruption_score(v, w, PAM)
output = query + '\n' + subject + '\n' + alignment + '\n' + str(disruptions) + '\n' + str(score)



print(output)
#with open('/Users/kevinroy/Desktop/test.txt', 'w') as outfile:
#    outfile.write(output)
#    outfile.close()

def shift_donor_containing_restriction_sites(genome_seq, chrom, chrom_length, DONOR_LENGTH, variant_start, variant_end, ref_seq, variant_seq, upstream_oligo_seq, downstream_oligo_seq, MINIMUM_HOMOLOGY, RESTRICTION_SITES_SCREENED_OUT, donors_with_restriction_sites_not_allowing_shifting, left_extension, right_extension):
    '''
    takes the donor DNA and the upstream and downstream oligo sequences
    checks if donor DNA contains restriction sites in RESTRICTION_SITES_SCREENED_OUT,
    and shift homologous arms to asymmetric lengths if necessary while maintaining MINIMUM_HOMOLOGY
    returns a shift integer (positive shifts donor to the right, negative to the left)
    '''       
    variant_length = len(variant_seq)
    left_side_donor_length = (DONOR_LENGTH - variant_length) // 2
    right_side_donor_length = DONOR_LENGTH - variant_length - left_side_donor_length
    end_of_chromosome_shift = 0
    at_end_of_chromosome = False
    at_beginning_of_chromosome = False
    if variant_end + right_side_donor_length > chrom_length:
        end_of_chromosome_shift = variant_end + right_side_donor_length - chrom_length
        right_side_donor_length -= end_of_chromosome_shift
        left_side_donor_length += end_of_chromosome_shift
        at_end_of_chromosome = True
    if variant_start - left_side_donor_length < 1:
        end_of_chromosome_shift = left_side_donor_length - variant_start
        left_side_donor_length -= end_of_chromosome_shift
        right_side_donor_length += end_of_chromosome_shift
        at_beginning_of_chromosome = True
    donor = genome_seq[chrom][variant_start - left_side_donor_length: variant_start] + variant_seq + genome_seq[chrom][variant_end: variant_end + right_side_donor_length] 
    target =  genome_seq[chrom][variant_start - left_side_donor_length: variant_start] + ' ' +  ref_seq + ' ' +  genome_seq[chrom][variant_end: variant_end + right_side_donor_length] 
    donor_with_flanking_seq = upstream_oligo_seq + donor + downstream_oligo_seq
    donor_shift = 0
    donor_contains_restriction_site_and_shift_will_not_allow_keeping_minimum_homology = False
    left_side_restriction_sites_present = []
    right_side_restriction_sites_present = []
    for restriction_site in RESTRICTION_SITES_SCREENED_OUT:
        if restriction_site in donor_with_flanking_seq:
            print('donor_harbors_restriction_site: ', restriction_site)
            restriction_site_idx = donor_with_flanking_seq.index(restriction_site)
            left_boundary = len(upstream_oligo_seq) + MINIMUM_HOMOLOGY
            right_boundary = len(donor_with_flanking_seq) - MINIMUM_HOMOLOGY - len(downstream_oligo_seq)
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
                    donors_with_restriction_sites_not_allowing_shifting += 1
                    donor_contains_restriction_site_and_shift_will_not_allow_keeping_minimum_homology = True
                    ## can include a subset of donors in front of the gRNA if the internal cloning site is in the middle of the donor
                    donor = genome_seq[chrom][variant_start - left_side_donor_length - left_extension: variant_start] + variant_seq + genome_seq[chrom][variant_end: variant_end + right_side_donor_length + right_extension] 
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
            donor = genome_seq[chrom][start_of_variant - left_side_donor_length: start_of_variant] + variant_seq + genome_seq[CHROM][end_of_variant: end_of_variant + right_side_donor_length]              
            target =  genome_seq[CHROM][start_of_variant - left_side_donor_length: start_of_variant] + ' ' +  ref_seq + ' ' +  genome_seq[CHROM][end_of_variant: end_of_variant + right_side_donor_length ]
        else:
            print('two restriction sites present on either side of donor', left_side_restriction_sites_present, right_side_restriction_sites_present)