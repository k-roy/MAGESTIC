'''
computation on bedgraph files
'''
import numpy, GTF_GFF_manipulation
import imp
GTF_GFF_manipulation = imp.reload(GTF_GFF_manipulation)

def load_genome(genome_file, Mito = False):
    '''
    input: genome in fasta format
    output: dictionary with chromosome names as keys, where the value is the chromosomal sequence as a string
    '''
    genome = open(genome_file, 'r')
    genome_dict = {}
    for line in genome:
        if 'm' in line or 'M' in line:  ## mito chromosomes are typically last
            break
        elif line[0] == '>':
            current_chromosome = line.strip()[1:]
            genome_dict[current_chromosome] = ''
        else:
            genome_dict[current_chromosome] += line.strip()
    genome.close()
    print((genome_file + " genome file loaded"))
    return genome_dict

def load_custom_genome(genome_file, Mito = False):
    '''
    input: genome in fasta format
    output: dictionary with chromosome names as keys, where the value is the chromosomal sequence as a string
    '''
    genome = open(genome_file, 'r')
    genome_dict = {}
    for line in genome:
        if line[0] == '>':
            current_chromosome = line.strip()[1:]
            genome_dict[current_chromosome] = ''
        else:
            genome_dict[current_chromosome] += line.strip().upper()
    genome.close()
    print((genome_file + " genome file loaded"))
    return genome_dict

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
    
def base_match(base1, base2, ambiguous_bases_allowed = True):
    '''
    takes two IUPAC characters,
    returns True if compatible, False otherwise
    '''
    compatible_bases = {'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT', 'M':'AC', 'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG', 'N':'ACGT'}
    if base1 == base2:
        return True
    elif ambiguous_bases_allowed:
        if base1 in compatible_bases and base2 in compatible_bases[base1]:
            return True
        elif base2 in compatible_bases and base1 in compatible_bases[base2]:
            return True
    else:
        return False

def hamming_distance(string1, string2, ambiguous_bases_allowed = True):
    '''
    input: two strings of equal length
    output: the number of mismatches
    '''
    mismatches = 0
    #if string1[3:] != 'CCCG':
    #    return 3
    for idx in range(len(string1)):
        if not base_match(string1[idx], string2[idx], ambiguous_bases_allowed):
            mismatches += 1
    return mismatches

def edit_distance(v, w):
    '''
      Input: Two nucleotide strings v and w, where v has length at most 1000 and w has length at most 100.
     Output: A highest-scoring fitting alignment between v and w. Use the simple scoring method in which
     matches count +1 and both the mismatch and indel penalties are 1.
    '''
    indel_penalty = 1
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
    edit_distance = max(len(v), len(w)) - alignment_score
    return edit_distance

def obtain_coordinates_for_sequence_in_motif_dict_format(query_sequence, mismatches_allowed, genome_dict, rc=True):
    '''
    takes a sequence and returns a dictionary of strands of chromosomes of lists of (start, stop) tuples for the sequence, and its rc if rc is set to True (True by default)
    motif dict format is: sample_dict[chromosome][ best_motif_coordinates ] = best_motif_strand, best_motif_match, lowest_mismatches 
    '''
    seq_length = len(query_sequence)
    annotations = generate_empty_nonstranded_bedgraph_dict()
    seq_occurence = 0
    rc_seq_occurences = 0
    for chromosome in genome_dict:
        print((chromosome, 'being analyzed for motif', query_sequence, 'with', mismatches_allowed, 'mismatches_allowed'))
        for idx in range( len(genome_dict[chromosome]) - seq_length ):
            subject_sequence = genome_dict[chromosome][idx: idx + seq_length]
            sense_HD = hamming_distance(query_sequence, subject_sequence)
            antisense_HD = hamming_distance(query_sequence, rev_comp(subject_sequence))
            if sense_HD <= mismatches_allowed:
                seq_occurence += 1
                annotations[chromosome][( idx, idx + seq_length )] = '+', subject_sequence, sense_HD
            if rc:
                if antisense_HD <= mismatches_allowed:
                    seq_occurence += 1
                    rc_seq_occurences += 1
                    annotations[chromosome][( idx, idx + seq_length )] = '-', subject_sequence, antisense_HD
        print(('total sense sequence occurrences = ', seq_occurence - rc_seq_occurences))
        print(('total rev comp sequence occurrences = ', rc_seq_occurences))
    print(('total sense sequence occurrences = ', seq_occurence - rc_seq_occurences))
    print(('total rev comp sequence occurrences = ', rc_seq_occurences))
    return annotations

def obtain_coordinates_for_sequence_and_write_outfile(query_sequence, mismatches_allowed, genome_dict, outfilename, rc=False):
    '''
    takes a sequence and returns a dictionary of strands of chromosomes of lists of (start, stop) tuples for the sequence, and its rc if rc is set to True (False by default)
    '''
    outfile = open(outfilename, 'w')
    seq_length = len(query_sequence)
    annotations = {}
    annotations['+'] = {}
    annotations['-'] = {}
    seq_occurence = 0
    rc_seq_occurences = 0
    consensus_sites = 0
    for chromosome in genome_dict:
        print((chromosome, 'being analyzed for motif', query_sequence, 'with', mismatches_allowed, 'mismatches_allowed'))
        if chromosome not in annotations['+']:
                annotations['+'][chromosome] = {}
                annotations['-'][chromosome] = {}
        for idx in range( len(genome_dict[chromosome]) - seq_length ):
            match_found = False
            subject_sequence = genome_dict[chromosome][idx: idx + seq_length]
            sense_HD = hamming_distance(query_sequence, subject_sequence)
            antisense_HD = hamming_distance(query_sequence, rev_comp(subject_sequence))
            if sense_HD <= mismatches_allowed:
                seq_occurence += 1
                annotations['+'][chromosome][seq_occurence] = (idx, idx + seq_length)
                match_found = True
            if rc:
                if antisense_HD <= mismatches_allowed:
                    seq_occurence += 1
                    rc_seq_occurences += 1
                    annotations['-'][chromosome][seq_occurence] = (idx, idx + seq_length)
                    match_found = True
            if match_found:
                if sense_HD < antisense_HD:
                    output = '\t'.join( [chromosome, str(idx), str(idx+ seq_length), str(sense_HD)]  ) + '\n'
                    outfile.write(output)
                else:
                    output = '\t'.join( [chromosome, str(idx), str(idx+ seq_length), str(antisense_HD)]  ) + '\n'
                    outfile.write(output)
                if min(sense_HD, antisense_HD) == 0:
                    consensus_sites += 1
    print('total sense sequence occurrences = ', seq_occurence - rc_seq_occurences)
    print('total rev comp sequence occurrences = ', rc_seq_occurences)
    print('total consensus sites', consensus_sites)
    return annotations
    
def get_range_of_ambiguous_adenosines(chromosome, coordinate, strand, genome_dict):
    '''
    takes a chromosomal coordinate and a genome dict
    returns a 0-based range of coordinates corresponding to the ambiguous adenosines DOWNSTREAM of the coordinate
    '''
    if strand == '+':
        end_coordinate = coordinate + 1
        while genome_dict[chromosome][end_coordinate] == 'A':
            end_coordinate += 1
        return coordinate, end_coordinate
    else:
        end_coordinate = coordinate - 1
        while genome_dict[chromosome][end_coordinate] == 'T':
            end_coordinate -= 1
        return end_coordinate, coordinate

def generate_empty_stranded_bedgraph_dict():
    '''
    takes no arguments
    returns a dictionary of strands of roman numeral-format chromosomes of empty dictionaries
    '''
    dictionary = {}
    for strand in '+-':
        dictionary[strand] = {}
        chromosome_numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI'.split(',')
        for chromosome_numeral in chromosome_numerals:
            dictionary[strand]['chr' + chromosome_numeral] = {}
    return dictionary
    
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
    
def generate_empty_motif_dict():
    '''
    takes no arguments
    returns a dictionary of chromosomes of empty lists
    (motif dictionaries appended with tuples of start, end, motif score)
    '''
    dictionary = {}
    chromosome_numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI'.split(',')
    for chromosome_numeral in chromosome_numerals:
        dictionary['chr' + chromosome_numeral] = []
    return dictionary

def bedgraph_to_dict(filename, chromosomes_to_include = 'ALL', include_mito = False  ):
    '''
    input: bedgraph, with four tab telimited columns: chromosome, start, end, signal
    output: strain dictionary of chromosome dictionaries of coordinate dictionaries with chromosome values
    '''
    infile = open(filename, 'r')
    sample_dict = {}
    for line in infile:
        info = line.strip('\n').split('\t')
        if len(info) == 4:
            chromosome, start, end, reads = info
            start, end = int(start), int(end)
            if chromosomes_to_include == 'ALL':
                if ('M' not in chromosome and 'm' not in chromosome) or include_mito:
                    if chromosome == 'Mito':
                        chromosome = 'chrMito'
                    if chromosome not in sample_dict:
                        sample_dict[chromosome] = {}
                    for idx in range(int(start), int(end)):
                        sample_dict[chromosome][idx] = float(reads)
                    
            elif chromosome in chromosomes_to_include:
                if ('M' not in chromosome and 'm' not in chromosome) or include_mito:
                    if chromosome == 'chrMito':
                        chromosome = 'Mito'
                    if chromosome not in sample_dict:
                        sample_dict[chromosome] = {}
                    for idx in range(int(start), int(end)):
                        sample_dict[chromosome][idx] = float(reads)
        else:
            chromosome, start, end, id_num, reads = info
            start, end = int(start), int(end)
            if chromosomes_to_include == 'ALL':
                if ('M' not in chromosome and 'm' not in chromosome) or include_mito:
                    if chromosome == 'Mito':
                        chromosome = 'chrMito'
                    if chromosome not in sample_dict:
                        sample_dict[chromosome] = {}
                    for idx in range(int(start), int(end)):
                        sample_dict[chromosome][idx] = float(reads)
                    
            elif chromosome in chromosomes_to_include:
                if ('M' not in chromosome and 'm' not in chromosome) or include_mito:
                    if chromosome == 'chrMito':
                        chromosome = 'Mito'
                    if chromosome not in sample_dict:
                        sample_dict[chromosome] = {}
                    for idx in range(int(start), int(end)):
                        sample_dict[chromosome][idx] = float(reads)
                
    infile.close()
    print('bedgraph loaded: ', filename)
    return sample_dict
    
def bedgraph_to_dict_excluding_regions(filename, strand, annotations, chromosomes_to_include = 'ALL', include_mito = False):
    '''
    input: bedgraph, with four tab telimited columns: chromosome, start, end, signal, and an annotations bedgraph with regions to exclude
    output: strain dictionary of chromosome dictionaries of coordinate dictionaries with chromosome values
    '''
    infile = open(filename, 'r')
    sample_dict = {}
    for line in infile:
        info = line.strip('\n').split('\t')
        chromosome, start, end, reads = info
        start, end = int(start), int(end)
        if chromosome not in sample_dict:
            sample_dict[chromosome] = {}
        if chromosomes_to_include == 'ALL':
            if ('M' not in chromosome and 'm' not in chromosome) or include_mito:
                if chromosome == 'Mito':
                    chromosome = 'chrMito'
                if chromosome not in sample_dict:
                    sample_dict[chromosome] = {}

                for idx in range(int(start), int(end)):
                    if chromosome in annotations[strand]:
                        if idx not in annotations[strand][chromosome]:
                            sample_dict[chromosome][idx] = float(reads) 
                    else:
                        sample_dict[chromosome][idx] = float(reads) 
        elif chromosome in chromosomes_to_include:
            if ('M' not in chromosome and 'm' not in chromosome) or include_mito:
                if chromosome == 'chrMito':
                    chromosome = 'Mito'
                if chromosome not in sample_dict:
                    sample_dict[chromosome] = {}
                for idx in range(int(start), int(end)):
                    if chromosome in annotations[strand]:
                        if idx not in annotations[strand][chromosome]:
                            sample_dict[chromosome][idx] = float(reads) 
                    else:
                        sample_dict[chromosome][idx] = float(reads) 
    return sample_dict
        
def generate_upregulated_bedgraphs_with_window(query_plus_strand_bedgraph, query_minus_strand_bedgraph, reference_plus_strand_bedgraph, reference_minus_strand_bedgraph, CLUSTER_WINDOW_SIZE, MIN_FOLD_UPREGULATION_OF_CLUSTER, MIN_READS_FOR_PEAK, min_read_value_in_place_of_zero = 0.01):
    '''
    takes query and reference bedgraphs (both strands), a window size to check upregulation, a minimum fold upgegulation, and if applicable, a normalization factor
    computes coordinates upregulated over the window size at least min fold, and returns upregulated bedgraphs (a tuple of both strands)
    '''
    query_total_reads = sum_total_bedgraph_reads(query_plus_strand_bedgraph, query_minus_strand_bedgraph)
    reference_total_reads = sum_total_bedgraph_reads(reference_plus_strand_bedgraph, reference_minus_strand_bedgraph)
    query_normalization = float(reference_total_reads) / query_total_reads 
    print('reference_total_reads', reference_total_reads, 'query_total_reads', query_total_reads, 'query_normalization', query_normalization)
    total_upregulated_sites = 0
    upregulated_plus_strand_bedgraph = {}
    for chromosome in query_plus_strand_bedgraph:
        if chromosome not in upregulated_plus_strand_bedgraph:
            upregulated_plus_strand_bedgraph[chromosome] = {}
        for coordinate in query_plus_strand_bedgraph[chromosome]:
            query_reads = float( compute_expression_value(chromosome, coordinate - CLUSTER_WINDOW_SIZE, coordinate + CLUSTER_WINDOW_SIZE, query_plus_strand_bedgraph) ) * query_normalization
            reference_reads = float( max( compute_expression_value(chromosome, coordinate - CLUSTER_WINDOW_SIZE, coordinate + CLUSTER_WINDOW_SIZE, reference_plus_strand_bedgraph), min_read_value_in_place_of_zero ) )
            if float( query_reads ) / reference_reads >= MIN_FOLD_UPREGULATION_OF_CLUSTER and query_plus_strand_bedgraph[chromosome][coordinate] >= MIN_READS_FOR_PEAK:
                upregulated_plus_strand_bedgraph[chromosome][coordinate] = query_plus_strand_bedgraph[chromosome][coordinate]
                total_upregulated_sites += 1
    upregulated_minus_strand_bedgraph = {}
    for chromosome in query_minus_strand_bedgraph:
        if chromosome not in upregulated_minus_strand_bedgraph:
            upregulated_minus_strand_bedgraph[chromosome] = {}
        for coordinate in query_minus_strand_bedgraph[chromosome]:
            query_reads = float( compute_expression_value(chromosome, coordinate - CLUSTER_WINDOW_SIZE, coordinate + CLUSTER_WINDOW_SIZE, query_minus_strand_bedgraph) ) * query_normalization
            reference_reads = float( max( compute_expression_value(chromosome, coordinate - CLUSTER_WINDOW_SIZE, coordinate + CLUSTER_WINDOW_SIZE, reference_minus_strand_bedgraph), min_read_value_in_place_of_zero ) )
            if float( query_reads ) / reference_reads >= MIN_FOLD_UPREGULATION_OF_CLUSTER and query_minus_strand_bedgraph[chromosome][coordinate] >= MIN_READS_FOR_PEAK:
                upregulated_minus_strand_bedgraph[chromosome][coordinate] = query_minus_strand_bedgraph[chromosome][coordinate]
                total_upregulated_sites += 1
    print('total_upregulated_sites =', total_upregulated_sites)
    return upregulated_plus_strand_bedgraph, upregulated_minus_strand_bedgraph
    
def generate_upregulated_bedgraphs(query_plus_strand_bedgraph,query_minus_strand_bedgraph, reference_plus_strand_bedgraph, reference_minus_strand_bedgraph, MIN_FOLD_UPREGULATION, MIN_READS_FOR_PEAK, min_read_value_in_place_of_zero = 0.01):
    '''
    takes query and reference bedgraphs (both strands), a minimum fold upgegulation, and if applicable, a normalization factor, and a value to use for the denominator in case the reference has zero reads
    computes coordinates upregulated at least min fold, and returns upregulated bedgraphs (a tuple of both strands)
    each coordinate assessed independently of those around it, i.e. no window required
    '''
    query_total_reads = sum_total_bedgraph_reads(query_plus_strand_bedgraph, query_minus_strand_bedgraph)
    reference_total_reads = sum_total_bedgraph_reads(reference_plus_strand_bedgraph, reference_minus_strand_bedgraph)
    query_normalization = query_total_reads / float(reference_total_reads)
    print('reference_total_reads', reference_total_reads, 'query_total_reads', query_total_reads, 'query_normalization', query_normalization)
    total_upregulated_sites = 0
    upregulated_plus_strand_bedgraph = {}
    for chromosome in query_plus_strand_bedgraph:
        if chromosome not in upregulated_plus_strand_bedgraph:
            upregulated_plus_strand_bedgraph[chromosome] = {}
        for coordinate in query_plus_strand_bedgraph[chromosome]:
            query_reads = float( query_plus_strand_bedgraph[chromosome][coordinate] ) * query_normalization
            reference_reads = float( reference_plus_strand_bedgraph[chromosome].get( coordinate, min_read_value_in_place_of_zero ) )
            if float( query_reads ) >= MIN_READS_FOR_PEAK and float( query_reads ) / reference_reads >= MIN_FOLD_UPREGULATION:
                upregulated_plus_strand_bedgraph[chromosome][coordinate] = query_plus_strand_bedgraph[chromosome][coordinate]
                total_upregulated_sites += 1
    upregulated_minus_strand_bedgraph = {}
    for chromosome in query_minus_strand_bedgraph:
        if chromosome not in upregulated_minus_strand_bedgraph:
            upregulated_minus_strand_bedgraph[chromosome] = {}
        for coordinate in query_minus_strand_bedgraph[chromosome]:
            query_reads = float( query_minus_strand_bedgraph[chromosome][coordinate] ) * query_normalization
            reference_reads = float( reference_minus_strand_bedgraph[chromosome].get( coordinate, min_read_value_in_place_of_zero ) )
            if float( query_reads ) >= MIN_READS_FOR_PEAK and float( query_reads ) / reference_reads >= MIN_FOLD_UPREGULATION:
                upregulated_minus_strand_bedgraph[chromosome][coordinate] = query_minus_strand_bedgraph[chromosome][coordinate]
                total_upregulated_sites += 1
                
    print('total_upregulated_sites =', total_upregulated_sites)
    return upregulated_plus_strand_bedgraph, upregulated_minus_strand_bedgraph
   
def intersect_bedgraphs(query1_plus_strand_bedgraph, query1_minus_strand_bedgraph, query2_plus_strand_bedgraph, query2_minus_strand_bedgraph, bp_buffer_for_peaks):
    '''
    takes two bedgraphs (both strands)
    returns a new bedgraph with the intersection of the two bedgraphs, using query1 for the read values for each coordiante
    i.e. how many peaks in one bedgraph are also in another, and allow some wiggle room for exact coordinates for peaks
    '''
    intersect_plus_strand, intersect_minus_strand = generate_empty_nonstranded_bedgraph_dict(), generate_empty_nonstranded_bedgraph_dict()
    total_query1_sites = 0
    overlap_sites = 0
    for chromosome in query1_plus_strand_bedgraph:
        for coordinate in query1_plus_strand_bedgraph[chromosome]:
            total_query1_sites += 1
            overlap = False
            for buffer_coord in range(coordinate - bp_buffer_for_peaks, coordinate + bp_buffer_for_peaks + 1):
                if coordinate in query2_plus_strand_bedgraph[chromosome]:
                    overlap = True
            if overlap:
                intersect_plus_strand[chromosome][coordinate] = float( query1_plus_strand_bedgraph[chromosome][coordinate] )
                overlap_sites += 1
    for chromosome in query1_minus_strand_bedgraph:
        for coordinate in query1_minus_strand_bedgraph[chromosome]:
            total_query1_sites += 1
            overlap = False
            for buffer_coord in range(coordinate - bp_buffer_for_peaks, coordinate + bp_buffer_for_peaks + 1):
                if coordinate in query2_minus_strand_bedgraph[chromosome]:
                    overlap = True
            if overlap:
                intersect_minus_strand[chromosome][coordinate] = float( query1_minus_strand_bedgraph[chromosome][coordinate] )
                overlap_sites += 1
    print('total_query1_sites =', total_query1_sites)
    print('overlap_sites =', overlap_sites)
    return intersect_plus_strand, intersect_minus_strand 
       

def compute_upregulated_and_not_upregulated_bedgraph_peaks(query_plus_strand_bedgraph,query_minus_strand_bedgraph, reference_plus_strand_bedgraph, reference_minus_strand_bedgraph, MIN_FOLD_UPREGULATION, MIN_READS_FOR_PEAK, parameters_outfile, min_read_value_in_place_of_zero = 0.01, MAX_FOLD_UPREGULATION = 1,):
    '''
    parses a query bedgraph into upregulated and not upregulated peaks relative to a reference bedgraph
    takes query and reference bedgraphs (both strands), a minimum fold upgegulation, and if applicable, a normalization factor, and a value to use for the denominator in case the reference has zero reads
    computes coordinates upregulated at least min fold, and returns upregulated bedgraphs (a tuple of both strands)
    '''
    query_total_reads = sum_total_bedgraph_reads(query_plus_strand_bedgraph, query_minus_strand_bedgraph)
    reference_total_reads = sum_total_bedgraph_reads(reference_plus_strand_bedgraph, reference_minus_strand_bedgraph)
    query_normalization = query_total_reads / float(reference_total_reads)
    print('reference_total_reads', reference_total_reads, 'query_total_reads', query_total_reads, 'query_normalization', query_normalization)
    total_sites = 0
    total_upregulated_sites = 0
    total_not_upregulated_sites = 0
    upregulated_plus_strand_bedgraph = {}
    not_upregulated_plus_strand_bedgraph = {}
    for chromosome in query_plus_strand_bedgraph:
        if chromosome not in upregulated_plus_strand_bedgraph:
            upregulated_plus_strand_bedgraph[chromosome] = {}
        if chromosome not in not_upregulated_plus_strand_bedgraph:
            not_upregulated_plus_strand_bedgraph[chromosome] = {}
        for coordinate in query_plus_strand_bedgraph[chromosome]:
            query_reads = float( query_plus_strand_bedgraph[chromosome][coordinate] ) * query_normalization
            reference_reads = float( reference_plus_strand_bedgraph[chromosome].get( coordinate, min_read_value_in_place_of_zero ) )
            if float( query_reads ) >= MIN_READS_FOR_PEAK:
                total_sites += 1
                if float( query_reads ) / reference_reads >= MIN_FOLD_UPREGULATION:
                    upregulated_plus_strand_bedgraph[chromosome][coordinate] = query_plus_strand_bedgraph[chromosome][coordinate]
                    total_upregulated_sites += 1
                elif float( query_reads ) / reference_reads <= MAX_FOLD_UPREGULATION:
                    not_upregulated_plus_strand_bedgraph[chromosome][coordinate] = query_plus_strand_bedgraph[chromosome][coordinate]
                    total_not_upregulated_sites += 1
    upregulated_minus_strand_bedgraph = {}
    not_upregulated_minus_strand_bedgraph = {}
    for chromosome in query_minus_strand_bedgraph:
        if chromosome not in upregulated_minus_strand_bedgraph:
            upregulated_minus_strand_bedgraph[chromosome] = {}
        if chromosome not in not_upregulated_minus_strand_bedgraph:
            not_upregulated_minus_strand_bedgraph[chromosome] = {}
        for coordinate in query_minus_strand_bedgraph[chromosome]:
            query_reads = float( query_minus_strand_bedgraph[chromosome][coordinate] ) * query_normalization
            reference_reads = float( reference_minus_strand_bedgraph[chromosome].get( coordinate, min_read_value_in_place_of_zero ) )
            if float( query_reads ) >= MIN_READS_FOR_PEAK:
                total_sites += 1
                if float( query_reads ) / reference_reads >= MIN_FOLD_UPREGULATION:
                    upregulated_minus_strand_bedgraph[chromosome][coordinate] = query_minus_strand_bedgraph[chromosome][coordinate]
                    total_upregulated_sites += 1
                elif float( query_reads ) / reference_reads <= MAX_FOLD_UPREGULATION:
                    not_upregulated_minus_strand_bedgraph[chromosome][coordinate] = query_minus_strand_bedgraph[chromosome][coordinate]
                    total_not_upregulated_sites += 1
    print('total_sites =', total_sites)
    print('total_upregulated_sites =', total_upregulated_sites)
    print('total_not_upregulated_sites =', total_not_upregulated_sites)
    parameters_outfile.write('total_sites = ' + str(total_sites) + '\n' + 'total_upregulated_sites = ' + str(total_upregulated_sites) + '\n' + 'total_not_upregulated_sites = ' + str(total_not_upregulated_sites) + '\n')
    return upregulated_plus_strand_bedgraph, upregulated_minus_strand_bedgraph, not_upregulated_plus_strand_bedgraph, not_upregulated_minus_strand_bedgraph
   
def sum_total_bedgraph_reads(bedgraph_dict_1, bedgraph_dict_2):
    '''
    input: a bedgraph dictionary (with 1 and 2 denoting plus and minus strands of the sample)
    output: the total signal in the dictionary
    '''
    total = 0
    for bedgraph_dict in bedgraph_dict_1, bedgraph_dict_2:
        for chromosome in bedgraph_dict:
            for idx in bedgraph_dict[chromosome]:
                total += bedgraph_dict[chromosome][idx]
    return total

def sum_nonstranded_bedgraph_reads(bedgraph_dict):
    '''
    input: a bedgraph dictionary (with 1 and 2 denoting plus and minus strands of the sample)
    output: the total signal in the dictionary
    '''
    total = 0
    for chromosome in bedgraph_dict:
        for idx in bedgraph_dict[chromosome]:
            total += bedgraph_dict[chromosome][idx]
    return total
    
def compute_expression_value(chromosome, start, end, bedgraph_dict):
    '''
    input: a query chromosome, start, and end coordinates, and a bedgraph_dict
    output: the sum of the expression value over the chromosomal coordinate range
    '''
    expression_value = 0
    if chromosome in bedgraph_dict:
        for idx in range(start, end):
            if idx in bedgraph_dict[chromosome]:
                expression_value += bedgraph_dict[chromosome][idx]
    return float(expression_value)

def return_expression_values(chromosome, start, end, bedgraph_dict, excluded_coordinates = [], pseudo_count_in_place_of_zero = 0):
    '''
    input: a query chromosome, start, and end coordinates, and a bedgraph_dict
    output: a list of the expression values over the chromosomal coordinate range
    '''
    expression_values = []
    for idx in range(start, end):
        if idx not in excluded_coordinates:
            expression_values.append(bedgraph_dict[chromosome].get(idx, pseudo_count_in_place_of_zero ))
        else:
            expression_values.append( pseudo_count_in_place_of_zero )
    return expression_values
    
def return_coordinate_with_max_reads_in_range(chromosome, start, end, bedgraph_dict):
    '''
    input: a query chromosome, start, and end coordinates, and a bedgraph_dict
    output: the coordinate with the max value over the chromosomal coordinate range
    '''
    coord_with_max_reads = None
    max_reads = 0
    for idx in range(start, end):
        reads = bedgraph_dict[chromosome].get(idx, 0)
        if reads > max_reads:
            max_reads = reads
            coord_with_max_reads = idx
    return coord_with_max_reads

   
def return_upregulated_expression_values(chromosome, start, end, bedgraph_dict_query, bedgraph_dict_subject, min_fold_upregulation ):
    '''
    input: a query chromosome, start, and end coordinates, and query bedgraph_dict, and a subject bedgraph dict
    output: a list of the expression values over the chromosomal coordinate range that are upregulated in the query,
    those that are not are given a zero value
    '''
    expression_values = []
    for idx in range(start, end):
        query = bedgraph_dict_query[chromosome].get(idx, 0)
        subject = bedgraph_dict_subject[chromosome].get(idx, 0)
        if subject == 0 or query / float(subject) >= min_fold_upregulation:
            expression_values.append(query)
        else:
            expression_values.append(0)
    return expression_values


def compute_differential_expression(chromosome, start, end, sample_A_replicates, sample_B_replicates, strand, gene, minimum_average_signal_per_base = 0):
    '''
    input: a query chromosome, start, and end coordinates, and a list of tuples of bedgraph dictionaries and normalization signals for sample_A_replicates and sample_B_replicates
    output: string representations of the difference in expression of sample_A over sample_B, mean signal per base across feature for A, mean signal per base across feature for B, the coefficient of variation for sample_A replicates, and the coefficient of variation  for sample_B replicates
    (the coefficient of variation is defined as the standard deviation divided by the mean, and is a relative measure of the reproducibility of the signal across replicates)
    USING READ COUNTS PER GENE RATHER THAN COVERAGE COULD BE BETTER    
    '''
    if chromosome == 'ChrMito':
        chromosome = 'chrMito'
    sample_A_expression_values = []
    sample_B_expression_values = []
    feature_length = end - start
#    print end, start, feature_length
    for replicate in sample_A_replicates:
        sample_A_expression_values.append( compute_expression_value(chromosome, start, end, replicate[0]) * replicate[1] )
    for replicate in sample_B_replicates:
        sample_B_expression_values.append( compute_expression_value(chromosome, start, end, replicate[0]) * replicate[1] )
    sample_A_mean = float ( sum ( sample_A_expression_values) ) /len (sample_A_replicates)
    sample_B_mean = float ( sum ( sample_B_expression_values) ) /len (sample_B_replicates)
    if sample_A_mean/feature_length < minimum_average_signal_per_base and sample_B_mean/feature_length < minimum_average_signal_per_base:
        return ['no sample B or A signal', 'NAN', 'NAN','NAN','NAN','NAN']
    sample_A_expression_values_array, sample_B_expression_values_array = numpy.array(sample_A_expression_values) , numpy.array(sample_B_expression_values)
    if sample_B_mean == 0:
        sample_A_coefficient_of_variation = numpy.std(sample_A_expression_values_array, axis = 0)/sample_A_mean
        return ['no sample B signal', str(sample_A_mean/feature_length), str(sample_A_mean/feature_length), 'NAN',  str(sample_A_coefficient_of_variation), 'NAN']
    elif sample_A_mean == 0:
        sample_B_coefficient_of_variation = numpy.std(sample_B_expression_values_array, axis = 0)/sample_B_mean
        return ['no sample A signal', str(sample_B_mean/feature_length), 'NAN', str(sample_B_mean/feature_length),  'NAN', str(sample_B_coefficient_of_variation)]
    sample_A_coefficient_of_variation = numpy.std(sample_A_expression_values_array, axis = 0)/sample_A_mean
    sample_B_coefficient_of_variation = numpy.std(sample_B_expression_values_array, axis = 0)/sample_B_mean
    differential_expression_for_A_over_B = sample_A_mean / sample_B_mean
    return (differential_expression_for_A_over_B, gene, chromosome, start, end, strand, (sample_A_mean/feature_length - sample_B_mean/feature_length ), sample_A_mean/feature_length, sample_B_mean/feature_length, sample_A_coefficient_of_variation, sample_B_coefficient_of_variation)

def write_bedgraph_dict_to_file(bedgraph_dict, outfile_name, normalization_factor = 1):
    '''
    input: bedgraph dictionary, name for outfile
    output: file version of bedgraph dicitonary
    '''
    outfile = open(outfile_name, 'w')
    for chromosome in bedgraph_dict:
        sorted_coordinates = sorted(bedgraph_dict[chromosome].keys())
        for coordinate in sorted_coordinates:
            reads = float(bedgraph_dict[chromosome][coordinate])
            #print chromosome, coordinate, reads
            output = chromosome + '\t' + str(coordinate) + '\t' + str(coordinate + 1) + '\t' + str(reads * normalization_factor) + '\n'
            #print output
            outfile.write(output)
    outfile.close()
    print(outfile_name + ' complete')
    return None
    
def calculate_center_of_mass(chromosome, start, end, bedgraph):
    '''
    takes a chromosomal range and a bedgraph
    returns the coordinate representing the center of mass of pA sites,
    no counts in the region returns a zero
    '''
    weighted_coordinate_sum = 0
    total_pA_reads = 0
    for coord in range(start, end):
        pA_reads = bedgraph[chromosome].get(coord, 0)
        total_pA_reads += pA_reads
        weighted_coordinate = pA_reads * coord
        weighted_coordinate_sum += weighted_coordinate
    if total_pA_reads > 0:
        return weighted_coordinate_sum / float(total_pA_reads)
    else:
        return 0

def check_fold_upregulation(chromosome, start, end, query_pA_sites, subject_pA_sites, MIN_REFERENCE_VALUE_TO_DIVIDE_BY):
    '''
    takes a chromosomal range and query and subject bedgraph
    returns the fold difference (float) in level of expression over the range of query over subject
    '''
    query_expression_level = compute_expression_value(chromosome, start, end, query_pA_sites)
    subject_expression_level = max( compute_expression_value(chromosome, start, end, subject_pA_sites), MIN_REFERENCE_VALUE_TO_DIVIDE_BY)
    ## print 'query_expression_level', query_expression_level, 'subject_expression_level', subject_expression_level
    return float(query_expression_level) / subject_expression_level
    
def extend_upregulated_cluster(chromosome, coordinate, query_pA_sites, subject_pA_sites, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_FOLD_UPREGULATION_OF_CLUSTER, MIN_REFERENCE_VALUE_TO_DIVIDE_BY, max_coordinate):
    '''
    takes a chromosomal coordinate in query pA sites, and two bedgraphs (query and subject pA sites), a coordinates examined list of integers
    computes an extended cluster, by agglomerating sites within FLANKING_NUCLEOTIDES (in both directions) as long as the region maintains MIN_FOLD_UPREGULATION_OF_CLUSTER
    if reference has no reads, gives an integer MIN_REFERENCE_VALUE_TO_DIVIDE_BY, and the max_coordindate encountered in the chromosome
    returns a tuple of cluster range, the distinct reads within the cluster, and a revised coordinates_examined list
    '''
    ## extend cluster to the left
    current_cluster_range = [coordinate, coordinate + 1]
    current_cluster_reads = [ query_pA_sites[chromosome][coordinate] ]
    previous_left_coord = coordinate
    coordinates_examined.append(coordinate)
    for left_coord in range(coordinate-1, 0, -1):
        if abs(previous_left_coord - left_coord) > FLANKING_NUCLEOTIDES:
            break
        if left_coord in query_pA_sites[chromosome]:
            if ( left_coord not in coordinates_examined ) and check_fold_upregulation(chromosome, left_coord, coordinate, query_pA_sites, subject_pA_sites, MIN_REFERENCE_VALUE_TO_DIVIDE_BY) >= MIN_FOLD_UPREGULATION_OF_CLUSTER:
                current_cluster_range[0] = left_coord
                previous_left_coord = left_coord
                coordinates_examined.append(left_coord)
                current_cluster_reads = [ query_pA_sites[chromosome][left_coord] ] + current_cluster_reads    
            else:
                break
    ## extend cluster to the right
    previous_right_coord = coordinate
    for right_coord in range(coordinate+1, max_coordinate + 1):
        if abs(right_coord - previous_right_coord) > FLANKING_NUCLEOTIDES:
            break
        if right_coord in query_pA_sites[chromosome]:
            if ( right_coord not in coordinates_examined ) and check_fold_upregulation(chromosome, coordinate, right_coord, query_pA_sites, subject_pA_sites, MIN_REFERENCE_VALUE_TO_DIVIDE_BY) >= MIN_FOLD_UPREGULATION_OF_CLUSTER:
                current_cluster_range[1] = right_coord
                previous_right_coord = right_coord
                coordinates_examined.append(right_coord)
                current_cluster_reads = current_cluster_reads + [ query_pA_sites[chromosome][right_coord] ]
            else:
                break
    ## print 'coordinate ', chromosome, coordinate, 'extended from ', previous_left_coord, 'to', previous_right_coord
    return current_cluster_range, current_cluster_reads, coordinates_examined   
    
def extend_cluster(chromosome, coordinate, query_pA_sites, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_READS_TO_EXTEND_CLUSTER, max_coordinate):
    '''
    takes a chromosomal coordinate in query pA sites, and two bedgraphs (query and subject pA sites)
    computes an extended cluster, by agglomerating sites within FLANKING_NUCLEOTIDES (in both directions) as long as the region maintains MIN_FOLD_UPREGULATION_OF_CLUSTER
    returns a tuple of cluster range, the distinct reads within the cluster, and a revised coordinates_examined list
    '''
    
    ## extend cluster to the left
    current_cluster_range = [coordinate, coordinate + 1]
    current_cluster_reads = [ query_pA_sites[chromosome][coordinate] ]
    previous_left_coord = coordinate
    coordinates_examined.append(coordinate)
    for left_coord in range(coordinate-1, 0, -1):
        if abs(previous_left_coord - left_coord) > FLANKING_NUCLEOTIDES:
            break
        if left_coord in query_pA_sites[chromosome]:
            if ( left_coord not in coordinates_examined ):
                current_cluster_range[0] = left_coord
                previous_left_coord = left_coord
                coordinates_examined.append(left_coord)
                current_cluster_reads = [ query_pA_sites[chromosome][left_coord] ] + current_cluster_reads    
            else:
                break
    ## extend cluster to the right
    previous_right_coord = coordinate
    for right_coord in range(coordinate+1, max_coordinate + 1):
        if abs(right_coord - previous_right_coord) > FLANKING_NUCLEOTIDES:
            break
        if right_coord in query_pA_sites[chromosome]:
            if ( right_coord not in coordinates_examined ):
                current_cluster_range[1] = right_coord
                previous_right_coord = right_coord
                coordinates_examined.append(right_coord)
                current_cluster_reads = current_cluster_reads + [ query_pA_sites[chromosome][right_coord] ]
            else:
                break
    ## print 'coordinate ', chromosome, coordinate, 'extended from ', previous_left_coord, 'to', previous_right_coord
    return current_cluster_range, current_cluster_reads, coordinates_examined
    
def determine_range_around_peak_accounting_for_signal(bedgraph, chromosome, start, end, current_cluster_reads, FRACTION_OF_TOTAL_CLUSTER_SIGNAL):
    '''
    takes a genomic window plus FLANKING_NUCLEOTIDES
    start with trimming either the left coord or right coord, whichever has lower reads
    stop when the window drops below PERCENTAGE_OF_TOTAL_CLUSTER_SIGNAL percentage of the total signal in that window
    returns how many of the coordinates around the peak are required to explain PERCENTAGE_OF_TOTAL_CLUSTER_SIGNAL percentage of the total signal in that window
    '''
    total_signal_in_range = float ( sum(current_cluster_reads) )
    ## find the max peak in the range
    most_reads = 0
    for coord in range(start, end + 1):
        if coord in bedgraph[chromosome]:
            reads = bedgraph[chromosome][coord]
            if reads > most_reads:
                most_reads = reads
                max_peak = coord
    right_side_running_total = 0
    left_side_running_total = 0
    right_coordinate = end
    left_coordinate = start
    distinct_peaks_accounting_for_signal = 1
    while ( right_side_running_total + left_side_running_total) / float(total_signal_in_range) <= ( 1 - FRACTION_OF_TOTAL_CLUSTER_SIGNAL ):
        right_coord_reads = bedgraph[chromosome].get(right_coordinate, 0)
        left_coord_reads = bedgraph[chromosome].get(left_coordinate, 0)
        while left_coord_reads == 0:
            left_coordinate += 1
            left_coord_reads = bedgraph[chromosome].get(left_coordinate, 0)
        while right_coord_reads == 0:
            right_coordinate -= 1
            right_coord_reads = bedgraph[chromosome].get(right_coordinate, 0)
        if right_coord_reads + right_side_running_total < left_coord_reads + left_side_running_total:
            right_coordinate -= 1
            right_side_running_total += right_coord_reads
        else:
            left_coordinate += 1
            left_side_running_total += left_coord_reads
            
        if (right_side_running_total + left_side_running_total) / total_signal_in_range >= (1-FRACTION_OF_TOTAL_CLUSTER_SIGNAL):
            reads = []
            for coordinate in range(left_coordinate, right_coordinate + 1):
                if coordinate in bedgraph[chromosome]:
                    reads.append(bedgraph[chromosome][coordinate])
                    distinct_peaks_accounting_for_signal += 1
    if distinct_peaks_accounting_for_signal == 0:
        print('0 distinct peaks in determine_range_around_peak_accounting_for_signal at peak ' + chromosome + ' ' + str(start) + ' to ' + str(end))
    return most_reads, max_peak, right_coordinate - left_coordinate, distinct_peaks_accounting_for_signal
    
def cluster_pA_sites(query_pA_sites_plus_strand, query_pA_sites_minus_strand, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types, outfilename, FLANKING_NUCLEOTIDES,  FRACTION_OF_TOTAL_CLUSTER_SIGNAL, MIN_READS_TO_EXTEND_CLUSTER = 1):
    '''
    takes two bedgraphs for the plus and minus strands, and clusters the sites
    writes each cluster to a file as the chromosome, start and end coordinates,  peak_coord, width_around_peak_accounting_for_signal, number_of_distinct_peaks, sequence 100 nt upstream, the sequence of the cluster, sequence 100 nt downstream
    returns the total number of clusters identified
    '''
    with open(outfilename, 'w') as outfile:
        output = '\t'.join( 'chromosome, max_peak_start, max_peak_end, strand, annotation_type, gene_name, strandedness, start, end, cluster_width, total_cluster_reads, max_peak_height,  width_around_peak_accounting_for_80%_of_total_signal, distinct_peaks_accounting_for_80%_of_total_signal'.split(', ') ) + '\n'
        outfile.write(output)
        print('cluster_pA_sites started')
        plus_strand_clusters = {}
        minus_strand_clusters = {}
        total_clusters = 0
        strand = '+'
        for chromosome in query_pA_sites_plus_strand:
            print(chromosome, strand)
            if chromosome not in plus_strand_clusters:
                plus_strand_clusters[chromosome] = {}
            coordinates_examined = []
            if query_pA_sites_plus_strand[chromosome] != {}:
                max_coordinate = max( query_pA_sites_plus_strand[chromosome].keys() )
                for coordinate in sorted(query_pA_sites_plus_strand[chromosome].keys() ):
                    if coordinate not in coordinates_examined:
                        total_clusters += 1
                        current_cluster_range, current_cluster_reads, coordinates_examined = extend_cluster(chromosome, coordinate, query_pA_sites_plus_strand, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_READS_TO_EXTEND_CLUSTER, max_coordinate)      
                        start, end = current_cluster_range
                        most_reads, max_peak, width_around_peak_accounting_for_signal, distinct_peaks_accounting_for_signal = determine_range_around_peak_accounting_for_signal(query_pA_sites_plus_strand, chromosome, start, end, current_cluster_reads, FRACTION_OF_TOTAL_CLUSTER_SIGNAL)             
                        plus_strand_clusters[chromosome][max_peak] =   float ( sum(current_cluster_reads) )  
                        annotation_type, gene_name, strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, max_peak, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types) 
                        ## insert the max_peak_coordinate into the clusters dictionary            
                        output = '\t'.join( [chromosome, str(max_peak), str(max_peak + 1), strand, annotation_type, gene_name, strandedness, str(start), str(end), str(end-start), str(sum(current_cluster_reads)), str(most_reads),  str(width_around_peak_accounting_for_signal), str(distinct_peaks_accounting_for_signal)  ] ) + '\n'
                        ## coordinates followed by the total range (in nt) of the cluster, the total reads in cluster, the number of reads for the peak coordinate in the cluster, the peak coordinate in the cluster,
                        ## the width around peak accounting for the signal, number of distinct peaks within the width accounting for the signal
                        outfile.write(output)
            else:
                print (chromosome, strand, 'empty')
        strand = '-'
        for chromosome in query_pA_sites_minus_strand:
            print(chromosome, strand)
            if chromosome not in minus_strand_clusters:
                minus_strand_clusters[chromosome] = {}
            coordinates_examined = []
            if query_pA_sites_minus_strand[chromosome] != {}:
                max_coordinate = max( query_pA_sites_minus_strand[chromosome].keys() )
                for coordinate in sorted(query_pA_sites_minus_strand[chromosome].keys() ):
                    ##print 'coordinate', coordinate, 'coordinates_examined', coordinates_examined
                    if coordinate not in coordinates_examined:
                        total_clusters += 1
                        current_cluster_range, current_cluster_reads, coordinates_examined = extend_cluster(chromosome, coordinate, query_pA_sites_minus_strand, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_READS_TO_EXTEND_CLUSTER, max_coordinate)      
                        start, end = current_cluster_range
                        most_reads, max_peak, width_around_peak_accounting_for_signal, distinct_peaks_accounting_for_signal = determine_range_around_peak_accounting_for_signal(query_pA_sites_minus_strand, chromosome, start, end, current_cluster_reads, FRACTION_OF_TOTAL_CLUSTER_SIGNAL)             
                        minus_strand_clusters[chromosome][max_peak] =   float ( sum(current_cluster_reads) )  
                        annotation_type, gene_name, strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, max_peak, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types) 
                        ## insert the max_peak_coordinate into the clusters dictionary            
                        output = '\t'.join( [chromosome, str(max_peak), str(max_peak + 1), strand, annotation_type, gene_name, strandedness, str(start), str(end), str(end-start), str(sum(current_cluster_reads)), str(most_reads),  str(width_around_peak_accounting_for_signal), str(distinct_peaks_accounting_for_signal)  ] ) + '\n'
                        ## coordinates followed by the total range (in nt) of the cluster, the total reads in cluster, the number of reads for the peak coordinate in the cluster, the peak coordinate in the cluster,
                        ## the width around peak accounting for the signal, number of distinct peaks within the width accounting for the signal
                        outfile.write(output)
            else:
                print (chromosome, strand, 'empty')
        print(total_clusters, ' total clusters for ', outfilename)
    outfile.close()
    return plus_strand_clusters, minus_strand_clusters

def cluster_upregulated_pA_sites(upregulated_pA_sites_plus_strand, query_pA_sites_plus_strand, subject_pA_sites_plus_strand, upregulated_pA_sites_minus_strand, query_pA_sites_minus_strand, subject_pA_sites_minus_strand, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types, outfilename, FLANKING_NUCLEOTIDES, MIN_FOLD_UPREGULATION_OF_CLUSTER, MIN_REFERENCE_VALUE_TO_DIVIDE_BY, FRACTION_OF_TOTAL_CLUSTER_SIGNAL):
    '''
    takes three bedgraphs, a query upregulated bedgraph, the query sites to cluster around the upregulated sites, and the subject pA sites to test for upregulation against
    writes each cluster to a file as the chromosome, start and end coordinates,  peak_coord, width_around_peak_accounting_for_signal, number_of_distinct_peaks, sequence 100 nt upstream, the sequence of the cluster, sequence 100 nt downstream
    returns the total number of clusters identified
    '''
    with open(outfilename, 'w') as outfile:
        header = '\t'.join( ['chromosome', 'peak coord', 'center of mass coordinate for cluster', 'strand', 'annotation type', 'gene name', 'gene ID', 'strandedness', 'cluster start', 'cluster end', 'cluster width', 'total cluster reads', 'peak height (reads)',  'width around peak accounting for ' + str(FRACTION_OF_TOTAL_CLUSTER_SIGNAL*100) + '% of signal', ' distinct number of three prime ends accounting for ' + str(FRACTION_OF_TOTAL_CLUSTER_SIGNAL*100) + '% of signal', 'fold upregulation' ] ) + '\n'
        outfile.write(header)        
        print('cluster_upregulated_pA_sites started')
        plus_strand_clusters = {}
        minus_strand_clusters = {}
        total_clusters = 0
        strand = '+'
        for chromosome in upregulated_pA_sites_plus_strand:
            print(chromosome, strand)
            if chromosome not in plus_strand_clusters:
                plus_strand_clusters[chromosome] = {}
    
            coordinates_examined = []
            max_coordinate = max( upregulated_pA_sites_plus_strand[chromosome].keys() )
            for coordinate in sorted(upregulated_pA_sites_plus_strand[chromosome].keys() ):
                if coordinate not in coordinates_examined:
                    
                    current_cluster_range, current_cluster_reads, coordinates_examined = extend_upregulated_cluster(chromosome, coordinate, query_pA_sites_plus_strand, subject_pA_sites_plus_strand, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_FOLD_UPREGULATION_OF_CLUSTER, MIN_REFERENCE_VALUE_TO_DIVIDE_BY, max_coordinate)               
                    start, end = current_cluster_range
                    center_of_mass = int( calculate_center_of_mass( chromosome, start, end, upregulated_pA_sites_plus_strand) )
                    most_reads, max_peak, width_around_peak_accounting_for_signal, distinct_peaks_accounting_for_signal = determine_range_around_peak_accounting_for_signal(query_pA_sites_plus_strand, chromosome, start, end, current_cluster_reads, FRACTION_OF_TOTAL_CLUSTER_SIGNAL)             
                    #if max_peak in    upregulated_pA_sites_plus_strand[chromosome]:      ##ensure that the max peak of the cluster is upregulated at least X-fold
                    total_clusters += 1
                    if max_peak in plus_strand_clusters[chromosome]:
                        print('peak conflict')
                    plus_strand_clusters[chromosome][max_peak] =   float ( sum(current_cluster_reads) )  
                    annotation_type, gene_name, gene_id, strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, max_peak, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types) 
                    if max_peak ==    822146  :
                        print (annotation_type, gene_name, gene_id, strandedness )
                    ## insert the max_peak_coordinate into the clusters dictionary
                    fold_upreglation = check_fold_upregulation(chromosome, start, end, query_pA_sites_plus_strand, subject_pA_sites_plus_strand, 1)  ## use 1 read for MIN_REFERENCE_VALUE_TO_DIVIDE_BY
                    output = '\t'.join( [chromosome, str(max_peak), str(center_of_mass), strand, annotation_type, gene_name, gene_id, strandedness, str(start), str(end), str(end-start), str(sum(current_cluster_reads)), str(most_reads),  str(width_around_peak_accounting_for_signal), str(distinct_peaks_accounting_for_signal), str(fold_upreglation)  ] ) + '\n'
                    ## coordinates followed by the total range (in nt) of the cluster, the total reads in cluster, the number of reads for the peak coordinate in the cluster, the peak coordinate in the cluster,
                    ## the width around peak accounting for the signal, number of distinct peaks within the width accounting for the signal
                    outfile.write(output)
        strand = '-'
        for chromosome in upregulated_pA_sites_minus_strand:
            print(chromosome, strand)
            if chromosome not in minus_strand_clusters:
                minus_strand_clusters[chromosome] = {}
            coordinates_examined = []
            max_coordinate = max( upregulated_pA_sites_minus_strand[chromosome].keys() )
            for coordinate in sorted(upregulated_pA_sites_minus_strand[chromosome].keys() ):
                ##print 'coordinate', coordinate, 'coordinates_examined', coordinates_examined
                if coordinate not in coordinates_examined:
                    
                    current_cluster_range, current_cluster_reads, coordinates_examined = extend_upregulated_cluster(chromosome, coordinate, query_pA_sites_minus_strand, subject_pA_sites_minus_strand, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_FOLD_UPREGULATION_OF_CLUSTER, MIN_REFERENCE_VALUE_TO_DIVIDE_BY, max_coordinate)               
                    start, end = current_cluster_range
                    center_of_mass = int( calculate_center_of_mass( chromosome, start, end, upregulated_pA_sites_minus_strand) )
                    most_reads, max_peak, width_around_peak_accounting_for_signal, distinct_peaks_accounting_for_signal = determine_range_around_peak_accounting_for_signal(query_pA_sites_minus_strand, chromosome, start, end, current_cluster_reads, FRACTION_OF_TOTAL_CLUSTER_SIGNAL)             
                    #if max_peak in upregulated_pA_sites_minus_strand[chromosome]:     ##ensure that the max peak of the cluster is upregulated at least X-fold
                    total_clusters += 1
                    minus_strand_clusters[chromosome][max_peak] =   float ( sum(current_cluster_reads) )  
                    annotation_type, gene_name, gene_id, strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, max_peak, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types) 
                    ## insert the max_peak_coordinate into the clusters dictionary  
                    fold_upreglation = check_fold_upregulation(chromosome, start, end, query_pA_sites_minus_strand, subject_pA_sites_minus_strand, 1)  ## use 1 read for MIN_REFERENCE_VALUE_TO_DIVIDE_BY
                    output = '\t'.join( [chromosome, str(max_peak), str(center_of_mass), strand, annotation_type, gene_name, gene_id, strandedness, str(start), str(end), str(end-start), str(sum(current_cluster_reads)), str(most_reads),  str(width_around_peak_accounting_for_signal), str(distinct_peaks_accounting_for_signal), str(fold_upreglation)  ] ) + '\n'
                    ## coordinates followed by the total range (in nt) of the cluster, the total reads in cluster, the number of reads for the peak coordinate in the cluster, the peak coordinate in the cluster,
                    ## the width around peak accounting for the signal, number of distinct peaks within the width accounting for the signal
                    outfile.write(output)
        print(total_clusters, ' total clusters for ', outfilename)
    outfile.close()
    return plus_strand_clusters, minus_strand_clusters
    
def cluster_upregulated_pA_sites_for_two_queries_and_two_references(upregulated_plus_strand_bedgraph, query_plus_strand_bedgraph, reference_plus_strand_bedgraph, query2_plus_strand_bedgraph, reference2_plus_strand_bedgraph, upregulated_minus_strand_bedgraph, query_minus_strand_bedgraph, reference_minus_strand_bedgraph, query2_minus_strand_bedgraph, reference2_minus_strand_bedgraph, annotations_type_dictionary, annotations_dictionary, PRIORITIZED_SEQUENCE_TYPES, outfilename, FLANKING_NUCLEOTIDES, MIN_FOLD_UPREGULATION_OF_CLUSTER, MIN_REFERENCE_VALUE_TO_DIVIDE_BY, FRACTION_OF_TOTAL_CLUSTER_SIGNAL, MIN_FOLD_UPREGULATION_OF_CLUSTER_QUERY2):
    '''
    takes three bedgraphs, a query upregulated bedgraph, the query sites to cluster around the upregulated sites, and the subject pA sites to test for upregulation against
    writes each cluster to a file as the chromosome, start and end coordinates,  peak_coord, width_around_peak_accounting_for_signal, number_of_distinct_peaks, sequence 100 nt upstream, the sequence of the cluster, sequence 100 nt downstream
    returns the total number of clusters identified
    '''
    with open(outfilename, 'w') as outfile:
        header = '\t'.join( ['chromosome', 'peak coord', 'center of mass coordinate for cluster', 'strand', 'annotation type', 'gene name', 'gene ID', 'strandedness', 'cluster start', 'cluster end', 'cluster width', 'total cluster reads', 'peak height (reads)',  'width around peak accounting for ' + str(FRACTION_OF_TOTAL_CLUSTER_SIGNAL*100) + '% of signal', ' distinct number of three prime ends accounting for ' + str(FRACTION_OF_TOTAL_CLUSTER_SIGNAL*100) + '% of signal', 'fold upregulation', 'fold upregulation in RRP6-FRB over RRP6/SEN1-FRB of the window spanning the peak with 20 bp upstream and downstream' ] ) + '\n'
        outfile.write(header)        
        print('cluster_upregulated_pA_sites started')
        plus_strand_clusters = {}
        minus_strand_clusters = {}
        total_clusters = 0
        strand = '+'
        for chromosome in upregulated_plus_strand_bedgraph:
            print(chromosome, strand)
            if chromosome not in plus_strand_clusters:
                plus_strand_clusters[chromosome] = {}
    
            coordinates_examined = []
            max_coordinate = max( upregulated_plus_strand_bedgraph[chromosome].keys() )
            for coordinate in sorted(upregulated_plus_strand_bedgraph[chromosome].keys() ):
                if coordinate not in coordinates_examined:
                    
                    current_cluster_range, current_cluster_reads, coordinates_examined = extend_upregulated_cluster(chromosome, coordinate, query_plus_strand_bedgraph, reference_plus_strand_bedgraph, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_FOLD_UPREGULATION_OF_CLUSTER, MIN_REFERENCE_VALUE_TO_DIVIDE_BY, max_coordinate)               
                    start, end = current_cluster_range
                    center_of_mass = int( calculate_center_of_mass( chromosome, start, end, upregulated_plus_strand_bedgraph) )
                    most_reads, max_peak, width_around_peak_accounting_for_signal, distinct_peaks_accounting_for_signal = determine_range_around_peak_accounting_for_signal(query_plus_strand_bedgraph, chromosome, start, end, current_cluster_reads, FRACTION_OF_TOTAL_CLUSTER_SIGNAL)             
                   # if max_peak in    upregulated_plus_strand_bedgraph[chromosome]:                     
                    
                    annotation_type, gene_name, gene_id, strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, max_peak, annotations_type_dictionary, annotations_dictionary, PRIORITIZED_SEQUENCE_TYPES) 
                    if max_peak ==    822146  :
                        print (annotation_type, gene_name, gene_id, strandedness )
                    ## insert the max_peak_coordinate into the clusters dictionary
                    fold_upreglation = check_fold_upregulation(chromosome, start, end, query_plus_strand_bedgraph, reference_plus_strand_bedgraph, 1)  ## use 1 read for MIN_REFERENCE_VALUE_TO_DIVIDE_BY
                    ## calculate fold upregulation for query2 and reference 2 for the region around the peak
                    
                    query2_fold_upreglation = check_fold_upregulation(chromosome,  max_peak - 20, max_peak + 20, query2_plus_strand_bedgraph, reference2_plus_strand_bedgraph, 0.1)
                  #  if query2_fold_upreglation >= MIN_FOLD_UPREGULATION_OF_CLUSTER_QUERY2:
                    plus_strand_clusters[chromosome][max_peak] =   float ( sum(current_cluster_reads) )  
                    total_clusters += 1
                    output = '\t'.join( [chromosome, str(max_peak), str(center_of_mass), strand, annotation_type, gene_name, gene_id, strandedness, str(start), str(end), str(end-start), str(sum(current_cluster_reads)), str(most_reads),  str(width_around_peak_accounting_for_signal), str(distinct_peaks_accounting_for_signal), str(fold_upreglation), str(query2_fold_upreglation)  ] ) + '\n'
                    ## coordinates followed by the total range (in nt) of the cluster, the total reads in cluster, the number of reads for the peak coordinate in the cluster, the peak coordinate in the cluster,
                    ## the width around peak accounting for the signal, number of distinct peaks within the width accounting for the signal
                    outfile.write(output)
        strand = '-'
        for chromosome in upregulated_minus_strand_bedgraph:
            print(chromosome, strand)
            if chromosome not in minus_strand_clusters:
                minus_strand_clusters[chromosome] = {}
            coordinates_examined = []
            max_coordinate = max( upregulated_minus_strand_bedgraph[chromosome].keys() )
            for coordinate in sorted(upregulated_minus_strand_bedgraph[chromosome].keys() ):
                ##print 'coordinate', coordinate, 'coordinates_examined', coordinates_examined
                if coordinate not in coordinates_examined:
                    current_cluster_range, current_cluster_reads, coordinates_examined = extend_upregulated_cluster(chromosome, coordinate, query_minus_strand_bedgraph, reference_minus_strand_bedgraph, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_FOLD_UPREGULATION_OF_CLUSTER, MIN_REFERENCE_VALUE_TO_DIVIDE_BY, max_coordinate)               
                    start, end = current_cluster_range
                    center_of_mass = int( calculate_center_of_mass( chromosome, start, end, upregulated_minus_strand_bedgraph) )
                    most_reads, max_peak, width_around_peak_accounting_for_signal, distinct_peaks_accounting_for_signal = determine_range_around_peak_accounting_for_signal(query_minus_strand_bedgraph, chromosome, start, end, current_cluster_reads, FRACTION_OF_TOTAL_CLUSTER_SIGNAL)             
                    # if max_peak in    upregulated_minus_strand_bedgraph[chromosome]:                      
                     
                    annotation_type, gene_name, gene_id, strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, max_peak, annotations_type_dictionary, annotations_dictionary, PRIORITIZED_SEQUENCE_TYPES) 
                    ## insert the max_peak_coordinate into the clusters dictionary  
                    fold_upreglation = check_fold_upregulation(chromosome, start, end, query_minus_strand_bedgraph, reference_minus_strand_bedgraph, 1)  ## use 1 read for MIN_REFERENCE_VALUE_TO_DIVIDE_BY
                    ## calculate fold upregulation for query2 and reference 2
                    query2_fold_upreglation = check_fold_upregulation(chromosome,  max_peak - 20, max_peak + 20,  query2_minus_strand_bedgraph, reference2_minus_strand_bedgraph, 0.1)
                  #  if query2_fold_upreglation >= MIN_FOLD_UPREGULATION_OF_CLUSTER_QUERY2:
                    minus_strand_clusters[chromosome][max_peak] =   float ( sum(current_cluster_reads) ) 
                    total_clusters += 1                        
                    output = '\t'.join( [chromosome, str(max_peak), str(center_of_mass), strand, annotation_type, gene_name, gene_id, strandedness, str(start), str(end), str(end-start), str(sum(current_cluster_reads)), str(most_reads),  str(width_around_peak_accounting_for_signal), str(distinct_peaks_accounting_for_signal), str(fold_upreglation), str(query2_fold_upreglation)  ] ) + '\n'
                    ## coordinates followed by the total range (in nt) of the cluster, the total reads in cluster, the number of reads for the peak coordinate in the cluster, the peak coordinate in the cluster,
                    ## the width around peak accounting for the signal, number of distinct peaks within the width accounting for the signal
                    outfile.write(output)
        print(total_clusters, ' total clusters for ', outfilename)
    outfile.close()
    return plus_strand_clusters, minus_strand_clusters    
    
def cluster_pA_sites_for_Rnt1p_targets(query_pA_sites_plus_strand, query_pA_sites_minus_strand, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types, outfilename, FLANKING_NUCLEOTIDES,  PERCENTAGE_OF_TOTAL_CLUSTER_SIGNAL, MIN_READS_TO_EXTEND_CLUSTER, predicted_plus_strand_tetraloops, predicted_minus_strand_tetraloops, WINDOW_TO_SEARCH_FOR_PREDICTED_RNT1_TETRALOOPS):
    '''
    takes two bedgraphs for the plus and minus strands, and clusters the sites, search for predicted tetraloops within 
    writes each cluster to a file as the chromosome, start and end coordinates,  peak_coord, width_around_peak_accounting_for_signal, number_of_distinct_peaks, sequence 100 nt upstream, the sequence of the cluster, sequence 100 nt downstream
    returns the total number of clusters identified
    '''
    with open(outfilename, 'w') as outfile:
        output = '\t'.join( 'chromosome, max_peak_start, max_peak_end, strand, annotation_type, gene_name, gene_id, strandedness, start, end, cluster_width, total_cluster_reads, max_peak_height,  width_around_peak_accounting_for_80%_of_total_signal, distinct_peaks_accounting_for_80%_of_total_signal, distance_to_nearest_NGNN_AAGU_tetraloop(10000=None)'.split(', ') ) + '\n'
        outfile.write(output)
        print('cluster_pA_sites started')
        plus_strand_clusters = {}
        minus_strand_clusters = {}
        total_clusters = 0
        strand = '+'
        for chromosome in query_pA_sites_plus_strand:
            print(chromosome, strand)
            if chromosome not in plus_strand_clusters:
                plus_strand_clusters[chromosome] = {}
            coordinates_examined = []
            if query_pA_sites_plus_strand[chromosome] != {}:
                max_coordinate = max( query_pA_sites_plus_strand[chromosome].keys() )
                for coordinate in sorted(query_pA_sites_plus_strand[chromosome].keys() ):
                    if coordinate not in coordinates_examined:
                        total_clusters += 1
                        current_cluster_range, current_cluster_reads, coordinates_examined = extend_cluster(chromosome, coordinate, query_pA_sites_plus_strand, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_READS_TO_EXTEND_CLUSTER, max_coordinate)      
                        start, end = current_cluster_range
                        most_reads, max_peak, width_around_peak_accounting_for_signal, distinct_peaks_accounting_for_signal = determine_range_around_peak_accounting_for_signal(query_pA_sites_plus_strand, chromosome, start, end, current_cluster_reads, PERCENTAGE_OF_TOTAL_CLUSTER_SIGNAL)             
                        plus_strand_clusters[chromosome][max_peak] =   float ( sum(current_cluster_reads) )
                        distance_from_tetraloop = 10000
                        for idx in range(max_peak - WINDOW_TO_SEARCH_FOR_PREDICTED_RNT1_TETRALOOPS, max_peak):
                            if idx in predicted_plus_strand_tetraloops[chromosome]:
                                distance_from_tetraloop = max_peak - idx 
                        for idx in range(max_peak, max_peak + 4 + WINDOW_TO_SEARCH_FOR_PREDICTED_RNT1_TETRALOOPS):
                            if idx in predicted_plus_strand_tetraloops[chromosome]:
                                distance_from_tetraloop = idx - (max_peak + 4)                     
                        annotation_type, gene_name, gene_id, strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, coordinate, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types) 
                        ## insert the max_peak_coordinate into the clusters dictionary            
                        output = '\t'.join( [chromosome, str(max_peak), str(max_peak + 1), strand, annotation_type, gene_name, gene_id, strandedness, str(start), str(end), str(end-start), str(sum(current_cluster_reads)), str(most_reads),  str(width_around_peak_accounting_for_signal), str(distinct_peaks_accounting_for_signal),  str(distance_from_tetraloop)  ] ) + '\n'
                        ## coordinates followed by the total range (in nt) of the cluster, the total reads in cluster, the number of reads for the peak coordinate in the cluster, the peak coordinate in the cluster,
                        ## the width around peak accounting for the signal, number of distinct peaks within the width accounting for the signal
                        outfile.write(output)
            else:
                print (chromosome, strand, 'empty')
        strand = '-'
        for chromosome in query_pA_sites_minus_strand:
            print(chromosome, strand)
            if chromosome not in minus_strand_clusters:
                minus_strand_clusters[chromosome] = {}
            coordinates_examined = []
            if query_pA_sites_minus_strand[chromosome] != {}:
                max_coordinate = max( query_pA_sites_minus_strand[chromosome].keys() )
                for coordinate in sorted(query_pA_sites_minus_strand[chromosome].keys() ):
                    ##print 'coordinate', coordinate, 'coordinates_examined', coordinates_examined
                    if coordinate not in coordinates_examined:
                        total_clusters += 1
                        current_cluster_range, current_cluster_reads, coordinates_examined = extend_cluster(chromosome, coordinate, query_pA_sites_minus_strand, coordinates_examined, FLANKING_NUCLEOTIDES, MIN_READS_TO_EXTEND_CLUSTER, max_coordinate)      
                        start, end = current_cluster_range
                        most_reads, max_peak, width_around_peak_accounting_for_signal, distinct_peaks_accounting_for_signal = determine_range_around_peak_accounting_for_signal(query_pA_sites_minus_strand, chromosome, start, end, current_cluster_reads, PERCENTAGE_OF_TOTAL_CLUSTER_SIGNAL)             
                        minus_strand_clusters[chromosome][max_peak] =   float ( sum(current_cluster_reads) )  
                        distance_from_tetraloop = 10000
                        for idx in range(max_peak - WINDOW_TO_SEARCH_FOR_PREDICTED_RNT1_TETRALOOPS, max_peak):
                            if idx in predicted_minus_strand_tetraloops[chromosome]:
                                distance_from_tetraloop = max_peak - idx 
                        for idx in range(max_peak, max_peak + 4 + WINDOW_TO_SEARCH_FOR_PREDICTED_RNT1_TETRALOOPS):
                            if idx in predicted_minus_strand_tetraloops[chromosome]:
                                distance_from_tetraloop = idx - (max_peak + 4)   
                        annotation_type, gene_name, gene_id, strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, coordinate, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types) 
                        ## insert the max_peak_coordinate into the clusters dictionary            
                        output = '\t'.join( [chromosome, str(max_peak), str(max_peak + 1), strand, annotation_type, gene_name, gene_id, strandedness, str(start), str(end), str(end-start), str(sum(current_cluster_reads)), str(most_reads),  str(width_around_peak_accounting_for_signal), str(distinct_peaks_accounting_for_signal),  str(distance_from_tetraloop)  ] ) + '\n'
                        ## coordinates followed by the total range (in nt) of the cluster, the total reads in cluster, the number of reads for the peak coordinate in the cluster, the peak coordinate in the cluster,
                        ## the width around peak accounting for the signal, number of distinct peaks within the width accounting for the signal
                        outfile.write(output)
            else:
                print (chromosome, strand, 'empty')
        print(total_clusters, ' total clusters for ', outfilename)
    outfile.close()
    return plus_strand_clusters, minus_strand_clusters

def bin_annotation_body_reads(BODY_SIZE, annotation_body_reads):
    '''
    takes a list of annotation gene body reads of variable length
    returns a binned list of reads of BODY_SIZE length
    '''
    if BODY_SIZE == 0:
        return []
    binned_reads = [0]*BODY_SIZE
    relative_size_difference = float(BODY_SIZE) / len(annotation_body_reads)
    for idx in range(len(annotation_body_reads)):
        binned_reads[min(int(idx*relative_size_difference), BODY_SIZE-1)] = annotation_body_reads[idx]
    return binned_reads
    
def aggregate_reads_for_genomic_feature(annotations, query_plus_strand, query_minus_strand, UPSTREAM_REGION_TO_INCLUDE, BODY_SIZE, DOWNSTREAM_REGION_TO_INCLUDE, MIN_READS_TO_INCLUDE_REGION = 1,  USE_ALL_ANNOTATIONS = False, SCALING_FACTOR = 1, annotations_coordinates_to_exclude = None, normalize_total_reads = True, normalize_both_strands = True, count_both_strands_as_separate_annotations = False):
    '''
    takes an annotations dictionary, query bedgraphs for the plus strand and 
    minus strand, the upstream and downstream region in bp in which to normalize 
    and aggregate reads, a body size for which to bin reads when annotations are of 
    different lengths, the minimum number of sense reads (in rpm) to include the region
    in the analysis (overriden if USE_ALL_ANNOTATIONS is True), a scaling factor for the 
    normalized expression values, and annotations_coordinates_to_exclude can be used
    to exclude 3' ends of abundant non-coding transcripts, which are contaminants in NET-seq data
    and which represent TRAMP-adenylation of mature products in nuclear exosome inactivated cells
    returns the annotations used in the anlysis, and the sense and antisense 
    normalized and aggregated reads as lists of floats
    at the end, sets the max of the absolute value of the sense and antisense to 100, and scales the data-set accordingly
    '''
    upstream_reads = [0]*UPSTREAM_REGION_TO_INCLUDE
    downstream_reads = [0]*DOWNSTREAM_REGION_TO_INCLUDE
    annotation_body_reads = [0]*BODY_SIZE
    aggregated_sense_reads = upstream_reads + annotation_body_reads + downstream_reads
    aggregated_antisense_reads = aggregated_sense_reads[:]
    total_genomic_features_included_in_this_anlysis = 0
    annotations_used = {}
    for strand in annotations:
        if strand not in annotations_used:
            annotations_used[strand] = {}
        for chromosome in annotations[strand]:
            if chromosome not in annotations_used[strand]:
                annotations_used[strand][chromosome] = {}
            for annotation in annotations[strand][chromosome]:
                ##print strand, chromosome, annotations[strand][chromosome][annotation], annotation
                if annotations_coordinates_to_exclude != None:                
                    excluded_coordinates = annotations_coordinates_to_exclude[strand][chromosome][annotation]
                else:
                    excluded_coordinates = []
                if strand == '+':
                    start = annotations[strand][chromosome][annotation][0]
                    end = annotations[strand][chromosome][annotation][1]

                    sense_upstream_reads = return_expression_values(chromosome, start - UPSTREAM_REGION_TO_INCLUDE, start, query_plus_strand, excluded_coordinates)
                    sense_downstream_reads = return_expression_values(chromosome, end, end + DOWNSTREAM_REGION_TO_INCLUDE, query_plus_strand, excluded_coordinates)
                    sense_annotation_reads = return_expression_values(chromosome, start, end, query_plus_strand, excluded_coordinates)
                    sense_binned_annotation_body_reads = bin_annotation_body_reads(BODY_SIZE, sense_annotation_reads)

                    antisense_upstream_reads = return_expression_values(chromosome, start - UPSTREAM_REGION_TO_INCLUDE, start, query_minus_strand)
                    antisense_downstream_reads = return_expression_values(chromosome, end, end + DOWNSTREAM_REGION_TO_INCLUDE, query_minus_strand)
                    antisense_annotation_reads = return_expression_values(chromosome, start, end, query_minus_strand)
                    antisense_binned_annotation_body_reads = bin_annotation_body_reads(BODY_SIZE, antisense_annotation_reads)

                    entire_sense_region_reads = sense_upstream_reads + sense_binned_annotation_body_reads + sense_downstream_reads
                    entire_antisense_region_reads = antisense_upstream_reads + antisense_binned_annotation_body_reads + antisense_downstream_reads
                elif strand == '-':
                    start = annotations[strand][chromosome][annotation][0]
                    end = annotations[strand][chromosome][annotation][1]

                    sense_upstream_reads = return_expression_values(chromosome, end, end + UPSTREAM_REGION_TO_INCLUDE, query_minus_strand, excluded_coordinates)[::-1]
                    sense_downstream_reads = return_expression_values(chromosome, start - DOWNSTREAM_REGION_TO_INCLUDE, start, query_minus_strand, excluded_coordinates)[::-1]
                    sense_annotation_reads = return_expression_values(chromosome, start, end, query_minus_strand, excluded_coordinates)[::-1]
                    sense_binned_annotation_body_reads = bin_annotation_body_reads(BODY_SIZE, sense_annotation_reads)

                    antisense_upstream_reads = return_expression_values(chromosome, end, end + UPSTREAM_REGION_TO_INCLUDE,  query_plus_strand)[::-1]
                    antisense_downstream_reads = return_expression_values(chromosome, start - DOWNSTREAM_REGION_TO_INCLUDE, start, query_plus_strand)[::-1]
                    antisense_annotation_reads = return_expression_values(chromosome, start, end, query_plus_strand)[::-1]
                    antisense_binned_annotation_body_reads = bin_annotation_body_reads(BODY_SIZE, antisense_annotation_reads)

                    entire_sense_region_reads = sense_upstream_reads + sense_binned_annotation_body_reads + sense_downstream_reads
                    entire_antisense_region_reads = antisense_upstream_reads + antisense_binned_annotation_body_reads + antisense_downstream_reads
                include_annotation = False
                if not normalize_both_strands:
                    if sum(entire_sense_region_reads) >= MIN_READS_TO_INCLUDE_REGION:
                        normalization_factor = SCALING_FACTOR /  max(entire_sense_region_reads) 
                        include_annotation = True
                        for idx in range(UPSTREAM_REGION_TO_INCLUDE + BODY_SIZE + DOWNSTREAM_REGION_TO_INCLUDE):
                            aggregated_sense_reads[idx] += entire_sense_region_reads[idx] *  normalization_factor
                    if sum(entire_antisense_region_reads) >= MIN_READS_TO_INCLUDE_REGION:
                        normalization_factor = SCALING_FACTOR  /  max(entire_antisense_region_reads) 
                        include_annotation = True
                        for idx in range(UPSTREAM_REGION_TO_INCLUDE + BODY_SIZE + DOWNSTREAM_REGION_TO_INCLUDE):
                            aggregated_antisense_reads[idx] += entire_antisense_region_reads[idx] * normalization_factor
    
                    if include_annotation:
                        annotations_used[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                        total_genomic_features_included_in_this_anlysis += 1
#                        if 350> entire_sense_region_reads.index( max(entire_sense_region_reads)) > 300:
#                            print(annotation)
                if normalize_both_strands:
                    if sum(entire_sense_region_reads + entire_antisense_region_reads) >= MIN_READS_TO_INCLUDE_REGION:
                        normalization_factor = SCALING_FACTOR / max(entire_sense_region_reads + entire_antisense_region_reads) 
                        include_annotation = True
                        for idx in range(UPSTREAM_REGION_TO_INCLUDE + BODY_SIZE + DOWNSTREAM_REGION_TO_INCLUDE):
                            aggregated_sense_reads[idx] += entire_sense_region_reads[idx] *  normalization_factor
                        for idx in range(UPSTREAM_REGION_TO_INCLUDE + BODY_SIZE + DOWNSTREAM_REGION_TO_INCLUDE):
                            aggregated_antisense_reads[idx] += entire_antisense_region_reads[idx] * normalization_factor
    
                    if include_annotation:
                        annotations_used[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                        total_genomic_features_included_in_this_anlysis += 1
#                        if 350> entire_sense_region_reads.index( max(entire_sense_region_reads)) > 300:
#                            print(annotation)
                #else:
                #    print annotation, 'has insufficient read coverage to be included in this analysis'
    if not normalize_total_reads:
        max_signal =  max(aggregated_sense_reads)
        normalization_value = SCALING_FACTOR / max_signal
        for idx in range(len(aggregated_sense_reads)):
            aggregated_sense_reads[idx] = aggregated_sense_reads[idx] * normalization_value
            
        max_signal =  max(aggregated_antisense_reads)
        normalization_value = SCALING_FACTOR / max_signal
        for idx in range(len(aggregated_antisense_reads)):       
            aggregated_antisense_reads[idx] = aggregated_antisense_reads[idx] * normalization_value
    if normalize_total_reads:
        total_signal =  max(aggregated_sense_reads)
        normalization_value = SCALING_FACTOR / total_signal
        for idx in range(len(aggregated_sense_reads)):
            aggregated_sense_reads[idx] = aggregated_sense_reads[idx] * normalization_value
            
        total_signal =  max(aggregated_antisense_reads)
        normalization_value = SCALING_FACTOR / total_signal
        for idx in range(len(aggregated_antisense_reads)):       
            aggregated_antisense_reads[idx] = aggregated_antisense_reads[idx] * normalization_value
#    if 0 in aggregated_sense_reads:
#        print(aggregated_sense_reads.index(0))
#    if 0 in aggregated_antisense_reads:
#        print(aggregated_antisense_reads.index(0))
    print(total_genomic_features_included_in_this_anlysis, 'total_genomic_features_included_in_this_anlysis')
    return annotations_used, total_genomic_features_included_in_this_anlysis, aggregated_sense_reads, aggregated_antisense_reads

def aggregate_nonstranded_reads_for_genomic_features(annotations, query_bedgraph, UPSTREAM_REGION_TO_INCLUDE, BODY_SIZE, DOWNSTREAM_REGION_TO_INCLUDE, MIN_READS_TO_INCLUDE_REGION, NORMALIZE = False, USE_ALL_ANNOTATIONS = False, SCALING_FACTOR = 1, max_stretch_of_zero_values_allowed = None, flanking_region_to_check_signal = 40):
    upstream_reads = [0]*UPSTREAM_REGION_TO_INCLUDE
    downstream_reads = [0]*DOWNSTREAM_REGION_TO_INCLUDE
    annotation_body_reads = [0]*BODY_SIZE
    aggregated_reads = upstream_reads + annotation_body_reads + downstream_reads
    total_genomic_features_included_in_this_anlysis = 0
    annotations_used = {}
    for strand in annotations:
        if strand not in annotations_used:
            annotations_used[strand] = {}
        for chromosome in annotations[strand]:
            if chromosome not in annotations_used[strand]:
                annotations_used[strand][chromosome] = {}
            for annotation in annotations[strand][chromosome]:
                start = annotations[strand][chromosome][annotation][0]
                end = annotations[strand][chromosome][annotation][1]
                if strand == '+':
                    upstream_reads = return_expression_values(chromosome, start - UPSTREAM_REGION_TO_INCLUDE, start, query_bedgraph)
                    downstream_reads = return_expression_values(chromosome, end, end + DOWNSTREAM_REGION_TO_INCLUDE, query_bedgraph)
                    annotation_reads = return_expression_values(chromosome, start, end, query_bedgraph)  
                elif strand == '-':
                    upstream_reads = return_expression_values(chromosome, end, end + UPSTREAM_REGION_TO_INCLUDE, query_bedgraph)[::-1]
                    downstream_reads = return_expression_values(chromosome, start - DOWNSTREAM_REGION_TO_INCLUDE, start, query_bedgraph)[::-1]
                    annotation_reads = return_expression_values(chromosome, start, end, query_bedgraph)[::-1]
                annotation_reads = bin_annotation_body_reads(BODY_SIZE, annotation_reads)
                entire_region_reads = upstream_reads + annotation_reads + downstream_reads
                reads_flanking_motif = upstream_reads[-flanking_region_to_check_signal:] + annotation_reads + downstream_reads[:flanking_region_to_check_signal]
                if max_stretch_of_zero_values_allowed == None or entire_region_reads.count(0) < max_stretch_of_zero_values_allowed:
                    # if sum(reads_flanking_motif) >= MIN_READS_TO_INCLUDE_REGION:
                    annotations_used[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                    total_genomic_features_included_in_this_anlysis += 1
                    normalization_factor = 1 * SCALING_FACTOR
                    if NORMALIZE:
                        if abs(sum(entire_region_reads)) > 0:
                            normalization_factor =  SCALING_FACTOR / abs(max(entire_region_reads))
                    for idx in range(len(entire_region_reads)):
                        if idx < len(aggregated_reads):
                            aggregated_reads[idx] += entire_region_reads[idx] * normalization_factor
#               else:
#                    print annotation, 'has insufficient read coverage to be included in this analysis'
    max_signal = max([abs(e) for e in aggregated_reads] )  ## check highest magnitude + or - value
    normalization_value = SCALING_FACTOR  / max_signal
    for idx in range(len(aggregated_reads)):
        aggregated_reads[idx] = aggregated_reads[idx] * normalization_value
    print(total_genomic_features_included_in_this_anlysis, 'total_genomic_features_included_in_this_anlysis')
    return annotations_used, total_genomic_features_included_in_this_anlysis, aggregated_reads

def aggregate_reads_for_genomic_features_on_multiple_data_sets(annotations, datasets, UPSTREAM_REGION_TO_INCLUDE, BODY_SIZE, DOWNSTREAM_REGION_TO_INCLUDE, MIN_TOTAL_READS_TO_INCLUDE_REGION, 
                                                               USE_ALL_ANNOTATIONS = False, pseudo_counts_to_add_in_place_of_zero = 0, SCALING_FACTOR = 100, annotations_coordinates_to_exclude = None, 
                                                               normalize_total_reads = True, normalize_both_strands = False, count_both_strands_as_separate_annotations = False, flanking_region_to_check_signal = 40):
    upstream_reads = [0]*UPSTREAM_REGION_TO_INCLUDE
    downstream_reads = [0]*DOWNSTREAM_REGION_TO_INCLUDE
    annotation_body_reads = [0]*BODY_SIZE
    aggregated_sense_reads = upstream_reads + annotation_body_reads + downstream_reads
    aggregated_antisense_reads = aggregated_sense_reads[:]
    for dataset in datasets:
        datasets[dataset][2] = aggregated_sense_reads[:]
        datasets[dataset][3] = aggregated_antisense_reads[:]
    total_genomic_features_included_in_this_anlysis = 0
    annotations_used = generate_empty_stranded_bedgraph_dict()
    max_pos_value_used = 0
    max_neg_value_used = 0
    total_features = 0
    for strand in annotations:
        for chromosome in annotations[strand]:
            for annotation in annotations[strand][chromosome]:
                total_features += 1
                start = annotations[strand][chromosome][annotation][0]
                end = annotations[strand][chromosome][annotation][1]
#                max_reads_in_flanking_window = 0
                if strand == '+':
                    for dataset in datasets:
                        query_plus_strand = datasets[dataset][0]
                        query_minus_strand = datasets[dataset][1]
                        sense_upstream_reads = return_expression_values(chromosome, start - UPSTREAM_REGION_TO_INCLUDE, start, query_plus_strand, [], pseudo_counts_to_add_in_place_of_zero)
                        sense_downstream_reads = return_expression_values(chromosome, end, end + DOWNSTREAM_REGION_TO_INCLUDE, query_plus_strand, [], pseudo_counts_to_add_in_place_of_zero)
                        sense_annotation_reads = return_expression_values(chromosome, start, end, query_plus_strand, [], pseudo_counts_to_add_in_place_of_zero)
                        sense_binned_annotation_body_reads = bin_annotation_body_reads(BODY_SIZE, sense_annotation_reads)
                        antisense_upstream_reads = return_expression_values(chromosome, start - UPSTREAM_REGION_TO_INCLUDE, start, query_minus_strand, [], pseudo_counts_to_add_in_place_of_zero)
                        antisense_downstream_reads = return_expression_values(chromosome, end, end + DOWNSTREAM_REGION_TO_INCLUDE, query_minus_strand, [], pseudo_counts_to_add_in_place_of_zero)
                        antisense_annotation_reads = return_expression_values(chromosome, start, end, query_minus_strand, [], pseudo_counts_to_add_in_place_of_zero)
                        antisense_binned_annotation_body_reads = bin_annotation_body_reads(BODY_SIZE, antisense_annotation_reads)
                        query_entire_sense_region_reads = sense_upstream_reads + sense_binned_annotation_body_reads + sense_downstream_reads
#                        reads_flanking_motif = upstream_reads[-flanking_region_to_check_signal:] + sense_binned_annotation_body_reads + downstream_reads[:flanking_region_to_check_signal]
#                        total_reads_flanking_motif = sum(reads_flanking_motif)
#                        if total_reads_flanking_motif > max_reads_in_flanking_window:
#                            max_reads_in_flanking_window = total_reads_flanking_motif
                        query_entire_antisense_region_reads = antisense_upstream_reads + antisense_binned_annotation_body_reads + antisense_downstream_reads
                        datasets[dataset][4] = query_entire_sense_region_reads
                        datasets[dataset][5] = query_entire_antisense_region_reads
                elif strand == '-':
                    for dataset in datasets:
                        query_plus_strand = datasets[dataset][0]
                        query_minus_strand = datasets[dataset][1]
                        sense_upstream_reads = return_expression_values(chromosome, end, end + UPSTREAM_REGION_TO_INCLUDE, query_minus_strand, [], pseudo_counts_to_add_in_place_of_zero)[::-1]
                        sense_downstream_reads = return_expression_values(chromosome, start - DOWNSTREAM_REGION_TO_INCLUDE, start, query_minus_strand, [], pseudo_counts_to_add_in_place_of_zero)[::-1]
                        sense_annotation_reads = return_expression_values(chromosome, start, end, query_minus_strand, [], pseudo_counts_to_add_in_place_of_zero)[::-1]
                        sense_binned_annotation_body_reads = bin_annotation_body_reads(BODY_SIZE, sense_annotation_reads)
                        antisense_upstream_reads = return_expression_values(chromosome, end, end + UPSTREAM_REGION_TO_INCLUDE,  query_plus_strand, [], pseudo_counts_to_add_in_place_of_zero)[::-1]
                        antisense_downstream_reads = return_expression_values(chromosome, start - DOWNSTREAM_REGION_TO_INCLUDE, start, query_plus_strand, [], pseudo_counts_to_add_in_place_of_zero)[::-1]
                        antisense_annotation_reads = return_expression_values(chromosome, start, end, query_plus_strand, [], pseudo_counts_to_add_in_place_of_zero)[::-1]
                        antisense_binned_annotation_body_reads = bin_annotation_body_reads(BODY_SIZE, antisense_annotation_reads)
                        query_entire_sense_region_reads = sense_upstream_reads + sense_binned_annotation_body_reads + sense_downstream_reads
#                        reads_flanking_motif = upstream_reads[-flanking_region_to_check_signal:] + sense_binned_annotation_body_reads + downstream_reads[:flanking_region_to_check_signal]
#                        total_reads_flanking_motif = sum(reads_flanking_motif)
#                        if total_reads_flanking_motif > max_reads_in_flanking_window:
#                            max_reads_in_flanking_window = total_reads_flanking_motif
                        query_entire_antisense_region_reads = antisense_upstream_reads + antisense_binned_annotation_body_reads + antisense_downstream_reads
                        datasets[dataset][4] = query_entire_sense_region_reads
                        datasets[dataset][5] = query_entire_antisense_region_reads

                total_sense_reads_for_this_genomic_region, total_antisense_reads_for_this_genomic_region = 0, 0
                for dataset in datasets:
                    total_sense_reads_for_this_genomic_region += sum(abs(e) for e in datasets[dataset][4] )
                    total_antisense_reads_for_this_genomic_region += sum(abs(e) for e in datasets[dataset][5] )
                                  
                if total_sense_reads_for_this_genomic_region > MIN_TOTAL_READS_TO_INCLUDE_REGION:
                    if normalize_both_strands:
                        total_reads = total_sense_reads_for_this_genomic_region + total_antisense_reads_for_this_genomic_region
#                        if total_reads >= MIN_TOTAL_READS_TO_INCLUDE_REGION:
                        annotations_used[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                        total_genomic_features_included_in_this_anlysis += 1
                        normalization_factor = 1.00
                        if normalize_total_reads:
                            normalization_factor = SCALING_FACTOR / total_reads  ## this is an important paramater that weights each locus for each dataset equally, as long as one of the datasets has the minimum reads for that locus
                        for dataset in datasets:
                            query_aggregated_sense_reads = datasets[dataset][2]
                            query_aggregated_antisense_reads  = datasets[dataset][3]
                            query_entire_sense_region_reads = datasets[dataset][4]
                            query_entire_antisense_region_reads = datasets[dataset][5] 
                            for idx in range(UPSTREAM_REGION_TO_INCLUDE + BODY_SIZE + DOWNSTREAM_REGION_TO_INCLUDE):
                                query_aggregated_sense_reads[idx] += query_entire_sense_region_reads[idx] * normalization_factor
                                query_aggregated_antisense_reads[idx] += query_entire_antisense_region_reads[idx] * normalization_factor
                            max_pos_value_for_this_dataset = max(abs(e) for e in query_aggregated_sense_reads)
                            max_neg_value_for_this_dataset = max(abs(e) for e in query_aggregated_antisense_reads)
                            if max_pos_value_for_this_dataset > max_pos_value_used:
                                max_pos_value_used = max_pos_value_for_this_dataset
                            if abs(max_neg_value_for_this_dataset) > max_neg_value_used:
                                max_neg_value_used = abs(max_neg_value_for_this_dataset)
                    #else:
                    #    print annotation, 'has insufficient read coverage to be included in this analysis'
                    else:
                        if total_sense_reads_for_this_genomic_region >= MIN_TOTAL_READS_TO_INCLUDE_REGION:
                            annotations_used[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                            total_genomic_features_included_in_this_anlysis += 1
                            normalization_factor = 1.00
                            if normalize_total_reads:
                                normalization_factor = SCALING_FACTOR / total_sense_reads_for_this_genomic_region
                            for dataset in datasets:
                                query_aggregated_sense_reads = datasets[dataset][2]
                                query_entire_sense_region_reads = datasets[dataset][4]
                                for idx in range(UPSTREAM_REGION_TO_INCLUDE + BODY_SIZE + DOWNSTREAM_REGION_TO_INCLUDE):
                                    query_aggregated_sense_reads[idx] += query_entire_sense_region_reads[idx] * normalization_factor
                                max_pos_value_for_this_dataset = max(query_aggregated_sense_reads)
                                if max_pos_value_for_this_dataset > max_pos_value_used:
                                    max_pos_value_used = max_pos_value_for_this_dataset
                                
                        if total_antisense_reads_for_this_genomic_region >= MIN_TOTAL_READS_TO_INCLUDE_REGION:
                            annotations_used[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                            # total_genomic_features_included_in_this_anlysis += 1
                            normalization_factor = 1.00
                            if normalize_total_reads:
                                normalization_factor = SCALING_FACTOR / total_antisense_reads_for_this_genomic_region
                            for dataset in datasets:
                                query_aggregated_antisense_reads  = datasets[dataset][3]
                                query_entire_antisense_region_reads = datasets[dataset][5] 
                                for idx in range(UPSTREAM_REGION_TO_INCLUDE + BODY_SIZE + DOWNSTREAM_REGION_TO_INCLUDE):
                                    query_aggregated_antisense_reads[idx] += query_entire_antisense_region_reads[idx] * normalization_factor
                                max_neg_value_for_this_dataset = max(abs(e) for e in query_aggregated_antisense_reads)
                                if abs(max_neg_value_for_this_dataset) > max_neg_value_used:
                                    max_neg_value_used = abs(max_neg_value_for_this_dataset)
    if normalize_both_strands:                       
        sense_datasets_scaling_factor = SCALING_FACTOR /  max(abs(max_pos_value_used)  , abs(max_neg_value_used)  )
        antisense_datasets_scaling_factor = SCALING_FACTOR /  max(abs(max_pos_value_used)  , abs(max_neg_value_used)  ) # max_pos_value_used #  ## picks the largest magnitude y-axis value for scaling the datasets together, does not change the relative magnitude for datasets
    else:
        sense_datasets_scaling_factor = SCALING_FACTOR / abs(max_pos_value_used)  
        antisense_datasets_scaling_factor = SCALING_FACTOR /  abs(max_neg_value_used)  
    print(max_pos_value_used, sense_datasets_scaling_factor, max_neg_value_used, antisense_datasets_scaling_factor)
    for dataset in datasets:
        temp_query_aggregated_sense_reads = [ sense_datasets_scaling_factor*e for e in datasets[dataset][2] ]
        temp_query_aggregated_antisense_reads  = [ antisense_datasets_scaling_factor*e for e in datasets[dataset][3] ]
        datasets[dataset][2] = temp_query_aggregated_sense_reads
        datasets[dataset][3] = temp_query_aggregated_antisense_reads
        
    print(total_genomic_features_included_in_this_anlysis, 'total_genomic_features_included_in_this_anlysis')
    print(total_features, 'total features considered')
    return annotations_used, total_genomic_features_included_in_this_anlysis, datasets
    

def smooth_line(reads, flanking_bp_smoothing_window):
    '''
    takes a list of reads and a smoothing window
    returns a new list with the read at each position averaged over the flanking_bp_smoothing_window upstream and downstream
    '''
    smoothed_reads = reads[:]
    for idx in range(flanking_bp_smoothing_window, len(reads) - flanking_bp_smoothing_window):
        smoothed_reads[idx] = sum( reads[idx - flanking_bp_smoothing_window: idx + flanking_bp_smoothing_window] ) / float(2* flanking_bp_smoothing_window)
    return smoothed_reads

def analyze_differential_expression_around_features(motif_dict, annotations_type_dict, annotations_dict, query_plus_strand, query_minus_strand, reference_plus_strand, reference_minus_strand, reference_RNA_seq_plus_strand, reference_RNA_seq_minus_strand, upstream_window, downstream_bp, min_rpm, min_fold_upregulation, min_reference_value_to_divide_by = 0.5, prioritized_sequence_types = ('snoRNA', 'snoRNA_300_nt_downstream', 'snoRNA_300_nt_upstream', 'tRNA', 'tRNA_100_nt_downstream', 'tRNA_100_nt_downstream', '3UTR', 'CDS', 'intron', '5UTR',  'CUTs',  'SUTs', 'XUTs', 'antisense_CDS', 'intergenic_region'), 
                                                    sides_of_feature = 'upstream,downstream', min_RNA_seq_rpm = 1, flanking_bp_to_examine_RNA_seq = 50):
    '''
    used to differentiate roadblocked from passive binding sites
    takes an upstream window as a tuple to examine signal (e.g., (8,20) to examine the signal from -20 to -8 for + strand, and +8 to +20 for minus strand)
    takes a second window to examine signal that bypasses the roadblock, (eg, 8 upstream motif to 20 downstream the motif)
    a minimum reads per million (rpm) in query to consider, a minimum fold upregulation, and the query and reference dicts
    RNA-seq bedgraph dicts are used to identify evidence of readthrough the binding site
    returns two dicts, one with a minimum signal within window_size upstream of the motif upregulated in query over reference, and the other not upregulated over reference
    '''
    ## dicts to track annotation types of upregulated / not_upregulated hits (with at least min_rpm)
    plus_strand_upregulated_sequence_type_hits, plus_strand_not_upregulated_sequence_type_hits, minus_strand_upregulated_sequence_type_hits, minus_strand_not_upregulated_sequence_type_hits = {}, {}, {}, {}
    ## establish the region to search for upregulated reads    
    proximal_idx, distal_idx = upstream_window
    upregulated_sites, not_upregulated_sites, no_pA_with_RNA_seq, no_pA_with_no_RNA_seq = {}, {}, {}, {}
    different_site_types_dicts = upregulated_sites, not_upregulated_sites, no_pA_with_RNA_seq, no_pA_with_no_RNA_seq
    for motif_dict_strand in motif_dict:
        for chromosome in motif_dict[motif_dict_strand]:
            if 'M' not in chromosome and 'm' not in chromosome: ## confirm that mtDNA is not analyzed
                for my_dict in different_site_types_dicts:
                    if chromosome not in my_dict:
                        my_dict[chromosome] = []
                for start, end in motif_dict[motif_dict_strand][chromosome]:
                    print(chromosome,start,end)
                    best_motif_match, lowest_mismatches = motif_dict[motif_dict_strand][chromosome][(start, end)]
                    
                    ## first check the plus strand
                    strand = '+'
                    query_expression_values = return_expression_values(chromosome, start - distal_idx, start - proximal_idx, query_plus_strand)
                    query_expression_values_with_downstream_region = return_expression_values(chromosome, start - distal_idx, end + downstream_bp, query_plus_strand)
                    reference_expression_values = return_expression_values(chromosome, start - distal_idx, start - proximal_idx, reference_plus_strand)
                    query_expression = sum( query_expression_values )
                    reference_expression = sum( reference_expression_values )
                    start_of_upstream_window = start - distal_idx
                    peak_coord = start_of_upstream_window + query_expression_values.index( max(query_expression_values) )
                    dist_to_peak_coord = peak_coord - start
                    plus_strand_fold_upregulation = query_expression / max ( float(reference_expression), min_reference_value_to_divide_by)
                    query_plus_strand_expression = query_expression           
                    plus_strand_hit = False
                    plus_strand_flanking_pA_hit = False   
                    plus_strand_RNA_seq_hit = False
                    annotation_type, gene_name, gene_id, gene_name_strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, peak_coord, annotations_type_dict, annotations_dict, prioritized_sequence_types)
                    annotation = GTF_GFF_manipulation.get_annotation_of_coordinate(strand, chromosome, start, annotations_dict)
                    upregulated_sites_plus_strand_output = (strand, start, end, motif_dict_strand, best_motif_match, lowest_mismatches, 'query_pA_signal_upstream_feature', plus_strand_fold_upregulation, query_expression, dist_to_peak_coord, annotation_type, gene_name, gene_name_strandedness, annotation )
                    if max(query_expression_values) >= min_rpm:  ## region is recorded if the query expression is above a minimum threshold                        
                        if plus_strand_fold_upregulation >= min_fold_upregulation and max(query_expression_values) >= max(query_expression_values_with_downstream_region):                     
                            plus_strand_hit = True                    
                            if annotation_type not in plus_strand_upregulated_sequence_type_hits:
                                plus_strand_upregulated_sequence_type_hits[annotation_type] = 1
                            else:
                                plus_strand_upregulated_sequence_type_hits[annotation_type] += 1
                    ## check if the downstream window contains sufficient signal in the reference to call a passive binding site
                    flanking_window_reference_signal = return_expression_values(chromosome, start - distal_idx, end + downstream_bp, reference_plus_strand)
                    plus_strand_flanking_window_reference_signal = sum(flanking_window_reference_signal)
                    start_of_window = start - distal_idx
                    peak_coord = start_of_window + flanking_window_reference_signal.index( max(flanking_window_reference_signal) )
                    dist_to_peak_coord = peak_coord - start
                    if plus_strand_flanking_window_reference_signal >= min_rpm:
                        plus_strand_flanking_pA_hit = True
                        not_upregulated_sites_plus_strand_output = (strand, start, end, motif_dict_strand, best_motif_match, lowest_mismatches, 'flanking_window_reference_pA_signal', flanking_window_reference_signal, query_expression, dist_to_peak_coord, annotation_type, gene_name, gene_name_strandedness, annotation )
                        if annotation_type not in plus_strand_not_upregulated_sequence_type_hits:
                            plus_strand_not_upregulated_sequence_type_hits[annotation_type] = 1
                        else:
                            plus_strand_not_upregulated_sequence_type_hits[annotation_type] += 1
                    ## if there are no pA sites nearby, check whether the region shows any evidence of readthrough transcription by examining RNA seq for the reference 
                    reference_RNA_seq_expression_values = return_expression_values(chromosome, start - flanking_bp_to_examine_RNA_seq, end + flanking_bp_to_examine_RNA_seq, reference_RNA_seq_plus_strand)
                    ## query_RNA_seq_expression_values = return_expression_values(chromosome, start - flanking_bp_to_examine_RNA_seq, end + flanking_bp_to_examine_RNA_seq, query_RNA_seq_plus_strand)
                    mean_reference_plus_strand_expression = sum(reference_RNA_seq_expression_values) / (float(2*flanking_bp_to_examine_RNA_seq) + (end-start))
                    plus_strand_output = (strand, start, end, motif_dict_strand, best_motif_match, lowest_mismatches, 'flanking_reference_RNA_seq_expression_values', mean_reference_plus_strand_expression, min( reference_RNA_seq_expression_values ), 'N/A', annotation_type, gene_name, gene_name_strandedness, annotation )
                    if min( reference_RNA_seq_expression_values ) >= min_RNA_seq_rpm:
                        plus_strand_RNA_seq_hit = True
                    ## check the minus strand
                    strand = '-'
                    query_expression_values = return_expression_values(chromosome, end + proximal_idx, end + distal_idx, query_minus_strand)
                    query_expression_values_with_downstream_region = return_expression_values(chromosome, start - downstream_bp, end + distal_idx, query_minus_strand)
                    reference_expression_values =  return_expression_values(chromosome,  end + proximal_idx, end + distal_idx, reference_minus_strand)            
                    query_expression = sum( query_expression_values )
                    query_minus_strand_expression = query_expression           
                    reference_expression = sum( reference_expression_values )
                    start_of_upstream_window = end + proximal_idx
                    peak_coord = start_of_upstream_window + query_expression_values.index( max(query_expression_values) )
                    dist_to_peak_coord = end - peak_coord
                    minus_strand_fold_upregulation = query_expression / max ( float(reference_expression), min_reference_value_to_divide_by)
                    minus_strand_hit = False
                    minus_strand_flanking_pA_hit = False
                    minus_strand_RNA_seq_hit = False
                    annotation_type, gene_name,gene_id, gene_name_strandedness = GTF_GFF_manipulation.get_annotation_type_and_gene_of_coordinate(strand, chromosome, peak_coord, annotations_type_dict, annotations_dict, prioritized_sequence_types)
                    annotation = GTF_GFF_manipulation.get_annotation_of_coordinate(strand, chromosome, start, annotations_dict)                    
                    upregulated_sites_minus_strand_output = (strand, start, end, motif_dict_strand, best_motif_match, lowest_mismatches, 'query_pA_signal_upstream_feature', minus_strand_fold_upregulation, query_expression, dist_to_peak_coord, annotation_type, gene_name, gene_name_strandedness, annotation )
                    if max(query_expression_values) >= min_rpm:
                        if minus_strand_fold_upregulation >= min_fold_upregulation and max(query_expression_values) >= max(query_expression_values_with_downstream_region):                     
                            minus_strand_hit = True                    
                            if annotation_type not in minus_strand_upregulated_sequence_type_hits:
                                minus_strand_upregulated_sequence_type_hits[annotation_type] = 1
                            else:
                                minus_strand_upregulated_sequence_type_hits[annotation_type] += 1
                    ## check if the downstream window contains sufficient signal in the reference to call a passive binding site
                    flanking_window_reference_signal = return_expression_values(chromosome, start - downstream_bp , end + distal_idx, reference_plus_strand)
                    minus_strand_flanking_window_reference_signal = sum(flanking_window_reference_signal)                    
                    peak_coord = start - downstream_bp + flanking_window_reference_signal.index( max(flanking_window_reference_signal) )
                    dist_to_peak_coord = end - peak_coord
                    if minus_strand_flanking_window_reference_signal >= min_rpm:
                        minus_strand_flanking_pA_hit = True
                        not_upregulated_sites_minus_strand_output = (strand, start, end, motif_dict_strand, best_motif_match, lowest_mismatches, 'flanking_window_reference_pA_signal', flanking_window_reference_signal, query_expression, dist_to_peak_coord, annotation_type, gene_name, gene_name_strandedness, annotation )
                        if annotation_type not in minus_strand_not_upregulated_sequence_type_hits:
                            minus_strand_not_upregulated_sequence_type_hits[annotation_type] = 1
                        else:
                            minus_strand_not_upregulated_sequence_type_hits[annotation_type] += 1
                    ## if there are no pA sites nearby, check whether the region shows any evidence of readthrough transcription by examining RNA seq for the reference 
                    reference_RNA_seq_expression_values = return_expression_values(chromosome, start - flanking_bp_to_examine_RNA_seq, end + flanking_bp_to_examine_RNA_seq, reference_RNA_seq_minus_strand)
                    ## query_RNA_seq_expression_values = return_expression_values(chromosome, start - flanking_bp_to_examine_RNA_seq, end + flanking_bp_to_examine_RNA_seq, query_RNA_seq_minus_strand)
                    mean_reference_minus_strand_expression = sum(reference_RNA_seq_expression_values) / (float(2*flanking_bp_to_examine_RNA_seq) + (end-start))
                    ## mean_query_plus_strand_expression = sum(query_RNA_seq_expression_values) / (float(2*flanking_bp_to_examine_RNA_seq) + (end-start))
                    minus_strand_output = (strand, start, end, motif_dict_strand, best_motif_match, lowest_mismatches, 'flanking_reference_RNA_seq_expression_values', mean_reference_plus_strand_expression, min( reference_RNA_seq_expression_values ), 'N/A', annotation_type, gene_name, gene_name_strandedness, annotation )
                    if min( reference_RNA_seq_expression_values ) >= min_RNA_seq_rpm:
                        minus_strand_RNA_seq_hit = True

                    if plus_strand_hit and minus_strand_hit:
                        if query_plus_strand_expression >= query_minus_strand_expression: ## use the strand with higher expression
                            upregulated_sites[chromosome].append(upregulated_sites_plus_strand_output)
                        else:
                            upregulated_sites[chromosome].append(upregulated_sites_minus_strand_output)
                    elif plus_strand_hit:
                        upregulated_sites[chromosome].append(upregulated_sites_plus_strand_output)
                    elif minus_strand_hit:
                        upregulated_sites[chromosome].append(upregulated_sites_minus_strand_output)
                    
                    elif plus_strand_flanking_pA_hit and minus_strand_flanking_pA_hit:
                        if plus_strand_flanking_window_reference_signal >= minus_strand_flanking_window_reference_signal:
                            not_upregulated_sites[chromosome].append(not_upregulated_sites_plus_strand_output)
                        else:
                            not_upregulated_sites[chromosome].append(not_upregulated_sites_minus_strand_output)
                    elif plus_strand_flanking_pA_hit:
                        not_upregulated_sites[chromosome].append(not_upregulated_sites_plus_strand_output)
                    elif minus_strand_flanking_pA_hit:
                        not_upregulated_sites[chromosome].append(not_upregulated_sites_minus_strand_output)
                        
                    else:
                        if plus_strand_RNA_seq_hit and minus_strand_RNA_seq_hit:
                            if mean_reference_plus_strand_expression >= mean_reference_minus_strand_expression: 
                                no_pA_with_RNA_seq[chromosome].append(plus_strand_output)
                            else:
                                no_pA_with_RNA_seq[chromosome].append(minus_strand_output)
                        elif plus_strand_RNA_seq_hit:
                            no_pA_with_RNA_seq[chromosome].append(plus_strand_output)
                        elif minus_strand_RNA_seq_hit:
                            no_pA_with_RNA_seq[chromosome].append(minus_strand_output)
                        else:
                            if mean_reference_plus_strand_expression >= mean_reference_minus_strand_expression: 
                                no_pA_with_no_RNA_seq[chromosome].append(plus_strand_output)
                            else:
                                no_pA_with_no_RNA_seq[chromosome].append(minus_strand_output)
                                
    sequence_type_dicts = plus_strand_upregulated_sequence_type_hits, plus_strand_not_upregulated_sequence_type_hits, minus_strand_upregulated_sequence_type_hits, minus_strand_not_upregulated_sequence_type_hits
    return sequence_type_dicts, different_site_types_dicts

def motif_bedgraph_to_dict(filename, R64_genome, motif, max_mismatches_allowed, bp_buffer_outside_footprint = 0):
    '''
    input: bedgraph, with four tab delimited columns: chromosome, start, end, signal
    output: motif dictionary of chromosome dictionaries of coordinate dictionaries of the overlapping motif with chromosome values
    if reverse complement is found to be a better match, the strand is denoted as '-'
    also keeps track of motif variants found and prints to std out
    '''
    binding_sites_found = 0
    infile = open(filename, 'r')
    sample_dict = {}
    motif_variants = {}
    for line in infile:
        info = line.strip().split('\t')
        chromosome, start, end, reads = info[:4]
        if 'M' not in chromosome and 'm' not in chromosome:  ## remove mtDNA from analysis
            start, end = int(start), int(end)
            if chromosome not in sample_dict:
                sample_dict[chromosome] = {}
            start = start - bp_buffer_outside_footprint
            end = end + bp_buffer_outside_footprint
            sequence = R64_genome[chromosome][start : end]
            lowest_mismatches = len(motif)
            best_motif_match = None
            best_motif_strand = None
            for idx in range(len(sequence) - len(motif)):
                ## check plus strand for motif
                query_motif = sequence[idx:idx+len(motif)]
                query_mismatches = hamming_distance( motif, query_motif )
                if query_mismatches < lowest_mismatches:
                    lowest_mismatches = query_mismatches
                    best_motif_match = query_motif
                    best_motif_strand = '+'
                    best_motif_coordinates = start + idx, start + idx + len(motif) 
                ## check minus strand for motif
                rc_query_motif = rev_comp(query_motif)
                rc_query_mismatches = hamming_distance( motif, rc_query_motif )
                if rc_query_mismatches < lowest_mismatches:
                    lowest_mismatches = rc_query_mismatches
                    best_motif_match = rc_query_motif
                    best_motif_strand = '-'
                    best_motif_coordinates = start + idx, start + idx + len(motif) 
            if best_motif_match != None and lowest_mismatches <= max_mismatches_allowed: 
                if best_motif_match not in motif_variants:
                    motif_variants[best_motif_match] = 1
                else:
                    motif_variants[best_motif_match] += 1
                sample_dict[chromosome][ best_motif_coordinates ] = best_motif_strand, best_motif_match, lowest_mismatches 
                binding_sites_found += 1
    for variant in motif_variants:
        print(variant, motif_variants[variant])
    print(str(binding_sites_found) + ' binding sites found with at most ', str(max_mismatches_allowed), 'mismatches for ', filename)
    infile.close()
    print('bedgraph processed for the motif: ', motif, filename)
    return sample_dict

def combine_motif_dicts(motif_dict_1, motif_dict_2, union = True):
    '''
    takes dictionaries of chromosomes with keys as tuples of start and stop coordinates and values as a tuple of best_motif_strand, best_motif_match, lowest_mismatches ( the direction of the motif is indicated by the strand )
    returns the union of the two dictionaries, with duplicate entries collapsed
    if union is False, returns the intersection
    '''
    combined_dict = {}
    motif_dict_1_motif_occurrences = 0
    motif_dict_2_motif_occurrences = 0
    intersection_occurrences = 0
    for chromosome in motif_dict_1:
        if chromosome not in combined_dict:
            combined_dict[chromosome] = {}
        for tpl in motif_dict_1[chromosome]:
            motif_dict_1_motif_occurrences += 1
            if union:
                combined_dict[chromosome][tpl] = motif_dict_1[chromosome][tpl]
            else:
                if tpl in motif_dict_2[chromosome]:
                    if motif_dict_1[chromosome][tpl] != motif_dict_2[chromosome][tpl]:
                        print('motif location found with conflicting entries', motif_dict_1[chromosome][tpl], motif_dict_2[chromosome][tpl])
                    combined_dict[chromosome][tpl] = motif_dict_1[chromosome][tpl]
                    intersection_occurrences += 1
    for chromosome in motif_dict_2:
        for tpl in motif_dict_2[chromosome]:
            motif_dict_2_motif_occurrences += 1
            if union:
                if tpl not in combined_dict[chromosome]:
                    combined_dict[chromosome][tpl] = motif_dict_2[chromosome][tpl]
                else:
                    intersection_occurrences += 1
    print('Motifs in motif_dict_1: ', motif_dict_1_motif_occurrences)
    print('Motifs in motif_dict_2: ', motif_dict_2_motif_occurrences)
    print('Motifs in both motif_dict_1 and motif_dict_2: ', intersection_occurrences)
    return combined_dict

def write_motif_bedgraph_dict_to_file(bedgraph_dict, outfile_name):
    '''
    input: bedgraph dictionary, name for outfile
    output: file version of bedgraph dicitonary
    '''
    outfile = open(outfile_name, 'w')
    for chromosome in bedgraph_dict:
        for coordinates in bedgraph_dict[chromosome]:
            start, end = coordinates
            best_motif_strand, best_motif_match, lowest_mismatches  = bedgraph_dict[chromosome][ (start, end) ]
            start, end = str(start), str(end)
            output = '\t'.join( [chromosome, start, end, best_motif_strand, best_motif_match, str(lowest_mismatches) ] ) + '\n'
            outfile.write(output)
    outfile.close()
    print(outfile_name + ' complete')
    return None
    
def motif_bedgraph_file_to_dict(bedgraph_filename):
    '''
    input: file version of bedgraph dicitonary
    output: bedgraph dictionary    
    '''
    motif_bedgraph_dict = generate_empty_stranded_bedgraph_dict()
    with open(bedgraph_filename, 'r') as infile:
        for line in infile:
            chromosome, start, end, best_motif_strand, best_motif_match, lowest_mismatches  = line.strip().split()
            start, end = int(start), int(end)
            motif_bedgraph_dict[best_motif_strand][chromosome][(start, end)] = best_motif_match, lowest_mismatches
    infile.close()
    print(bedgraph_filename + ' complete')
    return motif_bedgraph_dict
    
def write_upregulated_peaks_near_motif_bedgraph_dict_to_file(bedgraph_dict, outfile_name, dist_required_between_nearby_motifs = 2):
    '''
    input: bedgraph dictionary, name for outfile
    output: file version of bedgraph dicitonary
    '''
    motifs_used = generate_empty_nonstranded_bedgraph_dict()
    outfile = open(outfile_name, 'w')
    header = '\t'.join(['cluster peak strand', 'cluster peak chromosome', 'motif start','motif end', 'motif_strand', 'motif_sequence', 'mismatches from consensus', 'category', 'criteria for calling site', 'total cluster signal', 'distance from motif boundary to cluster peak (bp)', 'type of region', 'gene name', 'gene strand' ]) + '\n'
    outfile.write(header)
    for chromosome in bedgraph_dict:
        for entry in bedgraph_dict[chromosome]:
            strand, start, end, best_motif_strand, best_motif_match, lowest_mismatches, category, criteria_enabling_hit, query_expression, dist_to_peak_coord, annotation_type, gene_name, gene_name_strandedness, annotation = entry
            motif_site_nearby_used = False            
            for coord in range( int(start) - dist_required_between_nearby_motifs, int(end) + dist_required_between_nearby_motifs ):
                if coord in motifs_used[chromosome]:
                    motif_site_nearby_used = True
                else:
                    motifs_used[chromosome][coord] = 1  ## DUMMY VARIABLE
            if not motif_site_nearby_used:
                output = '\t'.join([strand, chromosome, str(start), str(end), best_motif_strand, best_motif_match, str(lowest_mismatches), category, str(criteria_enabling_hit), str(query_expression), str(dist_to_peak_coord), annotation_type, gene_name, gene_name_strandedness ]) + '\n'
                outfile.write(output)
    outfile.close()
    print(outfile_name + ' complete')
    return None
    
def write_sequence_from_bedgraph_file_to_fasta(bedgraph_infilename, fasta_outfilename, genome_dict, upstream_bp_to_include = 100, downstream_bp_to_include = 0 ):
    '''
    takes a begraph with the first four columns tab-separated (strand, chromosome, start, end, )
    returns the upstream sequence (determined by the strand)
    '''
    outfile = open(fasta_outfilename, 'w')
    with open(bedgraph_infilename, 'r') as infile:
        infile.readline()  ## read past the header
        for line in infile:
            info = line.strip().split()
            strand, chromosome, start, end  = info[:4]
            outfile.write('>' + '_'.join([ chromosome, start, end, strand]) + '\n')
            start, end = int(start), int(end)
            if strand == '+':
                sequence = genome_dict[chromosome][start - upstream_bp_to_include: end + downstream_bp_to_include]    
            else:
                sequence = rev_comp( genome_dict[chromosome][start - downstream_bp_to_include : end + upstream_bp_to_include]    )   
            outfile.write(sequence + '\n')
    outfile.close()
    print('completed fasta output for: ', fasta_outfilename)
    return None
    
def write_sequence_from_text_file_to_fasta(bedgraph_infilename, fasta_outfilename, genome_dict, upstream_bp_to_include = 100, downstream_bp_to_include = 0 ):
    '''
    takes a begraph with the first four columns tab-separated (chromosome, start, end, strand)
    returns the upstream sequence (determined by the strand)
    '''
    outfile = open(fasta_outfilename, 'w')
    with open(bedgraph_infilename, 'r') as infile:
        infile.readline()  ## read past the header
        for line in infile:
            info = line.strip().split()
            chromosome, start, end, strand  = info[:4]
            outfile.write('>' + '_'.join([ chromosome, start, end, strand]) + '\n')
            start, end = int(start), int(end)
            if strand == '+':
                sequence = genome_dict[chromosome][start - upstream_bp_to_include: end + downstream_bp_to_include]    
            else:
                sequence = rev_comp( genome_dict[chromosome][start - downstream_bp_to_include : end + upstream_bp_to_include]    )   
            outfile.write(sequence + '\n')
    outfile.close()
    print('completed fasta output for: ', fasta_outfilename)
    return None
    
def analyze_fasta_for_motifs(fasta_files, motif_list):
    '''
    takes a list of fasta filenames, and a list of motifs to search for
    returns the number of times each motif appears in all the fasta files
    '''
    
    total_bp_analyzed = 0
    total_sequences_analyzed = 0
    motif_counts = {}
    for motif in motif_list:
        motif_counts[motif] = 0
    for filename in fasta_files:
        print(filename)
        with open(filename, 'r') as infile:
            total_sequences = 0
            for line in infile:
                
                if line[0] != '>':
                    total_sequences += 1
                    total_sequences_analyzed += 1
                    line = line.strip()
                    total_bp_analyzed += len(line)
                    for motif in motif_list:
                        motif_counts[motif] += line.count(motif)
        print(filename, 'total sequences: ', total_sequences)
        infile.close()
    if total_sequences > 0:
        for motif in motif_list:
            print(motif, 'motif_counts[motif]', motif_counts[motif], 'total_sequences_analyzed', total_sequences_analyzed)
            motif_counts[motif] = motif_counts[motif] / float(total_sequences_analyzed ) ## * 4**( len (motif) )
            
    return motif_counts

def read_bedgraph_dict_from_file(filename, three_prime_ends_coordinates_to_exclude):
    '''
    input: filename for strand
    output: bedgraph dict with strands as first keys
    '''
    sites = {}
    num_lines = 0
    dummy_idx = 1
    cluster_peaks_near_three_prime_ends = 0
    with open(filename, 'r') as infile:
        for line in infile:
            chromosome, start, end, strand  = line.strip().split('\t')[:4]
            start, end = int(start), int(end)
            current_motif_length = end-start
            if start not in three_prime_ends_coordinates_to_exclude[strand][chromosome]:
                if strand not in sites:
                    sites[strand] = {}
                if chromosome not in sites[strand]:
                    sites[strand][chromosome] = {}
                sites[strand][chromosome][dummy_idx] =  (start,end)
                dummy_idx += 1
                num_lines += 1
            else:
                cluster_peaks_near_three_prime_ends += 1
    print('num_lines', num_lines)
    print('cluster_peaks_near_three_prime_ends', cluster_peaks_near_three_prime_ends)
    return sites, current_motif_length
      
def read_cluster_center_of_mass_from_file(filename, three_prime_ends_coordinates_to_exclude):
    '''
    input: filename for strand
    output: bedgraph dict with strands as first keys
    '''
    sites = {}
    num_lines = 0
    dummy_idx = 1
    cluster_peaks_near_three_prime_ends = 0
    with open(filename, 'r') as infile:
        infile.readline()
        for line in infile:
            chromosome, peak_coord, center_of_mass_coord, strand   = line.strip().split('\t')[:4] 
            num_sites_explaining_signal =   line.strip().split('\t')[-2]
            if int(num_sites_explaining_signal) > 0:
                peak_coord, center_of_mass_coord = int(peak_coord), int(float(center_of_mass_coord))
                current_motif_length = 1
                if center_of_mass_coord not in three_prime_ends_coordinates_to_exclude[strand][chromosome]:
                    if strand not in sites:
                        sites[strand] = {}
                    if chromosome not in sites[strand]:
                        sites[strand][chromosome] = {}
                    sites[strand][chromosome][dummy_idx] =  (center_of_mass_coord,center_of_mass_coord + 1)
                    dummy_idx += 1
                    num_lines += 1
                else:
                    cluster_peaks_near_three_prime_ends += 1
    print('num_lines', num_lines)
    print('cluster_peaks_near_three_prime_ends', cluster_peaks_near_three_prime_ends)
    return sites, current_motif_length
      
def check_peak_for_downstream_footprint(annotations, footprint_bedgraph, downstream_window):
    '''
    takes an annotation bedgraph (representing cluster peaks) and the DNase I footprints
    separates annotations into two dicts depending on whether the peak contains a
    footprint within the specified window downstream
    '''
    with_footprint = 0
    without_footprint = 0
    annotations_with_downstream_footprint = generate_empty_stranded_bedgraph_dict()
    annotations_without_downstream_footprint = generate_empty_stranded_bedgraph_dict()
    for strand in annotations: 
        for chromosome in annotations[strand]:
            for annotation in annotations[strand][chromosome]: 
                peak_coord = annotations[strand][chromosome][annotation][0]
                if strand == '+':
                    footprint_found = False
                    for coord in range(peak_coord, peak_coord + downstream_window):
                        if coord in footprint_bedgraph[chromosome]:
                            annotations_with_downstream_footprint[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                            footprint_found = True
                            with_footprint += 1
                            break
                    if not footprint_found:
                        without_footprint += 1
                        annotations_without_downstream_footprint[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                elif strand == '-':
                    footprint_found = False
                    for coord in range(peak_coord - downstream_window, peak_coord):
                        if coord in footprint_bedgraph[chromosome]:
                            annotations_with_downstream_footprint[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
                            footprint_found = True
                            with_footprint += 1
                            break
                    if not footprint_found:
                        without_footprint += 1
                        annotations_without_downstream_footprint[strand][chromosome][annotation] = annotations[strand][chromosome][annotation]
    print('with_footprint', with_footprint, 'without_footprint', without_footprint)
    return annotations_with_downstream_footprint, annotations_without_downstream_footprint
              
def write_annotation_dict_to_file(annotation_dict, outfile_name):
    '''
    input: bedgraph dictionary, name for outfile
    output: file version of bedgraph dicitonary
    '''
    outfile = open(outfile_name, 'w')
    for strand in annotation_dict:
        for chromosome in annotation_dict[strand]:
            for annotation in annotation_dict[strand][chromosome]:
                peak_coord = annotation_dict[strand][chromosome][annotation][0]
                output = '\t'.join([chromosome, str(peak_coord), str(peak_coord + 1), strand ]) + '\n'
                outfile.write(output)
    outfile.close()
    print(outfile_name + ' complete')
    return None
    