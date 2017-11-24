'''
set of tools to extract information from annotation files in GTF or GFF format
'''
import operator, bedgraph_computation
import imp
bedgraph_computation = imp.reload(bedgraph_computation)

HARD_DRIVE =  '/Users/kevinroy/Dropbox/' # '/Volumes/SANDISK128/'
INPUT_DIR = HARD_DRIVE + 'datasets/'

GENOME_FILEPATH=INPUT_DIR+'yeast_genome/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta'
GTF_FILEPATH=INPUT_DIR+'yeast_genome/R64_annotations.gtf'
GFF_FILEPATH=INPUT_DIR+'yeast_genome/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113_no_seq.gff'
TIF_SEQ_FILEPATH=INPUT_DIR+'S1_TIFs.txt'

def convert_chromosome_number_to_numeral(chromosome_number):
    chromosome_numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI'.split(',')
    chromosome_number_to_numeral = {}
    for idx in range(16):
        chromosome_number_to_numeral['chr' + str(idx+1)] = 'chr' + chromosome_numerals[idx]
    return chromosome_number_to_numeral[chromosome_number]

def get_parent_id(string):
    '''
    input: annotation in gff format
    output: parent_ID
    
    chrII	SGD	gene	2907	5009	.	-	.	ID=YBL111C;Name=YBL111C;Ontology_term=GO:0003674,GO:0005739,GO:0008150;Note=Helicase-like%20protein%20encoded%20within%20the%20telomeric%20Y%27%20element%3B%20relocalizes%20from%20mitochondrion%20to%20cytoplasm%20upon%20DNA%20replication%20stress;display=Helicase-like%20protein%20encoded%20within%20the%20telomeric%20Y%27%20element;dbxref=SGD:S000002151;orf_classification=Verified
    chrII	SGD	CDS	2907	4116	.	-	1	Parent=YBL111C_mRNA;Name=YBL111C_CDS;orf_classification=Verified
    chrII	SGD	CDS	4216	5009	.	-	0	Parent=YBL111C_mRNA;Name=YBL111C_CDS;orf_classification=Verified
    chrII	SGD	intron	4117	4215	.	-	.	Parent=YBL111C_mRNA;Name=YBL111C_intron;orf_classification=Verified
    chrII	SGD	mRNA	2907	5009	.	-	.	ID=YBL111C_mRNA;Name=YBL111C_mRNA;Parent=YBL111C
    '''
    dist = len('Parent=')
    if 'Parent=' not in string:
        raise ValueError
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'Parent=':
            id_start = idx+dist
            break
    for idx in range(id_start, len(string)):
        if string[idx] == ';':
            id_end = idx
            break
    return string[id_start:id_end]

def get_systematic_gene_name(string):
    '''
    input: annotation in gff format
    output: systematic gene name
    '''
    dist = len('ID=')
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'ID=':
            id_start = idx+dist
            break
    for idx in range(id_start, len(string)):
        if string[idx] == ';':
            id_end = idx
            break
    return string[id_start:id_end]

def get_common_gene_name(string):
    '''
    input: annotation in gff format
    output: common gene name, or systematic if common not present
    '''
    if 'gene=' not in string:
        dist = len('Name=')
        for idx in range(len(string)):
            if string[idx:idx+dist] == 'Name=':
                id_start = idx+dist
                break
        for idx in range(id_start, len(string)):
            if string[idx] == ';':
                id_end = idx
                break
        return string[id_start:id_end]
    dist = len('gene=')
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'gene=':
            id_start = idx+dist
            break
    for idx in range(id_start, len(string)):
        if string[idx] == ';':
            id_end = idx
            break
    return string[id_start:id_end]
    
def load_systematic_to_common_gene_name_dict(GFF_FILENAME):
    '''
    input: gff file with intron annotations
    output: dictionary of chromosomes of start and stop coordinate tuples of strand and gene name tuple
    '''
    systematic_to_common_gene_name_dict = {}
    with open(GFF_FILENAME,'r') as gff_file:
        for line in gff_file:
            info = line.strip().split('\t')
            chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
            if chromosome != 'chrMito' and chromosome != 'chrmt':  ## exlude mitochondrial introns
                start, end = int(start), int(end)
                if sequence_type == 'gene':
                    systematic_gene_name = get_systematic_gene_name(annotation)
                    common_gene_name = get_common_gene_name(annotation)
                    systematic_to_common_gene_name_dict[systematic_gene_name] = common_gene_name
    gff_file.close()
    return systematic_to_common_gene_name_dict

def load_common_to_systematic_gene_name_dict(GFF_FILENAME):
    '''
    input: gff file with intron annotations
    output: dictionary of chromosomes of start and stop coordinate tuples of strand and gene name tuple
    '''
    common_to_systematic_gene_name_dict = {}
    with open(GFF_FILENAME,'r') as gff_file:
        for line in gff_file:
            info = line.strip().split('\t')
            chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
            if chromosome != 'chrMito' and chromosome != 'chrmt':  ## exlude mitochondrial introns
                start, end = int(start), int(end)
                if sequence_type == 'gene':
                    systematic_gene_name = get_systematic_gene_name(annotation)
                    common_gene_name = get_common_gene_name(annotation)
                    common_to_systematic_gene_name_dict[common_gene_name] = systematic_gene_name
    gff_file.close()
    return common_to_systematic_gene_name_dict

def get_gene_name_from_gtf_annotation(string):
    '''
    input: annotation in gtf format
    output: gene_id
    '''
    dist = len('gene_name ')
    ID_found = False
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'gene_name ':
            ID_found = True
            id_start = idx+dist+1
            break
    for idx in range(id_start, len(string)):
        if string[idx] == '"':
            id_end = idx
            break
    if ID_found:
        return string[id_start:id_end]
    else:
        return None


def get_gene_id_from_gtf_annotation(string):
    '''
    input: annotation in gtf format
    output: gene_id
    '''
    dist = len('gene_id ')
    ID_found = False
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'gene_id ':
            ID_found = True
            id_start = idx+dist+1
            break
    for idx in range(id_start, len(string)):
        if string[idx] == '"':
            id_end = idx
            break
    if ID_found:
        return string[id_start:id_end]
    else:
        return None

def get_gene_id_from_gff_annotation(string):
    '''
    input: annotation in gtf format
    output: gene_id
    '''
    if 'ID=' not in string:
        return None
    dist = len('ID=')
    ID_found = False
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'ID=':
            ID_found = True
            id_start = idx+dist
            break
    if not ID_found:
        return None
    id_end = None
    for idx in range(id_start, len(string)):
        if string[idx] == ';':
            id_end = idx
            break
    if id_end == None:
        return None
    return string[id_start:id_end]


def get_parent_from_gff_annotation(string):
    '''
    input: annotation in gtf format
    output: Parent_id
    '''
    dist = len('Parent=')
    ID_found = False
    for idx in range(len(string)):
        if string[idx:idx+dist] == 'Parent=':
            ID_found = True
            id_start = idx+dist
            break
    if not ID_found:
        return None
    for idx in range(id_start, len(string)):
        if string[idx] == ';':
            id_end = idx
            break
    if ID_found:
        return string[id_start:id_end]
    else:
        return None

def get_gene_name_from_gff_annotation(string):
    '''
    input: annotation in gtf format
    output: Parent_id
    '''
    if 'Name=' not in string:
        return None
    dist = len('Name=')
    ID_found = False
    for idx in range(len(string)):

        if string[idx:idx+dist] == 'Name=':
            ID_found = True
            id_start = idx+dist
            break
    if not ID_found:
        return None
    for idx in range(id_start, len(string)):
        if string[idx] == ';':
            id_end = idx
            break
    if ID_found:
        return string[id_start:id_end]
    else:
        return None

def load_ORF_coordinates(gff_filename):
    '''
    takes in a gff file with CDS coordinates,
    loads a dictionary with dict[ORF_name] = [ (chrom, start, end, strand) ] 
    removes _CDS suffix
    additional exons get appended as separate tuples
    '''
    ORF_dict = {}
    with open(gff_filename, 'r') as gff:
        for line in gff:
            if line[0] != '#':
                try:
                    chrom, source, region, start, end, period, strand, frame, annotation =  line.split('\t')
                    start, end = int(start), int(end)
                except:
                    print(line)
                if region == 'CDS' and chrom != 'chrmt':
                    gene_name = get_gene_name_from_gff_annotation(annotation)[:-4]
                    if gene_name not in ORF_dict:
                        ORF_dict[gene_name] = [ (chrom, start, end, strand) ]
                    else:
                        ORF_dict[gene_name] += [ (chrom, start, end, strand) ]
    return ORF_dict

def get_ORF_seq_and_CDS_coords(ORF, ORF_dict, genome_seq):
    '''
    takes an ORF
    returns a list of codons for the ORF
    '''
    exons = ORF_dict[ORF]
    ORF_seq = []
    exon_coords = []
    strand = exons[0][3]
    frame_offset = 0
    stop_codon_indices = []
    stop_codons = 'TAA','TAG','TGA'
    if strand == '+':
        for exon in exons:
            chrom, start, end, strand = exon
            exon_seq = genome_seq[chrom][start-1:end]
#            if ORF == 'YPL198W_CDS':
#                print(ORF, exon, exon_seq)
            if frame_offset > 0:
                previous_partial_codon = ORF_seq[-1]
                rest_of_codon = exon_seq[:frame_offset]
                ORF_seq[-1] = previous_partial_codon + rest_of_codon
                if previous_partial_codon + rest_of_codon in stop_codons:
                    stop_codon_indices.append(len(ORF_seq))
                    if ORF == 'YPL198W_CDS' or ORF == 'YMR079W_CDS':
                        print(ORF, previous_partial_codon, rest_of_codon, ORF_seq)
                for idx in range( frame_offset ): 
                    exon_coords.append( start - 1 + idx )
            for idx in range(frame_offset, len(exon_seq), 3):
                frame_offset = 0 ## erase the previous frame offset, and allow the frame offset to be set at the end of the exon (i.e. last iteration through this loop)
                codon = exon_seq[idx:idx + 3]
                ORF_seq.append(codon)
                if codon in stop_codons:
                    stop_codon_indices.append(len(ORF_seq))  ## indicate the nucleotide location of stop codons in the ORF
                for exon_idx in range(len(codon)):
                    exon_coords.append( (start - 1 + exon_idx + idx) )
                if len(codon) < 3:  ## partial codons must be completed with the next exon
                    frame_offset = 3 - len(codon)
                    #print(ORF, ORF_seq, 'ORF_seq')
                    if ORF == 'KU42_CDS' or ORF == 'YMR079W_CDS':
                        print(ORF, ORF_seq, 'ORF_seq')
                        print('incomplete_codon in ', exon, 'from ', exons)
    else:
        for exon in exons[::-1]:
            chrom, start, end, strand = exon
            exon_seq = bedgraph_computation.rev_comp( genome_seq[chrom][start-1:end] )
#            if ORF == 'YDR424C_CDS':
#                print(ORF, exon, exon_seq)
            if frame_offset > 0:
                previous_partial_codon = ORF_seq[-1]
                rest_of_codon = exon_seq[:frame_offset]
                ORF_seq[-1] = previous_partial_codon + rest_of_codon
                if previous_partial_codon + rest_of_codon in stop_codons:
                    stop_codon_indices.append(len(ORF_seq))
                #print(ORF, previous_partial_codon, rest_of_codon, ORF_seq)
                for idx in range(frame_offset): 
                    exon_coords.append( end - idx )
            for idx in range(frame_offset, len(exon_seq), 3):
                frame_offset = 0  ## erase the previous frame offset, and allow the frame offset to be set at the end of the exon
                codon = exon_seq[idx:idx + 3]
                ORF_seq.append(codon)
                if codon in stop_codons:
                    stop_codon_indices.append(len(ORF_seq))
                for exon_idx in range(len(codon)):
                    exon_coords.append( (end - exon_idx - idx) )
                if len(codon) < 3:
                    # print('incomplete_codon in ', exon, 'from ', exons)
                    frame_offset = 3 - len(codon)
                    # print(ORF, ORF_seq, 'ORF_seq')
    stop_codons_found = ORF_seq.count('TAA') + ORF_seq.count('TGA') + ORF_seq.count('TAG')
    if stop_codons_found > 1:
        print(ORF + ' has ' + str(stop_codons_found) + ' stop codons at positions ' + str(stop_codon_indices) + ' in the ORF', exons)
    if stop_codons_found == 0:
        print(ORF + ' has no stop codon', exons)
    return chrom, strand, exon_coords, ORF_seq
    
def get_aa_num_to_codon_coords_aa(ORF_chrom, ORF_strand, ORF_exon_coords, ORF_seq, codon_to_aa):
    '''
    takes an ORF seq and the coordinates involved in the ORF
    returns a dictionary with  ORF_info[1] = 'ATG', (335,336,337), 'M' 
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
    ORF_info = {}
    for idx in range( len(ORF_seq) ):
        aa_num = idx + 1
        codon = ORF_seq[idx].upper()
#        print (codon_to_aa)
#        print (codon)
        aa = codon_to_aa[codon]
        aa_coords = tuple( ORF_exon_coords[idx*3: (idx+1)* 3] )
        ORF_info[aa_num] = ( codon, aa_coords, aa )
    return ORF_info, aa_num_to_exon

def obtain_coordinates_for_annotation_type_from_gtf(annotation_type, filter_annotation_info = '', exclude_Mito = True ):
    '''
    takes a annotation type and returns a dictionary of strands of chromosomes 
    of lists of (start, stop, annotation) tuples for the annotation type
    specify filter_annotation_info with a string to select only features that contain the string
    in the annotation column
    '''
    infile = open(GTF_FILEPATH, 'r')
    annotations = {}
    for line in infile:
        info = line.strip().split('\t')
        chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
        if exclude_Mito:
            if region == annotation_type and ( filter_annotation_info in annotation )  and 'M' not in chromosome:
                if strand == '.':
                    strand = '+'
                if strand not in annotations:
                    annotations[strand] = {}
                if chromosome not in annotations[strand]:
                    annotations[strand][chromosome] = {}
                start = int(start)
                end = int(end)
                gene_id = get_gene_id_from_gtf_annotation(annotation)
                if gene_id not in annotations[strand][chromosome]:
                   annotations[strand][chromosome][gene_id] = [start, end]
                else:
                    if end > annotations[strand][chromosome][gene_id][1]:
                        annotations[strand][chromosome][gene_id][1] = end
        else:
            if region == annotation_type and ( filter_annotation_info in annotation ):
                if strand == '.':
                    strand = '+'
                if strand not in annotations:
                    annotations[strand] = {}
                if chromosome not in annotations[strand]:
                    annotations[strand][chromosome] = {}
                start = int(start)
                end = int(end)
                gene_id = get_gene_id_from_gtf_annotation(annotation)
                if gene_id not in annotations[strand][chromosome]:
                   annotations[strand][chromosome][gene_id] = [start, end]
                else:
                    if end > annotations[strand][chromosome][gene_id][1]:
                        annotations[strand][chromosome][gene_id][1] = end
    return annotations

def obtain_coordinates_for_annotation_type_from_gff(annotation_type, filter_annotation_info = '' ):
    '''
    takes a annotation type and returns a dictionary of strands of chromosomes 
    of lists of (start, stop, annotation) tuples for the annotation type
    specify filter_annotation_info with a string to select only features that contain the string
    in the annotation column
    '''
    infile = open(GFF_FILEPATH, 'r')
    annotations = {}
    for line in infile:
        info = line.strip().split('\t')
        chromosome, SGD, sequence_type, start, end, period, strand, frame, annotation = info
        if sequence_type == annotation_type and ( filter_annotation_info in annotation ) and chromosome != 'chrmt':
            if strand == '.':
                strand = '+'
            if strand not in annotations:
                annotations[strand] = {}
            if chromosome not in annotations[strand]:
                annotations[strand][chromosome] = {}
            start = int(start)
            end = int(end)
            gene_id = get_gene_id_from_gff_annotation(annotation)
            if gene_id == None:
                gene_id = get_parent_from_gff_annotation(annotation)
            if gene_id not in annotations[strand][chromosome]:
               annotations[strand][chromosome][gene_id] = [start, end]
            else:
                if end > annotations[strand][chromosome][gene_id][1]:
                    annotations[strand][chromosome][gene_id][1] = end
        if sequence_type == 'snoRNA_gene':
            if 'C%2FD' not in annotation and 'H%2FACA' not in annotation:
                print('a snoRNA that is not box CD or HACA:', annotation)
    return annotations

def obtain_all_annotations_from_gtf():
    '''
    returns a dictionary of strands of chromosomes of start coordinates of 
    [ (stop coordinate, gene_name, sequence_type, region_type) ]
    '''
    infile = open(GTF_FILEPATH, 'r')
    longest_annotation_length = 0
    annotations = {}
    for strand in '+-':
        annotations[strand] = {}
    for line in infile:
        both_strands = False
        info = line.strip().split('\t')
        chromosome, region_type, sequence_type, start, end, period, strand, frame, annotation = info
        if sequence_type != 'chromosome':
            if strand == '.':
                strand = '+'
                both_strands = True
            if chromosome not in annotations['-']:
                annotations['-'][chromosome] = {}
            if chromosome not in annotations['+']:
                annotations['+'][chromosome] = {}
            start = int(start)
            end = int(end)
            annotation_length = end - start
            if annotation_length > longest_annotation_length:
                longest_annotation_length = annotation_length
            gene_name = get_gene_name_from_gtf_annotation(annotation)
            if gene_name == None:
                gene_name = get_gene_id_from_gtf_annotation(annotation)
            if start not in annotations[strand][chromosome]:
                annotations[strand][chromosome][start] = []
            annotations[strand][chromosome][start].append( (end, gene_name, sequence_type, region_type ) )
            if both_strands:
                strand = '-'
                if start not in annotations[strand][chromosome]:
                    annotations[strand][chromosome][start] = []
                annotations[strand][chromosome][start].append( (end, gene_name, sequence_type, region_type ) )
    print('GTF annotations loaded')
    infile.close()
    return annotations

def obtain_all_annotations_from_gff():
    '''
    returns a dictionary of strands of chromosomes of start coordinates of 
    [ (stop coordinate, gene_name, sequence_type, region_type) ]
    '''
    infile = open(GFF_FILEPATH, 'r')
    annotations = {}
    longest_annotation_length = 0
    for strand in '+-':
        annotations[strand] = {}
    for line in infile:
        both_strands = False
        info = line.strip().split('\t')
        chromosome, region_type, sequence_type, start, end, period, strand, frame, annotation = info
        if sequence_type != 'chromosome':
            if strand == '.':
                strand = '+'
                both_strands = True
            if chromosome not in annotations['-']:
                annotations['-'][chromosome] = {}
            if chromosome not in annotations['+']:
                annotations['+'][chromosome] = {}
            start = int(start)
            end = int(end)
            annotation_length = end - start
            if annotation_length > longest_annotation_length:
                longest_annotation_length = annotation_length
            gene_name = get_gene_from_gff_annotation(annotation)
            if gene_name == None:
                gene_name = get_gene_id_from_gff_annotation(annotation)
            if gene_name == None:
                gene_name = get_parent_from_gff_annotation(annotation)
            if gene_name == None:
                gene_name = 'intergenic'
            if start not in annotations[strand][chromosome]:
                annotations[strand][chromosome][start] = []
            annotations[strand][chromosome][start].append( (end, gene_name, sequence_type, region_type ) )
            if both_strands:
                strand = '-'
                if start not in annotations[strand][chromosome]:
                    annotations[strand][chromosome][start] = []
                annotations[strand][chromosome][start].append( (end, gene_name, sequence_type, region_type ) )
    infile.close()
    print('GFF annotations loaded')
    return annotations

def obtain_UTRs_from_TIF_seq(required_fraction_of_ORF_spanning_reads_for_UTR_calculation):
    '''
    returns the UTR termini as determined by TIF seq
    '''
    ## check sequence for downstream AG richness and screen these out because they likely represent internal poly(A) reverse transcription
    ## parameters for filtering
    DOWNSTREAM_NT_TO_CHECK = 10
    MAX_G_ALLOWED = 4 ## G can form a wobble pair
    MAX_CT_ALLOWED = 1
    
    genome_dict = bedgraph_computation.load_genome(GENOME_FILEPATH)
    chromosome_numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI'.split(',')
    chromosome_number_to_numeral = {}
    for idx in range(16):
        chromosome_number_to_numeral[str(idx+1)] = chromosome_numerals[idx]
    
    chromosome_lengths = {}
    for chromosome in genome_dict:
        chromosome_length = len(genome_dict[chromosome])
        chromosome_lengths[chromosome] = chromosome_length
    
    ORF_dictionary = {}
    for strand in '+','-':
        ORF_dictionary[strand] = {}
    flagged_isoforms = 0
    total_ORF_spanning_isoforms = 0
    with open(TIF_SEQ_FILEPATH, 'r') as infile:
        for line in infile:
            line = line.strip('\n')
            if "Covering one intact ORF" in line:
                total_ORF_spanning_isoforms += 1
                idx = line.find("Y", 0, -5)
                ORF = line[idx:]
                info = line.split(' ')
                chromosome, strand, t5, t3, YPD_reads = info[:5]
                chromosome = 'chr' + chromosome_number_to_numeral[chromosome]
                chromosome_length = chromosome_lengths[chromosome]
                t5 = int(t5)
                t3 = int(t3)
                A_count, G_count, CT_count = 0,0,0
                flagged = False
                if strand == '+':
                    for temp_coordinate in range(t3+1,min(chromosome_length,t3+DOWNSTREAM_NT_TO_CHECK)):
                        base = genome_dict[chromosome][temp_coordinate]
                        if base == 'A':
                            A_count += 1
                        elif base == 'G':
                            G_count += 1
                        else:
                            CT_count += 1
                    if ( G_count <= MAX_G_ALLOWED and CT_count <= MAX_CT_ALLOWED ):
                        flagged = True
                        flagged_isoforms += 1
                else:
                    A_count, G_count, CT_count = 0,0,0
                    flagged = False
                    for temp_coordinate in range(max(t3-DOWNSTREAM_NT_TO_CHECK,0),t3):
                        base = genome_dict[chromosome][temp_coordinate]
                        if base == 'T':
                            A_count += 1
                        elif base == 'C':
                            G_count += 1
                        else:
                            CT_count += 1
                    if ( G_count <= MAX_G_ALLOWED and CT_count <= MAX_CT_ALLOWED ):
                        flagged = True
                        flagged_isoforms += 1
                if not flagged:
                    YPD_reads = int(YPD_reads)
                    if YPD_reads > 0:
                        if chromosome not in ORF_dictionary[strand]:
                            ORF_dictionary[strand][chromosome] = {}
                        if ORF not in ORF_dictionary[strand][chromosome]:
                            ORF_dictionary[strand][chromosome][ORF] = {}
                            ORF_dictionary[strand][chromosome][ORF]['5UTR'] = []
                            ORF_dictionary[strand][chromosome][ORF]['3UTR'] = []
                        ORF_dictionary[strand][chromosome][ORF]['5UTR'].append( (t5, YPD_reads))
                        ORF_dictionary[strand][chromosome][ORF]['3UTR'].append( (t3, YPD_reads))
    infile.close()
    print('ORF dictionary loaded')
    print('total ORF spanning isoforms', total_ORF_spanning_isoforms)
    print('total number of flagged isoforms', flagged_isoforms)
    
    ## for each ORF, order the 5'UTR and 3'UTR coordinates in order of distance from the ORF
    ## compute the total reads for all 5'UTR coordinates, and all 3'UTR coordinates
    ## starting at the ORF-proximal UTR coordinate, walk along the coordinates until REQUIRED_FRACTION_FOR_TRANSCRIPT_TERMINI % of total UTR signal is accounted for
    ## this coordinate will be the start/end of >= REQUIRED_FRACTION_FOR_TRANSCRIPT_TERMINI % transcripts from this ORF    
    ## load the start and stop codons for each ORF
    ORF_start_stop_coordinates = {}
    infile = open(GTF_FILEPATH, 'r')
    
    ## create a dictionary for converting gene_id to gene_name
    gene_id_to_name = {}
    ORF_start_stop_coordinates = {}
    for strand in '+','-':
        ORF_start_stop_coordinates[strand] = {}
        ## load protein-coding annotations as a dictionary of strands of ORF systematic names with value as a list of lists of sequence feature types
        ## these can be of the type 3UTR, 5UTR, CDS, exon, start_codon, stop_codon followed by start and end coordinates
        ## e.g. ['5UTR', 123, 234]
    
    for line in infile:
        info = line.strip().split('\t')
        chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
        if 'M' not in chromosome:
            end = int(end)
            start = int(start)
            if region == 'protein_coding':
                gene_id = get_gene_id_from_gtf_annotation(annotation)
                gene_name = get_gene_name_from_gtf_annotation(annotation)
                gene_id_to_name[gene_id] = gene_name
                if gene_id not in ORF_start_stop_coordinates[strand]:
                    # the first entry in the ORF annotations is the chromosome
                    ORF_start_stop_coordinates[strand][gene_id] = [0,0]
                if strand == '+':
                    if sequence_type == 'start_codon':
                        ORF_start_stop_coordinates[strand][gene_id][0] = start
                    if sequence_type == 'stop_codon':
                        ORF_start_stop_coordinates[strand][gene_id][1] = end
                elif strand == '-':
                    if sequence_type == 'start_codon':
                        ORF_start_stop_coordinates[strand][gene_id][0] = end
                    if sequence_type == 'stop_codon':
                        ORF_start_stop_coordinates[strand][gene_id][1] = start
    infile.close()
    print('protein_coding annotations loaded')
    revised_termini = {}
    for strand in '+','-':
        revised_termini[strand] = {}
        for chromosome in ORF_dictionary[strand]:
            revised_termini[strand][chromosome] = {}
            for ORF in ORF_dictionary[strand][chromosome]:
                if strand == '+':
                    ORF_dictionary[strand][chromosome][ORF]['5UTR'] = sorted(ORF_dictionary[strand][chromosome][ORF]['5UTR'], reverse=True)
                    ORF_dictionary[strand][chromosome][ORF]['3UTR'] = sorted(ORF_dictionary[strand][chromosome][ORF]['3UTR'])
                else:
                    ORF_dictionary[strand][chromosome][ORF]['5UTR'] = sorted(ORF_dictionary[strand][chromosome][ORF]['5UTR'])
                    ORF_dictionary[strand][chromosome][ORF]['3UTR'] = sorted(ORF_dictionary[strand][chromosome][ORF]['3UTR'], reverse=True)
                fiveUTR_sum = 0
                threeUTR_sum = 0
                for coord in ORF_dictionary[strand][chromosome][ORF]['5UTR']:
                    fiveUTR_sum += coord[1]
                for coord in ORF_dictionary[strand][chromosome][ORF]['3UTR']:
                    threeUTR_sum += coord[1]
                fiveUTR_running_total = 0
                for coord in ORF_dictionary[strand][chromosome][ORF]['5UTR']:
                    fiveUTR_running_total += coord[1]
                    if fiveUTR_running_total / float(fiveUTR_sum) >= required_fraction_of_ORF_spanning_reads_for_UTR_calculation:
                        fiveUTR_terminus = coord[0]
                        break
                threeUTR_running_total = 0
                for coord in ORF_dictionary[strand][chromosome][ORF]['3UTR']:
                    threeUTR_running_total += coord[1]
                    if threeUTR_running_total / float(threeUTR_sum) >= required_fraction_of_ORF_spanning_reads_for_UTR_calculation:
                        threeUTR_terminus = coord[0]
                        break
                revised_termini[strand][ORF] = (fiveUTR_terminus, threeUTR_terminus)
    
    chromosome_letters = 'ABCDEFGHIJKLMNOP'
    chromosome_numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI'.split(',')
    chromosome_letter_to_number = {}
    for idx in range(len(chromosome_numerals)):
        numeral = chromosome_numerals[idx]
        letter = chromosome_letters[idx]
        chromosome_letter_to_number[letter] = numeral
    UTRs_filepath = INPUT_DIR +'UTR_termini_accounting_for_' + str(required_fraction_of_ORF_spanning_reads_for_UTR_calculation*100) + '%_of_transcripts_screened_for_internal_AG_mispriming.txt'
    UTRs = open(UTRs_filepath, 'w')
    for strand in '+','-':
        for ORF in revised_termini[strand]:
            if ORF[1] in chromosome_letter_to_number:
                chromosome = 'chr' + chromosome_letter_to_number[ORF[1]]
                ## output coordinates for five UTR start, start codon, stop codon, three UTR end
                ORF_start = str( ORF_start_stop_coordinates[strand][ORF][0] )
                ORF_end = str( ORF_start_stop_coordinates[strand][ORF][1] )
                UTR_info = '\t'.join([strand, chromosome,  ORF,  str(revised_termini[strand][ORF][0]), ORF_start, ORF_end, str(revised_termini[strand][ORF][1]) ])  + '\n'
                UTRs.write(UTR_info)
    UTRs.close()

    ORF_T_dict= {}
    num_of_ORF_Ts = 0
    with open(UTRs_filepath, 'r') as infile:
        for line in infile:
            strand, chromosome, gene_id, five_prime_UTR_start, start_codon_start, stop_codon_end, three_prime_UTR_end = line.strip().split()
            stop_codon_end, three_prime_UTR_end = int(stop_codon_end), int(three_prime_UTR_end)
            five_prime_UTR_start, start_codon_start = int(five_prime_UTR_start), int(start_codon_start)
            # three_UTR_length = abs(stop_codon_end - three_prime_UTR_end)
            if strand not in ORF_T_dict:
                ORF_T_dict[strand] = {}
            if chromosome not in ORF_T_dict[strand]:
                ORF_T_dict[strand][chromosome] = {}
            if strand == '+':
                ORF_T_dict[strand][chromosome][five_prime_UTR_start] = [(start_codon_start, gene_id_to_name[gene_id], gene_id, 'protein_coding_gene','5UTR' )]
                ## ORF_T_dict[strand][chromosome][start_codon_start] = [(stop_codon_end, ORF_name, 'CDS', 'protein_coding_gene')]
                ORF_T_dict[strand][chromosome][stop_codon_end] = [(three_prime_UTR_end, gene_id_to_name[gene_id], gene_id, 'protein_coding_gene', '3UTR' )]
            else:
                ORF_T_dict[strand][chromosome][start_codon_start] = [(five_prime_UTR_start, gene_id_to_name[gene_id], gene_id, 'protein_coding_gene', '5UTR')]
                ## ORF_T_dict[strand][chromosome][stop_codon_end ] = [(start_codon_start, ORF_name, 'CDS', 'protein_coding_gene')]
                ORF_T_dict[strand][chromosome][three_prime_UTR_end] = [(stop_codon_end, gene_id_to_name[gene_id], gene_id,  'protein_coding_gene', '3UTR')]
            num_of_ORF_Ts += 1
    infile.close()
    print('num_of_ORF_Ts', num_of_ORF_Ts)
    return ORF_T_dict, revised_termini

def obtain_all_annotations_from_gtf_UTR_from_TIF_seq_introns_from_gff(required_fraction_of_ORF_spanning_reads_for_UTR_calculation, CSXUT_FLANKING_NT = 50, TIF_seq_and_GTF_inclusive = True):
    '''
    returns a dictionary of strands of chromosomes of start coordinates of [ (stop coordinate, gene_name, sequence_type, region_type) ]
    '''

    annotations, revised_termini = obtain_UTRs_from_TIF_seq(required_fraction_of_ORF_spanning_reads_for_UTR_calculation)   
    infile = open(GFF_FILEPATH, 'r')
    for line in infile:
        info = line.strip().split('\t')
        chromosome, region_type, sequence_type, start, end, period, strand, frame, annotation = info
        if 'intron' in sequence_type and 'm' not in chromosome: ## exclude mitochondrial genes
            region_type = 'intron'
            start = int(start)
            end = int(end) + 1
            gene_name = get_gene_name_from_gff_annotation(annotation)
            gene_id = get_gene_id_from_gff_annotation(annotation)
            if gene_id == None:
                gene_id = 'intergenic'
            if gene_name == None:
                gene_name = gene_id
            if start not in annotations[strand][chromosome]:
                annotations[strand][chromosome][start] = []  ## create a list to allow for overlapping annotations
            annotations[strand][chromosome][start].append( (end, gene_name, gene_id, region_type, sequence_type ) )
    infile.close()
    print('GFF intron annotations loaded')
    
    infile = open(GTF_FILEPATH, 'r')
    for line in infile:
        info = line.strip().split('\t')
        chromosome, region_type, sequence_type, start, end, period, strand, frame, annotation = info
        if 'M' not in chromosome and 'exon' in sequence_type: ## exclude mitochondrial genes
            start = int(start)
            end = int(end) + 1
            if region_type in ('CUTs','SUTs','XUTs'):  ## termini for *UTs is ill-defined so add a buffer on either side
                end += CSXUT_FLANKING_NT
                start -= CSXUT_FLANKING_NT
            gene_name = get_gene_name_from_gtf_annotation(annotation)
            gene_id = get_gene_id_from_gtf_annotation(annotation)
            if gene_id == None:
                gene_id = 'intergenic'
            if gene_name == None:
                gene_name = gene_id
            if start not in annotations[strand][chromosome]:
                annotations[strand][chromosome][start] = []
            annotations[strand][chromosome][start].append( (end, gene_name, gene_id, region_type, sequence_type ) )
    print('GTF annotations loaded')
    infile.close()
    return annotations

def get_annotation_of_coordinate_for_splice_sites(strand, chromosome, coordinate, genome_dict, annotations_dict, pA_site_sequencing = False):
    '''
    takes a chromosomal coordinate and an annotations dictionary
    returns a list of annotations for the coordinate
    '''
    overlapping_annotations = ['unannotated splice site', 'intergenic' + '\t' + 'intergenic' + '\t' + 'intergenic' + '\t' + 'intergenic']
    annotation_found = False
    overlapping_ORF = False
    intron_found = False
    temp_coordinate = coordinate
    annotation_precedence_for_region_types = ('intron', 'CDS', '5UTR',  '3UTR')
    if pA_site_sequencing:
        if strand == '+':
            while genome_dict[chromosome][temp_coordinate + 1] == 'A':
                temp_coordinate += 1
        else:
            while genome_dict[chromosome][temp_coordinate - 1] == 'T':
                temp_coordinate -= 1
    for idx in range(temp_coordinate + 2, max(temp_coordinate - 10000, 0), -1):
        
        if idx in annotations_dict[strand][chromosome]:
            start = idx
            for region_type_priority in annotation_precedence_for_region_types:
                for entry in annotations_dict[strand][chromosome][start]:
                    region_type = entry[2]
                    if region_type_priority == region_type:
                        if region_type == 'intron':
                            intron_found = True
                        break
                        break
                        ## determine if coordinate of interest lies within an intron
            end, gene_name, sequence_type, region_type  = entry
            if 'intron' in sequence_type:
                if strand == '+':
                    if abs(temp_coordinate - start) < 5:
                        overlapping_annotations[0] = 'annotated splice site'
                else:
                    if abs(temp_coordinate - end) < 5:
                        overlapping_annotations[0] = 'annotated splice site'
            if temp_coordinate <= end + 2:
                annotation_found = True
                overlapping_annotations[1] = '\t'.join([region_type, sequence_type, gene_name, 'sense'])
                if sequence_type == 'CDS':
                    overlapping_ORF = True
                if intron_found:
                    break
            break
    if not overlapping_ORF:
        ## if the coordinate does not lie within an ORF, check whether the coordinate is antisense to an ORF
        if strand == '+':
            strand = '-'
        else:
            strand = '+'
        for idx in range(coordinate, max(coordinate - 10000, 0), -1):
            if idx in annotations_dict[strand][chromosome]:
                start = idx
                for entry in annotations_dict[strand][chromosome][start]:
                    if entry[2] in annotation_precedence_for_region_types:
                        end, gene_name, sequence_type, region_type  = entry
                        if coordinate <= end:
                            annotation_found = True
                            overlapping_annotations[1] = '\t'.join([region_type, sequence_type, gene_name, 'antisense'])
                            break
                break
#    if annotation_found:
#        print overlapping_annotations
#    else:
#        print strand, chromosome, coordinate, 'No annotation found!'
    return overlapping_annotations
    

def get_annotation_of_coordinate(strand, chromosome, coordinate, annotations_dict):
    '''
    takes a chromosomal coordinate and an annotations dictionary
    returns a list of annotations for the coordinate
    '''
    overlapping_annotations = [ 'intergenic', 'intergenic', 'intergenic', 'sense' ]
    priority_for_sequence_types = ('intron', 'exon', '5UTR',  '3UTR')
    for idx in range(coordinate, max(coordinate - 20000, 0), -1):
        if idx in annotations_dict[strand][chromosome]:
            start = idx
            for entry in annotations_dict[strand][chromosome][start]:
                end, gene_name, gene_id, region_type, sequence_type   = entry   ## region type is one of protein_coding, XUTs, CUTs, SUTs, tRNA, intergenic_region, rRNA, etc.
                                                                                ## sequence type is 5UTR, exon, 3UTR, intron, etc.
                for sequence_type_priority in priority_for_sequence_types:
                    if sequence_type_priority == sequence_type and coordinate <= end + 1:
                        # print(gene_name)
                        overlapping_annotations = [ sequence_type, gene_name, gene_id, 'sense' ]
                        break
                        break
                   # else:
                        # print sequence_type, sequence_type_priority
    if overlapping_annotations == [ 'intergenic', 'intergenic', 'intergenic', 'sense' ]:
        ## check whether the coordinate is antisense to an ORF
        if strand == '+':
            strand = '-'
        else:
            strand = '+'
        for idx in range(coordinate, max(coordinate - 20000, 0), -1):
            if idx in annotations_dict[strand][chromosome]:
                start = idx
                for entry in annotations_dict[strand][chromosome][start]:
                    end, gene_name, gene_id, region_type, sequence_type = entry
                    if sequence_type in priority_for_sequence_types and coordinate <= end + 1:
                        overlapping_annotations = [sequence_type, gene_name, gene_id, 'antisense' ] 
    return overlapping_annotations

def load_annotation_types(required_fraction_of_ORF_spanning_reads_for_UTR_calculation = .95, CSXUT_FLANKING_NT = 50, BP_TO_ADD_TO_UTRs = 50, FLANKING_BP_TO_ADD_TO_MISC_ANNOTATIONS = 100, BP_TO_ADD_DS_SNORNA = 300, BP_TO_ADD_US_SNORNA = 300, BP_TO_ADD_DS_TRNA = 100, BP_TO_ADD_US_TRNA = 100, TIF_seq_and_GTF_inclusive = True):
    '''
    loads an annotation for every chromosomal coordinate and returns the dictionary of strands of chromosomes of coordinates of annotations
    while slow to load and requiring much memory, this allows for very fast retrieval of chromosomal coordinate annotations for analyzing bedgraphs once loaded
    '''
    ## obtain 5' and 3' transcript termini that account for >= X % of ORF-Ts
    ## incorporate these data into a dictionary of strands of systematic ORF names with values as a two-item list of transcript termini (5'UTR,3'UTR)
    annotations, ORF_transcript_termini = obtain_UTRs_from_TIF_seq(required_fraction_of_ORF_spanning_reads_for_UTR_calculation)

    ## load annotations GTF file (converted from the pyCRAC annotations in EF2 to R64 coordinates)
    ## load protein-coding annotations as a dictionary of strands of ORF systematic names with value as a list of lists of sequence feature types
    ## these can be of the type 3UTR, 5UTR, CDS, exon, start_codon, stop_codon followed by start and end coordinates
    ## e.g. ['5UTR', 123, 234]
    ## For CUTs, SUTs, XUTs, allow for flexibility in start and stop coordinates with a CSXUT_FLANKING_NT value 
    infile = open(GTF_FILEPATH, 'r')
    annotations = {}
    ORF_annotations = {}
    for strand in '+','-':
        annotations[strand] = {}
        ORF_annotations[strand] = {}
    for line in infile:
        info = line.strip().split('\t')
        chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
        if chromosome == 'chrMito':
            chromosome = 'Mito'
        end = int(end) + 1
        start = int(start)
        if region == 'protein_coding':
            gene_id = get_gene_id_from_gtf_annotation(annotation)
            if gene_id not in ORF_annotations[strand]:
                # the first entry in the ORF annotations is the chromosome
                ORF_annotations[strand][gene_id] = [chromosome]
            ORF_annotations[strand][gene_id].append([sequence_type,start,end])
    infile.close()
    print('protein_coding annotations loaded')

    ## ORFs with more than one exon have an intervening intron that needs annotating
    for strand in '+','-':
        for ORF in ORF_annotations[strand]:
            exon_count = 0
            exons = []
            sequence_type_list = ORF_annotations[strand][ORF]
            for sequence_type in sequence_type_list:
                if 'exon' in sequence_type:
                    exon_count += 1
                    exons.append(sequence_type[1:])
            if exon_count > 1:
                exons = sorted(exons)
                for intron_num in range(len(exons)-1):
                    ORF_annotations[strand][ORF].append([ 'intron', exons[intron_num][1], exons[intron_num+1][0] ])
    print('introns added to protein_coding annotations')

    # correct pyCRAC R64 3'UTR coordinates with those obtained by Pelechano et al.
    # ORF dictionary are the coordinates from Pelechano et al.
    # ORF annotations are from pyCRAC annotations
    for strand in ['+','-']:
        for ORF in ORF_annotations[strand]:
            sequence_type_list = ORF_annotations[strand][ORF]
            three_prime_UTR_present = False
            five_prime_UTR_present = False
            if strand == '+':
                for sequence_type in sequence_type_list:
                    if 'start_codon' in sequence_type:
                        start_codon_begin = sequence_type[1]
                    elif 'stop_codon' in sequence_type:
                        stop_codon_end = sequence_type[2]
                for sequence_type in sequence_type_list:
                    if '3UTR' in sequence_type:
                        three_prime_UTR_present = True
                        end = sequence_type[2]
                        if ORF in ORF_transcript_termini[strand]:
                            TIF_end = ORF_transcript_termini[strand][ORF][1]
                            if TIF_seq_and_GTF_inclusive and TIF_end > end: 
                                sequence_type[2] = TIF_end
                            elif not TIF_seq_and_GTF_inclusive:
                                sequence_type[2] = TIF_end
                    if '5UTR' in sequence_type:
                        five_prime_UTR_present = True
                        start = sequence_type[1]
                        if ORF in ORF_transcript_termini[strand]:
                            TIF_start = ORF_transcript_termini[strand][ORF][0]
                            if TIF_seq_and_GTF_inclusive and TIF_start < start:     
                                sequence_type[1] = TIF_start
                            elif not TIF_seq_and_GTF_inclusive:
                                sequence_type[1] = TIF_start
                ## use Pelechano coordinates if pyCRAC UTR coordinates are unavailable
                if not three_prime_UTR_present:
                    if ORF in ORF_transcript_termini[strand]:
                        TIF_end = ORF_transcript_termini[strand][ORF][1]
                        ORF_annotations['+'][ORF].append( ['3UTR', stop_codon_end, TIF_end ] )
                if not five_prime_UTR_present:
                    if ORF in ORF_transcript_termini[strand]:
                        TIF_start = ORF_transcript_termini[strand][ORF][0]
                        ORF_annotations['+'][ORF].append( ['5UTR', TIF_start, start_codon_begin ] )
            elif strand == '-':
                for sequence_type in sequence_type_list:
                    if 'start_codon' in sequence_type:
                        start_codon_begin = sequence_type[2]
                    elif 'stop_codon' in sequence_type:
                        stop_codon_end = sequence_type[1]
                for sequence_type in sequence_type_list:
                    if '3UTR' in sequence_type:
                        three_prime_UTR_present = True
                        end = sequence_type[1]
                        if ORF in ORF_transcript_termini[strand]:
                            TIF_end = ORF_transcript_termini[strand][ORF][1]
                            if TIF_seq_and_GTF_inclusive and TIF_end < end: 
                                sequence_type[1] = TIF_end
                            elif not TIF_seq_and_GTF_inclusive:
                                sequence_type[1] = TIF_end
                    if '5UTR' in sequence_type:
                        five_prime_UTR_present = True
                        start = sequence_type[2]
                        if ORF in ORF_transcript_termini[strand]:
                            TIF_start = ORF_transcript_termini[strand][ORF][0]
                            if TIF_seq_and_GTF_inclusive and TIF_start > start:     
                                sequence_type[2] = TIF_start
                            elif not TIF_seq_and_GTF_inclusive:
                                sequence_type[2] = TIF_start
                ## use Pelechano coordinates if pyCRAC UTR coordinates are unavailable
                if not three_prime_UTR_present:
                    if ORF in ORF_transcript_termini[strand]:
                        TIF_end = ORF_transcript_termini[strand][ORF][1]
                        ORF_annotations['-'][ORF].append( ['3UTR', TIF_end, stop_codon_end ] )
                if not five_prime_UTR_present:
                    if ORF in ORF_transcript_termini[strand]:
                        TIF_start = ORF_transcript_termini[strand][ORF][0]
                        ORF_annotations['-'][ORF].append( ['5UTR', start_codon_begin, TIF_start  ] )
    print('5UTRs and 3UTRs from ORF annotations corrected with termini from TIF-seq')

    ## add snoRNA and tRNA annotations to the dictionary
    infile = open(GTF_FILEPATH, 'r')
    for line in infile:
        info = line.strip().split('\t')
        chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
        start = int(start)
        end = int(end) + 1
        if region in ('CUTs','SUTs','XUTs'):
            end += CSXUT_FLANKING_NT
            start -= CSXUT_FLANKING_NT   
        if chromosome not in annotations[strand]:
            annotations[strand][chromosome] = {}
        if region != 'protein_coding':
            for coordinate in range(start,end):
                if coordinate not in annotations[strand][chromosome]:
                    annotations[strand][chromosome][coordinate] = [region]
                else:
                    annotations[strand][chromosome][coordinate] += [region]
            region_to_be_added = False
            if region == 'snoRNA':
                region_to_be_added = True
                added_region = 'snoRNA_' + str(BP_TO_ADD_DS_SNORNA) + '_nt_downstream'
                if strand == '+':
                    added_start = end
                    added_end = end + BP_TO_ADD_DS_SNORNA
                elif strand == '-':
                    added_start = start - BP_TO_ADD_DS_SNORNA
                    added_end = start
            if region == 'tRNA':
                region_to_be_added = True
                added_region = 'tRNA_' + str(BP_TO_ADD_DS_TRNA) + '_nt_downstream'
                if strand == '+':
                    added_start = end
                    added_end = end + BP_TO_ADD_DS_TRNA
                elif strand == '-':
                    added_start = start - BP_TO_ADD_DS_TRNA
                    added_end = start
            if region_to_be_added:
                for coordinate in range(added_start,added_end):
                    if coordinate not in annotations[strand][chromosome]:
                        annotations[strand][chromosome][coordinate] = [added_region]
                    else:
                        annotations[strand][chromosome][coordinate] += [added_region]
            if region == 'snoRNA':
                region_to_be_added = True
                added_region = 'snoRNA_' + str(BP_TO_ADD_US_SNORNA) + '_nt_upstream'
                if strand == '+':
                    added_start = start - BP_TO_ADD_US_SNORNA
                    added_end = start
                elif strand == '-':
                    added_start = end
                    added_end = end + BP_TO_ADD_US_SNORNA
                for coordinate in range(added_start,added_end):
                    if coordinate not in annotations[strand][chromosome]:
                        annotations[strand][chromosome][coordinate] = [added_region]
                    else:
                        annotations[strand][chromosome][coordinate] += [added_region]
            if region == 'tRNA':
                region_to_be_added = True
                added_region = 'tRNA_' + str(BP_TO_ADD_US_TRNA) + '_nt_upstream'
                if strand == '+':
                    added_start = start - BP_TO_ADD_US_TRNA
                    added_end = start
                elif strand == '-':
                    added_start = end
                    added_end = end + BP_TO_ADD_US_TRNA
                for coordinate in range(added_start,added_end):
                    if coordinate not in annotations[strand][chromosome]:
                        annotations[strand][chromosome][coordinate] = [added_region]
                    else:
                        annotations[strand][chromosome][coordinate] += [added_region]
    infile.close()
    print('snoRNA and tRNA annotations by sequence position loaded')

    ## load DNA-based sequence types from the SGD .gff file
    misc_sequence_types = 'LTR_retrotransposon','binding_site','centromere','ARS','external_transcribed_spacer_region','internal_transcribed_spacer_region','long_terminal_repeat'
    infile = open(GFF_FILEPATH, 'r')
    for line in infile:
        info = line.strip().split('\t')
        chromosome, source, sequence_type, start, end, period, strand, frame, annotation = info
        end = int(end) + 1        
        if chromosome == 'chrMito':
            chromosome = 'Mito'
        if sequence_type in misc_sequence_types:
            if strand == '.':
                strand = '+'
                end = int(end) + FLANKING_BP_TO_ADD_TO_MISC_ANNOTATIONS
                start = int(start) - FLANKING_BP_TO_ADD_TO_MISC_ANNOTATIONS
                for coordinate in range(start,end):
                    if coordinate not in annotations[strand][chromosome]:
                        annotations[strand][chromosome][coordinate] = [sequence_type]
                    else:
                        annotations[strand][chromosome][coordinate] += [sequence_type]
                strand = '-'
                end = int(end) + FLANKING_BP_TO_ADD_TO_MISC_ANNOTATIONS
                start = int(start) - FLANKING_BP_TO_ADD_TO_MISC_ANNOTATIONS
                for coordinate in range(start,end):
                    if coordinate not in annotations[strand][chromosome]:
                        annotations[strand][chromosome][coordinate] = [sequence_type]
                    else:
                        annotations[strand][chromosome][coordinate] += [sequence_type]
            else:
                end = int(end) + FLANKING_BP_TO_ADD_TO_MISC_ANNOTATIONS
                start = int(start) - FLANKING_BP_TO_ADD_TO_MISC_ANNOTATIONS
                for coordinate in range(start,end):
                    if coordinate not in annotations[strand][chromosome]:
                        annotations[strand][chromosome][coordinate] = [sequence_type]
                    else:
                        annotations[strand][chromosome][coordinate] += [sequence_type]
    infile.close()
    print('Miscellaneous annotations by sequence position loaded')

    ## extend 3UTRs by a given window, and include antisense CDS as an annotation
    ## add the ORF annotations to the entire anntation dictionary
    for strand in ['+','-']:
        for ORF in ORF_annotations[strand]:
            chromosome = ORF_annotations[strand][ORF][0]
            if chromosome in annotations[strand]:
                for sequence_type in ORF_annotations[strand][ORF][1:]:
                    region = sequence_type[0]
                    start = sequence_type[1]
                    end = sequence_type[2]
                    if region == 'CDS':
                        ## include stop codon as CDS sequence
                        if strand == '+':
                            end += 3
                        else:
                            start -= 3
                    if region == '3UTR':
                        ## add N bases to 3'UTR regions to account for differences in sequencing/mapping methods between TIF-seq and DRS
                        if strand == '+':
                            end += BP_TO_ADD_TO_UTRs
                        else:
                            start -= BP_TO_ADD_TO_UTRs
                    for coordinate in range(start,end):
    
                        if coordinate not in annotations[strand][chromosome]:
                            annotations[strand][chromosome][coordinate] = [region]
                        else:
                            annotations[strand][chromosome][coordinate] += [region]
                    ## add antisense CDS as a region
                    if region == 'CDS':
                        if strand == '+':
                            antisense_strand = '-'
                        else:
                            antisense_strand = '+'
                        for coordinate in range(start,end):
                            if coordinate not in annotations[antisense_strand][chromosome]:
                                annotations[antisense_strand][chromosome][coordinate] = ['antisense_CDS']
                            else:
                                annotations[antisense_strand][chromosome][coordinate] += ['antisense_CDS']
            else:
                print(chromosome, 'not in annotations dict, triggered by:', ORF)
    print('all ORF annotations loaded')
    print()
    return annotations

def get_annotation_type_and_gene_of_coordinate(strand, chromosome, coordinate, annotations_type_dictionary, annotations_dictionary, prioritized_sequence_types):
    '''
    takes a coordinate and the two different types of annotations dictionaries, one for the sequence region types and one for the gene
    returns the annotation type, gene name (common name), gene id (systematic name), and strandedness relative to gene name as a tuple
    '''
    ## find the annotation type of the max peak
    sequence_types_of_coordinate = annotations_type_dictionary[strand][chromosome].get(coordinate, 'intergenic')
    annotation_found = False
    for sequence_type in prioritized_sequence_types:
        if sequence_type in sequence_types_of_coordinate:
            annotation_type = sequence_type
            annotation_found = True
            break
    if not annotation_found:
        annotation_type = 'intergenic'
        ## find the annotated gene in which the max peak resides
    annotation = get_annotation_of_coordinate(strand, chromosome, coordinate, annotations_dictionary) ## for ex., returns [ gene_name, gene_id, 'antisense' ] 
#    if annotation[0] != 'intergenic':
#        annotation_type = annotation[0]
	## [sequence_type, gene_name, gene_id, 'antisense' ] 
    if coordinate == 822146:
        print('TESTING')
        print(annotation_type, annotation[1], annotation[2], annotation[3])
    return annotation_type, annotation[1], annotation[2], annotation[3]

def global_pA_site_distribution_by_mass(plus_strand_bedgraph, minus_strand_bedgraph, annotations_dictionary, sequence_types_priority):
    '''
    takes pA site data in the form of two bedgraphs from both strands, the annotations dictionary, and the sequence types priority list to give preference to one sequence type in the case of coordinates with more than one annotation
    returns a sorted list of percentages for each sequence type
    '''
    sequence_type_counts = {}
    for sequence_type in sequence_types_priority:
        sequence_type_counts[sequence_type] = 0
    strand = '+'
    for chromosome in plus_strand_bedgraph:
        for coordinate in plus_strand_bedgraph[chromosome]:
            sequence_types_of_coordinate = annotations_dictionary[strand][chromosome].get(coordinate, ['intergenic'])
            annotation_found = False
            for sequence_type in sequence_types_priority:
                if sequence_type in sequence_types_of_coordinate:
                    sequence_type_counts[sequence_type] += plus_strand_bedgraph[chromosome][coordinate]
                    annotation_found = True
                    break
#            if not annotation_found:
#                print 'No annotation for', chromosome, coordinate, strand
    strand = '-'
    for chromosome in minus_strand_bedgraph:
        for coordinate in minus_strand_bedgraph[chromosome]:
            sequence_types_of_coordinate = annotations_dictionary[strand][chromosome].get(coordinate, ['intergenic'])
            annotation_found = False
            for sequence_type in sequence_types_priority:
                if sequence_type in sequence_types_of_coordinate:
                    sequence_type_counts[sequence_type] += minus_strand_bedgraph[chromosome][coordinate]
                    annotation_found = True
                    break
#            if not annotation_found:
#                print 'No annotation for', chromosome, coordinate, strand
    sorted_sequence_type_counts = sorted(list(sequence_type_counts.items()), key=operator.itemgetter(1))
    print(sorted_sequence_type_counts)
    return sorted_sequence_type_counts
    
def global_pA_site_distribution_by_coordinate(plus_strand_bedgraph, minus_strand_bedgraph, annotations_dictionary, sequence_types_priority):
    '''
    takes pA site data in the form of two bedgraphs from both strands, the annotations dictionary, and the sequence types priority list to give preference to one sequence type in the case of coordinates with more than one annotation
    returns a sorted list of percentages for each sequence type
    '''
    sequence_type_counts = {}
    for sequence_type in sequence_types_priority:
        sequence_type_counts[sequence_type] = 0
    strand = '+'
    for chromosome in plus_strand_bedgraph:
        for coordinate in plus_strand_bedgraph[chromosome]:
            sequence_types_of_coordinate = annotations_dictionary[strand][chromosome].get(coordinate, ['intergenic'])
            annotation_found = False
            for sequence_type in sequence_types_priority:
                if sequence_type in sequence_types_of_coordinate:
                    sequence_type_counts[sequence_type] += 1
                    annotation_found = True
                    break
            if not annotation_found:
                print('No annotation for', chromosome, coordinate, strand)
    strand = '-'
    for chromosome in minus_strand_bedgraph:
        for coordinate in minus_strand_bedgraph[chromosome]:
            sequence_types_of_coordinate = annotations_dictionary[strand][chromosome].get(coordinate, ['intergenic'])
            annotation_found = False
            for sequence_type in sequence_types_priority:
                if sequence_type in sequence_types_of_coordinate:
                    sequence_type_counts[sequence_type] += 1
                    annotation_found = True
                    break
            if not annotation_found:
                print('No annotation for', chromosome, coordinate, strand)
    sorted_sequence_type_counts = sorted(list(sequence_type_counts.items()), key=operator.itemgetter(1))
    return sorted_sequence_type_counts

def combine_annotations(sorted_sequence_type_counts, sequence_types_to_combine):
    '''
    takes the output of global_pA_site_distribution_by_mass, 
    which is a sorted list of tuples of (sequence type, counts) 
    and a tuple of tuples of sequence types to combine
    
    such as : 
    sorted_sequence_type_counts = {'snoRNA':100, 'snoRNA_300_nt_downstream':200, 'snoRNA_300_nt_upstream':50}
    sequence_types_to_combine = ( ('snoRNA', 'snoRNA_300_nt_downstream', 'snoRNA_300_nt_upstream'), ('tRNA', 'tRNA_100_nt_downstream', 'tRNA_100_nt_downstream'), ('3UTR',), ('CDS',), ('intron',), ('5UTR',), ('CUTs',), ('SUTs', 'XUTs'), ('antisense_CDS',), ('intergenic_region',) )  ## 'snRNA_300_nt_downstream', 'tRNA_100_nt_downstream', 'tRNA_100_nt_upstream'    
    
    returns a new sorted list of (combined_sequence types, counts)
    '''
    combined_sequence_type_dict = {}
#    print sequence_types_to_combine
#    print sorted_sequence_type_counts
    for sequence_types_tuple in sequence_types_to_combine:
#        print sequence_types_tuple
        combined_sequence_type_dict[sequence_types_tuple] = 0
    for sequence_type_and_counts in sorted_sequence_type_counts:
        sequence_type, counts = sequence_type_and_counts
#            print sequence_type, sorted_sequence_type_counts[sequence_type]
        for sequence_types_tuple in sequence_types_to_combine:
            for sequence_type_tpl in sequence_types_tuple:
                if sequence_type == sequence_type_tpl:
                    combined_sequence_type_dict[sequence_types_tuple] += counts
                    
#    print combined_sequence_type_dict
    sorted_combined_sequence_type_counts = sorted(list(combined_sequence_type_dict.items()), key=operator.itemgetter(1))
#    print sorted_combined_sequence_type_counts    
    return sorted_combined_sequence_type_counts

def get_three_prime_ends_with_flanking_region(lst, annotations, flanking_bp_to_include = 15):
    '''
    takes a list of annotation types, and the entire annotation list
    returns the chromosomal end coordinates for the snoRNAs in strand, chromosome, coord dictionary format
    '''
    selected_annotations = {}
    for strand in annotations:
        if strand not in selected_annotations:
            selected_annotations[strand] = {}
        for chromosome in annotations[strand]:
            if chromosome not in selected_annotations[strand]:
                selected_annotations[strand][chromosome] = []
            for start in annotations[strand][chromosome]:
                for entry in annotations[strand][chromosome][start]:
                    end, gene_name, sequence_type, region_type = entry
                    if region_type in lst:
                        if strand == '+':
                            for idx in range(end - flanking_bp_to_include, end + flanking_bp_to_include):
                                selected_annotations[strand][chromosome].append(idx) 
                        elif strand == '-':
                            for idx in range(start - flanking_bp_to_include, start + flanking_bp_to_include):
                                selected_annotations[strand][chromosome].append(idx)
    return selected_annotations
    