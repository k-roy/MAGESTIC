
'''
The goal of this script is to combine the sequencing replicates for each step 1 guide-donor-bc0 plasmid library.

This will require reverse complementing the KR1884-KR1953 products to have the same orientation as KR1952-KR1590, so that 4 datasets collapse into one.
Each PCR product had two replicates.

The different sequencing dates will be kept separate for now to allow for comparing Illumina 2 x 150 bp vs AVITI 2 x 150 bp vs 2 x 300 bp.

# inner_primers	for SpCas9 libraries
# KR1952-KR1590	   HNHNHtcggagctgcgattg gcaggcgcgccNNNNNNNNNNNNNNNNNNNNgtttgaagagcNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGATAGACGACTGGACAGCANNNNNNNNNNAGGAAAACAGAC AGTAACTCAGATTCVVDNDN
# KR1884-KR1953	   NHNHBBGAATCTGAGTTACT GTCTGTTTTCCTNNNNNNNNNNTGCTGTCCAGTCGTCTATCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNgctcttcaaacNNNNNNNNNNNNNNNNNNNNggcgcgcctgc caatcgcagctccgaDNDND
# KR1884-KR1953	rc HNHNHtcggagctgcgattg gcaggcgcgccNNNNNNNNNNNNNNNNNNNNgtttgaagagcNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGATAGACGACTGGACAGCANNNNNNNNNNAGGAAAACAGAC AGTAACTCAGATTCVVDNDN

# inner primers for LbCas12a libraries
# KR1952-KR1590	   HNHNHtcggagctgcgattg gcaggtttcaaagattaaataatttctactaagtgtagatNNNNNNNNNNNNNNNNNNNNtttcgaagagcNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNgttatgctggtcctaggtcgNNNNNNNNNNAGGAAAACAGAC AGTAACTCAGATTCVVDNDN
# KR1884-KR1953	   NHNHBBGAATCTGAGTTACT GTCTGTTTTCCTNNNNNNNNNNcgacctaggaccagcataacNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNgctcttcgaaaNNNNNNNNNNNNNNNNNNNNatctacacttagtagaaattatttaatctttgaaacctgc caatcgcagctccgaDNDND
# KR1884-KR1953 rc HNHNHtcggagctgcgattg gcaggtttcaaagattaaataatttctactaagtgtagatNNNNNNNNNNNNNNNNNNNNtttcgaagagcNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNgttatgctggtcctaggtcgNNNNNNNNNNAGGAAAACAGAC AGTAACTCAGATTCVVDNDN

input: collapsed fastq where the read name has the counts
@23
GCAGGCGCGCCGTGGGGTCAAATTGTTGCTGGTTTGAAGAGCCCATTACCTTCTCAGGGACCTGCCGCAGCGGCTCCCCCTGTACCTCAGCAACAATTCGATCCAACAGCACAGCTAAATTCTTTGATGAATATGCTTAACCAACAGCAGCAGCAACAACAACAAAGCTAACGATAGACGACTGGACAGCATAGTTTGGATAGGAAAACAGAC
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ

sample key with columns:
sample_number	library_ID	inner_primers	outer_primers	sequencing_dates
sample_1	V622	KR1952-KR1590	plate_5A	20250707,
sample_2	V622	KR1952-KR1590	plate_5A	20250707,
sample_3	V622	KR1884-KR1953	plate_5A	20250707,
sample_4	V622	KR1884-KR1953	plate_5A	20250707,

output: V plasmid tsv files with columns for 
for counts, leader, guide, RE_site, donor, SPS, bc0

The counts for sequencing replicates will be combined but kept separate by sequencing run and fwd counts (KR1952-KR1590) vs rev counts (KR1884-KR1953).
There will also be counts for total fwd and rev across all sequencing runs and replicates.
The counts columns are fwd_0607, rev_0607, fwd_0820, rev_0820, fwd_0821, rev_0821, total_fwd_counts, total_rev_counts.

sample sequence above split based on features:
GCAGGCGCGCC GTGGGGTCAAATTGTTGCTG GTTTGAAGAGC CCATTACCTTCTCAGGGACCTGCCGCAGCGGCTCCCCCTGTACCTCAGCAACAATTCGATCCAACAGCACAGCTAAATTCTTTGATGAATATGCTTAACCAACAGCAGCAGCAACAACAACAAAGCTAA CGATAGACGACTGGACAGCA TAGTTTGGAT AGGAAAACAGAC
'''

import pysam
import timeit
import os
from tqdm import tqdm

def rev_comp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdNn ', 'TGCAtgcaYRKMyrkmBVDHbvdhNn '))[::-1]

PROJECT_PATH = 'by_project/NNS/202250628_repeat_step_1_library_cloning/'

PROCESSED_DATA_DIR = '/path/to/processed_data/' + PROJECT_PATH
os.chdir(PROCESSED_DATA_DIR)
COLLAPSED_FASTQ_DIR = PROCESSED_DATA_DIR + 'collapsed_fastq/'
COUNTS_DIR = PROCESSED_DATA_DIR + 'guide_donor_bc0/'
for SUB_DIR in 'guide_donor_bc0/', 'log/':
    if not os.path.exists(PROCESSED_DATA_DIR + SUB_DIR):
        os.mkdir(PROCESSED_DATA_DIR + SUB_DIR)

KEYFILE_DIR = '/path/to/scripts_and_keyfiles/' + PROJECT_PATH + 'keyfiles/'
sample_keyfile = KEYFILE_DIR + 'step_1_sample_key.tsv'

max_reads_per_sample = 10000000000000000000
guide_length = 20
RE_site_length = len('GTTTGAAGAGC')
bc0_length= 10
SPS_length = 20
msd_fragment_length = 12
donor_fragment_length = 60
SPS_end_index = -msd_fragment_length - bc0_length
donor_end_index = -SPS_length - msd_fragment_length - bc0_length

sequencing_run_dates = '20250707',
project_name = '202250628_repeat_step_1_library_cloning'

def process_step_1_counts(step1_filename, inner_primers, sequencing_run_date, library_ID_to_bc0_info_to_counts):
    total_bc0_seqs_processed = 0
    total_counts_processed = 0
    start_time = timeit.default_timer()
    reads_processed_for_sample = 0
    with pysam.FastxFile(step1_filename) as fin:
        for entry in fin:
            reads_processed_for_sample += 1
            if reads_processed_for_sample > max_reads_per_sample:
                break
            total_bc0_seqs_processed += 1
            counts = int(entry.name)
            total_counts_processed += counts
            seq = entry.sequence
            qual = entry.quality
            seq_direction = 'fwd'
            if inner_primers in ('KR1884-KR1953', 'KR1882-KR1953'):
                seq = rev_comp(seq)
                qual = qual[::-1]
                if inner_primers == 'KR1882-KR1953':
                    seq += 'GAC'
                    qual += 'EEE'
                seq_direction = 'rev'
            if 'GTTTGAAGAGC' in seq:
                RE_site_found = True
                RE_site_index = seq.index('GTTTGAAGAGC')
                leader = seq[:RE_site_index-guide_length]
                guide = seq[RE_site_index-guide_length:RE_site_index]
                donor = seq[RE_site_index + RE_site_length:donor_end_index]
                SPS = seq[donor_end_index:SPS_end_index]
                bc0 = seq[SPS_end_index:-msd_fragment_length]
                msd = seq[-msd_fragment_length:]
            else:
                RE_site_found = False
                leader, guide, donor, SPS, bc0, msd = ['NA']*6
            donor_bc0_fragment = seq[-msd_fragment_length-donor_fragment_length:-msd_fragment_length]
            bc0_info = leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment, seq
            if bc0_info not in library_ID_to_bc0_info_to_counts[library_ID]:
                library_ID_to_bc0_info_to_counts[library_ID][bc0_info] = {}
                for seq_run_date in sequencing_run_dates:
                    for seq_dir in 'fwd', 'rev':
                        counts_name_initializer = seq_dir + '_' + seq_run_date[-4:]
                        library_ID_to_bc0_info_to_counts[library_ID][bc0_info][counts_name_initializer] = 0
            counts_name = seq_direction + '_' + sequencing_run_date[-4:]
            library_ID_to_bc0_info_to_counts[library_ID][bc0_info][counts_name] += counts
    fin.close()
    elapsed = timeit.default_timer() - start_time
    print(step1_filename + ' took ' + str(elapsed) + ' seconds to process\n')
    print(step1_filename + ' had ' + str(total_bc0_seqs_processed) + ' total_bc0_seqs_processed\n')
    print(step1_filename + ' had ' + str(total_counts_processed) + ' total_counts_processed\n')
    return(library_ID_to_bc0_info_to_counts)

# process files for each library_ID according to the sample key
library_ID_to_sample_info = {}
with open(sample_keyfile) as infile:
    header = infile.readline()
    for line in infile:
        sample_number, library_ID, inner_primers, outer_primers, dates_of_sequencing = line.strip().split('\t')
        if library_ID not in library_ID_to_sample_info:
            library_ID_to_sample_info[library_ID] = [(sample_number, inner_primers, outer_primers, dates_of_sequencing)]
        else:
            library_ID_to_sample_info[library_ID] += [(sample_number, inner_primers, outer_primers, dates_of_sequencing)]

for library_ID in tqdm(library_ID_to_sample_info, desc="Processing library_IDs"):
    library_ID_to_bc0_info_to_counts = {}
    logfilename = 'log/' + library_ID + '_step_1_guide_donor_bc0_counts_summary.log'
    logfile = open(logfilename, 'w')
    counts_outfilename = COUNTS_DIR + library_ID +'_guide_donor_bc0_counts.tsv'
    for sample_info in library_ID_to_sample_info[library_ID]:
        sample_number, inner_primers, outer_primers, dates_of_sequencing = sample_info
        for sequencing_run_date in dates_of_sequencing.split(','):
            if sequencing_run_date.strip() == '':
                continue
            sequencing_run_date = sequencing_run_date.strip()
            sample_name = sequencing_run_date + '_' + outer_primers + '_' + sample_number
            step1_filename = COLLAPSED_FASTQ_DIR + sample_name + '_collapsed.fastq'
            library_ID_to_bc0_info_to_counts = process_step_1_counts(step1_filename, inner_primers, sequencing_run_date, library_ID_to_bc0_info_to_counts)

    # tally total fwd and rev counts
    start_time = timeit.default_timer()
    for bc0_info in library_ID_to_bc0_info_to_counts[library_ID]:
        total_counts, total_fwd_counts, total_rev_counts = 0, 0, 0
        for counts_name in library_ID_to_bc0_info_to_counts[library_ID][bc0_info]:
            total_counts += library_ID_to_bc0_info_to_counts[library_ID][bc0_info][counts_name]
            if 'fwd' in counts_name:
                total_fwd_counts += library_ID_to_bc0_info_to_counts[library_ID][bc0_info][counts_name]
            else:
                total_rev_counts += library_ID_to_bc0_info_to_counts[library_ID][bc0_info][counts_name]
        library_ID_to_bc0_info_to_counts[library_ID][bc0_info]['total_counts'] = total_counts
        library_ID_to_bc0_info_to_counts[library_ID][bc0_info]['total_fwd_counts'] = total_fwd_counts
        library_ID_to_bc0_info_to_counts[library_ID][bc0_info]['total_rev_counts'] = total_rev_counts
        start_time = timeit.default_timer()

    # Dynamically determine all sequencing run date suffixes present in the data
    all_counts_names = set()
    for bc0_info in library_ID_to_bc0_info_to_counts[library_ID]:
        all_counts_names.update(library_ID_to_bc0_info_to_counts[library_ID][bc0_info].keys())
    # Separate out total counts and per-run counts
    total_counts_fields = ['total_counts', 'total_fwd_counts', 'total_rev_counts']
    per_run_counts_fields = sorted([name for name in all_counts_names if name not in total_counts_fields])
    header = ['library_ID'] + total_counts_fields + per_run_counts_fields + ['leader', 'guide', 'donor', 'SPS', 'bc0', 'msd', 'donor_bc0_fragment', 'seq']
    with open(counts_outfilename, 'w') as counts_outfile:
        counts_outfile.write('\t'.join(header) + '\n')
        for bc0_info in library_ID_to_bc0_info_to_counts[library_ID]:
            row = [library_ID]
            for field in total_counts_fields:
                row.append(str(library_ID_to_bc0_info_to_counts[library_ID][bc0_info].get(field, 0)))
            for field in per_run_counts_fields:
                row.append(str(library_ID_to_bc0_info_to_counts[library_ID][bc0_info].get(field, 0)))
            row += [str(e) for e in bc0_info]
            _ = counts_outfile.write('\t'.join(row) + '\n')
    elapsed = timeit.default_timer() - start_time
    print('Writing bc0_to_counts to file took ' + str(elapsed) + ' seconds to process\n')

logfile.close()



