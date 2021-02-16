# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 22:05:52 2019

@author: kevinroy
"""
import operator
from sys import argv

fastq_infilename = argv[1]
fastq_outfile = open(argv[2], 'w')
min_counts = int(argv[3])

seq_to_qual = {}
seq_to_counts = {}

with open(fastq_infilename, 'r') as infile:
    line_in_read = 0
    for line in infile:
        line_in_read += 1
        if line_in_read == 2:
            seq = line.strip()
        elif line_in_read == 4:
            qual = line.strip()
            line_in_read = 0
            if seq not in seq_to_counts:
                seq_to_counts[seq] = 0
            seq_to_counts[seq] += 1
            seq_to_qual[seq] = qual
infile.close()

sorted_seq = sorted(seq_to_counts.items(), key=operator.itemgetter(1), reverse=True)    

for seq_counts in sorted_seq:
    seq, counts = seq_counts
    if counts >= min_counts:
        	read_name = '@' + str(counts)
        	qual = seq_to_qual[seq]
        	fastq_outfile.write(read_name + '\n' + seq + '\n+\n' + qual + '\n')
fastq_outfile.close()   