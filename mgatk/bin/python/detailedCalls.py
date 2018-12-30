#!/usr/bin/python

###############################################################
# Given a bam and the fasta, write a file of variant attributes
###############################################################

import sys
import re
import pysam
import gzip

bamfile = sys.argv[1]
outfile = sys.argv[2]
mitoFasta = sys.argv[3]

def parse_fasta(filename):
	f = open(filename)
	sequences = {}
	for line in f:
		if line.startswith('>'):
			name = line[1:].strip()
			sequences[name] = ''
		else:
			sequences[name] = sequences[name] + line.strip()
	f.close()
	return sequences
	
fasta = parse_fasta(mitoFasta)
mito_genome, mito_seq = list(fasta.items())[0]
bam = pysam.AlignmentFile(bamfile, "rb")

out = open(outfile, 'w')
out.write(("\t".join(["ReadNumber", "PositionOnRead", "PositionOnChr", "ObservedAllele", "Quality"]) + "\n"))
line = 1
for entry in bam:
	
	seq = entry.query_sequence
	startbp = entry.reference_start
	n = len(seq)
	ref = mito_seq[startbp:(startbp+n)]
	ns = [i for i in range(n) if seq[i] != ref[i]]
	if(len(ns) < 5):
		for i in ns:
			out.write(("\t".join([str(line), str(i), str(i+startbp), seq[i], str(entry.query_qualities[i])]) + "\n"))
	line += 1
out.close()