#!/usr/bin/python

##########################
# Filter / clip .bam files
# Arguments: .bam, left clip, right clip, mito name
##########################

import sys
import pysam

bam = pysam.AlignmentFile(sys.argv[1], "rb")
out = pysam.AlignmentFile("-", "wb", template=bam)
for read in bam:
	q = read.qual
	read.seq = "N"*int(sys.argv[2]) + read.seq[int(sys.argv[2]):int(sys.argv[3])] + "N"*(abs(int(sys.argv[3])))
	read.qual = "!"*int(sys.argv[2]) + q[int(sys.argv[2]):int(sys.argv[3])] + "!"*(abs(int(sys.argv[3])))
	if(read.reference_name == sys.argv[4]):
		out.write(read)

