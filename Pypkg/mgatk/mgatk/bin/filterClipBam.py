#!/usr/bin/python

##########################
# Filter / clip .bam files
# Arguments: .bam, left clip, right clip, mito name
##########################

import sys
import pysam

bam = pysam.AlignmentFile(sys.argv[1], "rb")
out = pysam.AlignmentFile("-", "wb", template=bam)

if(int(sys.argv[2]) > 0 & int(sys.argv[3]) < 0):
	for read in bam:
		q = read.qual
		read.seq = "N"*int(sys.argv[2]) + read.seq[int(sys.argv[2]):int(sys.argv[3])] + "N"*(abs(int(sys.argv[3])))
		read.qual = "!"*int(sys.argv[2]) + q[int(sys.argv[2]):int(sys.argv[3])] + "!"*(abs(int(sys.argv[3])))
		if(read.reference_name == sys.argv[4]):
			out.write(read)
elif(int(sys.argv[3]) < 0):
	for read in bam:
		q = read.qual
		read.seq = read.seq[int(sys.argv[2]):int(sys.argv[3])] + "N"*(abs(int(sys.argv[3])))
		read.qual = q[int(sys.argv[2]):int(sys.argv[3])] + "!"*(abs(int(sys.argv[3])))
		if(read.reference_name == sys.argv[4]):
			out.write(read)
else:
		for read in bam:
			if(read.reference_name == sys.argv[4]):
				out.write(read)
				
