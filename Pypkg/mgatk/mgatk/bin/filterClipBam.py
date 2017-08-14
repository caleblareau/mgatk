#!/usr/bin/python

##########################
# Filter / clip .bam files
# Arguments: .bam, left clip, right clip, mito name
##########################

import sys
import re
import pysam

bamfile = sys.argv[1]
clipL = sys.argv[2]
clipR = sys.argv[3]
mtchr = sys.argv[4]
mtchr = sys.argv[4]

bam = pysam.AlignmentFile(bamfile, "rb")
out = pysam.AlignmentFile("-", "wb", template=bam)

# Clip Both
if(int(clipL) > 0 & int(sys.argv[3]) < 0):
	for read in bam:
		q = read.qual
		read.seq = "N"*int(clipL) + read.seq[int(clipL):int(clipR )] + "N"*(abs(int(clipR)))
		read.qual = "!"*int(clipL) + q[int(clipL):int(clipR )] + "!"*(abs(int(clipR)))
		if(read.reference_name == mtchr]):
			out.write(read)
			
# Clip right
elif(int(clipL) < 0):
	for read in bam:
		q = read.qual
		read.seq = read.seq[int(clipL):int(clipR)] + "N"*(abs(int(clipR)))
		read.qual = q[int(clipL):int(clipR)] + "!"*(abs(int(clipR)))
		if(read.reference_name == mtchr):
			out.write(read)

# No clipping, just filtering
else:
		for read in bam:
			if(read.reference_name == mtchr):
				out.write(read)
				
