#!/usr/bin/python

###############################################
# Summarizes the total number of reads per pos
##############################################

import sys
import re
import pysam

bamfile = sys.argv[1]
outpre = sys.argv[2]
mito_genome = sys.argv[3]
maxBP = sys.argv[4]
base_qual = sys.argv[5]
sample = sys.argv[6]

minBP = 0
bam = pysam.AlignmentFile(bamfile, "rb")
cc = bam.count_coverage(mito_genome, minBP, int(maxBP),quality_threshold=int(base_qual))
zipped_list = zip(list(cc[0]),list(cc[1]),list(cc[2]),list(cc[3]))
sums = [sum(item) for item in zipped_list]

def writeSparseMatrix(mid, vec):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(1,int(maxBP)-1):
			if(vec[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec[i])+"\n")
			
writeSparseMatrix("coverage", sums)
writeSparseMatrix("A", cc[0])
writeSparseMatrix("C", cc[1])
writeSparseMatrix("G", cc[2])
writeSparseMatrix("T", cc[3])
