#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position
###################################################

import sys
import re
import pysam

bamfile = sys.argv[1]
outpre = sys.argv[2]
mito_genome = sys.argv[3]
maxBP = sys.argv[4]
base_qual = sys.argv[5]
sample = sys.argv[6]

################
# Allele counts
################

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

################
# Quality scores
################

n = int(maxBP)

counts = [0.00000001] * n # initialize with a pseudo count to avoid dividing by zero
qualities = [0.0] * n

for read in bam:
	seq = read.seq
	quality = read.query_qualities
	for qpos, refpos in read.get_aligned_pairs(True):
		if qpos is not None and refpos is not None:
			counts[refpos] += 1
			qualities[refpos] += qpos

meanQual = [round(x/y,3) for x, y in zip(qualities, counts)]

qoutpre = outpre.replace("/temp/sparse_matrices/","/qc/quality/")
with open(qoutpre + ".quality.txt", 'w') as f:
	f.write(sample + "\n")
	for b in range(0,n):
		f.write(str(meanQual[b]) + "\n")

