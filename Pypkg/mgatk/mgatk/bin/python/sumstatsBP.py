#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position
###################################################

import sys
import re
import os
import pysam

bamfile = sys.argv[1]
outpre = sys.argv[2]
mito_genome = sys.argv[3]
maxBP = sys.argv[4]
base_qual = sys.argv[5]
sample = sys.argv[6]
fasta_file = sys.argv[7]

# Export Functions
def writeSparseMatrix(mid, vec):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(1,int(maxBP)-1):
			if(vec[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec[i])+"\n")


def writeSparseMatrix2(mid, vec1, vec2):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(1,int(maxBP)-1):
			if(vec1[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+"\n")

n = int(maxBP)
# BAQ
# initialize with a pseudo count to avoid dividing by zero
countsA = [0.00000001] * n 
countsC = [0.00000001] * n 
countsG = [0.00000001] * n 
countsT = [0.00000001] * n 

qualA = [0.0] * n
qualC = [0.0] * n
qualG = [0.0] * n
qualT = [0.0] * n

bam2 = pysam.AlignmentFile(bamfile, "rb")
for read in bam2:
	seq = read.seq
	quality = read.query_qualities
	for qpos, refpos in read.get_aligned_pairs(True):
		if qpos is not None and refpos is not None:
			if(seq[qpos] == "A"):
				qualA[refpos] += quality[qpos]
				countsA[refpos] += 1
			elif(seq[qpos] == "C"):
				qualC[refpos] += quality[qpos]
				countsC[refpos] += 1
			elif(seq[qpos] == "G"):
				qualG[refpos] += quality[qpos]
				countsG[refpos] += 1
			elif(seq[qpos] == "T"):
				qualT[refpos] += quality[qpos]
				countsT[refpos] += 1
			
meanQualA = [round(x/y,1) for x, y in zip(qualA, countsA)]
meanQualC = [round(x/y,1) for x, y in zip(qualC, countsC)]
meanQualG = [round(x/y,1) for x, y in zip(qualG, countsG)]
meanQualT = [round(x/y,1) for x, y in zip(qualT, countsT)]

# Allele Counts
minBP = 0
bam = pysam.AlignmentFile(bamfile, "rb")
cc = bam.count_coverage(mito_genome, minBP, int(maxBP),quality_threshold=int(base_qual))

writeSparseMatrix2("A", cc[0], meanQualA)
writeSparseMatrix2("C", cc[1], meanQualC)
writeSparseMatrix2("G", cc[2], meanQualG)
writeSparseMatrix2("T", cc[3], meanQualT)

zipped_list = zip(list(cc[0]),list(cc[1]),list(cc[2]),list(cc[3]))
sums = [sum(item) for item in zipped_list]
writeSparseMatrix("coverage", sums)
