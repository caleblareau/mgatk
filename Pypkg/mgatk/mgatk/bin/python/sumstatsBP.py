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
base_qual = float(sys.argv[5])
sample = sys.argv[6]
fasta_file = sys.argv[7]
alignment_quality = float(sys.argv[8])

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
	align_qual_read = read.mapping_quality
	for qpos, refpos in read.get_aligned_pairs(True):
		if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
			if(seq[qpos] == "A" and quality[qpos] > base_qual):
				qualA[refpos] += quality[qpos]
				countsA[refpos] += 1
			elif(seq[qpos] == "C" and quality[qpos] > base_qual):
				qualC[refpos] += quality[qpos]
				countsC[refpos] += 1
			elif(seq[qpos] == "G" and quality[qpos] > base_qual):
				qualG[refpos] += quality[qpos]
				countsG[refpos] += 1
			elif(seq[qpos] == "T" and quality[qpos] > base_qual):
				qualT[refpos] += quality[qpos]
				countsT[refpos] += 1
			
meanQualA = [round(x/y,1) for x, y in zip(qualA, countsA)]
meanQualC = [round(x/y,1) for x, y in zip(qualC, countsC)]
meanQualG = [round(x/y,1) for x, y in zip(qualG, countsG)]
meanQualT = [round(x/y,1) for x, y in zip(qualT, countsT)]

countsA = [ int(round(elem)) for elem in countsA ]
countsC = [ int(round(elem)) for elem in countsC ]
countsG = [ int(round(elem)) for elem in countsG ]
countsT = [ int(round(elem)) for elem in countsT ]


# Allele Counts
minBP = 0
bam = pysam.AlignmentFile(bamfile, "rb")

writeSparseMatrix2("A", countsA, meanQualA)
writeSparseMatrix2("C", countsC, meanQualC)
writeSparseMatrix2("G", countsG, meanQualG)
writeSparseMatrix2("T", countsT, meanQualT)

zipped_list = zip(list(countsA),list(countsC),list(countsG),list(countsT))
sums = [sum(item) for item in zipped_list]
writeSparseMatrix("coverage", sums)
