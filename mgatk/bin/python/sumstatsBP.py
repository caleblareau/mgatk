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
emit_base_qualities = sys.argv[9]

# Export Functions
def writeSparseMatrix(mid, vec):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec[i])+"\n")


def writeSparseMatrix2(mid, vec1, vec2):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec1[i] > 0 or vec2[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+"\n")

def writeSparseMatrix4(mid, vec1, vec2, vec3, vec4):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec1[i] > 0 or vec3[i] > 0):
				V.write(str(i+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+","+str(vec3[i])+","+str(vec4[i])+"\n")


n = int(maxBP)

# initialize with a pseudo count to avoid dividing by zero
countsA_fw = [0.00000001] * n 
countsC_fw = [0.00000001] * n 
countsG_fw = [0.00000001] * n 
countsT_fw = [0.00000001] * n 

qualA_fw = [0.0] * n
qualC_fw = [0.0] * n
qualG_fw = [0.0] * n
qualT_fw = [0.0] * n

countsA_rev = [0.00000001] * n 
countsC_rev = [0.00000001] * n 
countsG_rev = [0.00000001] * n 
countsT_rev = [0.00000001] * n 

qualA_rev = [0.0] * n
qualC_rev = [0.0] * n
qualG_rev = [0.0] * n
qualT_rev = [0.0] * n

bam2 = pysam.AlignmentFile(bamfile, "rb")
for read in bam2:
	seq = read.seq
	reverse = read.is_reverse
	quality = read.query_qualities
	align_qual_read = read.mapping_quality
	for qpos, refpos in read.get_aligned_pairs(True):
		if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
			if(seq[qpos] == "A" and quality[qpos] > base_qual):
				if(not reverse):
					qualA_fw[refpos] += quality[qpos]
					countsA_fw[refpos] += 1
				else:
					qualA_rev[refpos] += quality[qpos]
					countsA_rev[refpos] += 1
			elif(seq[qpos] == "C" and quality[qpos] > base_qual):
				if(not reverse):
					qualC_fw[refpos] += quality[qpos]
					countsC_fw[refpos] += 1
				else:
					qualC_rev[refpos] += quality[qpos]
					countsC_rev[refpos] += 1
			elif(seq[qpos] == "G" and quality[qpos] > base_qual):
				if(not reverse):
					qualG_fw[refpos] += quality[qpos]
					countsG_fw[refpos] += 1
				else:
					qualG_rev[refpos] += quality[qpos]
					countsG_rev[refpos] += 1
			elif(seq[qpos] == "T" and quality[qpos] > base_qual):
				if(not reverse):
					qualT_fw[refpos] += quality[qpos]
					countsT_fw[refpos] += 1
				else:
					qualT_rev[refpos] += quality[qpos]
					countsT_rev[refpos] += 1
			
meanQualA_fw = [round(x/y,1) for x, y in zip(qualA_fw, countsA_fw)]
meanQualC_fw = [round(x/y,1) for x, y in zip(qualC_fw, countsC_fw)]
meanQualG_fw = [round(x/y,1) for x, y in zip(qualG_fw, countsG_fw)]
meanQualT_fw = [round(x/y,1) for x, y in zip(qualT_fw, countsT_fw)]

countsA_fw = [ int(round(elem)) for elem in countsA_fw ]
countsC_fw = [ int(round(elem)) for elem in countsC_fw ]
countsG_fw = [ int(round(elem)) for elem in countsG_fw ]
countsT_fw = [ int(round(elem)) for elem in countsT_fw ]

meanQualA_rev = [round(x/y,1) for x, y in zip(qualA_rev, countsA_rev)]
meanQualC_rev = [round(x/y,1) for x, y in zip(qualC_rev, countsC_rev)]
meanQualG_rev = [round(x/y,1) for x, y in zip(qualG_rev, countsG_rev)]
meanQualT_rev = [round(x/y,1) for x, y in zip(qualT_rev, countsT_rev)]

countsA_rev = [ int(round(elem)) for elem in countsA_rev ]
countsC_rev = [ int(round(elem)) for elem in countsC_rev ]
countsG_rev = [ int(round(elem)) for elem in countsG_rev ]
countsT_rev = [ int(round(elem)) for elem in countsT_rev ]

# Allele Counts
bam = pysam.AlignmentFile(bamfile, "rb")

if(emit_base_qualities == "True"):
	writeSparseMatrix4("A", countsA_fw, meanQualA_fw, countsA_rev, meanQualA_rev)
	writeSparseMatrix4("C", countsC_fw, meanQualC_fw, countsC_rev, meanQualC_rev)
	writeSparseMatrix4("G", countsG_fw, meanQualG_fw, countsG_rev, meanQualG_rev)
	writeSparseMatrix4("T", countsT_fw, meanQualT_fw, countsT_rev, meanQualT_rev)
else:
	writeSparseMatrix2("A", countsA_fw, countsA_rev)
	writeSparseMatrix2("C", countsC_fw, countsC_rev)
	writeSparseMatrix2("G", countsG_fw, countsG_rev)
	writeSparseMatrix2("T", countsT_fw, countsT_rev)

zipped_list = zip(list(countsA_fw),list(countsC_fw),list(countsG_fw),list(countsT_fw), list(countsA_rev),list(countsC_rev),list(countsG_rev),list(countsT_rev))
sums = [sum(item) for item in zipped_list]
writeSparseMatrix("coverage", sums)
