#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position
###################################################

import sys
import re
import os
import pysam
import numpy as np
from collections import defaultdict

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

# organize reads into a dict where key is readname
bam2 = [x for x in pysam.AlignmentFile(bamfile, "rb")]
ordered_bam2 = defaultdict(list)
for read in bam2:
	ordered_bam2[read.query_name].append(read)

for read_name in ordered_bam2:
	# disregard singlets and multiplets
	if len(ordered_bam2[read_name]) != 2:
		continue
	
	# identify fwd and rev in a pair
	read0, read1 = ordered_bam2[read_name]
	if read0.is_reverse and not read1.is_reverse:
		fwd_read, rev_read = read1, read0
	elif not read0.is_reverse and read1.is_reverse:
		fwd_read, rev_read = read0, read1
	else:
		# disregard a pair if both are the same strand
		continue
	
	# gather what we need
	fwd_seq, rev_seq = fwd_read.query_sequence, rev_read.query_sequence
	fwd_quality, rev_quality = np.array(fwd_read.query_qualities), np.array(rev_read.query_qualities)
	fwd_align_qual_read, rev_align_qual_read = fwd_read.mapping_quality, rev_read.mapping_quality
	
	# check alignment quality
	if fwd_align_qual_read > alignment_quality and rev_align_qual_read > alignment_quality:
		# partition the pair into fwd-only, overlap, and rev-only
		overlap_length = fwd_read.get_overlap(rev_read.reference_start, rev_read.reference_end)
		if overlap_length == 0:
			fwd_use_idx = np.arange(len(fwd_seq))
			rev_use_idx = np.arange(len(rev_seq))
		else:
			# choose which strand to use in the overlap region based on quality score
			fwd_only_end = len(fwd_seq) - overlap_length
			rev_only_start = overlap_length
			fwd_overlap_quality = fwd_quality[fwd_only_end:]
			rev_overlap_quality = rev_quality[:rev_only_start]
			fwd_overlap_use_idx = np.where(fwd_overlap_quality > rev_overlap_quality)[0]
			rev_overlap_use_idx = np.where(fwd_overlap_quality < rev_overlap_quality)[0]
			
			# evenly assign bases with equal quality in overlap region
			equal_overlap_idx = np.where(fwd_overlap_quality == rev_overlap_quality)[0]
			equal_split = int(np.floor(len(equal_overlap_idx)/2))
			fwd_overlap_use_idx = np.concatenate([fwd_overlap_use_idx, equal_overlap_idx[:equal_split]])
			rev_overlap_use_idx = np.concatenate([rev_overlap_use_idx, equal_overlap_idx[equal_split:]])

			# merge the exclusive region and use idx in overlap region
			fwd_use_idx = np.concatenate([np.arange(fwd_only_end), fwd_overlap_use_idx + fwd_only_end])
			rev_use_idx = np.concatenate([rev_overlap_use_idx, np.arange(rev_only_start, len(rev_seq))])
	
	elif fwd_align_qual_read <= alignment_quality and rev_align_qual_read <= alignment_quality:
		# use none for either
		fwd_use_idx = np.array([])
		rev_use_idx = np.array([])
	
	elif fwd_align_qual_read > alignment_quality and rev_align_qual_read <= alignment_quality:
		# use none of rev and all of fwd
		fwd_use_idx = np.arange(len(fwd_seq))
		rev_use_idx = np.array([])
	
	elif fwd_align_qual_read <= alignment_quality and rev_align_qual_read > alignment_quality:
		# use all of rev and none of fwd
		fwd_use_idx = np.array([])
		rev_use_idx = np.arange(len(rev_seq))
	
	# handle fwd region
	fwd_aligned_pairs = fwd_read.get_aligned_pairs(True)
	fwd_region = [pair for pair in fwd_aligned_pairs if pair[0] in fwd_use_idx]
	for qpos, refpos in fwd_region:
		if refpos is not None and fwd_quality[qpos] > base_qual:
			if fwd_seq[qpos] == "A":
				qualA_fw[refpos] += fwd_quality[qpos]
				countsA_fw[refpos] += 1
			elif fwd_seq[qpos] == "C":
				qualC_fw[refpos] += fwd_quality[qpos]
				countsC_fw[refpos] += 1
			elif fwd_seq[qpos] == "G":
				qualG_fw[refpos] += fwd_quality[qpos]
				countsG_fw[refpos] += 1
			elif fwd_seq[qpos] == "T":
				qualT_fw[refpos] += fwd_quality[qpos]
				countsT_fw[refpos] += 1
	
	# handle rev region
	rev_aligned_pairs = rev_read.get_aligned_pairs(True)
	rev_region = [pair for pair in rev_aligned_pairs if pair[0] in rev_use_idx]
	for qpos, refpos in rev_region:
		if refpos is not None and rev_quality[qpos] > base_qual:
			if rev_seq[qpos] == "A":
				qualA_rev[refpos] += rev_quality[qpos]
				countsA_rev[refpos] += 1
			elif rev_seq[qpos] == "C":
				qualC_rev[refpos] += rev_quality[qpos]
				countsC_rev[refpos] += 1
			elif rev_seq[qpos] == "G":
				qualG_rev[refpos] += rev_quality[qpos]
				countsG_rev[refpos] += 1
			elif rev_seq[qpos] == "T":
				qualT_rev[refpos] += rev_quality[qpos]
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