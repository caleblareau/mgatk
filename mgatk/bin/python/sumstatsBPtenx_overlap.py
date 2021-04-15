#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position / strand
###################################################

import sys
import re
import os
import pysam
import numpy as np
from collections import defaultdict

bam_file = sys.argv[1]
barcodes_file = sys.argv[2]
out_pre = sys.argv[3]
max_bp = int(sys.argv[4])
base_qual = float(sys.argv[5])
fasta_file = sys.argv[6]
alignment_quality = float(sys.argv[7])
barcode_tag = sys.argv[8]

# Import barcodes
with open(barcodes_file) as barcode_file_handle:
    content = barcode_file_handle.readlines()
bcs = [x.strip() for x in content]
dna_letters = ['A','C','G','T']

def getBarcode(intags):
	'''
	Parse out the barcode per-read
	'''
	for tg in intags:
		if(barcode_tag == tg[0]):
			return(tg[1])
	return("NA")

# organize reads into a dict where key is [cellbc]$[readname]
bam_input = [x for x in pysam.AlignmentFile(bam_file, "rb")]
ordered_bam_input = defaultdict(list)
for read in bam_input:
	cell_barcode = getBarcode(read.tags)
	if cell_barcode != "NA":
		read_key = '$'.join([cell_barcode, read.query_name])
		ordered_bam_input[read_key].append(read)

# Dimension cell x position x letter x strand
# Coverage associated with the bases
ca =  np.zeros((len(bcs), max_bp, 4, 2), dtype=int)
for read_key in ordered_bam_input:
	# disregard singlets and multiplets
	if len(ordered_bam_input[read_key]) != 2:
		continue
	
	# identify fwd and rev in a pair
	read0, read1 = ordered_bam_input[read_key]
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

	cell_barcode = read_key.split('$')[0]
	c_idx = bcs.index(cell_barcode)

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
		if (refpos is not None) and (fwd_quality[qpos] > base_qual) and (fwd_seq[qpos] in dna_letters):
			l_idx = dna_letters.index(fwd_seq[qpos])
			ca[c_idx, refpos, l_idx, 0] += 1

	# handle rev region
	rev_aligned_pairs = rev_read.get_aligned_pairs(True)
	rev_region = [pair for pair in rev_aligned_pairs if pair[0] in rev_use_idx]
	for qpos, refpos in rev_region:
		if (refpos is not None) and (rev_quality[qpos] > base_qual) and (rev_seq[qpos] in dna_letters):
			l_idx = dna_letters.index(rev_seq[qpos])
			ca[c_idx, refpos, l_idx, 1] += 1


# Function to write the slice of the matrix that is associated with the 
def writeSparseMatrixLetter(letter, letter_idx):
	out_file_fn = out_pre + "."+letter+".txt"
	with open(out_file_fn,"w") as file_handle_fn:
		for cell_idx, cell_name in enumerate(bcs):
			
			# Pull out the stranded counts
			fw_vec = ca[cell_idx,:,letter_idx,0].ravel()
			rev_vec = ca[cell_idx,:,letter_idx,1].ravel()
			
			# Write each position
			for i in range(0,int(max_bp)):
				if(fw_vec[i] > 0 or rev_vec[i] > 0):
					file_handle_fn.write(str(i+1)+","+cell_name+","+str(fw_vec[i])+","+str(rev_vec[i])+"\n")

writeSparseMatrixLetter("A", 0)
writeSparseMatrixLetter("C", 1)
writeSparseMatrixLetter("G", 2)
writeSparseMatrixLetter("T", 3)

# Export the per-base coverage for the thrill of it and the depth
out_file_depth = out_pre.replace("/temp/sparse_matrices/", "/qc/depth/") + ".depth.txt"
out_file_coverage= out_pre + ".coverage.txt"
with open(out_file_coverage,"w") as file_handle_cov:
	with open(out_file_depth,"w") as file_handle_depth:
	
		# Loop over cells
		for cell_idx, cell_name in enumerate(bcs):
		
			# Pull out the summed counts per cell per position
			cov_vec = np.sum(ca[cell_idx,:,:,:], axis = (1,2)).tolist()
			depth = round(sum(cov_vec)/len(cov_vec),2)
			# Write each position
			for i in range(0,int(max_bp)):
				if(cov_vec[i] > 0):
					file_handle_cov.write(str(i+1)+","+cell_name+","+str(cov_vec[i])+"\n")
			
			# Now write the depth
			file_handle_depth.write(cell_name+"\t"+str(depth)+"\n")




