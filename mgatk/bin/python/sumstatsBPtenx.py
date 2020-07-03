#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position / strand
###################################################

import sys
import re
import os
import pysam
import numpy as np

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
bam_input = pysam.AlignmentFile(bam_file, "rb")
dna_letters = ['A','C','G','T']

def getBarcode(intags):
	'''
	Parse out the barcode per-read
	'''
	for tg in intags:
		if(barcode_tag == tg[0]):
			return(tg[1])
	return("NA")


# Dimension cell x position x letter x strand
# Coverage associated with the bases
ca =  np.zeros((len(bcs),max_bp,4,2), dtype=int)

for read in bam_input:
	
	if(read.is_reverse):
		s_idx = 1
	else:
		s_idx = 0
		
	# Get read attributes 
	seq = read.seq
	quality = read.query_qualities
	align_qual_read = read.mapping_quality
	cell_barcode = getBarcode(read.tags)
	if(cell_barcode != "NA"):
		c_idx = bcs.index(cell_barcode)
	
		for q_idx, p_idx in read.get_aligned_pairs(True):
			if q_idx is not None and p_idx is not None and align_qual_read > alignment_quality:
				if(quality[q_idx] > base_qual and seq[q_idx] in dna_letters):
					l_idx = dna_letters.index(seq[q_idx])
					ca[c_idx,p_idx,l_idx,s_idx] += 1


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




