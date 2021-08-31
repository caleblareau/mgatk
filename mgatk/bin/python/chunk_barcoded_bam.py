#!/usr/bin/python

import sys
import re
import pysam
import os
from collections import Counter

bamfile = sys.argv[1]
outfolder = sys.argv[2]
barcodeTag = sys.argv[3]
bcfile = sys.argv[4]
mtchr = sys.argv[5]
umitag = sys.argv[6]

base=os.path.basename(bamfile)
basename=os.path.basename(os.path.splitext(bcfile)[0])

def getBarcode(intags, tag_get):
	'''
	Parse out the barcode per-read
	'''
	for tg in intags:
		if(tag_get == tg[0]):
			return(tg[1])
	return("AA")

# Read in the barcodes
with open(bcfile) as barcode_file_handle:
    content = barcode_file_handle.readlines()
bc = [x.strip() for x in content] 

bam = pysam.AlignmentFile(bamfile, "rb")
outname = outfolder + "/" + basename + ".bam"
out = pysam.AlignmentFile(outname, "wb", template = bam)

# Make a DNA inspired additional barcode to account for potentially different sample indices
bases = "ACGT"
fauxdon = [a + b + c + d for a in bases for b in bases for c in bases for d in bases]

# Filter for reads that match the set of possible barcodes for this sample
try:
	Itr = bam.fetch(str(mtchr),multiple_iterators=False)
	for read in Itr:
		barcode_id = getBarcode(read.tags, barcodeTag)
		
		if(barcode_id in bc):
		
			# Now check for true UMI
			if(umitag != "XX"): 
				umi_id = getBarcode(read.tags, umitag)
			else:
				umi_id = ""
			
			# Make a fake UMI from 1) cell barcode + 2) captured umi + 3) experiment
			# all with just ACGTs so that picard doesn't bark at us. 
			# Only do this if the last string element is a number (i.e channel in 10x convention)
			if(barcode_id[-1].isnumeric()):
				split_two = barcode_id.split("-")
				faux_umi = split_two[0] + umi_id + fauxdon[(int(split_two[1]) - 1)]
			else:
				faux_umi = barcode_id + umi_id 
			read.tags = read.tags + [("MU", faux_umi)]
			out.write(read)
			

except OSError: # Truncated bam file from previous iteration handle
	print('Finished parsing bam')
	
bam.close()
out.close()

pysam.index(outname)

