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

base=os.path.basename(bamfile)
basename=os.path.basename(os.path.splitext(bcfile)[0])

def getBarcode(intags):
	'''
	Parse out the barcode per-read
	'''
	for tg in intags:
		if(barcodeTag == tg[0]):
			return(tg[1])
	return("NA")

# Read in the barcodes
with open(bcfile) as barcode_file_handle:
    content = barcode_file_handle.readlines()
bc = [x.strip() for x in content] 

bam = pysam.AlignmentFile(bamfile, "rb")
outname = outfolder + "/" + basename + ".bam"
out = pysam.AlignmentFile(outname, "wb", template = bam)

# Filter for reads that match the set of possible barcodes for this sample
try:
	for read in bam:
		barcode_id = getBarcode(read.tags)
		
		if(barcode_id in bc):
			out.write(read)
			
			
except OSError: # Truncated bam file from previous iteration handle
	print('Finished parsing bam')
	
bam.close()
out.close()

pysam.index(outname)

