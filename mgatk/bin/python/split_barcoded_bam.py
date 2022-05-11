#!/usr/bin/python

import sys
import re
import pysam
import os
from collections import Counter
from contextlib import contextmanager

bamfile = sys.argv[1]
outfolder = sys.argv[2]
barcodeTag = sys.argv[3]
bcfile = sys.argv[4]
mtchr = sys.argv[5]

base=os.path.basename(bamfile)
basename=os.path.splitext(base)[0]

def getBarcode(read):
	'''
	Parse out the barcode per-read
	'''
	if read.has_tag(barcodeTag):
		return read.get_tag(barcodeTag)
	else:
		return "NA"


def writePassingReads(bc_dict, mtchr):
	'''
	Write out reads to their corresponding files based on a barcode index
	'''
	bam = pysam.AlignmentFile(bamfile,'rb')
	Itr = bam.fetch(str(mtchr),multiple_iterators=False)
	for read in Itr:
		read_barcode = getBarcode(read)
		
		# If read barcode is in whitelist, then write it out
		if read_barcode in bc_dict:
			idx = bc_dict[read_barcode]
			file = fopen[idx]
			file.write(read)

# Read in the barcodes
with open(bcfile) as barcode_file_handle:
    content = barcode_file_handle.readlines()
bc = [x.strip() for x in content]

# Open up a bunch of files and write out reads for valid barcodes
@contextmanager
def multi_file_manager(files, mode='rt'):
	"""
	Open multiple files and make sure they all get closed.
	"""
	temp = pysam.AlignmentFile(bamfile, "rb")
	files = [pysam.AlignmentFile(file, "wb", template = temp) for file in files]
	temp.close()
	yield files
	for file in files:
		file.close()
		
# Final loop to write out passing reads
bambcfiles = [outfolder + "/" + bc1 + ".bam" for bc1 in bc]
bc_dict = {bc1: i  for i,bc1 in enumerate(bc)}
with multi_file_manager(bambcfiles) as fopen:
	writePassingReads(bc_dict, mtchr)



