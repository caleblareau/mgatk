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
min_barcodes = int(sys.argv[4])
mtchr = sys.argv[5]
quant_file = sys.argv[6]

base=os.path.basename(bamfile)
basename=os.path.splitext(base)[0]

def getBarcode(intags):
	'''
	Parse out the barcode per-read
	'''
	for tg in intags:
		if(barcodeTag == tg[0]):
			return(tg[1])
	return("NA")

def listBarcodes(mtchr):
	'''
	Make a giant list of observed barcodes at the mitochondrial chr
	'''
	barcodes_all = dict()
	bam = pysam.AlignmentFile(bamfile,'rb')
	Itr = bam.fetch(str(mtchr),multiple_iterators=False)
	
	for read in Itr:
		read_barcode = getBarcode(read.tags)
		barcodes_all[read_barcode] = barcodes_all.get(read_barcode, 0) + 1
	bam.close()
	return(barcodes_all)

def writePassingReads(bc, mtchr):
	'''
	Write out reads to their corresponding files based on a barcode index
	'''
	bam = pysam.AlignmentFile(bamfile,'rb')
	Itr = bam.fetch(str(mtchr),multiple_iterators=False)
	for read in Itr:
		read_barcode = getBarcode(read.tags)
		
		# If read barcode is in whitelist, then write it out
		if read_barcode in bc:
			idx = bc.index(read_barcode)
			file = fopen[idx]
			file.write(read)

# Quant barcodes and write it out
barcodes = listBarcodes(mtchr)
barcodes = {x : barcodes[x] for x in barcodes if barcodes[x] >= min_barcodes and x != "NA"}
bc = list(barcodes.keys())	
	
bcfile = open(quant_file, "w") 
for k, v in barcodes.items():
	bcfile.write(k +","+ str(v)+"\n")
bcfile.close() 

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
with multi_file_manager(bambcfiles) as fopen:
	writePassingReads(bc, mtchr)



