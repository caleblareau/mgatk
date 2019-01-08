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
passing_file = sys.argv[7]

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

def quantifyBarcodes(mtchr):
	'''
	Make a giant dictionary of observed barcodes at the mitochondrial chr
	'''
	barcodes_all = dict()
	bam = pysam.AlignmentFile(bamfile,'rb')
	Itr = bam.fetch(str(mtchr),multiple_iterators=False)
	
	for read in Itr:
		read_barcode = getBarcode(read.tags)
		barcodes_all[read_barcode] = barcodes_all.get(read_barcode, 0) + 1
	bam.close()
	return(barcodes_all)


# Quant barcodes and write it out
barcodes = quantifyBarcodes(mtchr)
barcodes = {x : barcodes[x] for x in barcodes if barcodes[x] >= min_barcodes and x != "NA"}
bc = list(barcodes.keys())	
	
quant_file_o = open(quant_file, "w") 
for k, v in barcodes.items():
	quant_file_o.write(k +","+ str(v)+"\n")
quant_file_o.close() 

passing_file_o = open(passing_file, "w") 
for k, v in barcodes.items():
	passing_file_o.write(k +"\n")
passing_file_o.close() 
