#!/usr/bin/python

###############################################
# Filter reads / clip .bam files
# Don't print() anything!!!! writing to STDOUT
##############################################

import sys
import re
import pysam

bamfile = sys.argv[1]
logfile = sys.argv[2]
mtchr = sys.argv[3]
proper_pair = sys.argv[4]
NHmax = sys.argv[5]
NMmax = sys.argv[6]

# https://github.com/pysam-developers/pysam/issues/509
bam = pysam.AlignmentFile(bamfile, "rb")
out = pysam.AlignmentFile("-", "wb", template = bam)

keepCount = 0
filtCount = 0

def filterReadTags(intags):
    '''
    Checks for aligner-specific read tags and filters
	'''
	
    for tg in intags:
    	if(('NH' == tg[0] and int(tg[1]) > int(NHmax)) or \
    		(('NM' == tg[0] or 'nM' == tg[0]) and int(tg[1]) > int(NMmax))):
        		return(False)
    return(True)

def pairing(read):
	'''
	Check if read is paired, properly paired, etc.
	'''
	
	if(proper_pair != "True"): # then user doesn't care to filter it
		return(True)
	else:
		return(read.is_proper_pair)

def processRead(read):
	global keepCount
	global filtCount
	if(filterReadTags(read.tags) and read.reference_name == mtchr and pairing(read)):
		keepCount += 1
		out.write(read)
	else:
		filtCount += 1

for read in bam:
	processRead(read)

with open(logfile , 'w') as outfile:
	outfile.write("Kept "+ str(keepCount) + "\n" + "Removed " + str(filtCount)+ "\n")

