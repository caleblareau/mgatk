#!/usr/bin/python

##########################
# Filter / clip .bam files
##########################

import sys
import re
import pysam

bamfile = sys.argv[1]
logfile = sys.argv[2]
clipL = sys.argv[3]
clipR = sys.argv[4]
mtchr = sys.argv[5]
NHmax = sys.argv[6]
NMmax = sys.argv[7]

# https://github.com/pysam-developers/pysam/issues/509
bam = pysam.AlignmentFile(bamfile, "rb")
out = pysam.AlignmentFile("-", "wb", template = bam)

keepCount = 0
filtCount = 0

def filterReadTags(intags):
    '''Checks for aligner-specific read tags and filters'''

    for tg in intags:
    	if(('NH' == tg[0] and int(tg[1]) > int(NHmax)) or \
    		(('NM' == tg[0] or 'nM' == tg[0]) and int(tg[1]) > int(NMmax))):
        		return(False)
    return(True)

def processRead(read):
	global keepCount
	global filtCount
	
	if(filterReadTags(read.tags)):
		keepCount += 1
		out.write(read)
	else:
		filtCount += 1

# Clip Both
if(int(clipL) > 0 and int(clipR) < 0):
	for read in bam:
		q = read.qual
		read.seq = "N"*int(clipL) + read.seq[int(clipL):int(clipR )] + "N"*(abs(int(clipR)))
		read.qual = "!"*int(clipL) + q[int(clipL):int(clipR )] + "!"*(abs(int(clipR)))
		processRead(read)
			
# Clip right
elif(int(clipR) < 0):
	for read in bam:
		q = read.qual
		read.seq = read.seq[:int(clipR)] + "N"*(abs(int(clipR)))
		read.qual = q[0:int(clipR)] + "!"*(abs(int(clipR)))
		processRead(read)
			
# Clip left
if(int(clipL) > 0):
	for read in bam:
		q = read.qual
		read.seq = "N"*int(clipL) + read.seq[int(clipL):]
		read.qual = "!"*int(clipL) + q[int(clipL):]
		processRead(read)
			
# No clipping, just filtering
else:
	for read in bam:
		processRead(read)

with open(logfile , 'w') as outfile:
	outfile.write("Kept "+ str(keepCount) + "\n" + "Removed " + str(filtCount)+ "\n")

