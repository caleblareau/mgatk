import sys
import re
import pysam
import gzip

#####################################
# Summarizes the quality per position
#####################################

import sys
import re
import pysam

bamfile = sys.argv[1]
outpre = sys.argv[2]
mito_genome = sys.argv[3]
maxBP = sys.argv[4]
sample = sys.argv[5]

minBP = 0

bam = pysam.AlignmentFile(bamfile, "rb")

n = int(maxBP)

counts = [0.00000001] * n # initialize with a pseudo count to avoid dividing by zero
qualities = [0.0] * n

for read in bam:
	seq = read.seq
	quality = read.query_qualities
	for qpos, refpos in read.get_aligned_pairs(True):
		if qpos is not None and refpos is not None:
			counts[refpos] += 1
			qualities[refpos] += qpos

meanQual = [round(x/y,3) for x, y in zip(qualities, counts)]

f = open(outpre + ".quality.txt", 'w')
f.write(sample + "\n")
for b in range(0,n):
	f.write(str(meanQual[b]) + "\n")
f.close()