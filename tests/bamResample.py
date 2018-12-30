#!/usr/bin/env python

import sys
import re
import pysam
import random
import os
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files and generate the XB (Cell Barcode) tag"
opts = OptionParser(usage=usage)
opts.add_option("--bam1", help="Filename of the first bam to be input")
opts.add_option("--bam2", help="Filename of the second bam to be input")
opts.add_option("--reads", help="Filename of the second bam to be input")
opts.add_option("--prop", help="Proportion of reads (between 0 and 1) to come from bam 1")

options, arguments = opts.parse_args()

bam1 = options.bam1
bam2 = options.bam2
nreads = float(options.reads)
prop = float(options.prop)


def bam_read_count_mito(bamfile):
	""" Return a tuple of the number of mapped and unmapped reads in a bam file """
	p = pysam.idxstats(bamfile)
	mapped = 0
	unmapped = 0
	for line in p.split("\n"):
		rname, rlen, nm, nu = line.rstrip().split("\t")
		if(rname == "chrM"):
			mapped += int(nm)
			unmapped += int(nu)
			return (mapped)

# Inter # of reads in bam file and how many are in the target
bam1counts = bam_read_count_mito(bam1)
bam2counts = bam_read_count_mito(bam2)

target1 = prop * nreads
target2 = (1-prop) * nreads

# Verify number of reads looks good
if(target1 > bam1counts):
	sys.exit("Not enough reads in bam 1 given configuration")
if(target2 > bam2counts):
	sys.exit("Not enough reads in bam 2 given configuration")

prop1hit = target1/bam1counts
prop2hit = target2/bam2counts

# Setup output file
mixy = "mix_" + str(prop) + "_" + os.path.basename(os.path.splitext(bam1)[0]) + "_" + str(1-prop) + "_" + os.path.basename(os.path.splitext(bam2)[0])
new_bam_name = mixy + ".bam"
coverageRatio = mixy + ".coverage.csv"

bam1io = pysam.AlignmentFile(bam1, "rb")
bam2io = pysam.AlignmentFile(bam2, "rb")

out = pysam.AlignmentFile("temp_" + new_bam_name, "wb", template = bam1io)

n = 16751
maxBP = 16751
countsA1 = [0] * n 
countsC1 = [0] * n 
countsG1 = [0] * n 
countsT1 = [0] * n 

countsA2 = [0] * n 
countsC2 = [0] * n 
countsG2 = [0] * n 
countsT2 = [0] * n 

# Parse bam data and write
n1 = 0
for read in bam1io:
	if(random.random() < prop1hit):
		seq = read.seq
		out.write(read)
		n1 = n1 + 1
		for qpos, refpos in read.get_aligned_pairs(True):
			if qpos is not None and refpos is not None:
				if(seq[qpos] == "A"):
					countsA1[refpos] += 1
				elif(seq[qpos] == "C"):
					countsC1[refpos] += 1
				elif(seq[qpos] == "G"):
					countsG1[refpos] += 1
				elif(seq[qpos] == "T"):
					countsT1[refpos] += 1
bam1io.close()

n2 = 0
for read in bam2io:
	if(random.random() < prop2hit):
		seq = read.seq
		out.write(read)
		n2 = n2 + 1
		for qpos, refpos in read.get_aligned_pairs(True):
			if qpos is not None and refpos is not None:
				if(seq[qpos] == "A"):
					countsA2[refpos] += 1
				elif(seq[qpos] == "C"):
					countsC2[refpos] += 1
				elif(seq[qpos] == "G"):
					countsG2[refpos] += 1
				elif(seq[qpos] == "T"):
					countsT2[refpos] += 1
bam2io.close()
out.close()

# Output the coverages
def writeSparseMatrixN(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11):
	with open(coverageRatio,"w") as V:
		V.write("BP,"+str(os.path.basename(os.path.splitext(bam1)[0]))+","+str(os.path.basename(os.path.splitext(bam2)[0])) + ",ratio,A1,C1,G1,T1,A2,C2,G2,T2\n")
		for i in range(0,int(maxBP)-1):
			V.write(str(i+1)+","+str(v1[i])+","+str(v2[i])+","+str(v3[i])+","+str(v4[i])+","+str(v5[i])+","+str(v6[i])+","+str(v7[i])+","+str(v8[i])+","+str(v9[i])+","+str(v10[i])+","+str(v11[i])+"\n")

zipped_list1 = zip(list(countsA1),list(countsC1),list(countsG1),list(countsT1))
sums1 = [sum(item) for item in zipped_list1]
zipped_list2 = zip(list(countsA2),list(countsC2),list(countsG2),list(countsT2))
sums2 = [sum(item) for item in zipped_list2]
ratio = [round((x/( x+ y + 0.00001)), 3) for x, y in zip(sums1, sums2)]
writeSparseMatrixN(sums1, sums2, ratio, countsA1, countsC1, countsG1, countsT1, countsA2, countsC2, countsG2, countsT2)

# Cleanup
pysam.sort("-o", new_bam_name, "temp_" + new_bam_name)
pysam.index(new_bam_name)
os.remove("temp_" + new_bam_name)


