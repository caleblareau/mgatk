#!/usr/bin/python

from os.path import join
import os
import subprocess
import sys
import shutil
import pysam
from ruamel import yaml

configFile = sys.argv[1]
inputbam = sys.argv[2]
outputbam = sys.argv[3]
sample = sys.argv[4]

with open(configFile, 'r') as stream:
	config = yaml.load(stream, Loader=yaml.Loader)

# Parse the configuration variables
indir = config["input_directory"]
outdir = config["output_directory"]
script_dir = config["script_dir"]

mito_genome = config["mito_genome"]
mito_length = str(config["mito_length"])
fasta_file = config["fasta_file"]

remove_duplicates = config["remove_duplicates"]
do_baq = config["baq"]
proper_paired = config["proper_paired"]

base_qual = str(config["base_qual"])
alignment_quality = config["alignment_quality"]
max_javamem  = config["max_javamem"]

NHmax = config["NHmax"]
NMmax = config["NMmax"]

# Software paths
java = "java"
python = "python"

# Script locations
filtclip_py = script_dir + "/bin/python/filterClipBam.py"
detailedcall_py = script_dir + "/bin/python/detailedCalls.py"
sumstatsBP_py = script_dir + "/bin/python/sumstatsBP.py"
MarkDuplicatesCall = java + " -Xmx"+max_javamem+"  -jar " + script_dir + "/bin/MarkDuplicates.jar"

# Prepare filepath locations
rmlog = outputbam.replace(".qc.bam", ".rmdups.log").replace("/temp/ready_bam/", "/logs/rmdupslogs/")
filtlog = outputbam.replace(".qc.bam", ".filter.log").replace("/temp/ready_bam/", "/logs/filterlogs/")
temp_bam0 = outputbam.replace(".qc.bam", ".temp0.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")
temp_bam1 = outputbam.replace(".qc.bam", ".temp1.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")
prefixSM = outdir + "/temp/sparse_matrices/" + sample
outputdepth = outdir + "/qc/depth/" + sample + ".depth.txt"

# 1) Filter bam files
pycall = " ".join([python, filtclip_py, inputbam, filtlog, mito_genome, proper_paired, NHmax, NMmax]) + " > " + temp_bam0
os.system(pycall)

# 2) Sort the filtered bam file
pysam.sort("-o", temp_bam1, temp_bam0)
pysam.index(temp_bam1)

# 2a) BAQ recalibration
if(do_baq == "True"):
	baq_bam = outputbam.replace(".qc.bam", ".baq.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")
	os.system("samtools calmd -bAr " + temp_bam1 + " " + fasta_file + " > " + baq_bam)
	os.system("mv " + baq_bam + " " + temp_bam1)
	os.system("rm " + temp_bam1 + ".bai")
	pysam.index(temp_bam1)

# 3) (Optional) Remove duplicates
if (remove_duplicates == "True"):
	mdc_long = MarkDuplicatesCall + " INPUT="+temp_bam1+" OUTPUT="+outputbam+" METRICS_FILE="+rmlog+" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT"
	os.system(mdc_long)
else: # just move the previous output
	os.system("mv " + temp_bam1 + " " + outputbam)
	os.system("rm " + temp_bam1 + ".bai")
pysam.index(outputbam)

# 4) Get allele counts per sample / base pair and per-base quality scores
alleleCountcall = " ".join([python, sumstatsBP_py, outputbam, prefixSM, mito_genome, mito_length, base_qual, sample, fasta_file, alignment_quality])
os.system(alleleCountcall)

# 5) Get depth from the coverage sparse matrix
with open(prefixSM + ".coverage.txt", 'r') as coverage:
	depth = 0
	for row in coverage:
		s = row.split(",")
		depth += int(s[2].strip())
with open(outputdepth, 'w') as d:
	d.write(sample + "\t" + str(round(float(depth)/float(mito_length),2)) + "\n")

