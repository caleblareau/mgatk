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
name = config["name"]
script_dir = config["script_dir"]

mito_genome = config["mito_genome"]
mito_length = str(config["mito_length"])
fasta_file = config["fasta_file"]

remove_duplicates = config["remove_duplicates"]
proper_paired = config["proper_paired"]
skip_indels = config["skip_indels"]

base_qual = str(config["base_qual"])
blacklist_percentile = config["blacklist_percentile"]
max_javamem  = config["max_javamem"]

clipL = config["clipl"]
clipR = config["clipr"]
NHmax = config["NHmax"]
NMmax = config["NMmax"]

detailed_calls = config["detailed_calls"]

# Software paths
java = "java"
Rscript = "Rscript"
samtools = "samtools"
python = "python"

# Script locations
filtclip_py = script_dir + "/bin/python/filterClipBam.py"
detailedcall_py = script_dir + "/bin/python/detailedCalls.py"
sumstatsBP_py = script_dir + "/bin/python/sumstatsBP.py"

depthTableQuery_R = script_dir + "/bin/R/depthTableQuery.R"
makeBlacklist_R = script_dir + "/bin/R/makeBlacklist.R"

MarkDuplicatesCall = java + " -Xmx"+max_javamem+"  -jar " + script_dir + "/bin/MarkDuplicates.jar"

# Prepare filepath locations
rmlog = outputbam.replace(".qc.bam", ".rmdups.log").replace("/temp/ready_bam/", "/logs/rmdupslogs/")
filtlog = outputbam.replace(".qc.bam", ".filter.log").replace("/temp/ready_bam/", "/logs/filterlogs/")
temp_bam0 = outputbam.replace(".qc.bam", ".temp0.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")
temp_bam1 = outputbam.replace(".qc.bam", ".temp1.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")
prefixSM = outdir + "/temp/sparse_matrices/" + sample
outputdepth = outdir + "/qc/depth/" + sample + ".depth.txt"

# 1) Filter bam files
pycall = " ".join([python, filtclip_py, inputbam, filtlog, clipL, "-"+clipR, mito_genome, NHmax, NMmax])
if(len(proper_paired) > 0):
	pycallout = pycall + proper_paried + " > " + temp_bam0 # samtools
else:
	pycallout = pycall+ " > " + temp_bam0
os.system(pycallout)

# 2) Sort the filtered bam file
pysam.sort("-o", temp_bam1, temp_bam0)
pysam.index(temp_bam1)

# 3) (Optional) Remove duplicates
if (remove_duplicates == "true"):
	mdc_long = MarkDuplicatesCall + " INPUT="+temp_bam1+" OUTPUT="+outputbam+" METRICS_FILE="+rmlog+" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT"
	os.system(mdc_long)
else: # just move the previous output
	os.system("mv " + temp_bam1 + " " + outputbam)
pysam.index(outputbam)

# 4) Get allele counts per sample / base pair and per-base quality scores
pycall = " ".join([python, sumstatsBP_py, outputbam, prefixSM, mito_genome, mito_length, base_qual, sample])
os.system(pycall)

# 5) Get depth from the coverage sparse matrix
with open(prefixSM + ".coverage.txt", 'r') as coverage:
	depth = 0
	for row in coverage:
		s = row.split(",")
		depth += int(s[2].strip())
with open(outputdepth, 'w') as d:
	d.write(sample + "\t" + str(round(float(depth)/float(mito_length),2)) + "\n")

