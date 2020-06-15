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

mito_genome = config["mito_chr"]
mito_length = str(config["mito_length"])
fasta_file = config["fasta_file"]

remove_duplicates = config["remove_duplicates"]
umi_barcode = config["umi_barcode"]
emit_base_qualities = config["emit_base_qualities"]

proper_paired = config["proper_paired"]
base_qual = str(config["base_qual"])
alignment_quality = config["alignment_quality"]
NHmax = config["NHmax"]
NMmax = config["NMmax"]

max_javamem  = config["max_javamem"]

# Software paths
java = "java"
python = "python"

# Script locations
filtclip_py = script_dir + "/bin/python/filterClipBam.py"
detailedcall_py = script_dir + "/bin/python/detailedCalls.py"
sumstatsBP_py = script_dir + "/bin/python/sumstatsBP.py"
picardCall = java + " -Xmx"+max_javamem+"  -jar " + script_dir + "/bin/picard.jar MarkDuplicates"

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

# See if we have UMIs
if(umi_barcode != "" and len(umi_barcode) == 2):
	umi_extra = " BARCODE_TAG=" + umi_barcode
else:
	umi_extra = "" 

# 3) (Optional) Remove duplicates
if (remove_duplicates == "True"):
	mdc_long = picardCall + " I="+temp_bam1+" O="+outputbam+" M="+rmlog+" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT QUIET=true VERBOSITY=ERROR USE_JDK_DEFLATER=true USE_JDK_INFLATER=true" + umi_extra 
	proc = subprocess.Popen(mdc_long, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
	out, err = proc.communicate()  # Read data from stdout and stderr
else: # just move the previous output
	os.system("mv " + temp_bam1 + " " + outputbam)
	os.system("rm " + temp_bam1 + ".bai")
pysam.index(outputbam)

# 4) Get allele counts per sample / base pair and per-base quality scores
alleleCountcall = " ".join([python, sumstatsBP_py, outputbam, prefixSM, mito_genome, mito_length, base_qual, sample, fasta_file, alignment_quality, emit_base_qualities])
os.system(alleleCountcall)

# 5) Get depth from the coverage sparse matrix
with open(prefixSM + ".coverage.txt", 'r') as coverage:
	depth = 0
	for row in coverage:
		s = row.split(",")
		depth += int(s[2].strip())
with open(outputdepth, 'w') as d:
	d.write(sample + "\t" + str(round(float(depth)/float(mito_length),2)) + "\n")

