import click
import os
import yaml
import os.path
import sys
import shutil
import random
import string
import itertools
import time
from pkg_resources import get_distribution
from subprocess import call, check_call
from .mgatkHelp import *

# ------------------------------
# Command line arguments
# ------------------------------	

@click.command()
@click.version_option()
@click.argument('mode')
@click.option('--input', default = ".", required=True, help='input directory; assumes .bam / .bam.bai files are present')
@click.option('--output', default="mgatk_out", required=True, help='Output directory for analysis')
@click.option('--mito-genome', default = "hg19", required=True, help='mitochondrial genome configuration. Choose hg19, mm10, or a custom .yaml file (see documentation)')
@click.option('--cluster-config', default = "", required=True, help='Cluster configuration for snakemake. See snakemake documentation for more details. Accepts .yaml and .json file formats.')
@click.option('--stingy', is_flag=True, help='Space-efficient analyses; remove non-vital intermediate files.')
@click.option('--atac-single', is_flag=True, help='Default parameters for ATAC-Seq single end read analyses.')
@click.option('--atac-paired', is_flag=True, help='Default parameters for ATAC-Seq paired end read analyses.')
@click.option('--rna-single', is_flag=True, help='Default parameters for RNA-Seq single end read analyses.')
@click.option('--rna-paired', is_flag=True, help='Default parameters for RNA-Seq paired end read analyses.')

@click.option('--filter-duplicates', is_flag=True, help='Filter marked (presumably PCR) duplicates')


@click.option('--keep-samples', default="ALL", help='Comma separated list of sample names to keep; ALL (special string) by default. Sample refers to basename of .bam file')
@click.option('--ignore-samples', default="NONE", help='Comma separated list of sample names to ignore; NONE (special string) by default. Sample refers to basename of .bam file')

@click.option('--skip-rds', is_flag=True, help='Generate plain-text only output. Otherwise, this generates a .rds obejct that can be immediately read into R')

def main(mode, input, output, mito_genome, cluster_config, stingy, atac_single, atac_paired, rna_single, rna_paired, keep_samples, ignore_samples, skip_rds):
	"""mgatk: Processing mitochondrial mutations."""
	__version__ = get_distribution('mgatk').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "mgatk v%s" % __version__)

	# -------------------------------
	# Verify dependencies
	# -------------------------------
	
	check_software_exists("R")
	check_software_exists("bcftools")
	check_software_exists("tabix")
	check_software_exists("samtools")
	check_software_exists("java")
	check_R_packages(['mgatk', 'ggplot2'])

	# -------------------------------
	# Determine samples for analysis
	# -------------------------------

	bams = os.popen('ls ' + input + '/*.bam').read().strip().split("\n")
	print(bams)
	
	def intersect(a, b):
		return list(set(a) & set(b))

	samples = bams
	
	if(keep_samples != "ALL"):
		keeplist = keep_samples.split(",")
		click.echo(gettime() + "Intersecting detected samples with user-retained ones: " + keep_samples)
		samples = intersect(samples, keeplist)
		
	if(ignore_samples != "NONE"):
		igslist = ignore_samples.split(",")
		for byesample in igslist:
			click.echo(gettime() + "Attempting to remove " + byesample + " from processing", logf)
			if byesample in samples: samples.remove(byesample)
	
	if not len(samples) > 0:
		sys.exit('ERROR: Could not import any samples from the user specification; check flags, logs and input configuration; QUITTING')


	# -------------------------------
	# Setup output folder
	# -------------------------------
	
	outfolder = output
	logfolder = outfolder + "/logs"
	internfolder = outfolder + "/.internal"
	parselfolder = internfolder + "/parseltongue"
	
	# Check if output directories exist; make if not
	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
	if not os.path.exists(logfolder):
		os.makedirs(logfolder)	
	if not os.path.exists(internfolder):
		os.makedirs(internfolder)
		with open(internfolder + "/README" , 'w') as outfile:
			outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")	
	if not os.path.exists(parselfolder):
		os.makedirs(parselfolder)
		with open(parselfolder + "/README" , 'w') as outfile:
			outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")	
	cwd = os.getcwd()
	logf = open(logfolder + "/base.mgatk.log", 'a')
		

	# -----------------------------------
	# Parse user-specified parameteres
	# -----------------------------------
	if(keep_indels):
		skip_indels = ""
	else:
		skip_indels = "--skip-indels "
	
	# -------------------
	# Process each sample
	# -------------------
	trimfolder = outfolder + "/01_raw"
	if not os.path.exists(trimfolder + "_process"):
		os.makedirs(trimfolder+ "_process")
				
	click.echo(gettime() + "First pass of samples", logf)
	
	snakedict1 = {'allsamples' : 'hey', 'input_directory' : input, 'output_directory' : output,
		'mitoQual' : mitoQual, 'skip_indels' : skip_indels}
	
	y1 = parselfolder + "/snake.scatter.yaml"
	with open(y1, 'w') as yaml_file:
		yaml.dump(snakedict1, yaml_file, default_flow_style=False)
		
	#snakecall1 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.Scatter --config cfp="' + y1 + '"'
	#os.system(snakecall1)
	click.echo(gettime() + "Sample scattering done.", logf)
	
	# Suspend logging
	logf.close()
	
	
