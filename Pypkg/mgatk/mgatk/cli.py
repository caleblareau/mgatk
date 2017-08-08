import click
import os
import os.path
import sys
import shutil
import yaml
import random
import string
import itertools
import time
from pkg_resources import get_distribution
from subprocess import call, check_call
from .proatacHelp import *
from .proatacProjectClass import *



# ------------------------------
# Command line arguments
# ------------------------------	

@click.command()
@click.version_option()
@click.option('--check', is_flag=True, help='[MODE] Check to see if all dependencies are properly configured.')
@click.option('--stingy', is_flag=True, help='Space-efficient analyses; remove non-vital intermediate files.')
@click.argument('manifest')


def main(manifest, check, stingy):
	"""Preprocessing ATAC and scATAC Data."""
	__version__ = get_distribution('proatac').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "Starting proatac pipeline v%s" % __version__)
	ymml = parse_manifest(manifest)

	
	# -------------------------------
	# Utility functions and variables
	# -------------------------------
	
	outfolder = os.path.abspath(ymml['project_dir']) 
	logfolder = outfolder + "/logs"
	internfolder = outfolder + "/.internal"
	parselfolder = internfolder + "/parseltongue"
	
	# Check if directories exist; make if not
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
	logf = open(logfolder + "/base.proatac.log", 'a')
		
	# Main Project Variable; output to .yaml so it's clear what's going on
	p = proatacProject(ymml, script_dir)
	projectdict = p.__dict__
	del projectdict["yaml"]
	
	# Copy original yaml
	cpyaml = 'cp ' + manifest + ' ' + logfolder +  "/supplied.yaml"
	os.system(cpyaml)
	
	y0 = logfolder + "/interpreted.params.yaml"
	with open(y0, 'w') as yaml_file:
		yaml.dump(projectdict, yaml_file, default_flow_style=False)
	
	# -------------------------------
	# Atypical analysis modes
	# -------------------------------	
	
	if check:
		sys.exit("Success! We're reasonably confident that all dependencies and files are good to go!")

	# -----------------------------------
	# ACTUAL PREPROCESSING / SNAKE MAKING
	# -----------------------------------
	
	click.echo(gettime() + "Project .yaml successfully loaded. ", logf)
	
	# -------------
	# Adapter Trim
	# -------------
	trimfolder = outfolder + "/01_trimmed"
	if not os.path.exists(trimfolder + "_reads"):
		os.makedirs(trimfolder+ "_reads")
				
	click.echo(gettime() + "Trimming samples", logf)
	
	snakedict1 = {'allsamples' : parselfolder + "/allsamples.csv", 'outdir' : outfolder,
		'scriptdir' : script_dir, 'PEAT' : p.peat_path, 'pigz' : p.pigz_path}
		
	y1 = parselfolder + "/snake.trim.yaml"
	with open(y1, 'w') as yaml_file:
		yaml.dump(snakedict1, yaml_file, default_flow_style=False)
		
	snakecall1 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.Trim --cores ' + p.max_cores + ' --config cfp="' + y1 + '"'
	os.system(snakecall1)
	click.echo(gettime() + "Sample trimming done.", logf)
	
	# -------------
	# Alignment
	# -------------
	if not os.path.exists(logfolder + "/bowtie2logs"):
		os.makedirs(logfolder+ "/bowtie2logs")
	if not os.path.exists(outfolder + "/02_aligned_reads"):
		os.makedirs(outfolder+ "/02_aligned_reads")
		
	click.echo(gettime() + "Aligning samples", logf)
	
	snakedict2 = {'bowtie2' : p.bowtie2_path, 'bowtie2index' : p.bowtie2_index,
		'outdir' : outfolder, 'samtools' : p.samtools_path}
	
	y2 = parselfolder + "/snake.align.yaml"
	with open(y2, 'w') as yaml_file:
		yaml.dump(snakedict2, yaml_file, default_flow_style=False)
	
	snakecall2 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.Align --cores ' + p.max_cores + ' --config cfp="' + y2 + '"'
	os.system(snakecall2)
	
	# -------------
	# bam process
	# -------------
	if not os.path.exists(outfolder + "/03_processed_reads"):
		os.makedirs(outfolder+ "/03_processed_reads")
	if not os.path.exists(outfolder + "/03_processed_reads/individual"):
		os.makedirs(outfolder+ "/03_processed_reads/individual")
	if not os.path.exists(logfolder + "/rmduplogs"):
		os.makedirs(logfolder+ "/rmduplogs")
		
	click.echo(gettime() + "Cleaning up .bam files", logf)
	
	# Determine chromosomes to keep / filter
	chrs = os.popen(p.samtools_path + " idxstats " +  outfolder + "/02_aligned_reads/* | cut -f1").read().strip().split("\n")
	rmchrlist = ["*", "chrY", "MT", "chrM"]
	keepchrs = [x for x in chrs if x not in rmchrlist and len(x) < int(p.chr_name_length)]
	
	snakedict3 = {'keepchrs' : keepchrs, 'read_quality' : p.read_quality, 'java' : p.java_path,
		'outdir' : outfolder, 'samtools' : p.samtools_path, 'project_name' : p.project_name, 'scriptdir' : script_dir}
		
	y3 = parselfolder + "/snake.bamprocess.yaml"
	with open(y3, 'w') as yaml_file:
		yaml.dump(snakedict3, yaml_file, default_flow_style=False)
	
	snakecall3 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.BamProcess --cores ' + p.max_cores + ' --config cfp="' + y3 + '"'
	os.system(snakecall3)
	
	# ---------------------
	# Mitochondria
	# ---------------------
	if(p.extract_mito):
		click.echo(gettime() + "Extracting mitochondrial reads", logf)
		
		if not os.path.exists(outfolder + "/03_processed_reads/mito"):
			os.makedirs(outfolder+ "/03_processed_reads/mito")
			
		# Determine mitochondrial chromosomes
		mitochrs = []
		for name in chrs:
			if 'MT' in name or 'chrM' in name:
				mitochrs.append(name)
				
		snakedictM = {'mitochrs' : mitochrs, 'outdir' : outfolder, 'samtools' : p.samtools_path, 
			'project_name' : p.project_name}
		yM = parselfolder + "/snake.mito.yaml"

		with open(yM, 'w') as yaml_file:
			yaml.dump(snakedictM, yaml_file, default_flow_style=False)
		
		snakecallM = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.MitoProcess --cores ' + p.max_cores + ' --config cfp="' + yM + '"'
		os.system(snakecallM)
		
		
	# ---------------------
	# Peaks
	# ---------------------
	if not os.path.exists(outfolder + "/04_qc"):
		os.makedirs(outfolder+ "/04_qc")
		
	snakedict4 = {
		'bedtools' : p.bedtools_path, 'blacklistFile' : p.blacklistFile, 'macs2' : p.macs2_path,
		'macs2_genome_size' : p.macs2_genome_size, 'n_peaks' : p.n_peaks, 'outdir' : outfolder,
		'peak_width': p.peak_width, 'project_name' : p.project_name, 'R' : p.R_path, 'script_dir' : script_dir
	}	
		
	y4 = parselfolder + "/snake.callpeaksone.yaml"
	with open(y4, 'w') as yaml_file:
		yaml.dump(snakedict4, yaml_file, default_flow_style=False)
	
	snakecall4 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.CallPeaksOne --cores ' + p.max_cores + ' --config cfp="' + y4 + '"'
	os.system(snakecall4)
	
	# ---------------------
	# Individual QC
	# ---------------------	
	if not os.path.exists(outfolder + "/04_qc/individualQC"):
		os.makedirs(outfolder+ "/04_qc/individualQC")
		
	snakedict5 = {
		'bedtools' : p.bedtools_path, 'outdir' : outfolder,
		'project_name' : p.project_name, 'R' : p.R_path, 'samtools' : p.samtools_path, 
		'script_dir' : script_dir
	}	
		
	y5 = parselfolder + "/snake.qcstats.yaml"
	with open(y5, 'w') as yaml_file:
		yaml.dump(snakedict5, yaml_file, default_flow_style=False)
	
	snakecall5 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.QCstats --cores ' + p.max_cores + ' --config cfp="' + y5 + '"'
	os.system(snakecall5)
		
	# Suspend logging
	logf.close()
	
	
