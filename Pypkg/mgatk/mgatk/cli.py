import click
import os
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs
import os.path
import sys
import shutil
import random
import string
import itertools
import time
import pysam
from pkg_resources import get_distribution
from subprocess import call, check_call
from .mgatkHelp import *

# ------------------------------
# Command line arguments
# ------------------------------	

@click.command()
@click.version_option()
@click.argument('mode', type=click.Choice(['call', 'check']))
@click.option('--input', '-i', required=True, help='input directory; assumes .bam files are present')
@click.option('--output', '-o', default="mgatk_out", help='Output directory for analysis')
@click.option('--name', '-n', default="mgatk",  help='Prefix for project name')

@click.option('--mito-genome', '-m', default = "hg19", required=True, help='mitochondrial genome configuration. Choose hg19, mm10, or a custom .fasta file; see documentation')
@click.option('--ncores', '-c', default = "detect", help='Number of cores to run job in parallel.')

@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see documentation.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')

@click.option('--atac-single', '-as', is_flag=True, help='Default parameters for ATAC-Seq single end read analyses; see documentation.')
@click.option('--atac-paired', '-ap',  is_flag=True, help='Default parameters for ATAC-Seq paired end read analyses; see documentation.')
@click.option('--rna-single', '-rs', is_flag=True, help='Default parameters for RNA-Seq single end read analyses; see documentation.')
@click.option('--rna-paired', '-rp', is_flag=True, help='Default parameters for RNA-Seq paired end read analyses; see documentation.')

@click.option('--NHmax', default = "1", help='Maximum number of read alignments allowed as governed by the NH flag.')
@click.option('--NMmax', default = "4", help='Maximum number of paired mismatches allowed represented by the NM/nM tags.')

@click.option('--remove-duplicates', '-rd', is_flag=True, help='Removed marked (presumably PCR) duplicates from Picard; not recommended for low-coverage RNA-Seq')
@click.option('--max-javamem', '-jm', default = "4000m", help='Maximum memory for java')

@click.option('--keep-indels', '-ki', is_flag=True, help='Keep marked indels for analysis; not recommended as this flag has not been well-tested')
@click.option('--proper-pairs', '-pp', is_flag=True, help='Require reads to be properly paired.')

@click.option('--base-qual', '-q', default = "20", help='Minimum base quality for deciding that a variant is real.')
@click.option('--blacklist-percentile', '-bp', default = "33", help='Samples with percentile depth below this number will be excluded when determining the blacklist.')

@click.option('--clipL', '-cl', default = "0", help='Number of variants to clip from left hand side of read.')
@click.option('--clipR', '-cr', default = "0", help='Number of variants to clip from right hand side of read.')

@click.option('--keep-samples', '-k', default="ALL", help='Comma separated list of sample names to keep; ALL (special string) by default. Sample refers to basename of .bam file')
@click.option('--ignore-samples', '-g', default="NONE", help='Comma separated list of sample names to ignore; NONE (special string) by default. Sample refers to basename of .bam file')

@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')
@click.option('--detailed-calls', '-dc', is_flag=True, help='Perform detailed variant calling; may be slow.')

@click.option('--skip-rds', '-s', is_flag=True, help='Generate plain-text only output. Otherwise, this generates a .rds obejct that can be immediately read into R')



def main(mode, input, output, name, mito_genome, ncores,
	cluster, jobs,
	atac_single, atac_paired, rna_single, rna_paired,
	nhmax, nmmax, 
	remove_duplicates, max_javamem, keep_indels, proper_pairs, blacklist_percentile,
	base_qual, clipl, clipr, keep_samples, ignore_samples,
	detailed_calls, keep_temp_files, skip_rds):
	
	"""
	mgatk: a mitochondrial genome analysis toolkit. \n
	MODE = ['call', 'one', 'check', 'gather'] \n
	See https://mgatk.readthedocs.io for more details.
	"""
	
	__version__ = get_distribution('mgatk').version
	script_dir = os.path.dirname(os.path.realpath(__file__))
	click.echo(gettime() + "mgatk v%s" % __version__)
	
	if(mode == "check"):
		click.echo(gettime() + "checking dependencies...")
		
	cwd = os.getcwd()

	# -------------------------------
	# Verify dependencies
	# -------------------------------

	check_software_exists("python")
	check_software_exists("samtools")
	
	if remove_duplicates:
		check_software_exists("java")

	check_software_exists("R")
	check_R_packages(['mgatk', 'ggplot2', "dtplyr", "dplyr"])
	
	# -------------------------------
	# Determine samples for analysis
	# -------------------------------
	
	bams = []
	bams = os.popen('ls ' + input + '/*.bam').read().strip().split("\n")

	if bams[0] == '':
		sys.exit('ERROR: Could not import any samples from the user specification; check flags, logs and input configuration; QUITTING')
	
	samples = []
	samplebams = []
	
	find = re.compile(r"^[^.]*")

	for bam in bams:
		#if(os.path.isfile(bam + ".bai")):
		samples.append(re.search(find, os.path.basename(bam)).group(0))
		samplebams.append(bam)
	
	if(keep_samples != "ALL"):
		keeplist = keep_samples.split(",")
		click.echo(gettime() + "Intersecting detected samples with user-retained ones: " + keep_samples)
		keepidx = findIdx(samples, keeplist)
		samples = [samples[i] for i in keepidx]
		samplebams = [samplebams[i] for i in keepidx]
		
	if(ignore_samples != "NONE"):
		iglist = ignore_samples.split(",")
		click.echo(gettime() + "Attempting to remove samples from processing:" + ignore_samples)
		rmidx = findIdx(samples, iglist)
		for index in sorted(rmidx, reverse=True):
			del samples[index]
			del samplebams[index]
    		
	if not len(samples) > 0:
		sys.exit('ERROR: Could not import any samples from the user specification; check flags, logs and input configuration; QUITTING')
	
	nsamplesNote = "The software will process " + str(len(samples)) + " samples"
	
	if(mode == "check"):
		sys.exit(gettime() + "mgatk check passed! "+nsamplesNote+" if same parameters are run in `call` mode")
	
	# -------------------------------
	# Configure the output setup
	# -------------------------------
	
	of = output
	tf = of + "/temp"
	qc = of + "/qc"
	
	folders = [of + "/logs", of + "/logs/filterlogs", of + "/fasta", of + "/.internal",
		 of + "/.internal/parseltongue", of + "/.internal/samples", of + "/final", 
		 tf, tf + "/ready_bam", tf + "/temp_bam", tf + "/sparse_matrices", tf + "/quality",
		 qc, qc + "/quality", qc + "/depth", qc + "/detailed"]

	mkfolderout = [make_folder(x) for x in folders]
	
	if (remove_duplicates):
			make_folder(of + "/logs/rmdupslogs")
	
	# Give the users some heads up
	with open(of + "/.internal" + "/README" , 'w') as outfile:
		outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")	
	with open(of + "/.internal" + "/parseltongue" + "/README" , 'w') as outfile:
		outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")
	with open(of + "/.internal" + "/samples" + "/README" , 'w') as outfile:
		outfile.write("This folder creates samples to be interpreted by Snakemake; don't modify it.\n\n")
	
	# Set up sample bam plain text file
	for i in range(len(samples)):
		with open(of + "/.internal/samples/" + samples[i] + ".bam.txt" , 'w') as outfile:
			outfile.write(samplebams[i])
	
	# Logging		
	logf = open(of + "/logs" + "/base.mgatk.log", 'a')
	click.echo(gettime() + "Starting analysis with mgatk", logf)
	click.echo(gettime() + nsamplesNote, logf)
	
	# -----------------------------------
	# Parse user-specified parameteres
	# -----------------------------------

	#-------------------------
	# Handle read filter flags
	#-------------------------
	
	if(rna_paired):
		click.echo(gettime() + "Adding read filtering flags for paired RNA.", logf)

	elif(atac_paired):
		click.echo(gettime() + "Adding read filtering flags for paired ATAC/WGS.", logf)

	elif(rna_single):
		click.echo(gettime() + "Adding read filtering flags for single-end RNA.", logf)

	elif(atac_single):
		click.echo(gettime() + "Adding read filtering flags for single-end ATAC/WGS.", logf)

	else:
		click.echo(gettime() + "No specific data type specified.", logf)

	#-----------------------------------
	# Potentially submit jobs to cluster
	#-----------------------------------
	
	snakeclust = ""
	njobs = int(jobs)
	if(njobs > 0 and cluster != ""):
		snakeclust = " --jobs " + jobs + " --cluster '" + cluster + "' "
	
	#-------------------
	# Handle .fasta file
	#-------------------
	
	supported_genomes = ['hg19', 'mm10']
	if any(mito_genome in s for s in supported_genomes):
		click.echo(gettime() + "Found designated mitochondrial genome: %s" % mito_genome, logf)
	fastaf, mito_genome, mito_seq, mito_length = handle_fasta(mito_genome, supported_genomes, script_dir, of, name)
		
	#-----------------------------
	# Other command line arguments
	#-----------------------------
	
	if(keep_indels):
		skip_indels = ""
	else:
		skip_indels = "--skip-indels "
		
	if(ncores == "detect"):
		ncores = str(available_cpu_count())
	else:
		ncores = str(ncores)
	
	proper_paired = ""
	if(proper_pairs):
		proper_paired = " | samtools view -f 0x2 -b - "
	
	click.echo(gettime() + "Processing .bams with "+ncores+" cores", logf)
	click.echo(gettime() + "Processing .bams with "+ncores+" cores")
	if(detailed_calls):
		click.echo(gettime() + "Also performing detailed variant calling.")
	
	#--------
	# Scatter
	#--------
	
	# add sqs to get .yaml to play friendly https://stackoverflow.com/questions/39262556/preserve-quotes-and-also-add-data-with-quotes-in-ruamel
	snakedict1 = {'input_directory' : sqs(input), 'output_directory' : sqs(output), 'script_dir' : sqs(script_dir),
		'fasta_file' : sqs(fastaf), 'mito_genome' : sqs(mito_genome), 'mito_length' : sqs(mito_length), 'name' : sqs(name),
		'base_qual' : sqs(base_qual), 'remove_duplicates' : sqs(remove_duplicates), 'blacklist_percentile' : sqs(blacklist_percentile), 
		'skip_indels' : sqs(skip_indels), 'clipl' : sqs(clipl), 'clipr' : sqs(clipr), 'proper_paired' : sqs(proper_paired),
		'NHmax' : sqs(nhmax), 'NMmax' : sqs(nmmax), 'detailed_calls' : sqs(detailed_calls), 'max_javamem' : sqs(max_javamem)}
	
	y_s = of + "/.internal/parseltongue/snake.scatter.yaml"
	with open(y_s, 'w') as yaml_file:
		yaml.dump(snakedict1, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
	
	snakecmd_scatter = 'snakemake'+snakeclust+' --snakefile ' + script_dir + '/bin/snake/Snakefile.Scatter --cores '+ncores+' --config cfp="' + y_s + '" -T'
	os.system(snakecmd_scatter)
	click.echo(gettime() + "mgatk successfully processed the supplied .bam files", logf)
	
	#-------
	# Gather
	#-------
	
	snakedict2 = {'mgatk_directory' : sqs(output), 'name' : sqs(name), 'script_dir' : sqs(script_dir)}
	
	y_g = of + "/.internal/parseltongue/snake.gather.yaml"
	with open(y_g, 'w') as yaml_file:
		yaml.dump(snakedict2, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
	
	snakecmd_gather = 'snakemake'+snakeclust+' --snakefile ' + script_dir + '/bin/snake/Snakefile.Gather --cores '+ncores+' --config cfp="' + y_g + '" -T'
	os.system(snakecmd_gather)
	click.echo(gettime() + "Successfully created final output files", logf)
	
	#--------
	# Cleanup
	#--------
	
	if keep_temp_files:
		click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified.", logf)
	else:
		shutil.rmtree(of + "/fasta")
		shutil.rmtree(of + "/.internal")
		shutil.rmtree(of + "/temp")
		if not detailed_calls:
			if os.path.exists(of + "/qc/detailed"):
				shutil.rmtree(of + "/qc/detailed")
		click.echo(gettime() + "Intermediate files successfully removed.", logf)
		
	# Suspend logging
	logf.close()
	
