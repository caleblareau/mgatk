import click
import os
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
from mgatk.mgatkHelp import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs
from multiprocessing import Pool

@click.command()
@click.version_option()
@click.option('--input', '-i', default = "", required=True, help='Input; either directory of singular .bam file; see wiki')
@click.option('--output', '-o', default="mgatk_out", help='Output directory for analysis.')
@click.option('--name', '-n', default="mgatk_del",  help='Prefix for project name')

@click.option('--mito-chromosome', '-mc', default = "chrM", required=True, help='Mitochondria chromosome name; see wiki')
@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')

@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see wiki.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')

@click.option('--left-coordinates', '-lc', default = "1", help='Comma separated values for right coordinate of deletions; see wiki')
@click.option('--right-coordinates', '-rc', default = "1000", help='Comma separated values for right coordinate of deletions; see wiki')

@click.option('--read-length', '-rl', default = "72", help='Expected length of a single read from the .bam file')
@click.option('--near-param', '-np', default = "9", help='Number of bases from the start of each read mates ignored for estimating heteroplasmy.')
@click.option('--far-param', '-fp', default = "24", help='Number of bases inside the insert of the read mates ignored for estimating heteroplasmy.')

@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')
@click.option('--snake-stdout', '-so', is_flag=True, help='Write snakemake log to sdout rather than a file.')

def main(input, output, name, mito_chromosome, ncores,
	cluster, jobs, 
	left_coordinates, right_coordinates,
	read_length, window_far, window_near,
	keep_temp_files, snake_stdout):
	
	"""
	mgatk-del: quantify deletion heteroplasmy in mtDNA. \n
	"""
	
	script_dir = os.path.dirname(os.path.realpath(__file__))
	cwd = os.getcwd()
	__version__ = get_distribution('mgatk').version
	click.echo(gettime() + "mgatk-del v%s" % __version__)
	
	# Determine cores
	if(ncores == "detect"):
		ncores = str(available_cpu_count())
	else:
		ncores = str(ncores)
	
	# Verify R dependencies
	check_software_exists("R")
	check_R_packages(["dplyr"])
	mito_chr = mito_chromosome
	
	# -------------------------------
	# Determine samples for analysis
	# -------------------------------
	bams = []
	bams = os.popen('ls ' + input + '/*.bam').read().strip().split("\n")

	if bams[0] == '':
		sys.exit('ERROR: Could not import any samples from the user specification; check --input parameter; QUITTING')

	samples = []
	samplebams = []

	# Loop over bam files
	for bam in bams:
		base=os.path.basename(bam)
		basename=os.path.splitext(base)[0]
		samples.append(basename)
		samplebams.append(bam)
		
	# parallel process to ensure that we have .bai files for each bam
	pool = Pool(processes=int(ncores))
	pm = pool.map(verify_bai, samplebams)
	pool.close()
	
	samples_fail = []
	for i in range(len(samples)):
		sample = samples[i]
		bam = samplebams[i]
		mito_length = -9
		if( not verify_sample_mitobam(bam, mito_chr, mito_length)):
			samples_fail.append(sample)
			
	if(len(samples_fail) > 0):
		click.echo(gettime() + "NOTE: the samples below either have 0 mtDNA reads at the specified chromosome or are mapped to an incorrectly specified reference mitochondrial genome")
		click.echo(gettime() + "Will remove samples from processing:" )
		rmidx = findIdx(samples, samples_fail)
		for index in sorted(rmidx, reverse=True):
			print("REMOVED: ", samples[index])
			del samples[index]
			del samplebams[index]
		
	if not len(samples) > 0:
		sys.exit('ERROR: Could not import any samples from the user specification. \nERROR: check flags, logs, and input configuration (including reference mitochondrial genome); \nQUITTING')
	
	# Make all of the output folders if necessary
	of = output; tf = of + "/temp"; logs = of + "/logs"
	folders = [of, tf, logs, of + "/.internal", tf + "/del",
		 of + "/.internal/parseltongue", of + "/.internal/samples", of + "/final"]

	mkfolderout = [make_folder(x) for x in folders]
		

	# Create internal README files 
	if not os.path.exists(of + "/.internal/README"):
		with open(of + "/.internal/README" , 'w') as outfile:
			outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")
	if not os.path.exists(of + "/.internal/parseltongue/README"):	
		with open(of + "/.internal/parseltongue/README" , 'w') as outfile:
			outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")
	if not os.path.exists(of + "/.internal/samples/README"):
		with open(of + "/.internal" + "/samples" + "/README" , 'w') as outfile:
			outfile.write("This folder creates samples to be interpreted by Snakemake; don't modify it.\n\n")
	
	nsamplesNote = "mgatk-del will process " + str(len(samples)) + " samples"
	click.echo(gettime() + "%s" % nsamplesNote)
	
	# Set up sample bam plain text file
	for i in range(len(samples)):
		with open(of + "/.internal/samples/" + samples[i] + ".bam.txt" , 'w') as outfile:
			outfile.write(samplebams[i])
	
	click.echo(gettime() + "Genotyping samples with "+ncores+" threads.")
	
	# add sqs to get .yaml to play friendly https://stackoverflow.com/questions/39262556/preserve-quotes-and-also-add-data-with-quotes-in-ruamel
	dict1 = {'input_directory' : sqs(input), 'output_directory' : sqs(output), 'script_dir' : sqs(script_dir),
		'name' : sqs(name),
		'left_coordinates' : sqs(left_coordinates), 'right_coordinates' : sqs(right_coordinates),
		'read_length' : sqs(read_length),
		'window_far' : sqs(window_far), 'window_near' : sqs(window_near)}
	
	# Potentially submit jobs to cluster
	snakeclust = ""
	njobs = int(jobs)
	if(njobs > 0 and cluster != ""):
		snakeclust = " --jobs " + jobs + " --cluster '" + cluster + "' "
		click.echo(gettime() + "Recognized flags to process jobs on a computing cluster.")
			
	y_s = of + "/.internal/parseltongue/snake.dels.yaml"
	with open(y_s, 'w') as yaml_file:
		yaml.dump(dict1, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
	
	cp_call = "cp " + y_s +  " " + logs + "/" + name + ".parameters_del.txt"
	os.system(cp_call)
	
	# Execute snakemake
	snake_stats = logs + "/" + name + ".snakemake_del.stats"
	snake_log = logs + "/" + name + ".snakemake_del.log"
	
	snake_log_out = ""
	if not snake_stdout:
		snake_log_out = ' &>' + snake_log
		
	snakecmd_del = 'snakemake'+snakeclust+' --snakefile ' + script_dir + '/singles_del/Snakefile.delSingles --cores '+ncores+' --config cfp="'  + y_s + '" --stats '+snake_stats + snake_log_out
	os.system(snakecmd_del)
	click.echo(gettime() + "mgatk-del successfully processed the supplied .bam files")
	click.echo(gettime() + "Successfully created final output files")
	
	if keep_temp_files:
		click.echo(gettime() + "Retaining temporary files.")
	else:
		byefolder = of
		shutil.rmtree(byefolder + "/.internal")
		shutil.rmtree(byefolder + "/temp")
		click.echo(gettime() + "Intermediate files successfully removed.")
	
