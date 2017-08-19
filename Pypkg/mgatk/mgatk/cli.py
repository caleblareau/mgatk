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
from .mgatkHelp import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs

@click.command()
@click.version_option()
@click.argument('mode', type=click.Choice(['call', 'check', 'one', 'gather', 'support']))
@click.option('--input', '-i', default = ".", required=True, help='Input; either directory of singular .bam file; see documentation')
@click.option('--output', '-o', default="mgatk_out", help='Output directory for analysis required for `call` and `one`; see documentation.')
@click.option('--name', '-n', default="mgatk",  help='Prefix for project name')

@click.option('--mito-genome', '-m', default = "hg19", required=True, help='mitochondrial genome configuration. Choose hg19, mm10, or a custom .fasta file; see documentation')
@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')

@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see documentation.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')

@click.option('--NHmax', default = "1", help='Maximum number of read alignments allowed as governed by the NH flag.')
@click.option('--NMmax', default = "4", help='Maximum number of paired mismatches allowed represented by the NM/nM tags.')

@click.option('--remove-duplicates', '-rd', is_flag=True, help='Removed marked (presumably PCR) duplicates from Picard; not recommended for low-coverage RNA-Seq')
@click.option('--max-javamem', '-jm', default = "4000m", help='Maximum memory for java')

@click.option('--proper-pairs', '-pp', is_flag=True, help='Require reads to be properly paired.')

@click.option('--base-qual', '-q', default = "20", help='Minimum base quality for deciding that a variant is real.')
@click.option('--blacklist-percentile', '-bp', default = "33", help='Samples with percentile depth below this number will be excluded when determining the blacklist.')

@click.option('--clipL', '-cl', default = "0", help='Number of variants to clip from left hand side of read.')
@click.option('--clipR', '-cr', default = "0", help='Number of variants to clip from right hand side of read.')

@click.option('--keep-samples', '-k', default="ALL", help='Comma separated list of sample names to keep; ALL (special string) by default. Sample refers to basename of .bam file')
@click.option('--ignore-samples', '-g', default="NONE", help='Comma separated list of sample names to ignore; NONE (special string) by default. Sample refers to basename of .bam file')

@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')
@click.option('--detailed-calls', '-dc', is_flag=True, help='Perform detailed variant calling; may be slow.')

@click.option('--skip-R', '-sr', is_flag=True, help='Generate plain-text only output. Otherwise, this generates a .rds obejct that can be immediately read into R for downstream analysis.')



def main(mode, input, output, name, mito_genome, ncores,
	cluster, jobs, nhmax, nmmax, 
	remove_duplicates, max_javamem, proper_pairs, blacklist_percentile,
	base_qual, clipl, clipr, keep_samples, ignore_samples,
	detailed_calls, keep_temp_files, skip_r):
	
	"""
	mgatk: a mitochondrial genome analysis toolkit. \n
	MODE = ['call', 'one', 'check', 'gather', 'support'] \n
	See https://mgatk.readthedocs.io for more details.
	"""
	
	script_dir = os.path.dirname(os.path.realpath(__file__))
	cwd = os.getcwd()
	__version__ = get_distribution('mgatk').version
	click.echo(gettime() + "mgatk v%s" % __version__)
	
	# Determine which genomes are available
	rawsg = os.popen('ls ' + script_dir + "/bin/anno/fasta/*.fasta").read().strip().split("\n")
	supported_genomes = [x.replace(script_dir + "/bin/anno/fasta/", "").replace(".fasta", "") for x in rawsg]  
	
	if(mode == "support"):
		click.echo(gettime() + "List of built-in genomes supported in mgatk:")
		click.echo(gettime() + str(supported_genomes))
		sys.exit(gettime() + 'Specify one of these genomes or provide your own .fasta file with the --mito-genome flag')
		
	if(mode == "check"):
		click.echo(gettime() + "checking dependencies...")
	

	# Verify dependencies
	check_software_exists("python")
	check_software_exists("samtools")
	
	if remove_duplicates:
		check_software_exists("java")
	
	if (mode == "call" or mode == "gather"):
		if not skip_r:	
			check_software_exists("R")
			check_R_packages(['mgatk', 'ggplot2', "dtplyr", "dplyr"])
	
	
	# -------------------------------
	# Determine samples for analysis
	# -------------------------------
	sampleregex = re.compile(r"^[^.]*")
	
	if(mode == "check" or mode == "call"):
	
		bams = []
		bams = os.popen('ls ' + input + '/*.bam').read().strip().split("\n")

		if bams[0] == '':
			sys.exit('ERROR: Could not import any samples from the user specification; check flags, logs and input configuration; QUITTING')
	
		samples = []
		samplebams = []
	
		for bam in bams:
			samples.append(re.search(sampleregex, os.path.basename(bam)).group(0))
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
	
		nsamplesNote = "mgatk will process " + str(len(samples)) + " samples"
		
	if(mode == "check"):
		# Exit gracefully
		sys.exit(gettime() + "mgatk check passed! "+nsamplesNote+" if same parameters are run in `call` mode")		
	elif(mode == "one"):
		# Input argument is assumed to be a .bam file
		filename, file_extension = os.path.splitext(input)
		if(file_extension != ".bam"):
			sys.exit('ERROR: in `one` mode, the input should be an individual .bam file.')
		if not os.path.exists(input):
			sys.exist('ERROR: No file found called "' + input + '"; please specify a valid .bam file')
		samples = [re.search(sampleregex, os.path.basename(input)).group(0)]
		samplebams = [input]
	

	if(mode == "call" or mode == "one"):
	
		# Make all of the output folders if necessary
		of = output; tf = of + "/temp"; qc = of + "/qc"
		folders = [of + "/logs", of + "/logs/filterlogs", of + "/fasta", of + "/.internal",
			 of + "/.internal/parseltongue", of + "/.internal/samples", of + "/final", 
			 tf, tf + "/ready_bam", tf + "/temp_bam", tf + "/sparse_matrices", tf + "/quality",
			 qc, qc + "/quality", qc + "/depth", qc + "/detailed"]

		mkfolderout = [make_folder(x) for x in folders]

		if(mode == "call"):
			# Logging		
			logf = open(output + "/logs" + "/base.mgatk.log", 'a')
			click.echo(gettime() + "Starting analysis with mgatk", logf)
			click.echo(gettime() + nsamplesNote, logf)

		if (remove_duplicates):
				make_folder(of + "/logs/rmdupslogs")
	
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
	
		# Set up sample bam plain text file
		for i in range(len(samples)):
			with open(of + "/.internal/samples/" + samples[i] + ".bam.txt" , 'w') as outfile:
				outfile.write(samplebams[i])
	
		#-------------------
		# Handle .fasta file
		#-------------------
	
		if any(mito_genome in s for s in supported_genomes):
			fastaf = script_dir + "/bin/anno/fasta/" + mito_genome + ".fasta"
			if(mode == "call"):
				click.echo(gettime() + "Found designated mitochondrial genome: %s" % mito_genome, logf)
		else:
			if os.path.exists(mito_genome):
				fastaf = mito_genome
			else:
				sys.exit('ERROR: Could not find file ' + mito_genome + '; QUITTING')
		fasta = parse_fasta(fastaf)	

		if(len(fasta.keys()) != 1):
			sys.exit('ERROR: .fasta file has multiple chromosomes; supply file with only 1; QUITTING')
		mito_genome, mito_seq = list(fasta.items())[0]
		mito_length = len(mito_seq)
	
		newfastaf = of + "/fasta/" + mito_genome + ".fasta"
		if not os.path.exists(newfastaf):
			shutil.copyfile(fastaf, newfastaf)
			fastaf = newfastaf
			pysam.faidx(fastaf)
	
			with open(of + "/final/" + mito_genome + "_refAllele.txt", 'w') as f:
				b = 1
				for base in mito_seq:
					f.write(str(b) + "\t" + base + "\n")
					b += 1
				f.close()
		
		# Other command line arguments	
		if(ncores == "detect"):
			ncores = str(available_cpu_count())
		else:
			ncores = str(ncores)


		click.echo(gettime() + "Processing .bams with "+ncores+" threads")
		if(detailed_calls):
			click.echo(gettime() + "Also performing detailed variant calling.")
	
		
		# add sqs to get .yaml to play friendly https://stackoverflow.com/questions/39262556/preserve-quotes-and-also-add-data-with-quotes-in-ruamel
		dict1 = {'input_directory' : sqs(input), 'output_directory' : sqs(output), 'script_dir' : sqs(script_dir),
			'fasta_file' : sqs(fastaf), 'mito_genome' : sqs(mito_genome), 'mito_length' : sqs(mito_length), 
			'base_qual' : sqs(base_qual), 'remove_duplicates' : sqs(remove_duplicates), 'blacklist_percentile' : sqs(blacklist_percentile), 
			'clipl' : sqs(clipl), 'clipr' : sqs(clipr), 'proper_paired' : sqs(proper_pairs),
			'NHmax' : sqs(nhmax), 'NMmax' : sqs(nmmax), 'detailed_calls' : sqs(detailed_calls), 'max_javamem' : sqs(max_javamem)}
		
		if(mode == "call"):
			
			# Potentially submit jobs to cluster
			snakeclust = ""
			njobs = int(jobs)
			if(njobs > 0 and cluster != ""):
				snakeclust = " --jobs " + jobs + " --cluster '" + cluster + "' "
				click.echo(gettime() + "Recognized flags to process jobs on a computing cluster.", logf)
				
			click.echo(gettime() + "Processing .bams with "+ncores+" threads", logf)
			
			y_s = of + "/.internal/parseltongue/snake.scatter.yaml"
			with open(y_s, 'w') as yaml_file:
				yaml.dump(dict1, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
			
			# Execute snakemake
			snakecmd_scatter = 'snakemake'+snakeclust+' --snakefile ' + script_dir + '/bin/snake/Snakefile.Scatter --cores '+ncores+' --config cfp="' + y_s + '" -T'
			os.system(snakecmd_scatter)
			click.echo(gettime() + "mgatk successfully processed the supplied .bam files", logf)
		
		if(mode == "one"):
		
			# Don't run this through snakemake as we may be trying to handle multiple at the same time
			sample = samples[0]
			inputbam = samplebams[0]
			outputbam = output + "/temp/ready_bam/"+sample+".qc.bam"
			
			y_s = of + "/.internal/samples/"+sample+".yaml"
			with open(y_s, 'w') as yaml_file:
				yaml.dump(dict1, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
			
			# Call the python script
			oneSample_py = script_dir + "/bin/python/oneSample.py"
			pycall = " ".join(['python', oneSample_py, y_s, inputbam, outputbam, sample])
			os.system(pycall)
	
	#-------
	# Gather
	#-------
	if(mode == "gather" or mode == "call"):
		
		if(mode == "call"):
			mgatk_directory = output
			
		elif(mode == "gather"): # in gather, the input argument specifies where things are
			logf = open(input + "/logs" + "/base.mgatk.log", 'a')
			click.echo(gettime() + "Gathering samples that were pre-called with `one`.", logf)
			mgatk_directory = input
			
		dict2 = {'mgatk_directory' : sqs(mgatk_directory), 'name' : sqs(name), 'script_dir' : sqs(script_dir)}
		y_g = mgatk_directory + "/.internal/parseltongue/snake.gather.yaml"
		with open(y_g, 'w') as yaml_file:
			yaml.dump(dict2, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)

		snakecmd_gather = 'snakemake'+snakeclust+' --snakefile ' + script_dir + '/bin/snake/Snakefile.Gather --config cfp="' + y_g + '" -T'
		os.system(snakecmd_gather)
		
		# Make .rds file from the output
		Rcall = "Rscript " + script_dir + "/bin/R/toRDS.R " + mgatk_directory + "/final " + name
		os.system(Rcall)
		
		click.echo(gettime() + "Successfully created final output files", logf)
	
	#--------
	# Cleanup
	#--------
	if(mode == "call" or mode == "gather"):
		if keep_temp_files:
			click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified.", logf)
		else:
			if(mode == "call"):
				byefolder = of
			if(mode == "gather"):
				byefolder = input
			
			shutil.rmtree(byefolder + "/fasta")
			shutil.rmtree(byefolder + "/.internal")
			shutil.rmtree(byefolder + "/temp")
			if not detailed_calls:
				if os.path.exists(byefolder + "/qc/detailed"):
					shutil.rmtree(byefolder + "/qc/detailed")
			click.echo(gettime() + "Intermediate files successfully removed.", logf)
		
		# Suspend logging
		logf.close()
	
