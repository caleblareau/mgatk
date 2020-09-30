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
import math

from pkg_resources import get_distribution
from subprocess import call, check_call
from .mgatkHelp import *
from ruamel import yaml
from itertools import repeat
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs
from multiprocessing import Pool

@click.command()
@click.version_option()
@click.argument('mode', type=click.Choice(['bcall', 'call', 'tenx', 'check','support']))
@click.option('--input', '-i', default = ".", required=True, help='Input; either directory of singular .bam file; see documentation. REQUIRED.')
@click.option('--output', '-o', default="mgatk_out", help='Output directory for analysis required for `call` and `bcall`. Default = mgatk_out')
@click.option('--name', '-n', default="mgatk",  help='Prefix for project name. Default = mgatk')

@click.option('--mito-genome', '-g', default = "rCRS", required=True, help='mitochondrial genome configuration. Choose hg19, hg38, mm10, (etc.) or a custom .fasta file; see documentation. Default = rCRS.')
@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')

@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see documentation.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')

@click.option('--barcode-tag', '-bt', default = "X",  help='Read tag (generally two letters) to separate single cells; valid and required only in `bcall` mode.')
@click.option('--barcodes', '-b', default = "",  help='File path to barcodes that will be extracted; useful only in `bcall` mode. If none supplied, mgatk will learn abundant barcodes from the bam file (threshold defined by the -mb tag).')
@click.option('--min-barcode-reads', '-mb', default = 1000,  help='Minimum number of mitochondrial reads for a barcode to be genotyped; useful only in `bcall` mode; will not overwrite the `--barcodes` logic. Default = 1000.')

@click.option('--NHmax', default = 1, help='Maximum number of read alignments allowed as governed by the NH flag. Default = 1.')
@click.option('--NMmax', default = 4, help='Maximum number of paired mismatches allowed represented by the NM/nM tags. Default = 4.')

@click.option('--keep-duplicates', '-kd', is_flag=True, help='Retained dupliate (presumably PCR) reads')
@click.option('--umi-barcode', '-ub', default = "",  help='Read tag (generally two letters) to specify the UMI tag when removing duplicates for genotyping.')

@click.option('--max-javamem', '-jm', default = "8000m", help='Maximum memory for java for running duplicate removal per core. Default = 8000m.')

@click.option('--proper-pairs', '-pp', is_flag=True, help='Require reads to be properly paired.')

@click.option('--base-qual', '-q', default = 0, help='Minimum base quality for inclusion in the genotype count. Default = 0.')
@click.option('--alignment-quality', '-aq', default = 0, help='Minimum alignment quality to include read in genotype. Default = 0.')
@click.option('--emit-base-qualities', '-eb', is_flag=True, help='Output mean base quality per alt allele as part of the final output.')

@click.option('--nsamples', '-ns', default = 7000, help='The number of samples / cells to be processed per iteration; Default = 7000. Supply 0 to try all.')

@click.option('--keep-samples', '-k', default="ALL", help='Comma separated list of sample names to keep; ALL (special string) by default. Sample refers to basename of .bam file')
@click.option('--ignore-samples', '-x', default="NONE", help='Comma separated list of sample names to ignore; NONE (special string) by default. Sample refers to basename of .bam file')

@click.option('--keep-temp-files', '-z', is_flag=True, help='Add this flag to keep all intermediate files.')
@click.option('--keep-qc-bams', '-qc', is_flag=True, help='Add this flag to keep the quality-controlled bams after processing.')

@click.option('--skip-R', '-sr', is_flag=True, help='Generate plain-text only output. Otherwise, this generates a .rds obejct that can be immediately read into R for downstream analysis.')
@click.option('--snake-stdout', '-so', is_flag=True, help='Write snakemake log to sdout rather than a file. May be necessary for certain HPC environments.')

def main(mode, input, output, name, mito_genome, ncores,
	cluster, jobs, barcode_tag, barcodes, min_barcode_reads,
	nhmax, nmmax, keep_duplicates, umi_barcode, max_javamem, 
	proper_pairs, base_qual, alignment_quality, emit_base_qualities,
	nsamples, keep_samples, ignore_samples,
	keep_temp_files, keep_qc_bams, skip_r, snake_stdout):
	
	"""
	mgatk: a mitochondrial genome analysis toolkit. \n
	MODE = ['bcall', 'call', 'tenx', 'check', 'support'] \n
	See https://github.com/caleblareau/mgatk/wiki for more details.
	"""
	
	script_dir = os.path.dirname(os.path.realpath(__file__))
	cwd = os.getcwd()
	__version__ = get_distribution('mgatk').version
	click.echo(gettime() + "mgatk v%s" % __version__)
	
	# Determine cores
	if(ncores == "detect"):
		ncores = str(available_cpu_count())
	else:
		ncores = str(ncores)
	
	# Now removing duplicates is the default, so recode variable names
	if(keep_duplicates):
		remove_duplicates = False
	else:
		remove_duplicates = True
	
	# Verify dependencies	
	if remove_duplicates:
		check_software_exists("java")
	
	if (mode == "call" or mode == "tenx" or mode == "bcall"):
		if not skip_r:
			check_software_exists("R")
			check_R_packages(["data.table", "SummarizedExperiment", "GenomicRanges", "Matrix"])
	
	
	# Determine which genomes are available
	rawsg = os.popen('ls ' + script_dir + "/bin/anno/fasta/*.fasta").read().strip().split("\n")
	supported_genomes = [x.replace(script_dir + "/bin/anno/fasta/", "").replace(".fasta", "") for x in rawsg]  
	
	if(mode == "support"):
		click.echo(gettime() + "List of built-in genomes supported in mgatk:")
		click.echo(gettime() + str(supported_genomes))
		sys.exit(gettime() + 'Specify one of these genomes or provide your own .fasta file with the --mito-genome flag')
		
	if(mode == "check"):
		click.echo(gettime() + "checking dependencies...")
	
	# Remember that I started off as bcall as this will become overwritten
	wasbcall = False
	if(mode == "bcall" or mode == "tenx"):
	
		if(barcode_tag == "X"):
			sys.exit('ERROR: in `'+mode+'` mode, must specify a valid read tag ID (generally two letters).')
			
		# Input argument is assumed to be a .bam file
		filename, file_extension = os.path.splitext(input)
		if(file_extension != ".bam"):
			sys.exit('ERROR: in `'+mode+'` mode, the input should be an individual .bam file.')
		if not os.path.exists(input):
			sys.exit('ERROR: No file found called "' + input + '"; please specify a valid .bam file.')
		if not os.path.exists(input + ".bai"):
			click.echo(gettime() + "Attempting to index: " + input)
			pysam.index(input)
			
			if not os.path.exists(input + ".bai"):
				sys.exit('ERROR: index your input .bam file for `bcall` mode.')
		click.echo(gettime() + "Found bam file: " + input + " for genotyping.")
		
		# Determine whether or not we have been supplied barcodes
		barcode_known = False
		if (os.path.exists(barcodes)) and (barcodes != ""):
			click.echo(gettime() + "Found file of barcodes to be parsed: " + barcodes)
			barcode_known = True
		else:
			if(mode == "tenx"):
				sys.exit(gettime() + 'Must specify a known barcode list with `tenx` mode')
			click.echo(gettime() + "Will determine barcodes with at least: " + str(min_barcode_reads) + " mitochondrial reads.")
			
		# Make temporary directory of inputs
		of = output; tf = of + "/temp"; bcbd = tf + "/barcoded_bams" # bcdb = barcoded bam directory
		folders = [of, tf, bcbd, of + "/final"]
		mkfolderout = [make_folder(x) for x in folders]
		
		# Handle fasta requirements
		fastaf, mito_chr, mito_length = handle_fasta_inference(mito_genome, supported_genomes, script_dir, mode, of)
		idxs = pysam.idxstats(input).split("\n")
		
		# Handle common mtDNA reference genome errors
		bam_length = 0
		for i in idxs:
			if(i.split("\t")[0] == mito_chr):
				bam_length = int(i.split("\t")[1])
		
		if(mito_length == bam_length):
			click.echo(gettime() + "User specified mitochondrial genome matches .bam file")
		elif(bam_length == 16569):
			click.echo(gettime() + "User specified mitochondrial genome does NOT match .bam file; using rCRS instead (length == 16569)")
			fastaf, mito_chr, mito_length = handle_fasta_inference("rCRS", supported_genomes, script_dir, mode, of)
		elif(bam_length == 16571):
			click.echo(gettime() + "User specified mitochondrial genome does NOT match .bam file; using hg19 instead (length == 16571)")
			fastaf, mito_chr, mito_length = handle_fasta_inference("hg19", supported_genomes, script_dir, mode, of)
		else:
			click.echo(gettime() + "User specified mitochondrial genome does NOT match .bam file; correctly specify reference genome or .fasta file")
			quit()
		
		# Actually call the external script based on user input
		if(not barcode_known):
			barc_quant_file = of + "/final/barcodeQuants.tsv"
			passing_barcode_file = of + "/final/passingBarcodes.tsv"
			find_barcodes_py = script_dir + "/bin/python/find_barcodes.py"
			
			pycall = " ".join(['python', find_barcodes_py, input, bcbd, barcode_tag, str(min_barcode_reads), mito_chr, barc_quant_file, passing_barcode_file])
			os.system(pycall)
			barcodes = passing_barcode_file

		# Potentially split the valid barcodes into smaller files if we need to
		if(mode == "bcall"):
			barcode_files = split_barcodes_file(barcodes, nsamples, output)
			split_barcoded_bam_py = script_dir + "/bin/python/split_barcoded_bam.py"
			
			# Loop over the split sample files
			for i in range(len(barcode_files)):
				one_barcode_file = barcode_files[i]
				pycall = " ".join(['python', split_barcoded_bam_py, input, bcbd, barcode_tag, one_barcode_file, mito_chr])
				os.system(pycall)
				
			# Update everything to appear like we've just set `call` on the set of bams
			mode = "call"
			input = bcbd 
			wasbcall = True
			
		if(mode == "tenx"):
			barcode_files = split_barcodes_file(barcodes, math.ceil(file_len(barcodes)/int(ncores)), output)
			samples = [os.path.basename(os.path.splitext(sample)[0]) for sample in barcode_files] 
			samplebams = [of + "/temp/barcoded_bams/" + sample + ".bam" for sample in samples]
			
			if(umi_barcode == ""):
				umi_barcode = "XX"
			
			# Enact the split in a parallel manner
			pool = Pool(processes=int(ncores))
			pmblah = pool.starmap(split_chunk_file, zip(barcode_files, repeat(script_dir), repeat(input), repeat(bcbd), repeat(barcode_tag), repeat(mito_chr), repeat(umi_barcode)))
			pool.close()
			umi_barcode = "MU"
		
	
		click.echo(gettime() + "Finished determining/splitting barcodes for genotyping.")
		
	# -------------------------------
	# Determine samples for analysis
	# -------------------------------
	if(mode == "check" or mode == "call"):
	
		bams = []
		bams = os.popen('ls ' + input + '/*.bam').read().strip().split("\n")

		if bams[0] == '':
			sys.exit('ERROR: Could not import any samples from the user specification; check flags, logs and input configuration; QUITTING')
	
		samples = []
		samplebams = []
		
		if(not wasbcall):
			fastaf, mito_chr, mito_length = handle_fasta_inference(mito_genome, supported_genomes, script_dir, mode, output, write_files = False)
		
		# Loop over bam files
		for bam in bams:
			base=os.path.basename(bam)
			basename=os.path.splitext(base)[0]
			samples.append(basename)
			samplebams.append(bam)
			
		# parallel process to ensure we have .bai files for each bam
		pool = Pool(processes=int(ncores))
		pm = pool.map(verify_bai, samplebams)
		pool.close()
		
		samples_fail = []
		for i in range(len(samples)):
			sample = samples[i]
			bam = samplebams[i]
			if( not verify_sample_mitobam(bam, mito_chr, mito_length)):
				samples_fail.append(sample)
				
		if(keep_samples != "ALL"):
			keeplist = keep_samples.split(",")
			click.echo(gettime() + "Intersecting detected samples with user-retained ones: " + keep_samples)
			keepidx = findIdx(samples, keeplist)
			samples = [samples[i] for i in keepidx]
			samplebams = [samplebams[i] for i in keepidx]
		
		if(ignore_samples != "NONE"):
			iglist = ignore_samples.split(",")
			click.echo(gettime() + "Will remove samples from processing:" + ignore_samples)
			rmidx = findIdx(samples, iglist)
			for index in sorted(rmidx, reverse=True):
				del samples[index]
				del samplebams[index]
				
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
	
		nsamplesNote = "mgatk will process " + str(len(samples)) + " samples"
		
	if(mode == "check"):
		# Exit gracefully
		sys.exit(gettime() + "mgatk check passed! "+nsamplesNote+" if same parameters are run in `call` mode")
			

	if(mode == "call" or mode == "tenx"):
	
		# Make all of the output folders if necessary
		of = output; tf = of + "/temp"; qc = of + "/qc"; logs = of + "/logs"
		folders = [logs, of + "/logs/filterlogs", of + "/fasta", of + "/.internal",
			 of + "/.internal/parseltongue", of + "/.internal/samples", of + "/final", 
			 tf, tf + "/ready_bam", tf + "/temp_bam", tf + "/sparse_matrices", tf + "/quality",
			 qc, qc + "/quality", qc + "/depth"]

		mkfolderout = [make_folder(x) for x in folders]
		
		#-------------------
		# Handle .fasta file
		#-------------------
		if((mode == "call" and wasbcall == False)):
			fastaf, mito_genmito_chrome, mito_length = handle_fasta_inference(mito_genome, supported_genomes, script_dir, mode, of)
			print(gettime() + "Found designated mitochondrial chromosome: %s" % mito_chr)
			
		if(mode == "call" or mode == "tenx"):
			# Logging		
			logf = open(output + "/logs" + "/base.mgatk.log", 'a')
			click.echo(gettime() + "Starting analysis with mgatk", logf)
			
			if(mode == "call"):
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
		
		click.echo(gettime() + "Genotyping samples with "+ncores+" threads")
		
		# add sqs to get .yaml to play friendly https://stackoverflow.com/questions/39262556/preserve-quotes-and-also-add-data-with-quotes-in-ruamel
		dict1 = {'input_directory' : sqs(input), 'output_directory' : sqs(output), 'script_dir' : sqs(script_dir),
			'fasta_file' : sqs(fastaf), 'mito_chr' : sqs(mito_chr), 'mito_length' : sqs(mito_length), 'name' : sqs(name),
			'base_qual' : sqs(base_qual), 'remove_duplicates' : sqs(remove_duplicates),
			'barcode_tag' : sqs(barcode_tag), 'umi_barcode' : sqs(umi_barcode),
			'alignment_quality' : sqs(alignment_quality), 'emit_base_qualities' : sqs(emit_base_qualities),
			'proper_paired' : sqs(proper_pairs),
			'NHmax' : sqs(nhmax), 'NMmax' : sqs(nmmax), 'max_javamem' : sqs(max_javamem)}
		
		# Potentially submit jobs to cluster
		snakeclust = ""
		njobs = int(jobs)
		if(njobs > 0 and cluster != ""):
			snakeclust = " --jobs " + jobs + " --cluster '" + cluster + "' "
			click.echo(gettime() + "Recognized flags to process jobs on a computing cluster.", logf)
			
		click.echo(gettime() + "Processing samples with "+ncores+" threads", logf)
		
		y_s = of + "/.internal/parseltongue/snake.scatter.yaml"
		with open(y_s, 'w') as yaml_file:
			yaml.dump(dict1, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
		
		cp_call = "cp " + y_s +  " " + logs + "/" + name + ".parameters.txt"
		os.system(cp_call)
	
		if(mode == "call"):
			
			# Execute snakemake
			snake_stats = logs + "/" + name + ".snakemake_scatter.stats"
			snake_log = logs + "/" + name + ".snakemake_scatter.log"
			
			snake_log_out = ""
			if not snake_stdout:
				snake_log_out = ' &>' + snake_log
				
			snakecmd_scatter = 'snakemake'+snakeclust+' --snakefile ' + script_dir + '/bin/snake/Snakefile.Scatter --cores '+ncores+' --config cfp="'  + y_s + '" --stats '+snake_stats + snake_log_out
			os.system(snakecmd_scatter)
			
		elif(mode == "tenx"):
			
			# Execute snakemake
			snake_stats = logs + "/" + name + ".snakemake_tenx.stats"
			snake_log = logs + "/" + name + ".snakemake_tenx.log"
			
			snake_log_out = ""
			if not snake_stdout:
				snake_log_out = ' &>' + snake_log
				
			snakecmd_tenx = 'snakemake'+snakeclust+' --snakefile ' + script_dir + '/bin/snake/Snakefile.tenx --cores '+ncores+' --config cfp="'  + y_s + '" --stats '+snake_stats + snake_log_out
			os.system(snakecmd_tenx)
		
		click.echo(gettime() + "mgatk successfully processed the supplied .bam files", logf)

	
	#-------
	# Gather
	#-------
	if(mode == "call"):
		
		mgatk_directory = output
		dict2 = {'mgatk_directory' : sqs(mgatk_directory), 'name' : sqs(name),
			'script_dir' : sqs(script_dir)}
		y_g = mgatk_directory + "/.internal/parseltongue/snake.gather.yaml"
		with open(y_g, 'w') as yaml_file:
			yaml.dump(dict2, yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
		
		# Snakemake gather
		snake_stats = logs + "/" + name + ".snakemake_gather.stats"
		snake_log = logs + "/" + name + ".snakemake_gather.log"
		
		snake_log_out = ""
		if not snake_stdout:
			snake_log_out = ' &>' + snake_log 
			
		snakecmd_gather = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.Gather --cores 1 --config cfp="' + y_g + '" --stats '+snake_stats + snake_log_out
		os.system(snakecmd_gather)
	
	if(mode == "call" or mode == "tenx"):
	
		# Make .rds file from the output
		Rcall = "Rscript " + script_dir + "/bin/R/toRDS.R " + output + "/final " + name
		os.system(Rcall)
		click.echo(gettime() + "Successfully created final output files", logf)
	
	#--------
	# Cleanup
	#--------
	if(mode == "call" or mode == "tenx"):
		if keep_qc_bams:
			click.echo(gettime() + "Final bams retained since --keep-qc-bams was specified.", logf)
			dest = shutil.move(of + "/temp/ready_bam", of + "/qc_bam")  

		if keep_temp_files:
			click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified.", logf)
		else:
			shutil.rmtree(of+ "/fasta")
			shutil.rmtree(of + "/.internal")
			shutil.rmtree(of + "/temp")
			click.echo(gettime() + "Intermediate files successfully removed.", logf)
		
		# Suspend logging
		logf.close()
	
