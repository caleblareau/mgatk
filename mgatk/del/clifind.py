import click
import os
import os.path
import sys
import shutil
import random
import string
import time
import pysam
import numpy as np
import re
import csv

from pkg_resources import get_distribution


@click.command()
@click.version_option()
@click.option('--input', '-i', default = ".", required=True, help='Input; a single .bam file of reads to be processed.')
@click.option('--mito-chromosome', '-mc', default = "chrM", required=True, help='Name of mtDNA chromosome in bam file (e.g. chrM or MT)')
@click.option('--output', '-o', default = "mgatkdel_find", required=True, help='Name of output files prefix')

def main(input, mito_chromosome, output):
	
	"""
	mgatk-del-find: detect possible deletion junctions from bam files. \n
	See: `mgatk-del-find --help`
	"""
	
	
	def gettime(): 
		"""
		Matches `date` in Linux
		"""
		return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + 
			time.strftime("%Z ") + time.strftime("%Y")+ ": ")

	
	script_dir = os.path.dirname(os.path.realpath(__file__))
	cwd = os.getcwd()
	__version__ = get_distribution('mgatk').version
	click.echo(gettime() + "mgatk-del-find v%s" % __version__)
	R_plot_script = script_dir + "/bulk_del/plot_deletion_breaks_bulk.R"
	
	bam_in = pysam.AlignmentFile(input, "rb")
	click.echo(gettime() + "Processing .bam file.")
	
	def process_cigar_for_clip_position(cigar, tuple):

		pos = None

		# Case 1/2: start of read
		if(cigar[1] == "H" or cigar[2] == "H" or
		   cigar[1] == "S" or cigar[2] == "S"):
			pos = tuple[0][1]

		# Case 2: end of read
		if(cigar[-1] == "H" or cigar[-1] == "S"):
			pos = tuple[-1][1]

		return(pos)
		
	# function to get tags
	def getTag(intags, tag):
		for tg in intags:
			if(tag == tg[0]):
				return(tg[1])
		return(None)
		
	def SA_pos(SA_tag):
		position = int(SA_tag.split(',')[1]) - 1 # 0-based
		if SA_tag.split(',')[3][1] == "M" or SA_tag.split(',')[3][2] == "M":
			position += (int(SA_tag.split(',')[3].split('M')[0]) - 1)
		return(position)
		
	#iniate empty lists for storing counts
	# loop and extract per position clipping counts
	clip_pos_count = [0] * 16569
	SA_count = [0] * 16569
	
	for read in bam_in.fetch(mito_chromosome):
		seq = read.seq
		reverse = read.is_reverse
		#quality = read.query_qualities
		#align_qual_read = read.mapping_quality
		cigar_string = read.cigarstring
		positions = read.get_reference_positions()
		tuple = read.get_aligned_pairs(True)
		if(positions and cigar_string):
			clip_pos = process_cigar_for_clip_position(cigar_string, tuple)
			if clip_pos is not None:
				clip_pos_count[clip_pos] += 1
			SA_tag=getTag(read.tags,'SA')
			if SA_tag is not None:
				if str(SA_tag.split(',')[0]) == str(mito_chromosome):
					position = SA_pos(SA_tag)
					#print(str(positions[0]) + " " + str(positions[-1]) + " " + str(position) + " " + str(SA_tag))
					SA_count[position] +=1	
					
	clip_pos_count = np.array(clip_pos_count)
	SA_count = np.array(SA_count)
	click.echo(gettime() + "Finished processing bam file.")
	
	# get per base coverage; collapse to N bases not ber nucleotide
	# and convert to list

	cov =  bam_in.count_coverage(mito_chromosome, quality_threshold = 0, read_callback = "nofilter")
	cov_out =  np.array(np.add(np.add(cov[0],cov[1]), np.add(cov[2],cov[3])).tolist())

	outfile_clip = output + ".clip.tsv"
	with open(outfile_clip, 'w') as f:
		writer = csv.writer(f, delimiter='\t')
		idx = np.argsort(-clip_pos_count)
		writer.writerow(["position", "coverage", "clip_count", "SA"])
		writer.writerows(zip(np.array(range(1,16570))[idx],cov_out[idx],clip_pos_count[idx], SA_count[idx]))
	
	# Make the R call
	click.echo(gettime() + "Visualizing results.")
	Rcall = "Rscript " + R_plot_script + " " + outfile_clip
	os.system(Rcall)
		