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
@click.option('--input', '-i', default = ".", required=True, help='Input folder; should contain all bams (1 per cell) that will be assessed for deletion heteroplasmy.')
@click.option('--mito-chromosome', '-c', default = "chrM", required=True, help='Name of mtDNA chromosome in bam file (e.g. chrM or MT)')


def main(input, mito_chromosome):
	
	"""
	mgatk-del: detect possible deletion junctions from bam files. \n
	See: `mgatk-del-fin --help`
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
	click.echo(gettime() + "mgatk-del v%s" % __version__)
	
		
	