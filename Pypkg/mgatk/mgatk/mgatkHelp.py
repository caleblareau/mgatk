import yaml
import itertools
import time
import re
import os
import sys
import csv

def string_hamming_distance(str1, str2):
    """
    Fast hamming distance over 2 strings known to be of same length.
    In information theory, the Hamming distance between two strings of equal 
    length is the number of positions at which the corresponding symbols 
    are different.
    eg "karolin" and "kathrin" is 3.
    """
    return sum(itertools.imap(operator.ne, str1, str2))


def rev_comp(seq):
    """
    Fast Reverse Compliment
    """  
    tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(tbl[s] for s in seq[::-1])

def parse_manifest(manifest):
	"""
	Basic function to parse yaml/yml
	"""
	samples = []
	if manifest.endswith(('.yaml', '.yml')):
		with open(manifest, 'r') as f: 
			m = yaml.load(f)
		return m
	else:
		click.echo(gettime() + "Please specify a valid .yaml file for analysis")

def gettime(): 
	"""
	Matches `date` in Linux
	"""
	return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + 
		time.strftime("%Z ") + time.strftime("%Y")+ ": ")
		

def findIdx(list1, list2):
	"""
	Return the indices of list1 in list2
	"""
	return [i for i, x in enumerate(list1) if x in list2]

def check_R_packages(required_packages, R_path):
	"""
	Determines whether or not R packages are properly installed
	"""
	installed_packages = os.popen(R_path + ''' -e "installed.packages()" | awk '{print $1}' | sort | uniq''').read().strip().split("\n")
	if(not set(required_packages) < set(installed_packages)):
		sys.exit("ERROR: cannot find the following R package: " + str(set(required_packages) - set(installed_packages)) + "\n" + 
			"Install it in your R console and then try rerunning proatac (but there may be other missing dependencies).")
		

def process_seq_dir(d, logf, listAllSamples):
	"""
	Function that takes a dictionary parsed from the main .yaml
	and returns something more coherent to be processed downstream
	Could parameterize this function further given more versions of scATAC
	that would be worth processing.
	"""
	
	name = d['name']
	version = d['version']
	directory =  d['dir']
	fastq_path = d['fastq_path']
	
	# Handle whether the reads are R1/R2 (default) or something else
	# Some attempt on the author's part to be savvy here. 
	if 'read' in d:
		read = d['read']
	else:
		read = ["R1", "R2"]
	if(len(read) != 2):
		sys.exit("ERROR: invalid input-- " + read + " -- for 'read' in .yaml; need length 2 python vector")
	
	# Go get the fastq files
	noSID = re.sub("{sample_id}", "", fastq_path)
	read1 = re.sub("{read}", read[0], noSID)
	read2 = re.sub("{read}", read[1], noSID)
	
	reads1 = os.popen("ls " + directory + "/*" + read1).read().strip().split("\n")
	reads2 = os.popen("ls " + directory + "/*" + read2).read().strip().split("\n")
	
	# Parse out the sample names
	temp1 = [re.sub(read1, "", x) for x in reads1]
	samples = [re.sub(directory + "/", "", x) for x in temp1]
	
	# Remove / filter samples as requested if applicable; also remove those from the file listing
	if 'keep_samples' in d:
		keep_samples = d['keep_samples']
	else:
		keep_samples = ""
		
	if(keep_samples != ""):
		keeplist = findIdx(samples, keep_samples) # indicies of samples to be kept
		samples = [samples[i] for i in keeplist]
		reads1 = [reads1[i] for i in keeplist]
		reads2 = [reads2[i] for i in keeplist]
	
	if 'remove_samples' in d:
		remove_samples = d['remove_samples']
	else:
		remove_samples = ""
		
	if(remove_samples != ""):
		rmlist = findIdx(samples, remove_samples) # indices of samples that should be removed
		allidx = range(len(samples)) # 1:n
		keeplist2 = [x for x in allidx if x not in rmlist] # indices in 1:n not in list of indices to be removed
		samples = [samples[i] for i in keeplist2]
		reads1 = [reads1[i] for i in keeplist2]
		reads2 = [reads2[i] for i in keeplist2]
	
	with open(listAllSamples, 'a') as f:
		writer = csv.writer(f)
		rows = zip(samples,reads1,reads2)
		for row in rows:
			writer.writerow(row)
	
	