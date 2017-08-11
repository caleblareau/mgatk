import itertools
import time
import shutil
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

def check_R_packages(required_packages):
	"""
	Determines whether or not R packages are properly installed
	"""
	R_path = shutil.which("R")
	installed_packages = os.popen(R_path + ''' -e "installed.packages()" | awk '{print $1}' | sort | uniq''').read().strip().split("\n")
	if(not set(required_packages) < set(installed_packages)):
		sys.exit("ERROR: cannot find the following R package: " + str(set(required_packages) - set(installed_packages)) + "\n" + 
			"Install it in your R console and then try rerunning proatac (but there may be other missing dependencies).")

def check_software_exists(tool):
	tool_path = shutil.which(tool)
	if(str(tool_path) == "None"):
		sys.exit("ERROR: cannot find "+tool+" in environment; add it to user PATH environment")
