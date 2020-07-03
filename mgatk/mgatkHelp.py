import itertools
import time
import shutil
import re
import os
import sys
import subprocess
import pysam
import filecmp
import math
import pysam

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
			"Install it in your R console and then try rerunning mgatk (but there may be other missing dependencies).")

def check_software_exists(tool):
	tool_path = shutil.which(tool)
	if(str(tool_path) == "None"):
		sys.exit("ERROR: cannot find "+tool+" in environment; add it to user PATH environment")


def parse_fasta(filename):
	"""
	Imports specified .fasta file
	"""
	f = open(filename)
	sequences = {}
	for line in f:
		if line.startswith('>'):
			name = line[1:].strip()
			sequences[name] = ''
		else:
			sequences[name] = sequences[name] + line.strip()
	f.close()
	return sequences

def verify_bai(bamfile):
	'''
	Function that indexes bam file from input if missing
	'''
	bai_file = bamfile + ".bai"
	if(not os.path.exists(bai_file)):
		pysam.index(bamfile)

# Helper function for chunking the bam
def split_chunk_file(one_barcode_file, script_dir, input, bcbd, barcode_tag, mito_chr, umi_barcode):
	chunk_bam_py = script_dir + "/bin/python/chunk_barcoded_bam.py"
	pycall = " ".join(['python', chunk_bam_py, input, bcbd, barcode_tag, one_barcode_file, mito_chr, umi_barcode])
	os.system(pycall)

def verify_sample_mitobam(bam, mito_chr, mito_length):
	idxs = pysam.idxstats(bam).split("\n")
	nReads = 0
	bam_length = 0
	
	# Pull out essentials from idxstats
	for i in idxs:
		if(i.split("\t")[0] == mito_chr):
			bam_length = int(i.split("\t")[1])
			nReads = int(i.split("\t")[2])
			
	# If missing for some reason and don't want to check
	if(mito_length == -9):
		mito_length = bam_length
	return(bam_length == mito_length and nReads > 0)


def handle_fasta_inference(mito_genome, supported_genomes, script_dir, mode, of, write_files = True):
	"""
	Determines what's going on with the mitochondrial genome
	based on user input / existing data
	"""
	if any(mito_genome in s for s in supported_genomes):
		fastaf = script_dir + "/bin/anno/fasta/" + mito_genome + ".fasta"
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
	
	if(write_files):
		make_folder(of + "/fasta/")
		make_folder(of + "/final/")
	
	newfastaf = of + "/fasta/" + mito_genome + ".fasta"
	
	# Need a special logic to potentially over-write the existing 
	# mito fasta file if they wind up not being the same from guess/checking hg19
	writeFA = False
	if os.path.exists(newfastaf):
		if filecmp.cmp(fastaf, newfastaf,shallow=False):
			pass
		else:
			writeFA = True
	else:
		writeFA = True
		
	if writeFA and write_files:
			
		shutil.copyfile(fastaf, newfastaf)
		fastaf = newfastaf
		pysam.faidx(fastaf)
	
		with open(of + "/final/" + mito_genome + "_refAllele.txt", 'w') as f:
			b = 1
			for base in mito_seq:
				f.write(str(b) + "\t" + base + "\n")
				b += 1
			f.close()
	return(fastaf, mito_genome, mito_length)
	

def make_folder(folder):
	"""
	Function to only make a given folder if it does not already exist
	"""
	if not os.path.exists(folder):
		os.makedirs(folder)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def split_barcodes_file(barcode_file, nsamples, output):
	"""
	Function to only make a given folder if it does not already exist
	"""
	n_samples_observed = file_len(barcode_file)
	
	# See if we need to do anything... if user didn't specify or we have relatively few samples, just return the barcode file back
	if(n_samples_observed < nsamples or nsamples == 0):
		return([barcode_file])
	else:
		# Need to split files into a maximum of n samples
		total_files = math.ceil(n_samples_observed / nsamples)
		lines_per_file = nsamples
		
		# Set up a temporary output folder to route files to
		full_output_folder = output + "/temp" + "/barcode_files"
		make_folder(full_output_folder)
		smallfile = None
		counter = 0
		with open(barcode_file) as bigfile:
			for lineno, line in enumerate(bigfile):
				if lineno % lines_per_file == 0:
					if smallfile:
						smallfile.close()
					counter = counter + 1
					small_filename = full_output_folder + "/barcodes." + str(counter) + ".txt" 
					smallfile = open(small_filename, "w")
				smallfile.write(line)
		if smallfile:
			smallfile.close()
		barcodes_files = [full_output_folder + "/barcodes." + str(x) + ".txt" for x in list(range(1, total_files + 1))]
		return(barcodes_files)



# https://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-using-python	
def available_cpu_count():
    """
    Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program
	"""
    # cpuset
    # cpuset may restrict the number of *available* processors
    try:
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
                      open('/proc/self/status').read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    # http://code.google.com/p/psutil/
    try:
        import psutil
        return psutil.cpu_count()   # psutil.NUM_CPUS on old versions
    except (ImportError, AttributeError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])

        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0:
            return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                  stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)

        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        res = 0
        for pd in pseudoDevices:
            if re.match(r'^cpuid@[0-9]+$', pd):
                res += 1

        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0:
            return res
    except OSError:
        pass

    raise Exception('Can not determine number of CPUs on this system')
    
