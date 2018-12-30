#!/usr/bin/python

import sys
import pysam
from itertools import groupby

def ranges(lst):
    '''
    Yield range of consecutive numbers
    '''
    pos = (j - i for i, j in enumerate(lst))
    t = 0
    for i, els in groupby(pos):
        l = len(list(els))
        el = lst[t]
        t += l
        yield [el, el+l]

def strand(read):
	'''
	Yield strand for a given read
	'''
	if(read.is_reverse):
		return "-"
	else:
		return("+")

def assemble_insertions(positions, ref_pos, seq, strand):
	'''
	Assembles insert positions based on index position and read + genomic positions
	'''
	
	# Handle if at the left end of read
	if(positions[0][0] == 0):
		del positions[0]
	
	# Return nothing if we find nothing
	if(len(positions) == 0):
		return([])
	
	# Reference start position of insert
	starts = [ref_pos[x[0]-1] for i,x in enumerate(positions)]
	
	# Length of insert
	lens = [x[1] - x[0] for i,x in enumerate(positions) if x[1] < len(ref_pos)] # guards against far right end matching
	
	# Index starting position of sequence + the full thing
	idxs = [ref_pos.index(x) for i,x in enumerate(starts)]
	inserts = [str(starts[i]) + "_" + seq[(idxs[i]-1):(idxs[i] + lens[i])] + "_" + strand for i,x in enumerate(lens)]
	return(inserts)

bamfile = sys.argv[1]
bamFP = pysam.Samfile(bamfile, "rb")

def missing_elements(L):
	'''
	Function that given a sequence of integars plus the start / end
	values, returns the values that are missing
	'''
	start, end = L[0], L[-1]
	return sorted(set(range(start, end + 1)).difference(L))

insertions = []
deletions = []
for read in bamFP:
	if( not( read.is_unmapped ) ):   #if it's mapped
		cigarLine=read.cigar
		for (cigarType,cigarLength) in cigarLine:
			if(cigarType == 0):
				pass
			elif(cigarType == 1): #insertions
			
				# get read bp where there is no reference alignment
				ref = read.get_reference_positions(full_length=True)
				g = list(ranges([i for i,x in enumerate(ref) if x == None]))
				insertions.append(assemble_insertions(g, ref, read.query_alignment_sequence, strand(read)))
				#print(g)
			elif(cigarType == 2): #deletion
				# get reference genome positions where missing
				g = list(ranges(missing_elements(read.get_reference_positions(full_length=True))))
				#print(g)
			else:
				pass
bamFP.close()
insertions = [item for sublist in insertions for item in sublist]
print(insertions)