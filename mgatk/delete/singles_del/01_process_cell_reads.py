import re
import os
import pysam
from collections import Counter
import sys
from optparse import OptionParser
# Parse out data
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process 1 bam file with specific parameters defined"
opts = OptionParser(usage=usage)
opts.add_option("-i", "--input", help="input single cell bam file")
opts.add_option("-o", "--output", help="output filepath")
options, arguments = opts.parse_args()
inbam = options.input
output_file = options.output
bam_in = pysam.AlignmentFile(inbam, "rb")
def process_cigar_for_clip_position(cigar, tuple):
    pos = 0
    # Case 1/2: start of read
    if(cigar[1] == "H" or cigar[2] == "H" or
       cigar[1] == "S" or cigar[2] == "S"):
        pos = tuple[0][1] + 1 #offset 1-base and differentiate 0 pos/noclip
    # Case 2: end of read
    if(cigar[-1] == "H" or cigar[-1] == "S"):
        pos = tuple[-1][1] + 1 #offset 1-base and differentiate 0 pos/noclip
    return(pos)
def left_clip(cigar):
    if(cigar[1] == "H" or cigar[2] == "H" or
       cigar[1] == "S" or cigar[2] == "S"):
        return(1)
    else:
        return(0)
def right_clip(cigar):
    if(cigar[-1] == "S" or cigar[-1] == "H"):
        return(1)
    else:
        return(0)
# fuzzy match lev distance
def lev(a, b):
    if not a: return len(b)
    if not b: return len(a)
    return min(lev(a[1:], b[1:])+(a[0] != b[0]), lev(a[1:], b)+1, lev(a, b[1:])+1)
def get_clipped_string(cigar, tuple, seq):
    barcode_5mer="CTGTC"
    barcode_5mer_rc="GACAG"
    # Case 1/2: start of read
    if(cigar[1] == "H" or cigar[2] == "H" or
       cigar[1] == "S" or cigar[2] == "S"):
        pos = tuple[0][0]
        seq_out = seq[:pos][-4:]
        return(lev(barcode_5mer, seq_out) <=1 or
               lev(barcode_5mer_rc, seq_out) <=1 )
    # Case 2: end of read
    if(cigar[-1] == "H" or cigar[-1] == "S"):
        pos = tuple[-1][0]
        seq_out = seq[pos + 1:][:4]
        return(lev(barcode_5mer, seq_out) <=1 or
               lev(barcode_5mer_rc, seq_out) <=1 )
    return(False)
def get_n_clipped(cigar):
    out = "0"
    # Case 1/2: start of read
    if(cigar[1] == "H" or cigar[2] == "H" or
       cigar[1] == "S" or cigar[2] == "S"):
         out = re.split('[A-Z]', cigar)[0]
    # Case 2: end of read
    if(cigar[-1] == "H" or cigar[-1] == "S"):
        out = re.split('[A-Z]', cigar)[-2]
    return(out)

# Open all the output files and spit out the filtered data
# Based on where the matching value originates
clip_pos_count = Counter()
outfile_handle = open(output_file, 'w')
for read in bam_in:
    seq = read.seq
    reverse = read.is_reverse
    cigar_string = read.cigarstring
    positions = read.get_reference_positions()
    tuple = read.get_aligned_pairs(True)
    if(positions and cigar_string and seq):
        start = str(positions[0] + 1) #offset to 1 base
        end = str(positions[-1] + 1) #offset to 1 base
        rc = str(right_clip(cigar_string))
        lc = str(left_clip(cigar_string))
        clip_pos = str(process_cigar_for_clip_position(cigar_string, tuple))
        is_barcode = str(get_clipped_string(cigar_string, tuple, seq))
        clip_pos_count[clip_pos] += 1
        read_name = str(read.query_name)
        n_clipped = get_n_clipped(cigar_string)
        list_of_outs = [start, end, lc, rc, clip_pos, read_name, is_barcode, n_clipped]
        outfile_handle.write("\t".join(list_of_outs) + "\n")
bam_in.close()