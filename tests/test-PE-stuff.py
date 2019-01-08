from collections import defaultdict
import pysam

filename="humanbam/MGH60-P6-A11.mito.bam"

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    fail_read = 0
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

bam = pysam.AlignmentFile(filename, 'rb')
for read1, read2 in read_pair_generator(bam):
	if(read1 != None and read2 != None):
		print(read1.query_name, read1.reference_start, read2.reference_start, read1.is_reverse, read2.is_reverse)
