import sys
import re
import os
import pysam
import gzip

bamfile = sys.argv[1]
outfile = bamfile.rsplit( ".", 1 )[ 0 ] + ".alleleDetails.csv.gz"

samfile = pysam.AlignmentFile(bamfile, "rb" )
out = gzip.open(outfile, 'wb')
for pileupcolumn in samfile.pileup("chrM"):
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
            out.write(('%s, %s, %s, %s, %s\n' %
                  (pileupread.alignment.query_name,
                   pileupcolumn.pos +1,
                   pileupread.alignment.query_sequence[pileupread.query_position],
                   pileupread.alignment.query_qualities[pileupread.query_position],
                   pileupread.alignment.is_reverse)).encode('utf-8'))
samfile.close()
out.close()
