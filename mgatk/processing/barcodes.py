import os
from contextlib import contextmanager

import pysam


def _get_barcode(intags, barcode_tag):
    for tg in intags:
        if barcode_tag == tg[0]:
            return tg[1]
    return "NA"


def quantify_barcodes(bamfile, barcode_tag, min_barcodes, mtchr, quant_file, passing_file):
    """
    Count reads per barcode at the mitochondrial locus and write out
    1) a full quantification table and 2) the barcodes passing min_barcodes.
    """
    barcodes_all = dict()
    bam = pysam.AlignmentFile(bamfile, "rb")
    for read in bam.fetch(str(mtchr), multiple_iterators=False):
        read_barcode = _get_barcode(read.tags, barcode_tag)
        barcodes_all[read_barcode] = barcodes_all.get(read_barcode, 0) + 1
    bam.close()

    barcodes = {
        bc: count for bc, count in barcodes_all.items()
        if count >= min_barcodes and bc != "NA"
    }

    with open(quant_file, "w") as quant_file_o:
        for bc, count in barcodes.items():
            quant_file_o.write(bc + "," + str(count) + "\n")

    with open(passing_file, "w") as passing_file_o:
        for bc in barcodes:
            passing_file_o.write(bc + "\n")

    return barcodes


@contextmanager
def _multi_file_manager(bamfile, files):
    """Open multiple output bams (templated off bamfile) and ensure they all get closed."""
    temp = pysam.AlignmentFile(bamfile, "rb")
    handles = [pysam.AlignmentFile(f, "wb", template=temp) for f in files]
    temp.close()
    try:
        yield handles
    finally:
        for handle in handles:
            handle.close()


def split_barcoded_bam(bamfile, outfolder, barcode_tag, bcfile, mtchr):
    """Split a bam into one file per barcode listed in bcfile."""
    with open(bcfile) as barcode_file_handle:
        bc = [x.strip() for x in barcode_file_handle.readlines()]

    bambcfiles = [os.path.join(outfolder, bc1 + ".bam") for bc1 in bc]
    bc_dict = {bc1: i for i, bc1 in enumerate(bc)}

    with _multi_file_manager(bamfile, bambcfiles) as fopen:
        bam = pysam.AlignmentFile(bamfile, "rb")
        for read in bam.fetch(str(mtchr), multiple_iterators=False):
            read_barcode = _get_barcode(read.tags, barcode_tag)
            if read_barcode in bc_dict:
                fopen[bc_dict[read_barcode]].write(read)
        bam.close()


def chunk_barcoded_bam(bamfile, outfolder, barcode_tag, bcfile, mtchr, umitag):
    """
    Write reads matching a whitelist of barcodes into a single chunked bam,
    tagging each read with a synthetic all-ACGT "MU" barcode (cell barcode +
    UMI + a fake per-channel suffix) so downstream dedup tools like Picard
    are happy with the tag contents.
    """
    basename = os.path.splitext(os.path.basename(bcfile))[0]

    with open(bcfile) as barcode_file_handle:
        bc = set(x.strip() for x in barcode_file_handle.readlines())

    bam = pysam.AlignmentFile(bamfile, "rb")
    outname = os.path.join(outfolder, basename + ".bam")
    out = pysam.AlignmentFile(outname, "wb", template=bam)

    bases = "ACGT"
    fauxdon = [a + b + c + d for a in bases for b in bases for c in bases for d in bases]

    def _get_tag_or_fallback(read, tag):
        try:
            return read.get_tag(tag)
        except Exception:
            return "AA"

    try:
        for read in bam.fetch(str(mtchr), multiple_iterators=False):
            barcode_id = _get_tag_or_fallback(read, barcode_tag)

            if barcode_id in bc:
                if umitag != "XX":
                    umi_id = _get_tag_or_fallback(read, umitag)
                else:
                    umi_id = ""

                # Make a fake UMI from 1) cell barcode + 2) captured umi + 3) experiment
                # all with just ACGTs so that picard doesn't bark at us.
                # Only do this if the last string element is a number (i.e channel in 10x convention)
                if barcode_id[-1].isnumeric():
                    split_two = barcode_id.split("-")
                    faux_umi = split_two[0] + umi_id + fauxdon[int(split_two[1]) - 1]
                else:
                    faux_umi = barcode_id + umi_id
                read.tags = read.tags + [("MU", faux_umi)]
                out.write(read)
    except OSError:  # Truncated bam file from previous iteration handle
        pass

    bam.close()
    out.close()
    pysam.index(outname)
