import re

import pysam


def _process_cigar_for_clip_position(cigar, aligned_pairs):
    pos = 0
    # Case 1/2: start of read
    if cigar[1] == "H" or cigar[2] == "H" or cigar[1] == "S" or cigar[2] == "S":
        pos = aligned_pairs[0][1] + 1  # offset 1-base and differentiate 0 pos/noclip
    # Case 2: end of read
    if cigar[-1] == "H" or cigar[-1] == "S":
        pos = aligned_pairs[-1][1] + 1  # offset 1-base and differentiate 0 pos/noclip
    return pos


def _left_clip(cigar):
    if cigar[1] == "H" or cigar[2] == "H" or cigar[1] == "S" or cigar[2] == "S":
        return 1
    return 0


def _right_clip(cigar):
    if cigar[-1] == "S" or cigar[-1] == "H":
        return 1
    return 0


def _lev(a, b):
    """Fuzzy match (Levenshtein) distance."""
    if not a:
        return len(b)
    if not b:
        return len(a)
    return min(_lev(a[1:], b[1:]) + (a[0] != b[0]), _lev(a[1:], b) + 1, _lev(a, b[1:]) + 1)


def _get_clipped_string(cigar, aligned_pairs, seq):
    barcode_5mer = "CTGTC"
    barcode_5mer_rc = "GACAG"
    # Case 1/2: start of read
    if cigar[1] == "H" or cigar[2] == "H" or cigar[1] == "S" or cigar[2] == "S":
        pos = aligned_pairs[0][0]
        seq_out = seq[:pos][-4:]
        return _lev(barcode_5mer, seq_out) <= 1 or _lev(barcode_5mer_rc, seq_out) <= 1
    # Case 2: end of read
    if cigar[-1] == "H" or cigar[-1] == "S":
        pos = aligned_pairs[-1][0]
        seq_out = seq[pos + 1:][:4]
        return _lev(barcode_5mer, seq_out) <= 1 or _lev(barcode_5mer_rc, seq_out) <= 1
    return False


def _get_n_clipped(cigar):
    out = "0"
    # Case 1/2: start of read
    if cigar[1] == "H" or cigar[2] == "H" or cigar[1] == "S" or cigar[2] == "S":
        out = re.split("[A-Z]", cigar)[0]
    # Case 2: end of read
    if cigar[-1] == "H" or cigar[-1] == "S":
        out = re.split("[A-Z]", cigar)[-2]
    return out


def process_cell_reads(inbam, output_file):
    """
    Scan a single-cell bam for soft/hard-clipped reads and record, per read,
    the clip position/side plus a fuzzy check for a nearby deletion-junction
    barcode motif. Used by mgatk-del to build per-cell deletion evidence.
    """
    bam_in = pysam.AlignmentFile(inbam, "rb")
    with open(output_file, "w") as outfile_handle:
        for read in bam_in:
            seq = read.seq
            cigar_string = read.cigarstring
            positions = read.get_reference_positions()
            aligned_pairs = read.get_aligned_pairs(True)
            if positions and cigar_string and seq:
                start = str(positions[0] + 1)  # offset to 1 base
                end = str(positions[-1] + 1)  # offset to 1 base
                rc = str(_right_clip(cigar_string))
                lc = str(_left_clip(cigar_string))
                clip_pos = str(_process_cigar_for_clip_position(cigar_string, aligned_pairs))
                is_barcode = str(_get_clipped_string(cigar_string, aligned_pairs, seq))
                read_name = str(read.query_name)
                n_clipped = _get_n_clipped(cigar_string)
                list_of_outs = [start, end, lc, rc, clip_pos, read_name, is_barcode, n_clipped]
                outfile_handle.write("\t".join(list_of_outs) + "\n")
    bam_in.close()
