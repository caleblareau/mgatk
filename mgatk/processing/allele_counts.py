import pysam


def _write_sparse_matrix(outpre, mid, sample, max_bp, vec):
    with open(outpre + "." + mid + ".txt", "w") as V:
        for i in range(0, int(max_bp)):
            if vec[i] > 0:
                V.write(str(i + 1) + "," + sample + "," + str(vec[i]) + "\n")


def _write_sparse_matrix2(outpre, mid, sample, max_bp, vec1, vec2):
    with open(outpre + "." + mid + ".txt", "w") as V:
        for i in range(0, int(max_bp)):
            if vec1[i] > 0 or vec2[i] > 0:
                V.write(str(i + 1) + "," + sample + "," + str(vec1[i]) + "," + str(vec2[i]) + "\n")


def _write_sparse_matrix4(outpre, mid, sample, max_bp, vec1, vec2, vec3, vec4):
    with open(outpre + "." + mid + ".txt", "w") as V:
        for i in range(0, int(max_bp)):
            if vec1[i] > 0 or vec3[i] > 0:
                V.write(
                    str(i + 1) + "," + sample + "," + str(vec1[i]) + "," + str(vec2[i]) + ","
                    + str(vec3[i]) + "," + str(vec4[i]) + "\n"
                )


def sumstats_bp(bamfile, outpre, mito_genome, max_bp, base_qual, sample, fasta_file,
                alignment_quality, emit_base_qualities):
    """Per-position, per-strand allele counts (and optional mean base qualities) for a bulk/bcall sample."""
    max_bp = int(max_bp)
    base_qual = float(base_qual)
    alignment_quality = float(alignment_quality)

    n = max_bp
    counts_a_fw = [0.00000001] * n
    counts_c_fw = [0.00000001] * n
    counts_g_fw = [0.00000001] * n
    counts_t_fw = [0.00000001] * n

    qual_a_fw = [0.0] * n
    qual_c_fw = [0.0] * n
    qual_g_fw = [0.0] * n
    qual_t_fw = [0.0] * n

    counts_a_rev = [0.00000001] * n
    counts_c_rev = [0.00000001] * n
    counts_g_rev = [0.00000001] * n
    counts_t_rev = [0.00000001] * n

    qual_a_rev = [0.0] * n
    qual_c_rev = [0.0] * n
    qual_g_rev = [0.0] * n
    qual_t_rev = [0.0] * n

    bam2 = pysam.AlignmentFile(bamfile, "rb")
    for read in bam2:
        seq = read.seq
        reverse = read.is_reverse
        quality = read.query_qualities
        align_qual_read = read.mapping_quality
        for qpos, refpos in read.get_aligned_pairs(True):
            if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
                if seq[qpos] == "A" and quality[qpos] > base_qual:
                    if not reverse:
                        qual_a_fw[refpos] += quality[qpos]
                        counts_a_fw[refpos] += 1
                    else:
                        qual_a_rev[refpos] += quality[qpos]
                        counts_a_rev[refpos] += 1
                elif seq[qpos] == "C" and quality[qpos] > base_qual:
                    if not reverse:
                        qual_c_fw[refpos] += quality[qpos]
                        counts_c_fw[refpos] += 1
                    else:
                        qual_c_rev[refpos] += quality[qpos]
                        counts_c_rev[refpos] += 1
                elif seq[qpos] == "G" and quality[qpos] > base_qual:
                    if not reverse:
                        qual_g_fw[refpos] += quality[qpos]
                        counts_g_fw[refpos] += 1
                    else:
                        qual_g_rev[refpos] += quality[qpos]
                        counts_g_rev[refpos] += 1
                elif seq[qpos] == "T" and quality[qpos] > base_qual:
                    if not reverse:
                        qual_t_fw[refpos] += quality[qpos]
                        counts_t_fw[refpos] += 1
                    else:
                        qual_t_rev[refpos] += quality[qpos]
                        counts_t_rev[refpos] += 1

    mean_qual_a_fw = [round(x / y, 1) for x, y in zip(qual_a_fw, counts_a_fw)]
    mean_qual_c_fw = [round(x / y, 1) for x, y in zip(qual_c_fw, counts_c_fw)]
    mean_qual_g_fw = [round(x / y, 1) for x, y in zip(qual_g_fw, counts_g_fw)]
    mean_qual_t_fw = [round(x / y, 1) for x, y in zip(qual_t_fw, counts_t_fw)]

    counts_a_fw = [int(round(elem)) for elem in counts_a_fw]
    counts_c_fw = [int(round(elem)) for elem in counts_c_fw]
    counts_g_fw = [int(round(elem)) for elem in counts_g_fw]
    counts_t_fw = [int(round(elem)) for elem in counts_t_fw]

    mean_qual_a_rev = [round(x / y, 1) for x, y in zip(qual_a_rev, counts_a_rev)]
    mean_qual_c_rev = [round(x / y, 1) for x, y in zip(qual_c_rev, counts_c_rev)]
    mean_qual_g_rev = [round(x / y, 1) for x, y in zip(qual_g_rev, counts_g_rev)]
    mean_qual_t_rev = [round(x / y, 1) for x, y in zip(qual_t_rev, counts_t_rev)]

    counts_a_rev = [int(round(elem)) for elem in counts_a_rev]
    counts_c_rev = [int(round(elem)) for elem in counts_c_rev]
    counts_g_rev = [int(round(elem)) for elem in counts_g_rev]
    counts_t_rev = [int(round(elem)) for elem in counts_t_rev]

    if emit_base_qualities == "True":
        _write_sparse_matrix4(outpre, "A", sample, max_bp, counts_a_fw, mean_qual_a_fw, counts_a_rev, mean_qual_a_rev)
        _write_sparse_matrix4(outpre, "C", sample, max_bp, counts_c_fw, mean_qual_c_fw, counts_c_rev, mean_qual_c_rev)
        _write_sparse_matrix4(outpre, "G", sample, max_bp, counts_g_fw, mean_qual_g_fw, counts_g_rev, mean_qual_g_rev)
        _write_sparse_matrix4(outpre, "T", sample, max_bp, counts_t_fw, mean_qual_t_fw, counts_t_rev, mean_qual_t_rev)
    else:
        _write_sparse_matrix2(outpre, "A", sample, max_bp, counts_a_fw, counts_a_rev)
        _write_sparse_matrix2(outpre, "C", sample, max_bp, counts_c_fw, counts_c_rev)
        _write_sparse_matrix2(outpre, "G", sample, max_bp, counts_g_fw, counts_g_rev)
        _write_sparse_matrix2(outpre, "T", sample, max_bp, counts_t_fw, counts_t_rev)

    zipped_list = zip(
        list(counts_a_fw), list(counts_c_fw), list(counts_g_fw), list(counts_t_fw),
        list(counts_a_rev), list(counts_c_rev), list(counts_g_rev), list(counts_t_rev),
    )
    sums = [sum(item) for item in zipped_list]
    _write_sparse_matrix(outpre, "coverage", sample, max_bp, sums)


def sumstats_bp_overlap(bamfile, outpre, mito_genome, max_bp, base_qual, sample, fasta_file,
                        alignment_quality, emit_base_qualities):
    """Same as sumstats_bp, but only counts each base once in overlapping mate-pair regions."""
    import numpy as np
    from collections import defaultdict

    max_bp = int(max_bp)
    base_qual = float(base_qual)
    alignment_quality = float(alignment_quality)

    n = max_bp
    counts_a_fw = [0.00000001] * n
    counts_c_fw = [0.00000001] * n
    counts_g_fw = [0.00000001] * n
    counts_t_fw = [0.00000001] * n

    qual_a_fw = [0.0] * n
    qual_c_fw = [0.0] * n
    qual_g_fw = [0.0] * n
    qual_t_fw = [0.0] * n

    counts_a_rev = [0.00000001] * n
    counts_c_rev = [0.00000001] * n
    counts_g_rev = [0.00000001] * n
    counts_t_rev = [0.00000001] * n

    qual_a_rev = [0.0] * n
    qual_c_rev = [0.0] * n
    qual_g_rev = [0.0] * n
    qual_t_rev = [0.0] * n

    # organize reads into a dict where key is readname
    bam2 = [x for x in pysam.AlignmentFile(bamfile, "rb")]
    ordered_bam2 = defaultdict(list)
    for read in bam2:
        ordered_bam2[read.query_name].append(read)

    for read_name in ordered_bam2:
        # disregard singlets and multiplets
        if len(ordered_bam2[read_name]) != 2:
            continue

        # identify fwd and rev in a pair
        read0, read1 = ordered_bam2[read_name]
        if read0.is_reverse and not read1.is_reverse:
            fwd_read, rev_read = read1, read0
        elif not read0.is_reverse and read1.is_reverse:
            fwd_read, rev_read = read0, read1
        else:
            # disregard a pair if both are the same strand
            continue

        # gather what we need
        fwd_seq, rev_seq = fwd_read.query_sequence, rev_read.query_sequence
        fwd_quality, rev_quality = np.array(fwd_read.query_qualities), np.array(rev_read.query_qualities)
        fwd_align_qual_read, rev_align_qual_read = fwd_read.mapping_quality, rev_read.mapping_quality

        # check alignment quality
        if fwd_align_qual_read > alignment_quality and rev_align_qual_read > alignment_quality:
            # partition the pair into fwd-only, overlap, and rev-only
            overlap_length = fwd_read.get_overlap(rev_read.reference_start, rev_read.reference_end)
            if overlap_length == 0:
                fwd_use_idx = np.arange(len(fwd_seq))
                rev_use_idx = np.arange(len(rev_seq))
            else:
                # choose which strand to use in the overlap region based on quality score
                fwd_only_end = len(fwd_seq) - overlap_length
                rev_only_start = overlap_length
                fwd_overlap_quality = fwd_quality[fwd_only_end:]
                rev_overlap_quality = rev_quality[:rev_only_start]
                fwd_overlap_use_idx = np.where(fwd_overlap_quality > rev_overlap_quality)[0]
                rev_overlap_use_idx = np.where(fwd_overlap_quality < rev_overlap_quality)[0]

                # evenly assign bases with equal quality in overlap region
                equal_overlap_idx = np.where(fwd_overlap_quality == rev_overlap_quality)[0]
                equal_split = int(np.floor(len(equal_overlap_idx) / 2))
                fwd_overlap_use_idx = np.concatenate([fwd_overlap_use_idx, equal_overlap_idx[:equal_split]])
                rev_overlap_use_idx = np.concatenate([rev_overlap_use_idx, equal_overlap_idx[equal_split:]])

                # merge the exclusive region and use idx in overlap region
                fwd_use_idx = np.concatenate([np.arange(fwd_only_end), fwd_overlap_use_idx + fwd_only_end])
                rev_use_idx = np.concatenate([rev_overlap_use_idx, np.arange(rev_only_start, len(rev_seq))])

        elif fwd_align_qual_read <= alignment_quality and rev_align_qual_read <= alignment_quality:
            # use none for either
            fwd_use_idx = np.array([])
            rev_use_idx = np.array([])

        elif fwd_align_qual_read > alignment_quality and rev_align_qual_read <= alignment_quality:
            # use none of rev and all of fwd
            fwd_use_idx = np.arange(len(fwd_seq))
            rev_use_idx = np.array([])

        elif fwd_align_qual_read <= alignment_quality and rev_align_qual_read > alignment_quality:
            # use all of rev and none of fwd
            fwd_use_idx = np.array([])
            rev_use_idx = np.arange(len(rev_seq))

        # handle fwd region
        fwd_aligned_pairs = fwd_read.get_aligned_pairs(True)
        fwd_region = [pair for pair in fwd_aligned_pairs if pair[0] in fwd_use_idx]
        for qpos, refpos in fwd_region:
            if refpos is not None and fwd_quality[qpos] > base_qual:
                if fwd_seq[qpos] == "A":
                    qual_a_fw[refpos] += fwd_quality[qpos]
                    counts_a_fw[refpos] += 1
                elif fwd_seq[qpos] == "C":
                    qual_c_fw[refpos] += fwd_quality[qpos]
                    counts_c_fw[refpos] += 1
                elif fwd_seq[qpos] == "G":
                    qual_g_fw[refpos] += fwd_quality[qpos]
                    counts_g_fw[refpos] += 1
                elif fwd_seq[qpos] == "T":
                    qual_t_fw[refpos] += fwd_quality[qpos]
                    counts_t_fw[refpos] += 1

        # handle rev region
        rev_aligned_pairs = rev_read.get_aligned_pairs(True)
        rev_region = [pair for pair in rev_aligned_pairs if pair[0] in rev_use_idx]
        for qpos, refpos in rev_region:
            if refpos is not None and rev_quality[qpos] > base_qual:
                if rev_seq[qpos] == "A":
                    qual_a_rev[refpos] += rev_quality[qpos]
                    counts_a_rev[refpos] += 1
                elif rev_seq[qpos] == "C":
                    qual_c_rev[refpos] += rev_quality[qpos]
                    counts_c_rev[refpos] += 1
                elif rev_seq[qpos] == "G":
                    qual_g_rev[refpos] += rev_quality[qpos]
                    counts_g_rev[refpos] += 1
                elif rev_seq[qpos] == "T":
                    qual_t_rev[refpos] += rev_quality[qpos]
                    counts_t_rev[refpos] += 1

    mean_qual_a_fw = [round(x / y, 1) for x, y in zip(qual_a_fw, counts_a_fw)]
    mean_qual_c_fw = [round(x / y, 1) for x, y in zip(qual_c_fw, counts_c_fw)]
    mean_qual_g_fw = [round(x / y, 1) for x, y in zip(qual_g_fw, counts_g_fw)]
    mean_qual_t_fw = [round(x / y, 1) for x, y in zip(qual_t_fw, counts_t_fw)]

    counts_a_fw = [int(round(elem)) for elem in counts_a_fw]
    counts_c_fw = [int(round(elem)) for elem in counts_c_fw]
    counts_g_fw = [int(round(elem)) for elem in counts_g_fw]
    counts_t_fw = [int(round(elem)) for elem in counts_t_fw]

    mean_qual_a_rev = [round(x / y, 1) for x, y in zip(qual_a_rev, counts_a_rev)]
    mean_qual_c_rev = [round(x / y, 1) for x, y in zip(qual_c_rev, counts_c_rev)]
    mean_qual_g_rev = [round(x / y, 1) for x, y in zip(qual_g_rev, counts_g_rev)]
    mean_qual_t_rev = [round(x / y, 1) for x, y in zip(qual_t_rev, counts_t_rev)]

    counts_a_rev = [int(round(elem)) for elem in counts_a_rev]
    counts_c_rev = [int(round(elem)) for elem in counts_c_rev]
    counts_g_rev = [int(round(elem)) for elem in counts_g_rev]
    counts_t_rev = [int(round(elem)) for elem in counts_t_rev]

    if emit_base_qualities == "True":
        _write_sparse_matrix4(outpre, "A", sample, max_bp, counts_a_fw, mean_qual_a_fw, counts_a_rev, mean_qual_a_rev)
        _write_sparse_matrix4(outpre, "C", sample, max_bp, counts_c_fw, mean_qual_c_fw, counts_c_rev, mean_qual_c_rev)
        _write_sparse_matrix4(outpre, "G", sample, max_bp, counts_g_fw, mean_qual_g_fw, counts_g_rev, mean_qual_g_rev)
        _write_sparse_matrix4(outpre, "T", sample, max_bp, counts_t_fw, mean_qual_t_fw, counts_t_rev, mean_qual_t_rev)
    else:
        _write_sparse_matrix2(outpre, "A", sample, max_bp, counts_a_fw, counts_a_rev)
        _write_sparse_matrix2(outpre, "C", sample, max_bp, counts_c_fw, counts_c_rev)
        _write_sparse_matrix2(outpre, "G", sample, max_bp, counts_g_fw, counts_g_rev)
        _write_sparse_matrix2(outpre, "T", sample, max_bp, counts_t_fw, counts_t_rev)

    zipped_list = zip(
        list(counts_a_fw), list(counts_c_fw), list(counts_g_fw), list(counts_t_fw),
        list(counts_a_rev), list(counts_c_rev), list(counts_g_rev), list(counts_t_rev),
    )
    sums = [sum(item) for item in zipped_list]
    _write_sparse_matrix(outpre, "coverage", sample, max_bp, sums)
