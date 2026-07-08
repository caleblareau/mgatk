import pysam


def _filter_read_tags(intags, nh_max, nm_max):
    """Reject reads whose NH/NM(nM) tags exceed the configured maxima."""
    for tg in intags:
        if ("NH" == tg[0] and int(tg[1]) > int(nh_max)) or (
            ("NM" == tg[0] or "nM" == tg[0]) and int(tg[1]) > int(nm_max)
        ):
            return False
    return True


def filter_clip_bam(bamfile, out_bam, logfile, mtchr, proper_pair, nh_max, nm_max):
    """
    Filter a bam down to reads on mtchr passing NH/NM tag and (optionally)
    proper-pairing thresholds, writing the result to out_bam and a
    keep/remove count summary to logfile.
    """
    bam = pysam.AlignmentFile(bamfile, "rb")

    # Modify the header to account for CRA v2 nonsense
    # https://github.com/pysam-developers/pysam/issues/509
    new_header = str(bam.header).replace("@HD\tSO:coordinate", "@HD\tVN:1.5\tSO:coordinate")
    out = pysam.AlignmentFile(out_bam, "wb", text=new_header)

    def pairing(read):
        if proper_pair != "True":  # then user doesn't care to filter it
            return True
        return read.is_proper_pair

    keep_count = 0
    filt_count = 0
    for read in bam:
        if _filter_read_tags(read.tags, nh_max, nm_max) and read.reference_name == mtchr and pairing(read):
            keep_count += 1
            out.write(read)
        else:
            filt_count += 1

    bam.close()
    out.close()

    with open(logfile, "w") as outfile:
        outfile.write("Kept " + str(keep_count) + "\n" + "Removed " + str(filt_count) + "\n")
