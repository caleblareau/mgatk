import os

import pysam
from dedup import picardlike


def dedup_bam(temp_bam, output_bam, rmlog, remove_duplicates, barcode_tag=""):
    """
    Either deduplicate temp_bam -> output_bam (duplicate records dropped,
    temp_bam assumed already sorted) or just promote temp_bam to output_bam
    unchanged, then index output_bam.

    barcode_tag, if non-empty, groups duplicate detection by that read tag's
    value (e.g. a cell barcode), so duplicates are only merged within reads
    sharing the same tag value -- see dedup.picardlike.mark_duplicates for why
    this is a deliberate design choice rather than Picard's own BARCODE_TAG.
    """
    if remove_duplicates == "True":
        dup_idx, _ = picardlike.deduplicate(
            temp_bam, output_bam, barcode_tag=barcode_tag if barcode_tag else None,
        )
        with pysam.AlignmentFile(temp_bam, "rb") as bam:
            total = sum(1 for _ in bam.fetch(until_eof=True))
        with open(rmlog, "w") as f:
            f.write("total_reads\t{}\n".format(total))
            f.write("duplicate_reads\t{}\n".format(len(dup_idx)))
    else:
        os.replace(temp_bam, output_bam)
        os.remove(temp_bam + ".bai")
    pysam.index(output_bam)
