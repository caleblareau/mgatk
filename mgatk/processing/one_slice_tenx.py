import pysam

from .allele_counts_tenx import sumstats_bp_tenx, sumstats_bp_tenx_overlap
from .dedup import dedup_bam
from .filtering import filter_clip_bam


def process_one_slice_tenx(config, input_bam, sample):
    """
    tenx per-slice pipeline: filter -> sort -> (optional) dedup -> per-cell
    allele counting. `config` is the parsed mgatk tenx config (the same dict
    Snakemake loads via `configfile:` for `tenx` runs).
    """
    outdir = config["output_directory"]

    mito_genome = config["mito_chr"]
    mito_length = str(config["mito_length"])
    fasta_file = config["fasta_file"]

    remove_duplicates = config["remove_duplicates"]
    barcode_tag = config["barcode_tag"]
    umi_barcode = config["umi_barcode"]

    handle_overlap = config["handle_overlap"]
    proper_paired = config["proper_paired"]
    base_qual = str(config["base_qual"])
    alignment_quality = config["alignment_quality"]
    nh_max = config["NHmax"]
    nm_max = config["NMmax"]

    barcodes_file = outdir + "/temp/barcode_files/" + sample + ".txt"
    out_pre = outdir + "/temp/sparse_matrices/" + sample

    output_bam = outdir + "/temp/ready_bam/" + sample + ".qc.bam"
    rmlog = output_bam.replace(".qc.bam", ".rmdups.log").replace("/temp/ready_bam/", "/logs/rmdupslogs/")
    filtlog = output_bam.replace(".qc.bam", ".filter.log").replace("/temp/ready_bam/", "/logs/filterlogs/")
    temp_bam0 = output_bam.replace(".qc.bam", ".temp0.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")
    temp_bam1 = output_bam.replace(".qc.bam", ".temp1.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")

    # 1) Filter bam file down to mtDNA reads passing NH/NM/pairing thresholds
    filter_clip_bam(input_bam, temp_bam0, filtlog, mito_genome, proper_paired, nh_max, nm_max)

    # 2) Sort the filtered bam file
    pysam.sort("-o", temp_bam1, temp_bam0)
    pysam.index(temp_bam1)

    # 3) (Optional) Remove duplicates; tenx mode always dedups within the (cell/UMI) barcode tag
    dedup_bam(temp_bam1, output_bam, rmlog, remove_duplicates, umi_barcode)

    # 4) Get per-cell allele counts
    if handle_overlap == "True":
        sumstats_bp_tenx_overlap(output_bam, barcodes_file, out_pre, mito_length, base_qual,
                                  fasta_file, alignment_quality, barcode_tag)
    else:
        sumstats_bp_tenx(output_bam, barcodes_file, out_pre, mito_length, base_qual,
                         fasta_file, alignment_quality, barcode_tag)
