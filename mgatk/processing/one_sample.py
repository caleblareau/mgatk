import pysam

from .allele_counts import sumstats_bp, sumstats_bp_overlap
from .dedup import dedup_bam
from .filtering import filter_clip_bam


def process_one_sample(config, inputbam, outputbam, sample):
    """
    Bulk/bcall per-sample pipeline: filter -> sort -> (optional) dedup ->
    allele counting -> depth summary. `config` is the parsed mgatk scatter
    config (the same dict Snakemake loads via `configfile:` for `call`/`bcall`
    runs).
    """
    outdir = config["output_directory"]

    mito_genome = config["mito_chr"]
    mito_length = str(config["mito_length"])
    fasta_file = config["fasta_file"]

    remove_duplicates = config["remove_duplicates"]
    umi_barcode = config["umi_barcode"]
    emit_base_qualities = config["emit_base_qualities"]

    handle_overlap = config["handle_overlap"]
    proper_paired = config["proper_paired"]
    base_qual = str(config["base_qual"])
    alignment_quality = config["alignment_quality"]
    nh_max = config["NHmax"]
    nm_max = config["NMmax"]

    rmlog = outputbam.replace(".qc.bam", ".rmdups.log").replace("/temp/ready_bam/", "/logs/rmdupslogs/")
    filtlog = outputbam.replace(".qc.bam", ".filter.log").replace("/temp/ready_bam/", "/logs/filterlogs/")
    temp_bam0 = outputbam.replace(".qc.bam", ".temp0.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")
    temp_bam1 = outputbam.replace(".qc.bam", ".temp1.bam").replace("/temp/ready_bam/", "/temp/temp_bam/")
    prefix_sm = outdir + "/temp/sparse_matrices/" + sample
    output_depth = outdir + "/qc/depth/" + sample + ".depth.txt"

    # 1) Filter bam file down to mtDNA reads passing NH/NM/pairing thresholds
    filter_clip_bam(inputbam, temp_bam0, filtlog, mito_genome, proper_paired, nh_max, nm_max)

    # 2) Sort the filtered bam file
    pysam.sort("-o", temp_bam1, temp_bam0)
    pysam.index(temp_bam1)

    # 3) (Optional) Remove duplicates; barcode-aware dedup only if a real 2-letter UMI tag is set
    umi_tag = umi_barcode if (umi_barcode != "" and len(umi_barcode) == 2) else ""
    dedup_bam(temp_bam1, outputbam, rmlog, remove_duplicates, umi_tag)

    # 4) Get allele counts per sample / base pair and per-base quality scores
    if handle_overlap == "True":
        sumstats_bp_overlap(outputbam, prefix_sm, mito_genome, mito_length, base_qual, sample,
                             fasta_file, alignment_quality, emit_base_qualities)
    else:
        sumstats_bp(outputbam, prefix_sm, mito_genome, mito_length, base_qual, sample,
                    fasta_file, alignment_quality, emit_base_qualities)

    # 5) Get depth from the coverage sparse matrix
    with open(prefix_sm + ".coverage.txt", "r") as coverage:
        depth = 0
        for row in coverage:
            s = row.split(",")
            depth += int(s[2].strip())
    with open(output_depth, "w") as d:
        d.write(sample + "\t" + str(round(float(depth) / float(mito_length), 2)) + "\n")
