import glob
import gzip

import matplotlib
matplotlib.use("Agg")  # this runs headless as part of the pipeline; never needs a GUI backend
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

LETTERS = list("ATCG")


def _load_mgatk_output(output_dir, mito_length):
    # assuming mgatk output naming convention
    base_files = [glob.glob(output_dir + "*.{}.txt.gz".format(nt))[0] for nt in "ATCG"]

    base_coverage_dict = dict()
    for i, nt in enumerate("ATCG"):
        cur_base_data = pd.read_csv(gzip.open(base_files[i]), header=None)

        # gather coverage for forward strand
        fwd_base_df = cur_base_data[[0, 1, 2]].pivot_table(index=1, columns=0)
        fwd_base_df.columns = [x[1] for x in fwd_base_df.columns.values]  # flatten weird multiindex after pivot
        fwd_base_df.index.name = None
        missing_pos = [x for x in range(1, mito_length + 1) if x not in fwd_base_df.columns]
        fwd_base_df[missing_pos] = 0  # fill in missing positions
        fwd_base_df = fwd_base_df.fillna(0).sort_index(axis=1)  # assume all nan are true zeroes

        # gather coverage for forward strand
        rev_base_df = cur_base_data[[0, 1, 3]].pivot_table(index=1, columns=0)
        rev_base_df.columns = [x[1] for x in rev_base_df.columns.values]
        rev_base_df.index.name = None
        missing_pos = [x for x in range(1, mito_length + 1) if x not in rev_base_df.columns]
        rev_base_df[missing_pos] = 0
        rev_base_df = rev_base_df.fillna(0).sort_index(axis=1)

        # organize base data into a dict
        base_coverage_dict[nt] = (fwd_base_df, rev_base_df)

    return base_coverage_dict


def _gather_possible_variants(base_coverage_dict, reference_file, mito_length):
    # sum across cells and strands for each base and position
    aggregated_genotype = pd.DataFrame(np.zeros((4, mito_length)), index=list("ATCG"), columns=np.arange(1, mito_length + 1))
    for nt in base_coverage_dict:
        # sum across cells for each strand separately
        fwd_base_df, rev_base_df = base_coverage_dict[nt]
        fwd_base_sum, rev_base_sum = fwd_base_df.sum(), rev_base_df.sum()

        # sequencing artifact if a base/position is only nonzero for one strand across cells, ignore them
        masking = ~((fwd_base_sum > 0) & (rev_base_sum > 0))  # True if position not >0 for both strands
        fwd_base_sum[masking], rev_base_sum[masking] = 0, 0

        # sum across strands
        aggregated_genotype.loc[nt, :] = fwd_base_sum + rev_base_sum

    # make a reference set of tuples (pos, ref_base)
    ref_set = [x.strip().split() for x in open(reference_file, "r").readlines()]
    ref_n_positions = [int(x[0]) for x in ref_set if x[1].upper() not in LETTERS]
    ref_set = set([(int(x[0]), x[1].upper()) for x in ref_set if x[1].upper() in LETTERS])
    ref_dict = dict(ref_set)

    # make an observed set of tuples which are nonzero
    non_zero_idx = np.where(aggregated_genotype > 0)
    non_zero_bases = [LETTERS[i] for i in non_zero_idx[0]]
    non_zero_pos = [int(i + 1) for i in non_zero_idx[1]]
    observed_set = list(zip(non_zero_pos, non_zero_bases))
    observed_set = set([x for x in observed_set if x[0] not in ref_n_positions])  # disregard positions in ref with N

    # take difference between observed and reference
    variant_set = observed_set - ref_set
    # sort by (position, observed base) rather than position alone: multiple variants can
    # share a position, and sorting on position only left those ties broken by set iteration
    # order, which is hash-randomized per process and therefore not reproducible across runs
    variants = sorted([(x[0], ref_dict[x[0]], x[1]) for x in list(variant_set)], key=lambda x: (x[0], x[2]))  # (pos, ref base, obs base)

    return variants


def call_variants(mgatk_out_dir, sample_prefix, mito_length, low_coverage_threshold, mito_genome):
    """
    Aggregate per-cell allele-count matrices into candidate variants, compute
    heteroplasmy/VMR/strand-correlation summary stats, and write the
    variant_stats/cell_heteroplasmic_df tables plus the VMR-vs-strand-correlation plot.
    """
    mito_length = int(mito_length)
    low_coverage_threshold = int(low_coverage_threshold)

    base_coverage_dict = _load_mgatk_output(mgatk_out_dir, mito_length)
    cell_barcodes = base_coverage_dict["A"][0].index

    # total coverage per position per cell
    total_coverage = pd.DataFrame(np.zeros((len(cell_barcodes), mito_length)), index=cell_barcodes, columns=np.arange(1, mito_length + 1))
    for nt in base_coverage_dict:
        total_coverage += base_coverage_dict[nt][0]
        total_coverage += base_coverage_dict[nt][1]

    # exclude low coverage cells from variant calling
    cell_barcodes = total_coverage.index[total_coverage.mean(axis=1) > low_coverage_threshold]
    for nt in base_coverage_dict:
        base_coverage_dict[nt] = (base_coverage_dict[nt][0].loc[cell_barcodes, :], base_coverage_dict[nt][1].loc[cell_barcodes, :])
    total_coverage = total_coverage.loc[cell_barcodes, :]

    # call potential variants
    variants = _gather_possible_variants(base_coverage_dict, mgatk_out_dir + mito_genome + "_refAllele.txt", mito_length)
    variant_names = ["{}{}>{}".format(x[0], x[1], x[2]) for x in variants]

    # build two <cell by variant tables>, one for each strand
    total_coverage_variant_df = []
    fwd_cell_variant_df, rev_cell_variant_df = [], []
    for i, var in enumerate(variants):
        pos, base = var[0], var[2]
        total_coverage_variant_df.append(total_coverage[pos])
        fwd_cell_variant_df.append(base_coverage_dict[base][0][pos].values)
        rev_cell_variant_df.append(base_coverage_dict[base][1][pos].values)
    total_coverage_variant_df = pd.DataFrame(np.array(total_coverage_variant_df).T, index=cell_barcodes, columns=variant_names)
    fwd_cell_variant_df = pd.DataFrame(np.array(fwd_cell_variant_df).T, index=cell_barcodes, columns=variant_names)
    rev_cell_variant_df = pd.DataFrame(np.array(rev_cell_variant_df).T, index=cell_barcodes, columns=variant_names)
    all_cell_variant_df = fwd_cell_variant_df + rev_cell_variant_df

    # heteroplasmic ratio
    heteroplasmic_df = all_cell_variant_df / total_coverage_variant_df

    # strand correlation
    mask_idx = (fwd_cell_variant_df + rev_cell_variant_df) == 0  # set 0 on both strands to nan to exclude from correlation calculation
    fwd_cell_variant_df[mask_idx] = np.nan
    rev_cell_variant_df[mask_idx] = np.nan
    variant_strand_corr = fwd_cell_variant_df.corrwith(rev_cell_variant_df).round(3)

    # vmr
    variant_mean = all_cell_variant_df.sum() / total_coverage_variant_df.sum()
    variant_var = heteroplasmic_df.var()
    variant_vmr = variant_var / (variant_mean + 0.00000000001)

    # compute other summary stats
    variant_position = [x[0] for x in variants]
    variant_nucleotide = ["{}>{}".format(x[1], x[2]) for x in variants]
    variant_n_cells_conf_detected = ((fwd_cell_variant_df >= 2) & (rev_cell_variant_df >= 2)).sum()
    variant_n_cells_over_5 = (heteroplasmic_df >= 0.05).sum()
    variant_n_cells_over_10 = (heteroplasmic_df >= 0.1).sum()
    variant_n_cells_over_20 = (heteroplasmic_df >= 0.2).sum()
    variant_n_cells_over_95 = (heteroplasmic_df >= 0.95).sum()
    max_heteroplasmy = heteroplasmic_df.max()
    variant_mean_coverage = total_coverage_variant_df.mean()

    # pack summary stats
    variant_output = pd.DataFrame([variant_position, variant_nucleotide, variant_names,
                                   variant_vmr, variant_mean, variant_var,
                                   variant_n_cells_conf_detected, variant_n_cells_over_5,
                                   variant_n_cells_over_10, variant_n_cells_over_20, variant_n_cells_over_95,
                                   max_heteroplasmy, variant_strand_corr, variant_mean_coverage]).T
    variant_output.columns = ["position", "nucleotide", "variant", "vmr", "mean", "variance",
                              "n_cells_conf_detected", "n_cells_over_5",
                              "n_cells_over_10", "n_cells_over_20", "n_cells_over_95",
                              "max_heteroplasmy", "strand_correlation", "mean_coverage"]
    variant_output[["vmr", "mean", "variance", "strand_correlation", "mean_coverage", "max_heteroplasmy"]] = \
        variant_output[["vmr", "mean", "variance", "strand_correlation", "mean_coverage", "max_heteroplasmy"]].astype(float)

    # exclude variants with less than three cells
    multi_cell_variants = variant_output[variant_output["n_cells_conf_detected"] >= 3]["variant"]
    heteroplasmic_df = heteroplasmic_df[multi_cell_variants]

    # generate caleb plot
    plt.figure(figsize=(10, 8))
    plt.scatter(variant_output[variant_output["variant"].isin(multi_cell_variants)]["strand_correlation"],
                np.log10(variant_output[variant_output["variant"].isin(multi_cell_variants)]["vmr"]), s=5)
    plt.axhline(np.log10(0.01), color="red", alpha=0.4, linestyle=":")
    plt.axvline(0.65, color="red", alpha=0.4, linestyle=":")
    plt.xlabel("strand correlation", fontsize=20)
    plt.ylabel("log10(VMR)", fontsize=20)

    # save results
    plt.savefig(mgatk_out_dir + sample_prefix + ".vmr_strand_plot.png")
    plt.close()
    variant_output.to_csv(mgatk_out_dir + sample_prefix + ".variant_stats.tsv.gz", sep="\t", compression="gzip", index=False)
    heteroplasmic_df.to_csv(mgatk_out_dir + sample_prefix + ".cell_heteroplasmic_df.tsv.gz", sep="\t", compression="gzip")
