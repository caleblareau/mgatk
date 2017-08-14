##!/usr/bin/env Rscript

# Reads depth table and computes quantiles
# before determining which samples to keep

args <- commandArgs(trailingOnly = TRUE)
depthfile <- args[1]
outfile <- args[2]
percentile_keep <- args[3]

df <- read.table(depthfile, header = FALSE)
perc.rank <- function(x) trunc(rank(x))/length(x)
df$perc <- perc.rank(df[,2])
write.table(df[df$perc*100 > as.numeric(percentile_keep),1], file = outfile,
	row.names = FALSE, col.names = FALSE, quote = FALSE)
