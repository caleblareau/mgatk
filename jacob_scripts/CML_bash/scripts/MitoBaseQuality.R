#!/usr/bin/env Rscript

# input: single sample mpileup cleaned txt output for BQ and BAQ
# output: QC pdf, high quality alt allele positions 
# author: Jacob C Ulirsch 25 June 2017

# Load libraries
library(ggplot2)
library(BuenColors)

# First argument is path to input file without file type
args <- commandArgs(trailingOnly=TRUE)

# Read in BAQ input
BAQ_tbl <- read.table(paste0(args[1], ".BAQ.txt"))

# Calculate mean BAQ for each position
BAQ_tbl$refBAQ <- BAQ_tbl$V6 / (BAQ_tbl$V2 + BAQ_tbl$V3)
BAQ_tbl$altBAQ <- BAQ_tbl$V8 / (BAQ_tbl$V4 + BAQ_tbl$V5)
BAQ_tbl$totalBAQ <- (BAQ_tbl$V6 + BAQ_tbl$V8) / (BAQ_tbl$V2 + BAQ_tbl$V3 + BAQ_tbl$V4 + BAQ_tbl$V5)
BAQ_tbl$refdepth <- BAQ_tbl$V2 + BAQ_tbl$V3
BAQ_tbl$altdepth <- BAQ_tbl$V4 + BAQ_tbl$V5
BAQ_tbl$totaldepth <- BAQ_tbl$V2 + BAQ_tbl$V3 + BAQ_tbl$V4 + BAQ_tbl$V5

# Calculate alt allele ratio (0-1)
BAQ_tbl$ratio <- (BAQ_tbl$V4 + BAQ_tbl$V5) / (BAQ_tbl$V2 + BAQ_tbl$V3 + BAQ_tbl$V4 + BAQ_tbl$V5)

# Read in BQ input
BQ_tbl <- read.table(paste0(args[1], ".BQ.txt"))

# Calculate mean BQ for each position
BAQ_tbl$refBQ <- BQ_tbl$V6 / (BQ_tbl$V2 + BQ_tbl$V3)
BAQ_tbl$altBQ <- BQ_tbl$V8 / (BQ_tbl$V4 + BQ_tbl$V5)
BAQ_tbl$totalBQ <- BQ_tbl$V8 / (BQ_tbl$V2 + BQ_tbl$V3 + BQ_tbl$V4 + BQ_tbl$V5)

# Create pdf of BAQ/BQ/ratio plot
pdf(paste0(args[1], ".BAQplot.pdf"), width=7, height=5)
p <- ggplot(BAQ_tbl, aes(x=altBAQ, y=log10(ratio*100+1), color=altBQ)) + geom_point() + geom_point(data=BAQ_tbl[BAQ_tbl$altBAQ>20 & BAQ_tbl$ratio > 0.005, ]) +theme_bw() + pretty_plot() + labs(x = "Alternate BAQ", y = "AF") +  scale_color_gradientn(colors = jdb_palette("brewer_violet", type = c("continuous"))) + geom_hline(yintercept = log10(1.5)) +  geom_vline(xintercept = 20)
p
dev.off()

# Write txt file of high quality alternative allele positions
write.table(BAQ_tbl[BAQ_tbl$altBAQ>20 & BAQ_tbl$ratio > 0.02, ]$V1, paste0(args[1], ".QCpos.txt"), quote=F, row.names=F, col.names=F)


