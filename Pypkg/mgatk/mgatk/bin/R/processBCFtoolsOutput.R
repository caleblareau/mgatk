#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(dtplyr)))
suppressMessages(suppressWarnings(require(dplyr)))

# i/o
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
altfile <- args[2]
covfile <- args[3]

# Import bcftools output and make coverage and alt alleles file
dt <- fread(infile, sep = "\t")
altdt <- dt[dt[[5]]!=0,c(1,2,3,5)]
covdt <- dt[,c(1,3,4)]
remove(dt)
setkey(covdt,NULL)
covdt <- unique(covdt)
write.table(altdt, file = altfile, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(covdt, file = covfile, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
