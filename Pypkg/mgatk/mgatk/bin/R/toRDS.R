#!/usr/bin/env Rscript

# Script that takes raw output from the mgatk python
# package and preapres an R object for analysis

args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]
name <- args[2]

SE <- mgatk::importMito(folder)
saveRDS(SE, file = paste0(folder, "/", name, ".rds"))
