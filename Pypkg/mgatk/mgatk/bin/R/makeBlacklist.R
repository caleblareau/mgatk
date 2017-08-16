#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(dtplyr)))
suppressMessages(suppressWarnings(require(dplyr)))

# i/o
args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
mitoLength <- as.numeric(args[3])

# More summary statistics
BAQ_tbl <- fread(paste0("zcat < ", infile))
BAQ_tbl$altBAQ <- BAQ_tbl$V8 / (BAQ_tbl$V4 + BAQ_tbl$V5)

# Compute min and mean BAQ
minMax_tbl <- BAQ_tbl %>% na.omit %>% group_by(V1) %>% summarize(altmean=mean(na.omit(altBAQ)),altmin=quantile(na.omit(altBAQ),0.2))

# Make sure each variant is represented
dummy <- data.table(V1 = 1:mitoLength, altmean = -100, altmin = -100)

# Build output table
out <- rbind(minMax_tbl,dummy)[, .(aMin = max(altmin), aMean = max(altmean)), by = V1]
out <- out[order(out$V1, decreasing = FALSE),]
colnames(out) <- c("Position" , "aMin", "aMean")
out$aMin <- round(out$aMin,1)
out$aMean <- round(out$aMean,1)
write.table(data.frame(out), file = outfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

