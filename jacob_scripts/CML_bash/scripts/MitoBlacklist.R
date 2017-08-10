#!/usr/bin/env Rscript

# input: multiple sample mpileup cleaned txt output for BAQ
# output: QC pdf, multiple sample blacklist positions
# author: Jacob C Ulirsch 25 June 2017

# Load libraries
library(ggplot2)
library(BuenColors)
library(dplyr)

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
BAQ_tbl$ratio <- BAQ_tbl$altdepth / BAQ_tbl$totaldepth

# Make blacklist
BAQ_tbl.plot <- BAQ_tbl %>% group_by(V1) %>% summarize(altmean=mean(na.omit(altBAQ)),altmin=quantile(na.omit(altBAQ),0.2),altmax=quantile(na.omit(altBAQ),0.8))
BAQ_tbl.plot.bl <- BAQ_tbl.plot %>% filter(altmin < mean(altmin) - 1.645 * sd(altmin) | altmean < mean(altmean) - 1.645 * sd(altmean))
bl <- BAQ_tbl.plot.bl$V1
bl.df <- data.frame(blacklist=bl)

# Write txt file of high quality alternative allele positions
write.table(bl.df,paste0(args[1],".BLpos.txt"),quote=F,row.names=F,col.names = F)

# Create pdf of BAQ plot with blacklist
pdf(paste0(args[1],".BLplot.pdf"),width=7,height=5)
p <- ggplot(BAQ_tbl.plot, aes(x = V1)) +
  geom_ribbon(aes(ymin=altmin,ymax=altmax), colour = jdb_palette("brewer_green")[3]) +
  geom_line(aes(y=altmean),size=0.2, colour = jdb_palette("brewer_green")[9]) +
  geom_point(data=bl.df, aes(x=blacklist, y=5),size=0.1, colour = jdb_palette("brewer_red")[7]) +
  #geom_point(aes(x=V1, y=state),size=0.1, colour = jdb_palette("brewer_red")[7]) +
  pretty_plot() +
  labs(x = "Mitochondrial Genome Position", y = "BAQ")
p
dev.off()

# In progress HMM 
#set.seed(1)
#hmm <- depmix(altmin ~ 1, nstates = 3, data=BAQ_tbl.plot)
#hmm.pars <- c(0.01,0.985,0.005,0.95,0.05,0.00,0.03,0.095,0.02,0.00,0.05,0.95,-4.7,2,0,1,6.8,5.2)
#fixed <- c(0,0,0,0,0,1,0,0,0,1,0,0,1,1,1,1,1,1)
#names(hmm.pars) <- names(c(unlist(getpars(hmm))))
#hmm <- setpars(hmm,hmm.pars)
#hmmfit <- fit(hmm, verbose = TRUE,fixed=fixed)
#summary(hmmfit)
#post_probs <- posterior(hmmfit)
#BAQ_tbl.plot$state <- post_probs$state


