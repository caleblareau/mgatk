#' ---
#' title: "Perform lineage tracing using michondrial DNA mutations: CML"
#' author: "Jacob C Ulirsch"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' Load libraries
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
library(ff)
library(ffbase)
library(ffbase2)
library(readr)
library(dplyr)
library(gplots)
library(reshape)
library(ape)
library(ggplot2)
library(Rtsne)
library(BuenColors)
library(preprocessCore)
library(ineq)
library(RColorBrewer)
library(cellscape)
library(ggtree)
source("distanceMetrics.R")
"%ni%" <- Negate("%in%")
library(iheatmapr)
library(dtplyr)

#' Read in SNV allele counts as an ffdf
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
options(fftempdir = "../ff/")
mito <- read.table.ffdf(file="/home/julirsch/Desktop/sankaranlab/Mito/16_CML/calls/All.SNPs.MNP.clean.txt", colClasses=c("integer","factor","factor","integer"))
names(mito) <- c("pos", "allele", "cell", "reads")
mito$pos_allele <- as.ff(as.factor(paste0(as.character(mito[,]$pos), "_", as.character(mito[,]$allele))))
dim(mito)

#' Load in reference of mitochondrial genome
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mtDNA.ref <- read_delim("../data/hg19_mtDNA.txt", delim = "\t", col_names = F)
names(mtDNA.ref) <- c("pos", "allele")
mtDNA.ref$pos_allele <- paste0(mtDNA.ref$pos, "_", mtDNA.ref$allele)
ref <- (mtDNA.ref %>% filter())$pos_allele

#' Load in blacklist
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mtDNA.bl <- read_delim("../data/19_All.SNPs.Q0.BLpos.txt", delim = "\t", col_names = F)
names(mtDNA.bl) <- c("pos")

#' Calculate read depth at each base for each sample
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
#Calculate total read depth
mito.cov <- ffdfdply(x=mito, split=as.character.ff(mito$pos), BATCHBYTES=500000000, FUN=function(x) { require(doBy); summaryBy(reads ~ pos + cell, data=x, keep.names=TRUE, FUN=sum)})
names(mito.cov) <- c("pos", "cell", "depth")
#Calculate read depth per cell
mito.cov.cell <- ffdfdply(x=mito.cov, split=as.character.ff(mito.cov$cell), BATCHBYTES=1000000000, FUN=function(x) { require(doBy); summaryBy(log10(depth) ~ cell, data=x, keep.names=TRUE, FUN=function(x) {c(median = median(x), sum <- sum(10^x))})})
mito.cov.cell <- as.data.frame(mito.cov.cell)
names(mito.cov.cell) <- c("cell", "depth.median","depth.sum")
#Calculate read depth per pos
mito.cov.pos <- ffdfdply(x=mito.cov, split=as.character.ff(mito.cov$pos), BATCHBYTES=500000000, FUN=function(x) { require(doBy); result <- summaryBy(log10(depth) ~ pos, data=x, keep.names=TRUE, FUN=function(x) {c(med = median(x), me <- mean(x), s = sd(x), qlow=quantile(x,0.1), qhigh=quantile(x,0.9))})})
names(mito.cov.pos) <- c("pos", "depth.median", "depth.mean", "depth.sd", "depth.qlow", "depth.qhigh")
mito.cov.pos <- as.data.frame(mito.cov.pos)

# Make pass filters for cells and positions
pass.cell <- mito.cov.cell[mito.cov.cell$depth.median>2,]$cell
pass.pos <- mito.cov.pos[mito.cov.pos$depth.median>2,]$pos

#' QC mitochorondrial reads 
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Merge coverage
mito.2 <- merge(mito, mito.cov, by = c("cell", "pos"))
# Filter non-ATCG bases and reference bases
mito.2 <-
  mito.2 %>%
  filter(allele != "<*>") %>%
  filter(!pos_allele %in% ref) %>%
  filter(cell %in% pass.cell) %>%
  as.data.frame()
# Calculate alt allele ratio
mito.2$ratio <- mito.2$reads/mito.2$depth
mito.2$ratio <- round(mito.2$ratio,4)
# QC is pos not pos_allele specific so remove lower ratio or count pos_allele at each pos
mito.filt <- mito.2 %>%
  group_by(pos) %>%
  filter(reads == max(reads)) %>%
  #filter(ratio == max(ratio)) %>%
  .$pos_allele %>%
  unique() %>%
  as.vector()
# Also remove blacklist positions
mito.2 <- mito.2 %>%
  filter(pos_allele %in% mito.filt) %>%
  filter(pos_allele %ni% mtDNA.bl$pos)

#' Subset samples by heteroplasmic variants
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mito.3 <- mito.2 %>%
  #dplyr::filter(cell %in% c("mito.G10_S34","mito.C9_S26","mito.D6_S29")) %>%
  #dplyr::filter(cell %in% cml1266.srr) %>%
  group_by(pos_allele) %>%
  dplyr::summarize(maxratio = max(na.omit(ratio)), minratio = min(na.omit(ratio)), maxdeltaratio = max(dist(na.omit(ratio),method="maximum")), gini=ineq(ratio,type="Gini",na.rm=T), var=var(na.omit(ratio)), minaltcounts=min(na.omit(depth*ratio)), maxaltcounts=max(na.omit(depth*ratio)))
het_var <- mito.3 %>%
  filter(maxratio > 0.1, maxaltcounts > 10) %>%
  filter(pos %in% pass.pos) %>%
  .$pos_allele
mito.plot <- mito.2 %>%
  filter(pos_allele %in% het_var) %>%
  #filter(grepl("BM",cell)) %>%
  #dplyr::filter(cell %in% cml1266.srr) %>%
  dplyr::select(cell,pos_allele,pos,allele,ratio,depth)
# Create ratio matrix
mito.ratio <- cast(mito.plot, cell ~ pos + allele, value = "ratio")
row.names(mito.ratio) <- mito.ratio$cell
mito.ratio <- mito.ratio[,-1]
mito.ratio.mat <- sqrt(as.matrix(mito.ratio))
row.names(mito.ratio.mat) <- row.names(mito.ratio)
colnames(mito.ratio.mat) <- colnames(mito.ratio)
mito.ratio.mat.t <- t(mito.ratio.mat)
# Make coverage matrix
mito.depth <- cast(mito.plot, cell ~ pos + allele, value = "depth")
row.names(mito.depth) <- mito.depth$cell
mito.depth <- mito.depth[,-1]
mito.depth.mat <- as.matrix(mito.depth)
row.names(mito.depth.mat) <- row.names(mito.depth)
colnames(mito.depth.mat) <- colnames(mito.depth)
mito.depth.mat.t <- t(mito.depth.mat)
# Annotate donors
annot <- read.table("../data/bestMetaData_CML.txt",header=T,sep="\t")
mito.donor <- row.names(mito.ratio.mat)
mito.donor <- annot[row.names(mito.ratio.mat),]$Patient

#' Unsupervised reconstruction of lineages
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=9, fig.height=7, fig.align='center'
# Make distance matrix (can choose euc, abs, or sqrt -- sqrt seems to work best)
dist.mat <- dist_mito(df.ratio=mito.ratio.mat.t, df.coverage = log10(mito.depth.mat.t), dist.type = "abs")
# Do neighbor joining
mito.nj <- nj(dist.mat)
info <- split(row.names(dist.mat),mito.donor)
mito.nj <- groupOTU(mito.nj, info)
# Plot NJ results a few different ways
#ggtree(mito.nj,aes(col=group), layout="circular", branch.length="none") + theme(legend.position="bottom") + scale_color_manual(values = jdb_palette("lawhoops")[c(20,1,10,15,3,13,18)]) + geom_tiplab(aes(angle=angle),size=0.8)
#ggtree(mito.nj,aes(col=group), layout="circular") + theme(legend.position="bottom") + scale_color_manual(values = jdb_palette("lawhoops")) + geom_tiplab(aes(angle=angle),size=0.8)
# Plot NJ with heatmap of mutations
mito.order.df <- subset(fortify(mito.nj),isTip)
mito.order <- mito.order.df$label[order(mito.order.df$y, decreasing=TRUE)]
mito.ratio.t <- mito.ratio.mat[mito.order,]
png("/dev/null")
hm <- heatmap.2(data.matrix(mito.ratio.mat.t), trace = "none", cexCol = 0.2, cexRow = 0.25, dendrogram = "col", Rowv=FALSE, col = colorRampPalette(rev(jdb_palette("brewer_violet")[c(9,5,1)]))(length(breaks.val) - 1), breaks = breaks.val)
dev.off()
p <- ggtree(mito.nj,aes(col=group)) + scale_color_manual(values = jdb_palette("lawhoops")) + geom_tiplab(align=TRUE,size=1.5,linesize=0.1) + theme_tree2() + scale_x_continuous(expand=c(0, 0.1))
p#mito.ratio.t[mito.ratio.t > 0.12] <- 0.12
gheatmap(p,mito.ratio.mat[,labels(hm$colDendrogram)], offset = 0.1,width=2,colnames=F) + scale_fill_gradientn(colors = jdb_palette("brewer_violet",type="continuous")) + theme(legend.position="bottom")

p <- ggtree(mito.nj,aes(col=group)) + scale_color_manual(values = jdb_palette("lawhoops")) + geom_tiplab(align=TRUE,size=0,linesize=0.1) + theme_tree2() + scale_x_continuous(expand=c(0, 0.01))
#mito.ratio.t[mito.ratio.t > 0.12] <- 0.12
gheatmap(p,mito.ratio.mat[,labels(hm$colDendrogram)], offset = 0,width=5,colnames=F) + scale_fill_gradientn(colors = jdb_palette("brewer_violet",type="continuous")) + theme(legend.position="bottom")


#' Unsupervised tSNE of lineages
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=9, fig.height=7, fig.align='center'
# Calculate tSNE
set.seed(12345)
mito.ratio.mat[is.na(mito.ratio.mat)] <- 0 # Set missing to 0 for PCA
mito.tsne <- Rtsne(mito.ratio.mat,is_distance=F, check_duplicates = FALSE, pca = TRUE, perplexity=30, theta=0.0, dims=2)
mito.tsne.plot <- as.data.frame(mito.tsne$Y)
names(mito.tsne.plot) <- c("dim1","dim2")
mito.tsne.plot$Class <- mito.donor
# Plot tSNE
ggplot(mito.tsne.plot, aes(x=dim1, y=dim2,color=Class)) +
  geom_point(size=0.1) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE 2D Embedding of Cells by mtDNA Mutations") +
  theme_light(base_size=12) +
  theme(strip.background = element_blank(),
        strip.text.x     = element_blank(),
        axis.text.x      = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.border     = element_blank()) +
  #  scale_color_manual(values = jdb_palette("lawhoops")) +
  theme(legend.position="none", legend.box = "horizontal", legend.text = element_text(size=8)) #+ 
#guides(shape=guide_legend(override.aes=list(color=0.2,size=0.5)))
