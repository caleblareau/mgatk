#' ---
#' title: "Perform lineage tracing using michondrial DNA mutations: analysis of TF1 and CD34+ clonal colonies"
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

#' Read in SNV allele counts as an ffdf
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
options(fftempdir = "../ff/")
mito <- read.table.ffdf(file="/home/julirsch/Desktop/sankaranlab/Mito/19_colonies4/calls/All.SNPs.MNP.clean.txt")
names(mito) <- c("pos", "allele", "cell", "reads")
mito$pos_allele <- as.ff(as.factor(paste0(as.character(mito[,]$pos), "_", as.character(mito[,]$allele))))
mito.names <- read.table("../data/19_colonynames5.txt",header=T,sep="\t")
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
#ffsave(mito.cov,file="../processed/19_mito.cov.rds")
#ffload("../processed/19_mito.cov_20170325.rds")
#Calculate read depth per cell
mito.cov.cell <- ffdfdply(x=mito.cov, split=as.character.ff(mito.cov$cell), BATCHBYTES=1000000000, FUN=function(x) { require(doBy); summaryBy(log10(depth) ~ cell, data=x, keep.names=TRUE, FUN=function(x) {c(median = median(x), sum <- sum(10^x))})})
mito.cov.cell <- as.data.frame(mito.cov.cell)
names(mito.cov.cell) <- c("cell", "depth.median","depth.sum")
#Calculate read depth per pos
mito.cov.pos <- ffdfdply(x=mito.cov, split=as.character.ff(mito.cov$pos), BATCHBYTES=500000000, FUN=function(x) { require(doBy); result <- summaryBy(log10(depth) ~ pos, data=x, keep.names=TRUE, FUN=function(x) {c(med = median(x), me <- mean(x), s = sd(x), qlow=quantile(x,0.1), qhigh=quantile(x,0.9))})})
names(mito.cov.pos) <- c("pos", "depth.median", "depth.mean", "depth.sd", "depth.qlow", "depth.qhigh")
mito.cov.pos <- as.data.frame(mito.cov.pos)

#' QC mitochorondrial reads 
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Merge coverage
mito.2 <- merge(mito, mito.cov, by = c("cell", "pos"))
# Filter non-ATCG bases and reference bases
mito.2 <-
  mito.2 %>%
  filter(allele != "<*>") %>%
  filter(!pos_allele %in% ref) %>%
  as.data.frame()
# Calculate alt allele ratio
mito.2$ratio <- mito.2$reads/mito.2$depth
mito.2$ratio <- round(mito.2$ratio,4)
# QC is pos not pos_allele specific so remove lower ratio pos_allele at each pos
mito.filt <- mito.2 %>%
  group_by(pos) %>%
  filter(ratio == max(ratio)) %>%
  .$pos_allele %>%
  unique() %>%
  as.vector()
# Also remove blacklist positions
mito.2 <- mito.2 %>%
  filter(pos_allele %in% mito.filt) %>%
  filter(pos_allele %ni% mtDNA.bl$pos)

#' Subset based upon feature (pos_allele) characteristics 
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Calculate feature statistics
mito.3 <- mito.2 %>%
  as.data.frame() %>%
  group_by(pos_allele) %>%
  summarize(maxratio = max(ratio), minratio = min(ratio), maxdeltaratio = max(dist(ratio,method="maximum")), gini=ineq(ratio,type="Gini",na.rm=T), var=var(ratio))
# Filter and include only variants above detection threshold with some minimum differences between samples
het_var <- mito.3 %>%
  filter(maxratio > 0.005, maxdeltaratio > 0.005, minratio < 0.9) %>%
  .$pos_allele
mito.plot <- mito.2 %>%
  filter(pos_allele %in% het_var) %>%
  dplyr::select(cell,pos_allele,pos,allele,ratio,depth)

#' Hierarchical clustering of samples/cells and variants
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=9, fig.height=7, fig.align='center'
# Transform from long to wide
mito.plot.w <- cast(mito.plot, cell ~ pos + allele, value = "ratio")
row.names(mito.plot.w) <- mito.plot.w$cell
mito.plot.w <- mito.plot.w[,-1]
# Sqrt of ratios makes patterns more obvious since many interesting variants are lower frequency
mito.matrix <- sqrt(as.matrix(mito.plot.w))
# Setup and plot heatmap
row.names(mito.matrix) <- row.names(mito.plot.w)
colnames(mito.matrix) <- colnames(mito.plot.w)
breaks.val = seq(0, 1, by = 0.01)
hm <- heatmap.2(mito.matrix, trace = "none", cexCol = 0.2, cexRow = 0.25, dendrogram = "both", col = colorRampPalette(rev(jdb_palette("brewer_violet")[c(9,5,1)]))(length(breaks.val) - 1), breaks = breaks.val)

#' Interactive investigation of variants at known lineages
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Create dataframe for cellscape
mito.mutmat <- merge(mito.plot,mito.names,by.x="cell",by.y="samples")
names(mito.mutmat) <- c("cell","pos_allele","coord","variant_base","VAF","depth","single_cell_id")
mito.mutmat$chr <- "MT"
# Need to edit the code so that we can pick cell colors, scale colors, and have a better scale
mito.mutmat$VAF <- sqrt(mito.mutmat$VAF) / sqrt(0.12); mito.mutmat[mito.mutmat$VAF > 1,]$VAF <- 1
mito.edges <- read.table("../data/19_mito_tree_edges.txt",header=T)
mito.annot <- read.table("../data/19_mito_annot.txt",header=T)
# Need to edit code to work on order
mito.mut_order <- paste0("MT:",gsub("_.*","",labels(hm$colDendrogram)))
cellscape(mut_data=mito.mutmat, tree_edges=mito.edges, sc_annot=mito.annot, mut_order=mito.mut_order)

#' Unsupervised reconstruction of lineages
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=9, fig.height=7, fig.align='center'
# Create dataframe for ggtree
mito.mutmat2 <- merge(mito.plot,mito.names,by.x="cell",by.y="samples")
mito.coverage <- cast(mito.mutmat2, pos_allele ~ names, value = "depth")[,-1]
mito.ratio  <- cast(mito.mutmat2, pos_allele ~ names, value = "ratio")
row.names(mito.ratio) <- mito.ratio[,1]
mito.ratio <- mito.ratio[,-1]
mito.ratio.t  <- cast(mito.mutmat2, names ~ pos_allele, value = "ratio")
row.names(mito.ratio.t) <- mito.ratio.t[,1]
mito.ratio.t <- mito.ratio.t[,-1]
# Make distance matrix (can choose euc, abs, or sqrt -- sqrt seems to work best)
dist.mat <- dist_mito(df.ratio=sqrt(mito.ratio), df.coverage=mito.coverage, dist.type = "abs")
# Do neighbor joining
mito.nj <- nj(dist.mat)
info <- split(as.character(mito.annot$single_cell_id),mito.annot$genotype)
mito.nj <- groupOTU(mito.nj, info)
# Plot NJ results a few different ways
ggtree(mito.nj,aes(col=group), layout="circular", branch.length="none") + theme(legend.position="bottom") + scale_color_manual(values = jdb_palette("lawhoops")) + geom_tiplab(aes(angle=angle),size=0.8)
ggtree(mito.nj,aes(col=group), layout="circular") + theme(legend.position="bottom") + scale_color_manual(values = jdb_palette("lawhoops")) + geom_tiplab(aes(angle=angle),size=0.8)
# Plot NJ with heatmap of mutations
mito.order.df <- subset(fortify(mito.nj),isTip)
mito.order <- mito.order.df$label[order(mito.order.df$y, decreasing=TRUE)]
mito.ratio.t <- mito.ratio.t[mito.order,]
png("/dev/null")
hm <- heatmap.2(data.matrix(sqrt(mito.ratio.t)), trace = "none", cexCol = 0.2, cexRow = 0.25, dendrogram = "col", Rowv=FALSE, col = colorRampPalette(rev(jdb_palette("brewer_violet")[c(9,5,1)]))(length(breaks.val) - 1), breaks = breaks.val)
dev.off()
p <- ggtree(mito.nj,aes(col=group)) + scale_color_manual(values = jdb_palette("lawhoops")) + geom_tiplab(align=TRUE,size=1.5,linesize=0.1) + theme_tree2() + scale_x_continuous(expand=c(0, 30))
mito.ratio.t[mito.ratio.t > 0.12] <- 0.12
gheatmap(p,sqrt(mito.ratio.t[,labels(hm$colDendrogram)]),offset = 50,width=2,colnames=F) + scale_fill_gradientn(colors = jdb_palette("brewer_violet",type="continuous")) + theme(legend.position="bottom")

#' Unsupervised tSNE of lineages
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=9, fig.height=7, fig.align='center'
# Calculate tSNE
mito.annot2 <- merge(mito.names,mito.annot,by.x="names",by.y="single_cell_id")
mito.annot2 <- mito.annot2[order(as.vector(mito.annot2$names)),]
mito.ratio.t <- mito.ratio.t[sort(as.vector(mito.annot2$names)),]
set.seed(1)
mito.tsne <- Rtsne(mito.ratio.t,is_distance=F, check_duplicates = FALSE, pca = TRUE, perplexity=10, theta=0.0, dims=2)
mito.tsne.plot <- as.data.frame(mito.tsne$Y)
names(mito.tsne.plot) <- c("dim1","dim2")
mito.tsne.plot$Class <- mito.annot2$genotype
# Plot tSNE
ggplot(mito.tsne.plot, aes(x=dim1, y=dim2, color=Class)) +
  geom_point(size=1.25) +
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
  scale_color_manual(values = jdb_palette("lawhoops"))


# Estimate clonality of bulk population
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=5, fig.height=5, fig.align='center'
# Calculate tSNE
library(glmnet)
# Set up Xs and Ys
bulk1 <- "20170323-TF1-bulk1"
bulk2 <- "20170323-TF1-bulk2"
mix1 <- "20170514-TF1-mix1"
mix2 <- "20170514-TF1-mix2"
gen1 <- c("20170323-TF1-A9", "20170514-TF1-B11-E8", "20170323-TF1-C7", "20170323-TF1-D3", "20170514-TF1-F4-E9", "20170323-TF1-G11", "20170323-TF1-B3", "20170514-TF1-B5-C3", "20170323-TF1-B9", "20170323-TF1-C4", "20170323-TF1-C10", "20170323-TF1-D2")
mito.ratio.bulk1 <- mito.ratio %>% dplyr::select(bulk1)
mito.ratio.bulk1.vec <- sqrt(as.vector(mito.ratio.bulk1[,1]))
mito.ratio.bulk2 <- mito.ratio %>% dplyr::select(bulk2)
mito.ratio.bulk2.vec <- sqrt(as.vector(mito.ratio.bulk2[,1]))
mito.ratio.mix1 <- mito.ratio %>% dplyr::select(mix1)
mito.ratio.mix1.vec <- sqrt(as.vector(mito.ratio.mix1[,1]))
mito.ratio.mix2 <- mito.ratio %>% dplyr::select(mix2)
mito.ratio.mix2.vec <- sqrt(as.vector(mito.ratio.mix2[,1]))
mito.ratio.gen1 <- mito.ratio %>% dplyr::select(gen1)
mito.ratio.gen1.mat <- sqrt(as.matrix(as.data.frame(mito.ratio.gen1)))
colnames(mito.ratio.gen1.mat) <- colnames(mito.ratio.gen1)

# Make truth
mito.ratio.mix1.truth <- sqrt(1/3*mito.ratio[,"20170514-TF1-B11-E8"] + 1/3*mito.ratio[,"20170514-TF1-B5-C3"] + 1/3*mito.ratio[,"20170514-TF1-F4-E9"])
mito.ratio.mix2.truth <- sqrt(4/15*mito.ratio[,"20170514-TF1-B11-E8"] + 1/15*mito.ratio[,"20170514-TF1-B5-C3"] + 10/15*mito.ratio[,"20170514-TF1-F4-E9"])

# Set up constraints on the coefficients
zero.vec <- rep(0,dim(mito.ratio.gen1)[2])
one.vec <- rep(1,dim(mito.ratio.gen1)[2])

# Choice of lambda
lambda.val <- "lambda.1se"
#lambda.val <- "lambda.min"

# Choice of CV fold
fold <- dim(mito.ratio.gen1)[1] # LOO

#Fit lasso model and cross validate (mix1)
fit.cv <- cv.glmnet(mito.ratio.gen1.mat,mito.ratio.mix1.vec,nfolds=fold,intercept=F,alpha=1,lower.limits=zero.vec,upper.limits=one.vec,standardize=TRUE)
plot(fit.cv)
fit.pred <- predict(fit.cv, as.matrix(mito.ratio.gen1.mat, s = lambda.val))
cor(cbind(fit.pred, mito.ratio.mix1.vec))[1, 2]^2
plot(fit.pred, mito.ratio.mix1.vec, cex=0.2)
plot(mito.ratio.mix1.truth, fit.pred, cex=0.2)
plot(mito.ratio.mix1.truth,mito.ratio.mix1.vec, cex=0.2,ylim=c(0,0.4),xlim=c(0,0.4),ylab="AF observed in mixture sample",xlab="AF from known combination of clones")
abline(lm(mito.ratio.mix1.vec ~ fit.pred), col = "red", lwd = 1)
abline(a=0,b=1, col = "red", lwd = 1, lty=3)
fit.coef <- coef(fit.cv, s=lambda.val); print(sum(fit.coef))
if(sum(fit.coef) > 1) {
  fit.coef <- coef(fit.cv, s=lambda.val) / sum(coef(fit.cv, s=lambda.val))
}
fit.coef.df <- as.data.frame(as.matrix(fit.coef))
names(fit.coef.df) <- "value"
fit.coef.df$cell <- gsub(".*-TF1-","",row.names(fit.coef.df))
fit.coef.df$cell <- gsub("-.*","",fit.coef.df$cell)
fit.coef.df$type <- "estimate"
fit.truth <- data.frame(value=c(1/3,1/3,1/3), cell=c("B11","B5","F4"), type=c("truth","truth","truth"))
fit.coef.df <- rbind(fit.coef.df,fit.truth)
ggplot(fit.coef.df, aes(x=type, y=value, fill=cell, color=cell)) + geom_col() + theme_bw() + pretty_plot() + scale_color_manual(values = jdb_palette("lawhoops")) + scale_fill_manual(values = jdb_palette("lawhoops")) + ylim(0,1)

#Fit lasso model and cross validate (mix2)
fit.cv <- cv.glmnet(mito.ratio.gen1.mat,mito.ratio.mix2.vec,nfolds=fold,intercept=F,alpha=1,lower.limits=zero.vec,upper.limits=one.vec,standardize=TRUE)
plot(fit.cv)
fit.pred <- predict(fit.cv, as.matrix(mito.ratio.gen1.mat, s = lambda.val))
cor(cbind(fit.pred, mito.ratio.mix2.vec))[1, 2]^2
#plot(fit.pred, mito.ratio.mix2.vec, cex = 0.2)
plot(mito.ratio.mix2.truth,mito.ratio.mix2.vec, cex=0.2,ylim=c(0,0.4),xlim=c(0,0.4),ylab="AF observed in mixture sample",xlab="AF from known combination of clones")
#plot(mito.ratio.mix2.truth,fit.pred, cex=0.2)
#abline(lm(mito.ratio.mix2.vec ~ fit.pred), col = "red", lwd = 1)
abline(a=0,b=1, col = "red", lwd = 1, lty=3)
fit.coef <- coef(fit.cv, s=lambda.val); print(sum(fit.coef))
if(sum(fit.coef) > 1) {
  fit.coef <- coef(fit.cv, s=lambda.val) / sum(coef(fit.cv, s=lambda.val))
}
fit.coef.df <- as.data.frame(as.matrix(fit.coef))
names(fit.coef.df) <- "value"
fit.coef.df$cell <- gsub(".*-TF1-","",row.names(fit.coef.df))
fit.coef.df$cell <- gsub("-.*","",fit.coef.df$cell)
fit.coef.df$type <- "estimate"
fit.truth <- data.frame(value=c(4/15,1/15,10/15), cell=c("B11","B5","F4"), type=c("truth","truth","truth"))
fit.coef.df <- rbind(fit.coef.df,fit.truth)
ggplot(fit.coef.df, aes(x=type, y=value, fill=cell, color=cell)) + geom_col() + theme_bw() + pretty_plot() + scale_color_manual(values = jdb_palette("lawhoops")) + scale_fill_manual(values = jdb_palette("lawhoops")) + ylim(0,1)

#Fit lasso model and cross validate (bulk1)
fit.cv <- cv.glmnet(mito.ratio.gen1.mat,mito.ratio.bulk1.vec,nfolds=fold,intercept=F,alpha=1,lower.limits=zero.vec,upper.limits=one.vec,standardize=TRUE)
plot(fit.cv)
fit.pred <- predict(fit.cv, as.matrix(mito.ratio.gen1.mat, s = lambda.val))
cor(cbind(fit.pred, mito.ratio.bulk1.vec))[1, 2]^2
plot(fit.pred, mito.ratio.bulk1.vec, cex = 0.2)
abline(lm(mito.ratio.bulk1.vec ~ fit.pred), col = "red", lwd = 1)
abline(a=0,b=1, col = "red", lwd = 1, lty=3)
fit.coef <- coef(fit.cv, s=lambda.val); print(sum(fit.coef))
if(sum(fit.coef) > 1) {
  fit.coef <- coef(fit.cv, s=lambda.val) / sum(coef(fit.cv, s=lambda.val))
}
fit.coef.df <- as.data.frame(as.matrix(fit.coef))
names(fit.coef.df) <- "value"
fit.coef.df$cell <- row.names(fit.coef.df)
ggplot(fit.coef.df, aes(x="estimate", y=value, color=cell,fill=cell)) + geom_col() + theme_bw() + pretty_plot() + scale_color_manual(values = jdb_palette("lawhoops")) + scale_fill_manual(values = jdb_palette("lawhoops")) + ylim(0,1)

#Fit lasso model and cross validate (bulk2)
fit.cv <- cv.glmnet(mito.ratio.gen1.mat,mito.ratio.bulk2.vec,nfolds=fold,intercept=F,alpha=1,lower.limits=zero.vec,upper.limits=one.vec,standardize=TRUE)
plot(fit.cv)
fit.pred <- predict(fit.cv, as.matrix(mito.ratio.gen1.mat, s = lambda.val))
cor(cbind(fit.pred, mito.ratio.bulk2.vec))[1, 2]^2
plot(fit.pred, mito.ratio.bulk2.vec, cex = 0.2, ylab="AF observed in bulk TF1 cell line", xlab="Estimated AF")
abline(lm(mito.ratio.bulk2.vec ~ fit.pred), col = "red", lwd = 1)
abline(a=0,b=1, col = "red", lwd = 1, lty=3)
fit.coef <- coef(fit.cv, s=lambda.val); print(sum(fit.coef)) 
if(sum(fit.coef) > 1) {
  fit.coef <- coef(fit.cv, s=lambda.val) / sum(coef(fit.cv, s=lambda.val))
}
fit.coef.df <- as.data.frame(as.matrix(fit.coef))
names(fit.coef.df) <- "value"
fit.coef.df$cell <- gsub(".*-TF1-","",row.names(fit.coef.df))
fit.coef.df$cell <- gsub("-.*","",fit.coef.df$cell)
ggplot(fit.coef.df, aes(x="estimate", y=value, color=cell,fill=cell)) + geom_col() + theme_bw() + pretty_plot() + scale_color_manual(values = jdb_palette("lawhoops")) + scale_fill_manual(values = jdb_palette("lawhoops")) + ylim(0,1)

#' Unsupervised reconstruction of lineages + TF zscores
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE, fig.width=9, fig.height=7, fig.align='center'
# Read in z-score data.frame
all <- readRDS("../processed/TF1_calledPeaks.27June17.deviationScores.rds")
TF1mat <- all@assays$data$z
TF1mat.t <- data.frame(t(TF1mat))
row.names(TF1mat.t) <- gsub("\\.","-",row.names(TF1mat.t))
TF1mat.t2 <- merge(mito.names,TF1mat.t,by.x="samples",by.y="row.names")
row.names(TF1mat.t2) <- TF1mat.t2$names
TF1mat.t2 <- TF1mat.t2[,-1:-2]
# Plot NJ with heatmap of z scores
mito.order.df <- subset(fortify(mito.nj),isTip)
mito.order <- mito.order.df$label[order(mito.order.df$y, decreasing=TRUE)]
TF1mat.t2 <- TF1mat.t2[mito.order,]
TF1mat.t2.sub <- TF1mat.t2 %>% dplyr::select(vdfsort[1:20,]$TF)
png("/dev/null")
hm <- heatmap.2(data.matrix(TF1mat.t2.sub), trace = "none", cexCol = 0.2, cexRow = 0.25, dendrogram = "col", Rowv=FALSE, col = colorRampPalette(rev(jdb_palette("brewer_violet")[c(9,5,1)]))(length(breaks.val) - 1), breaks = breaks.val)
TF1mat.t2.sub <- apply(TF1mat.t2.sub,2,minmaxscale)
dev.off()
p <- ggtree(mito.nj,aes(col=group)) + scale_color_manual(values = jdb_palette("lawhoops")) + geom_tiplab(align=TRUE,size=1.5,linesize=0.1) + theme_tree2() + scale_x_continuous(expand=c(0, 30))
gheatmap(p,TF1mat.t2.sub[,labels(hm$colDendrogram)],colnames=T,font.size=0.2,offset = 50,width=2) + scale_fill_gradientn(colors = jdb_palette("solar_extra",type="continuous")) + theme(legend.position="bottom")

TF1counts <- read.table("../data/25June17_allSamples_Combined_250bp.counts.txt",header=T)
TF1counts.t <- t(TF1counts)
row.names(TF1counts.t) <- gsub("\\.","-",row.names(TF1counts.t))
TF1counts.t2 <- merge(mito.names,TF1counts.t,by.x="samples",by.y="row.names")
row.names(TF1counts.t2) <- TF1counts.t2$names
TF1counts.t2 <- TF1counts.t2[,-1:-2]
TF1counts.2 <- t(TF1counts.t2)
TFcounts.cor <- cor(TF1counts.2)
TFcounts.dist <- as.dist(1-cor(TFcounts.cor))
TFcounts.dist <- dist(TF1counts.2)
mito.nj <- nj(TFcounts.dist)
