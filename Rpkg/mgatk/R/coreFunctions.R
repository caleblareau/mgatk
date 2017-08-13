#' @include mgatk.R
NULL

#' Compute weighted distance function.
#'
#' \code{computeWeightedDistance} takes a matrix of features x samples.
#' Essentially computing the pairwise distance difference between
#'
#' @param x Data matrix of ratios (for mito: minor allele frequencies)
#' @param y Data matrix of weights (for mito: coverage)
#' @param method = "abs" The distance metric for the data.
#'
#' @return One or more outputs. See above
#' @import methods
#' @import Rcpp
#' @examples
#'
#' x <- matrix(runif(1000), nrow = 20) # 50 samples
#' y <- matrix(rpois(1000, 10), nrow = 20) # 50 samples
#' d <- computeWeightedDistance(x,y)
#'
#' @export
setGeneric(name = "computeWeightedDistance", def = function(x, y = NULL, method = "euc")
  standardGeneric("computeWeightedDistance"))

#' @rdname computeWeightedDistance
setMethod("computeWeightedDistance", signature("ANY", "ANY", "ANY"),
          definition = function(x, y = NULL, method = "euc") {

            stopifnot(method %in% c("abs", "euc", "sqrt"))

            if (is.null(y)) y <- matrix(1, nrow = dim(x)[1], ncol = dim(x)[2])

            stopifnot(dim(x)[1] == dim(y)[1])
            stopifnot(dim(x)[2] == dim(y)[2])

            stopifnot(is.numeric(x))
            stopifnot(is.numeric(y))

            # Remove NAs from the data
            x[is.na(x)] <- 0
            y[is.na(y)] <- 1
            y[y < 0] <- 0.1 # charity count

            if(method == "abs"){
              mat <- calcWdist_abs(x, y)
            } else if (method == "euc"){
              mat <- calcWdist_euclidean(x, y)
            } else if (method == "sqrt"){
              mat <- calcWdist_sqrt(x, y)
            }
            return(mat)
          })

#' Call mitochondrial heteroplasmic variants
#'
#' \code{callVariants} takes a RangedSummarizedExperiment
#' and infers which allele, the frequence per sample,
#' and other meta data is approproriate for down-stream
#' analysis in mgatk.
#'
#' @param mitoSE A RangedSummarizeExperiment initialized by
#' mgatk that will have variants called.
#' @param filterRefAllele Default = TRUE If the reference allele
#' is included in the mcols slot of the rowRanges of the mitoSE
#' object, then filter these when doing variant calling.
#'
#' @return A RangedSummarizedExperiment object with variants
#' called. This adds a 'variantFreq' slot to assays with the
#' per-sample variant allele frequency as well as meta data
#' to the rowRanges that summarizes the variants
#'
#' @examples
#'
#' datafile <- paste0(system.file('extdata',package='mgatk'), "/mgatk_CML.txt.gz")
#' mitoSE <- importMito.txt(datafile)
#' vSE <- callVariants(mitoSE, filterRefAlle = FALSE)
#'
#' @export
setGeneric(name = "callVariants", def = function(mitoSE, filterRefAllele = TRUE)
  standardGeneric("callVariants"))

#' @rdname callVariants
setMethod("callVariants", signature("SummarizedExperiment", "ANY"),
          definition = function(mitoSE, filterRefAllele = TRUE){
  DNAplus <- c('A', 'C', 'G', 'T', 'coverage')
  stopifnot(all(DNAplus %in% names(assays(mitoSE))))

  d <- matACTG(mitoSE)
  dmax <- apply(d/rowSums(d),1,max)


  prevalenceOrder <- orderACTG(mitoSE)
  df <- data.frame(mcols(rowRanges(mitoSE)))

  # Make vector of reference alleles; non-sense if
  # user doesn't want to filter or it's not available
  if(filterRefAllele & "refAllele" %in% colnames(df)){
    refAllele <- df$refAllele
  }else{
    refAllele <- rep("N", dim(mitoSE)[1])
  }

  # Loop over each variant position
  loopout <- lapply(1:dim(mitoSE)[1], function(i){
      if(is.na(dmax[i]) | dmax[i] == 1){
        keep <- FALSE
        base <- "N"
        freq <- 0
        sampleFreq <- rep(0, dim(mitoSE)[2])
      } else {
        keep <- TRUE
        bases <- strsplit(prevalenceOrder, split = ",")[[1]]
        base <- ifelse(bases[1] == refAllele[i], bases[2], bases[1])
        freq <- dmax[i]
        sampleFreq <- assays(mitoSE)[[base]][i,]/assays(mitoSE)[["coverage"]][i,]
      }
      list(data.frame(keep, base, freq, stringsAsFactors = FALSE), sampleFreq)
  })

  # Build qc data frame/table
  getpos <- function(listList,pos) listList[[pos]]
  qcDF <- data.table::rbindlist(lapply(loopout, getpos, 1))

  # With variants called, build new metadata
  SE <- mitoSE[qcDF$keep,]
  mcols(rowRanges(SE)) <- cbind(mcols(rowRanges(mitoSE))[qcDF$keep,],
    data.frame(
      base = qcDF[qcDF$keep,][["base"]],
      totalFreq = qcDF[qcDF$keep,][["freq"]],
      stringsAsFactors = FALSE
    ))
  assays(SE)$variantFreq <- t(sapply(loopout, getpos, 2))[qcDF$keep,]
  return(SE)

})

####################
# Internal Functions
####################

# Function returns a data.matrix
# of the total allele counts per nucleotide / position
matACTG <- function(mitoSE){
  d <- data.matrix(data.frame(
    A = rowSums((assays(mitoSE)[["A"]])),
    C = rowSums((assays(mitoSE)[["C"]])),
    G = rowSums((assays(mitoSE)[["G"]])),
    T = rowSums((assays(mitoSE)[["T"]]))
  ))
  return(d)
}


# Function returns a vector of a comma separated list of DNA
# nucleotides by decreasing prevalence in a Summarized Experiment
orderACTG <- function(mitoSE){
  DNA <- c('A', 'C', 'G', 'T')
  stopifnot(all(DNA %in% names(assays(mitoSE))))
  d <- matACTG(mitoSE)

  orderACTGvec <- function(vec){
    names(vec) <- DNA
    paste(names(vec)[order(vec, decreasing = TRUE)], collapse = ",")
  }
  acgtVec <- apply(d, 1, orderACTGvec)
  return(acgtVec)
}
