#' @include mgatk.R
NULL

#' Get the mitochondrial chromosome from a MAE object
#'
#' \code{getMitoChr} takes a MultiAssayExperiment Object
#' and returns a character of the mitochondrial chromosome
#' name (e.g. "chrM" or "MT")
#'
#' @param mitoMAE A MultiAssayExperiment initilized by mgatk
#'
#' @return Character associated with the mitochondrial
#' chromosomes used in mgatk
#'
#' @examples
#'
#' folder <-paste0(system.file('extdata',package='mgatk'),"/glioma/final")
#' mitoMAE <- importMito(folder)
#' getMitoChr(mitoMAE)
#'
#' @export
setGeneric(name = "getMitoChr", def = function(mitoMAE)
  standardGeneric("getMitoChr"))

#' @rdname getMitoChr
setMethod("getMitoChr", signature("MultiAssayExperiment"),
          definition = function(mitoMAE) {
  chr <- levels(seqnames(mitoMAE@ExperimentList[["coverage"]]@rowRanges))[1]
  return(chr)
})

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

