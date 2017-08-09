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
