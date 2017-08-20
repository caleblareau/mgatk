#' @include blacklist.R
NULL

#' Get a set of reduced features for plotting
#'
#' \code{mitoFeatureReduction} takes a RangedSummarizedExperiment
#' used in mgatk and returns a data.frame of reduced features
#' that can be used for plotting.
#'
#' The data frame also includes all meta data that is
#' included in the SummarizedExperiment's colData.
#'
#' @param mitoSE A RangedSummarizedExperiment intialized by mgatk
#' @param nPC Default = 10; Number of principal components
#' to be returned and
#' @param nTSNE Default = 2; Number of tSNE coordinates to
#' be returned
#' @param perplexity Default = 30; Perplexity for tSNE computation
#'
#' @return A data.frame containing all of the
#' colData in the original SummarizedExperiment
#' as well as PCs and tSNE coordiantes
#'
#' @importFrom irlba irlba
#' @importFrom Rtsne Rtsne
#' @import SummarizedExperiment
#'
#' @examples
#' # Import an object
#' folder <-paste0(system.file('extdata',package='mgatk'),"/glioma/final")
#' mitoSE <- importMito(folder)
#'
#' # Filter blacklist
#' mitoSE2 <- filterKnownBlacklist(mitoSE, "hg19_TF1")
#'
#' #df <- mitoFeatureReduction(mitoSE2, nPC = 3, nTSNE = 2, perplexity = 30)
#'
#' @export
setGeneric(name = "mitoFeatureReduction", def = function(mitoSE, nPC = 10, nTSNE = 2, perplexity = 30)

  standardGeneric("mitoFeatureReduction"))

#' @rdname mitoFeatureReduction
setMethod("mitoFeatureReduction", signature("SummarizedExperiment", "ANY", "ANY", "ANY"),
          definition = function(mitoSE, nPC = 10, nTSNE = 2, perplexity = 30){

  mat <- assays(mitoSE)[["frequency"]]

  pcs <- irlba(mat, nv = nPC)
  pcdf <- pcs$v
  colnames(pcdf) <- paste0("PC", as.character(1:nPC))

  tSNE <- Rtsne(pcdf, pca = FALSE, perplexity = perplexity)
  tSNEdf <- tSNE$Y[,1:nTSNE]
  colnames(tSNEdf) <- paste0("tSNE", as.character(1:nTSNE))

  df <- data.frame(colData(mitoSE), tSNEdf, pcdf)
  return(df)
})
