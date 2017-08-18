#' @include io.R
NULL


#' Filter variants based on pre-defined blacklist files
#'
#' \code{filterKnownBlacklist} takes a RangedSummarizedExperiment
#' used in mgatk and returns a subsetted object where mitochondrial
#' variant positions are removed based on pre-defined.
#'
#' The name of the filter includes the genome build and the
#' experimental data source that was used to determine the
#' blacklist.
#'
#' Current filters that are valid include:
#'
#' 'hg19_TF1'
#'
#' @param mitoSE A RangedSummarizedExperiment intialized by mgatk
#' @param filter A
#'
#' @param ... Additional parameters to pass to the
#' importMito.explict function
#'
#' @return Another mgatk object that is an S4 class
#' RangedSummarizedExperiment that removes variants
#' that are included in the specified blacklist.
#'
#' @seealso filterKnownBlacklist
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @importFrom utils read.table
#'
#' @examples
#' # Import an object
#' folder <-paste0(system.file('extdata',package='mgatk'),"/glioma/final")
#' mitoSE <- importMito(folder)
#' dim(mitoSE)
#'
#' # Filter blacklist
#' mitoSE2 <- filterKnownBlacklist(mitoSE, "hg19_TF1")
#' dim(mitoSE)
#'
#' @export
setGeneric(name = "filterKnownBlacklist", def = function(mitoSE, filter)

  standardGeneric("filterKnownBlacklist"))

#' @rdname filterKnownBlacklist
setMethod("filterKnownBlacklist", signature("SummarizedExperiment", "character"),
          definition = function(mitoSE, filter){

  stopifnot(filter %in% c("hg19_TF1"))
  "%ni%" <- Negate("%in%")

  blfile <- paste0(system.file('extdata',package='mgatk'),"/blacklist/", filter,".txt")
  blacklist <- read.table(blfile)
  return(mitoSE[start(rowRanges(mitoSE)) %ni% blacklist[,1],])

  })
