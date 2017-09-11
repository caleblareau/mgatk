#' @include io.R
NULL

#' Make a data.frame of reduced features for plotting
#'
#' \code{filterKnownBlacklist} takes a MultiAssayExperiment
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
#' @param mitoMAE A MultiAssayExperiment intialized by mgatk
#' @param filter A character signifing the position file to
#' filter mitochondrial variants from.
#'
#' @return Another mgatk object that is an S4 class
#' MultiAssayExperiment that removes variants
#' that are included in the specified blacklist.
#'
#' @seealso filterKnownBlacklist
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import MultiAssayExperiment
#' @importFrom utils read.table
#'
#' @examples
#' # Import an object
#' folder <-paste0(system.file('extdata',package='mgatk'),"/glioma/final")
#' mitoMAE <- importMito(folder)
#' dim(mitoMAE@ExperimentList[["coverage"]])
#'
#' # Filter blacklist
#' mitoMAE2 <- filterKnownBlacklist(mitoMAE, "hg19_TF1")
#' dim(mitoMAE2@ExperimentList[["coverage"]])
#'
#' @export
setGeneric(name = "filterKnownBlacklist", def = function(mitoMAE, filter)

  standardGeneric("filterKnownBlacklist"))

#' @rdname filterKnownBlacklist
setMethod("filterKnownBlacklist", signature("MultiAssayExperiment", "character"),
          definition = function(mitoMAE, filter){

  stopifnot(filter %in% c("hg19_TF1"))
  "%ni%" <- Negate("%in%")

  blfile <- paste0(system.file('extdata',package='mgatk'),"/blacklist/", filter,".txt")
  blacklist <- read.table(blfile)
  grblack <- GenomicRanges::GRanges(seqnames = getMitoChr(mitoMAE),
                   IRanges::IRanges(blacklist[,1], width = 1))
  subsetted <- subsetByRow(mitoMAE, grblack)
  return(subsetted)
})
