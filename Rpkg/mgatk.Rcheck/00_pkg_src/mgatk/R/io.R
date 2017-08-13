#' @include coreFunctions.R
NULL

#' Import data from plain text file into
#'
#' \code{importMito.txt} takes a long file of position, allele,
#' sample, and coverage (non-zero) in some order.
#'
#' @param datafile The filepath of a plain text or gzipped file
#' that contains data in a long matrix format of the position,
#' allele, sample, and coverage. This will be imported and a
#' data object will be rendered that enables downstream
#' mgatk analyses
#' @param pos Default = 1 column index of the genomic position
#' @param allele Default = 2 column index of the allele (A/C/G/T)
#' @param sample Default = 3 column index of the sample name
#' @param coverage Default = 4 column index that shows the number
#' of reads for the sample
#' @param mito Default = "chrM" The `seqname` of the mitochondrial
#' genome to be initialized in the
#'
#' @return An initialized mgatk object that is an S4 class
#' RangedSummarizedExperiment
#'
#' @import Matrix
#' @importFrom data.table fread
#' @importFrom tools file_ext
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @examples
#'
#' datafile <- paste0(system.file('extdata',package='mgatk'), "/mgatk_CML.txt.gz")
#' mitoSE <- importMito.txt(datafile)
#'
#' @export
setGeneric(name = "importMito.txt", def = function(datafile, pos = 1, allele = 2, sample = 3, coverage = 4, mito = "chrM")
  standardGeneric("importMito.txt"))

#' @rdname importMito.txt
setMethod("importMito.txt", signature("character", "ANY", "ANY", "ANY", "ANY", "ANY"),
          definition = function(datafile, pos = 1, allele = 2, sample = 3, coverage = 4, mito = "chrM"){

  stopifnot(length(file) == 1)

  # fread datafile in
  if(tools::file_ext(datafile) == "gz"){
    dt <- fread(paste0("zcat < ", datafile), stringsAsFactors = TRUE)
  } else if(tools::file_ext(datafile) %in% c("txt", "csv", "tsv")){
    dt <- fread(paste0(datafile), stringsAsFactors = TRUE)
  } else{
    stop("Provide a valid file format (.gz, .txt, .csv, or .tsv)")
  }

  stopifnot(all(dim(dt)[2] >= c(pos, allele, sample, coverage)))

  # Handle column naming based on user input
  ct <- paste0("X", 1:dim(dt)[2])
  ct[pos] <- "pos"; ct[allele] <- "allele"
  ct[sample] <- "sample"; ct[coverage] <- "coverage"
  colnames(dt) <- ct

  # Set up downstream processing
  samples <- levels(dt[[sample]])
  maxpos <- max(dt$pos)
  sdt <- lapply(split(1:nrow(dt), dt[[allele]]), function(x) dt[x])

  # cast the split data tables into sparse matrices
  # extra element appended to all vectors is for correct
  # dimension defintions
  sparseMatrixMake <- function(letter){
    Matrix::sparseMatrix(
      i = c(sdt[[letter]][["pos"]],maxpos),
      j = c(as.numeric(sdt[[letter]][["sample"]]), 1),
      x = c(sdt[[letter]][["coverage"]],0)
    )
  }

  SMlist <- lapply(c("A","C","G","T"), sparseMatrixMake)
  if(length(unique(lapply(SMlist, dim))) > 1) stop("Error importing")
  names(SMlist) <- c("A","C","G","T")

  # Make summary matrices
  SMlist$coverage <- SMlist[["A"]] + SMlist[["C"]] + SMlist[["G"]]+ SMlist[["T"]]

  # Make GRanges
  row_g <- GRanges(seqnames = mito, IRanges(1:maxpos, width = 1))

  # Make a summarized experiment
  SE <- SummarizedExperiment(
    assays = SMlist,
    colData = DataFrame(samples = samples
                        ),
    rowData = row_g
  )

  # Annotate with prevalence
  mcols(rowRanges(SE)) <- DataFrame(
    rowRanges(SE),
    prevalenceOrder = orderACTG(SE)
  )

  return(SE)
})
