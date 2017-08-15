#' @include coreFunctions.R
NULL

#' Import data from plain text file into
#'
#' \code{importMito.txt} takes a long file of position, allele,
#' sample, and count (non-zero) in some order as variant calls
#' and a separate file of the position and sample coverage and
#' produces a RangedSummarizedExperiment object that serves
#' as the backbone of the R interface for mgatk.
#'
#' There are seven total v.* and c.* parameters that are the specified
#' index associated with the feature of interest. If the data was
#' pre-processed with the mgatk python package, you shouldn't need
#' to modify these at all.
#'
#' @param variantCallFile The filepath of a plain text or gzipped file
#' that contains data in a long matrix format of the position,
#' allele, sample, and coverage. This will be imported and a
#' data object will be rendered that enables downstream
#' mgatk analyses
#' @param coverageFile The filepath of a plain text or gzipped file
#' that contains data in a long matrix format of the position,
#' sample, and coverage. This will be imported and a
#' data object will be rendered that enables downstream
#' mgatk analyses
#'
#' @param mito Default = "chrM" The `seqname` of the mitochondrial
#' genome to be initialized in the
#' @param keepACGT Default = FALSE Keep assays with sparse matrix
#' counts of the ACGT alleles per base / individual.
#'
#' @param v.pos Default = 1 column index of the genomic position
#' in the variant calls file
#' @param v.allele Default = 2 column index of the allele (A/C/G/T)
#' in the variant calls file
#' @param v.sample Default = 3 column index of the sample name
#' in the variant calls file
#' @param v.count Default = 4 column index that shows the number
#' of reads for the sample in the variant calls file
#' @param c.pos Default = 1 column index of the genomic position
#' in the coverage file
#' @param c.sample Default = 2 column index of the sample name
#' in the coverage file
#' @param c.coverage Default = 3 column index that shows the number
#' of reads covering a sample / position
#'
#' @return An initialized mgatk object that is an S4 class
#' RangedSummarizedExperiment that includes variant
#' frequency and coverage as assay slots.
#'
#' @import Matrix
#' @importFrom data.table fread dcast.data.table
#' @importFrom tools file_ext
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @examples
#'
#' path <-paste0(system.file('extdata',package='mgatk'),"/mouse_tcr/")
#' variantCallFile <-paste0(path, "mouse_tcr.variantCalls.txt")
#' coverageFile <- paste0(path, "mouse_tcr.coverage.txt")
#' mitoSE <- importMito.txt(variantCallFile = variantCallFile, coverageFile = coverageFile)
#' dim(mitoSE)
#'
#' @export
setGeneric(name = "importMito.txt",
           def = function(variantCallFile, coverageFile, mito = "chrM", keepACGT = FALSE,
                          v.pos = 1, v.allele = 2, v.sample = 3, v.count = 4,
                          c.pos = 1, c.sample = 2, c.coverage = 3)
  standardGeneric("importMito.txt"))

#' @rdname importMito.txt
setMethod("importMito.txt", signature("character", "character", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY", "ANY"),
          definition = function(variantCallFile, coverageFile, mito = "chrM", keepACGT = FALSE,
                                v.pos = 1, v.allele = 2, v.sample = 3, v.count = 4,
                                c.pos = 1, c.sample = 2, c.coverage = 3){

  stopifnot(length(variantCallFile) == 1)
  stopifnot(length(coverageFile) == 1)

  # fread variantCallFile in
  if(tools::file_ext(variantCallFile) == "gz"){
    dt <- fread(paste0("zcat < ", variantCallFile), stringsAsFactors = TRUE)
  } else if(tools::file_ext(variantCallFile) %in% c("txt", "csv", "tsv")){
    dt <- fread(paste0(variantCallFile), stringsAsFactors = TRUE)
  } else{
    stop("Provide a valid file format for the variant call file (.gz, .txt, .csv, or .tsv)")
  }

  stopifnot(all(dim(dt)[2] >= c(v.pos, v.allele, v.sample, v.count)))

  # fread coverageFile in
  if(tools::file_ext(coverageFile) == "gz"){
    cov <- fread(paste0("zcat < ", coverageFile), stringsAsFactors = TRUE)
  } else if(tools::file_ext(coverageFile) %in% c("txt", "csv", "tsv")){
    cov <- fread(paste0(coverageFile), stringsAsFactors = TRUE)
  } else{
    stop("Provide a valid file format for the coverage call file (.gz, .txt, .csv, or .tsv)")
  }

  stopifnot(all(dim(cov)[2] >= c(c.pos, c.sample, c.coverage)))

  # Handle column naming based on user input
  ct <- paste0("X", 1:dim(dt)[2])
  ct[v.pos] <- "pos"; ct[v.allele] <- "allele"
  ct[v.sample] <- "sample"; ct[v.count] <- "count"
  colnames(dt) <- ct

  ct <- paste0("X", 1:dim(cov)[2])
  ct[c.pos] <- "pos";
  ct[c.sample] <- "sample"; ct[c.coverage] <- "coverage"
  colnames(cov) <- ct

  # Set up downstream processing including robust ordering
  # The coverage file could have slightly more variants / 
  # individual samples depending on the calls, so base it
  # of of them 
  
  samples <- levels(cov[[c.sample]])
  dt$sample <- factor(dt$sample, levels = samples)
  cov$sample <- factor(cov$sample, levels = samples)
  maxpos <- max(cov$pos)
  maxsamples <- length(samples)

  # determine allele with greatest count
  aldt <- dcast.data.table(dt, pos ~ allele, max, value.var = "count")
  allele <- colnames(aldt)[2:5][apply(aldt[,2:5], 1, which.max)]

  # Make vector of alleles for late
  allelesAll <- unname(allele[as.character(1:maxpos)])

  # Extract columns with variants per sample with counts
  names(allele) <- as.character(aldt[["pos"]])
  dt$chosenAllele <- allele[as.character(dt[["pos"]])]
  d <- dcast.data.table(dt[dt$allele == dt$chosenAllele], pos + sample ~ ., max, value.var = "count")
  colnames(d) <- c("pos", "sample", "count")
  d$sample <- factor(d$sample, levels = samples)

  # Add the frequency, coverage sensitive to the order
  # extra element appended to all vectors is for correct dim
  SMlist <- list()

  SMlist$coverage <- Matrix::sparseMatrix(
    i = c(cov[["pos"]],maxpos),
    j = c(as.numeric(cov[["sample"]]), maxsamples),
    x = c(cov[["coverage"]],0)
  )

  SMlist$freq <- Matrix::sparseMatrix(
    i = c(d[["pos"]],maxpos),
    j = c(as.numeric(d[["sample"]]), maxsamples),
    x = c(d[["count"]],0)
  ) / (SMlist$coverage + 0.001) # charity count preserves zeros

  if(keepACGT){

    # cast the split data tables into sparse matrices
    sdt <- lapply(split(1:nrow(dt), dt[[v.allele]]), function(x) dt[x])

    sparseMatrixMake <- function(letter){
      Matrix::sparseMatrix(
        i = c(sdt[[letter]][["pos"]],maxpos),
        j = c(as.numeric(sdt[[letter]][["sample"]]), maxsamples),
        x = c(sdt[[letter]][["count"]],0)
      )
    }

    ACGTlist <- lapply(c("A","C","G","T"), sparseMatrixMake)
    if(length(unique(lapply(ACGTlist, dim))) > 1) stop("Error importing")
    names(ACGTlist) <- c("A","C","G","T")
    SMlist <- c(SMlist, ACGTlist)
    remove(sdt)
  }
  remove(dt)

  # Make GRanges and include the allele chosen for analysis
  row_g <- GRanges(seqnames = mito, IRanges(1:maxpos, width = 1),
                   mcols = DataFrame(allele = allelesAll))

  # Make a summarized experiment
  SE <- SummarizedExperiment(
    assays = SMlist,
    colData = DataFrame(samples = samples, row.names = samples),
    rowData = row_g
  )

  # Remove only positions that we got coverage for
  return(SE[!is.na(allelesAll),])
})
