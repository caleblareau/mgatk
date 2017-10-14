#' @include coreFunctions.R
NULL

#' Import data from plain text files explicitly point
#' to the requisite input files and render a RangedSummarized
#' Experiment.
#'
#' \code{importMito.explicit} takes a sparse matrix file of position,
#' sample, and count (non-zero) in some order as variant calls
#' (one file per allele) as well as
#' a separate file of the position and sample coverage and
#' produces a MultiAssayExperiment object that serves
#' as the backbone of the R interface for mgatk.
#'
#' The filepath of a plain text or gzipped file
#' that contains data in a long matrix format of the position,
#' allele, sample, and coverage. This will be imported and a
#' data object will be rendered that enables downstream
#' mgatk analyses
#'
#' @param Afile Contains the position, sample, count of A alleles
#' @param Cfile Contains the position, sample, count of C alleles
#' @param Gfile Contains the position, sample, count of G alleles
#' @param Tfile Contains the position, sample, count of T alleles
#'
#' @param coverageFile The filepath of a plain text or gzipped file
#' that contains data in a sparse matrix format of the position,
#' sample, and coverage. This will be imported and a
#' data object will be rendered that enables downstream
#' mgatk analyses.
#'
#' @param depthFile The filepath of a plain text or gzipped file
#' that contains two columns indicating the sample name and the
#' second the mean coverage of the sample about the mitochondrial
#' genome
#'
#' @param referenceAlleleFile The filepath of a plain text or gzipped file
#' that contains two columns indicating the position and designated
#' reference allele.
#'
#' @param mitoChr Default = "chrM" The `seqname` of the mitochondrial
#' genome to be initialized in the RangedSummarizedExperiment object
#'
#' @return An initialized mgatk object that is an S4 class
#' MultiAssayExperiment that includes variant
#' count and coverage as assay slots.
#'
#' @import Matrix
#' @importFrom data.table fread dcast.data.table
#' @importFrom tools file_ext
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#' @import GenomicRanges
#' @import MultiAssayExperiment
#' @examples
#'
#' path <-paste0(system.file('extdata',package='mgatk'),"/glioma/final/")
#'
#' Afile <-paste0(path, "glio.A.txt")
#' Cfile <-paste0(path, "glio.C.txt")
#' Gfile <-paste0(path, "glio.G.txt")
#' Tfile <-paste0(path, "glio.T.txt")
#'
#' coverageFile <- paste0(path, "glio.coverage.txt")
#' depthFile <- paste0(path, "glio.depthTable.txt")
#' referenceAlleleFile <- paste0(path, "chrM_refAllele.txt")
#'
#' mitoSE <- importMito.explicit(Afile, Cfile, Gfile, Tfile,
#'   coverageFile, depthFile, referenceAlleleFile)
#' dim(mitoSE)
#'
#' @export
setGeneric(name = "importMito.explicit",
           def = function(Afile, Cfile, Gfile, Tfile,
                          coverageFile, depthFile, referenceAlleleFile,
                          mitoChr = "chrM")

  standardGeneric("importMito.explicit"))

#' @rdname importMito.explicit
setMethod("importMito.explicit", signature("character", "character", "character", "character",
                                      "character", "character", "character", "ANY"),
          definition = function(Afile, Cfile, Gfile, Tfile,
                          coverageFile, depthFile, referenceAlleleFile,
                          mitoChr = "chrM"){

  variantFiles <- list(Afile, Cfile, Gfile, Tfile)
  metaFiles <- list(coverageFile, depthFile, referenceAlleleFile)

  nullout <- lapply(c(variantFiles, metaFiles), function(file){
    stopifnot(length(file) == 1)
  })

  # Set up downstream processing including robust ordering
  # The coverage file could have slightly more variants /
  # individual samples depending on the calls, so base it
  # of of them
  importDT <- function(file){
    if(tools::file_ext(file) == "gz"){
      cov <- data.table::fread(paste0("zcat < ", file), stringsAsFactors = TRUE)
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      cov <- data.table::fread(paste0(file), stringsAsFactors = TRUE)
    } else{
      stop("Provide a valid file format for the  file (.gz, .txt, .csv, or .tsv)")
    }
  }

  cov <- importDT(coverageFile)

  samplesOrder <- levels(cov[[2]])
  maxpos <- max(cov[[1]])
  maxsamples <- length(samplesOrder)

  # make coverage a sparse matrix
  covmat <- Matrix::sparseMatrix(
    i = cov[[1]],
    j = as.numeric(cov[[2]]),
    x = cov[[3]]
  )
  remove(cov)

  # Import Counts and BAQ
  importSMs <- function(file){
    # fread the individual variant calls in
    if(tools::file_ext(file) == "gz"){
      dt <-  data.table::fread(paste0("zcat < ", file), stringsAsFactors = TRUE)
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      dt <-  data.table::fread(paste0(file), stringsAsFactors = TRUE)
    } else{
      stop("Provide a valid file format for the variant call file (.gz, .txt, .csv, or .tsv)")
    }

    dt$sample <- factor(dt$sample, levels = samplesOrder)

    counts <- Matrix::sparseMatrix(
      i = c(dt[[1]],maxpos),
      j = c(as.numeric(dt[[2]]), maxsamples),
      x = c(dt[[3]],0)
    )

    BAQ <- Matrix::sparseMatrix(
      i = c(dt[[1]],maxpos),
      j = c(as.numeric(dt[[2]]), maxsamples),
      x = c(dt[[4]],0)
    )
    remove(dt)
    return(list("counts" = counts, "BAQ" = BAQ))
  }

  ACGT <- lapply(variantFiles, importSMs)
  names(ACGT) <- c("A", "C", "G", "T")

  # Make a long matrix of BAQ and Counts for non-reference alleles
  ref <- importDT(referenceAlleleFile)
  whichA <- which(ref[["V2"]][1:maxpos] != "A")
  whichC <- which(ref[["V2"]][1:maxpos] != "C")
  whichG <- which(ref[["V2"]][1:maxpos] != "G")
  whichT <- which(ref[["V2"]][1:maxpos] != "T")

  longBAQ <- rbind(
    ACGT[["A"]][["BAQ"]][whichA,],
    ACGT[["C"]][["BAQ"]][whichC,],
    ACGT[["G"]][["BAQ"]][whichG,],
    ACGT[["T"]][["BAQ"]][whichT,]
  )

  longCounts <- rbind(
    ACGT[["A"]][["counts"]][whichA,],
    ACGT[["C"]][["counts"]][whichC,],
    ACGT[["G"]][["counts"]][whichG,],
    ACGT[["T"]][["counts"]][whichT,]
  )
  letterz <- c(rep("A", length(whichA)), rep("C", length(whichC)), rep("G", length(whichG)), rep("T", length(whichT)))
  remove(ACGT)

  # Create colData
  depth <- data.frame(importDT(depthFile))
  sdf <- merge(data.frame(sample = samplesOrder), depth, by.x = "sample", by.y = "V1")
  rownames(sdf) <- samplesOrder
  colnames(sdf) <- c("sample", "depth")

  # Make row Ranges for each object
  row_g_cov <- GenomicRanges::GRanges(seqnames = mitoChr,
                   IRanges::IRanges(1:maxpos, width = 1))
  GenomicRanges::mcols(row_g_cov) <- data.frame(refAllele = ref[[2]][1:maxpos])

  row_g_allele <- GenomicRanges::GRanges(seqnames = mitoChr,
                   IRanges::IRanges(1:maxpos, width = 1))[c(whichA, whichC, whichG, whichT)]

  GenomicRanges::mcols(row_g_allele) <- data.frame(refAllele = (ref[[2]][1:maxpos])[c(whichA, whichC, whichG, whichT)],
                                                   altAllele = letterz)

  # Make summarized experiments and
  coverage <- SummarizedExperiment::SummarizedExperiment(
    assays = list("coverage" = covmat),
    colData = S4Vectors::DataFrame(sdf),
    rowData = row_g_cov
  )

  alleles <- SummarizedExperiment::SummarizedExperiment(
    assays = list("BAQ" = longBAQ, "counts" = longCounts),
    colData = S4Vectors::DataFrame(sdf),
    rowData = row_g_allele
  )

  remove(longBAQ)
  remove(longCounts)

  MAE <- MultiAssayExperiment::MultiAssayExperiment(
    list("alleles" = alleles, "coverage" = coverage),
    colData = S4Vectors::DataFrame(sdf)
  )
  return(MAE)
})


#' Import data from folder output of mgatk python run.
#'
#' \code{importMito} takes a path to a folder that contains
#' all the requisitive input files for setting up an analysis
#' in mgatk. This is typically the `final` folder when running
#' the python mgatk utility.
#'
#' @param folder Filepath to folder (likely ending in 'final')
#' that contains all the vital input / output files for mgatk
#'
#' @param ... Additional parameters to pass to the
#' importMito.explict function
#'
#' @return An initialized mgatk object that is an S4 class
#' MultiAssayExperiment that includes variant
#' frequency and coverage as assay slots.
#'
#' @seealso importMito.explicit
#'
#' @examples
#' folder <-paste0(system.file('extdata',package='mgatk'),"/glioma/final")
#'
#' mitoSE <- importMito(folder)
#' dim(mitoSE)
#'
#' @export
setGeneric(name = "importMito", def = function(folder, ...) standardGeneric("importMito"))

#' @rdname importMito
setMethod("importMito", signature("character"), definition = function(folder, ...){

  files <- list.files(folder, full.names = TRUE)

  checkGrep <- function(hit){
    if(length(hit) != 1){
      stop("Improper folder specification; file missing / extra file present. See documentation")
    } else {
      return(hit)
    }
  }

  # Set up file paths
  Afile <- files[checkGrep(grep(".A.txt", files))]
  Cfile <- files[checkGrep(grep(".C.txt", files))]
  Gfile <- files[checkGrep(grep(".G.txt", files))]
  Tfile <- files[checkGrep(grep(".T.txt", files))]
  coverageFile <- files[checkGrep(grep(".coverage.txt", files))]
  depthFile <- files[checkGrep(grep(".depthTable.txt", files))]
  referenceAlleleFile <- files[checkGrep(grep("refAllele.txt", files))]

  # Parse out the mitochondrial genome name from the file name
  sv <- strsplit(gsub("_refAllele.txt", "", basename(referenceAlleleFile)), split = "[.]")[[1]]
  mitoChr <- sv[length(sv)]

  SE <- importMito.explicit(Afile, Cfile, Gfile, Tfile,
                      coverageFile, depthFile, referenceAlleleFile, mitoChr, ...)
  return(SE)
})
