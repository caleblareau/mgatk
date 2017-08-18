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
#' produces a RangedSummarizedExperiment object that serves
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
#' @param maxFreq Default = 0.5 The maximum allele frequency for
#' a called variant
#'
#'
#' @param minCoverage Default = 10 The maximum number of total
#' reads for the variant to be carried forward
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
#' path <-paste0(system.file('extdata',package='mgatk'),"/glioma/final/")
#'
#' Afile <-paste0(path, "mgatk.A.txt")
#' Cfile <-paste0(path, "mgatk.C.txt")
#' Gfile <-paste0(path, "mgatk.G.txt")
#' Tfile <-paste0(path, "mgatk.T.txt")
#'
#' coverageFile <- paste0(path, "mgatk.coverage.txt")
#' depthFile <- paste0(path, "mgatk.depthTable.txt")
#' referenceAlleleFile <- paste0(path, "mgatk.chrM_refAllele.txt")
#'
#' mitoSE <- importMito.explicit(Afile, Cfile, Gfile, Tfile,
#'   coverageFile, depthFile, referenceAlleleFile)
#' dim(mitoSE)
#'
#' @export
setGeneric(name = "importMito.explicit",
           def = function(Afile, Cfile, Gfile, Tfile,
                          coverageFile, depthFile, referenceAlleleFile,
                          mitoChr = "chrM", maxFreq = 0.5, minCoverage = 10)

  standardGeneric("importMito.explicit"))

#' @rdname importMito.explicit
setMethod("importMito.explicit", signature("character", "character", "character", "character",
                                      "character", "character", "character", "ANY", "ANY", "ANY"),
          definition = function(Afile, Cfile, Gfile, Tfile,
                          coverageFile, depthFile, referenceAlleleFile,
                          mitoChr = "chrM", maxFreq = 0.5, minCoverage = 10){

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
      cov <- fread(paste0("zcat < ", file), stringsAsFactors = TRUE)
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      cov <- fread(paste0(file), stringsAsFactors = TRUE)
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

  importSM <- function(file){
    # fread the individual variant calls in
    if(tools::file_ext(file) == "gz"){
      dt <- fread(paste0("zcat < ", file), stringsAsFactors = TRUE)
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      dt <- fread(paste0(file), stringsAsFactors = TRUE)
    } else{
      stop("Provide a valid file format for the variant call file (.gz, .txt, .csv, or .tsv)")
    }

    dt$sample <- factor(dt$sample, levels = samplesOrder)

    mat <- Matrix::sparseMatrix(
      i = c(dt[[1]],maxpos),
      j = c(as.numeric(dt[[2]]), maxsamples),
      x = c(dt[[3]],0)
    ) / (covmat + 0.00000000001)
    remove(dt)
    return(round(mat,4))
  }

  ACGT <- lapply(variantFiles, importSM)
  names(ACGT) <- c("A", "C", "G", "T")

  # Call variants
  freqMat <- sapply(ACGT, Matrix::rowMeans)
  ref <- importDT(referenceAlleleFile)

  # Make a matrix with 0s at the reference allele
  refMat <- Matrix::sparseMatrix(
    i = (ref[[1]])[1:maxpos],
    j = as.numeric(factor(ref[[2]], levels = c("A", "C", "G", "T"))[1:maxpos]),
    x = -1
  ) + 1
  colnames(refMat) <- c("A", "C", "G", "T")

  # Zeros where reference allele or > maxFreq
  altAllele <-  c("A", "C", "G", "T")[apply(freqMat*refMat*(freqMat < maxFreq), 1, which.max)]
  As <- which(altAllele == "A"); Cs <- which(altAllele == "C")
  Gs <- which(altAllele == "G"); Ts <- which(altAllele == "T")

  # Get allele-specific frequences
  freq <- rbind(ACGT[["A"]][As,], ACGT[["C"]][Cs,],
                ACGT[["G"]][Gs,], ACGT[["T"]][Ts,])
  remove(ACGT)
  allorder <- c(As, Cs, Gs, Ts)
  freq <- freq[order(allorder),]

  # Add column meta data
  depth <- data.frame(importDT(depthFile))
  sdf <- merge(data.frame(sample = samplesOrder), depth, by.x = "sample", by.y = "V1")
  rownames(sdf) <- samplesOrder
  colnames(sdf) <- c("sample", "depth")

  # Make GRanges and include the allele chosen for analysis
  row_g <- GenomicRanges::GRanges(seqnames = mitoChr,
                   IRanges::IRanges(1:maxpos, width = 1),
                   mcols = DataFrame(refAllele = ref[[2]][1:maxpos], altAllele = altAllele))

  # Make a summarized experiment
  SE <- SummarizedExperiment::SummarizedExperiment(
    assays = list("frequency" = freq, "coverage" = covmat),
    colData = DataFrame(sdf),
    rowData = row_g
  )

  # Remove only positions that we got coverage for
  return(SE[Matrix::rowSums(covmat) > 10, ])
})


#' Import data from folder output of mgatk python run.
#'
#' \code{importMito} takes a path to a folder that contains
#' all the requisitive input files for setting up an analysis
#' in mgatk. This is typically the `final` folder when running
#' the python mgatk utility.
#'
#' @param minCoverage Default = 10 The maximum number of total
#' reads for the variant to be carried forward
#'
#' @param ... Additional parameters to pass to the
#' importMito.explict function
#'
#' @return An initialized mgatk object that is an S4 class
#' RangedSummarizedExperiment that includes variant
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
setGeneric(name = "importMito", def = function(folder, ...)

  standardGeneric("importMito"))

#' @rdname importMito
setMethod("importMito", signature("character"),
          definition = function(folder){

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
                      coverageFile, depthFile, referenceAlleleFile, mitoChr...)
  return(SE)
})
