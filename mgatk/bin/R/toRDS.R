#!/usr/bin/env Rscript
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(SummarizedExperiment)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(data.table)))

options(warn=-1)

if(FALSE){
  
}

# Explicit import of mgatk output files
importMito.explicit <- function(Afile, Cfile, Gfile, Tfile,
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
      cov <- suppressMessages(data.table::fread(paste0("zcat < ", file), stringsAsFactors = TRUE))
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      cov <- suppressMessages(data.table::fread(paste0(file), stringsAsFactors = TRUE))
    } else{
      stop("Provide a valid file format for the  file (.gz, .txt, .csv, or .tsv)")
    }
  }
  
  cov <- importDT(coverageFile)
  
  # Make a long matrix of bq and Counts for non-reference alleles
  ref <- importDT(referenceAlleleFile)
  maxpos <- max(ref[[1]])
  
  samplesOrder <- levels(cov[[2]])
  maxsamples <- length(samplesOrder)
  
  # make coverage a sparse matrix
  covmat <- Matrix::sparseMatrix(
    i = c(cov[[1]], maxpos), 
    j = c(as.numeric(cov[[2]]), maxsamples), 
    x = c(cov[[3]], 0)
  )
  remove(cov)
  
  # Import Counts and qualities
  importSMs <- function(file){
    # fread the individual variant calls in
    if(tools::file_ext(file) == "gz"){
      dt <-  suppressMessages(data.table::fread(paste0("zcat < ", file), stringsAsFactors = TRUE))
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      dt <-  suppressMessages(data.table::fread(paste0(file), stringsAsFactors = TRUE))
    } else{
      stop("Provide a valid file format for the alt allele abundances (.gz, .txt, .csv, or .tsv)")
    }
    
    # Encore same ordering as coverage for the sample factor 
    dt[[2]] <- factor(as.character(dt[[2]]), levels = samplesOrder)
    
    if(dim(dt)[2] == 6){
      # including base qualities
      
      qual_fw <- Matrix::sparseMatrix(
        i = c(dt[[1]],maxpos),
        j = c(as.numeric(dt[[2]]), maxsamples),
        x = c(dt[[4]],0)
      )
      
      qual_rev <- Matrix::sparseMatrix(
        i = c(dt[[1]],maxpos),
        j = c(as.numeric(dt[[2]]), maxsamples),
        x = c(dt[[6]],0)
      )
      
      count_fw_idx <- 3
      count_rev_idx <- 5
      base_quals = TRUE
      
    } else {
      # No base qualities are emitted
      count_fw_idx <- 3
      count_rev_idx <- 4
      base_quals = FALSE
    }
    
    counts_fw <- Matrix::sparseMatrix(
      i = c(dt[[1]],maxpos),
      j = c(as.numeric(dt[[2]]), maxsamples),
      x = c(dt[[count_fw_idx]],0)
    )
    
    counts_rev <- Matrix::sparseMatrix(
      i = c(dt[[1]],maxpos),
      j = c(as.numeric(dt[[2]]), maxsamples),
      x = c(dt[[count_rev_idx]],0)
    )
    
    remove(dt)
    
    # Return more substantial list only if there are base qualities involved
    if(base_quals){
      return(list("counts_fw" = counts_fw, "qual_fw" = qual_fw,
                  "counts_rev" = counts_rev, "qual_rev" = qual_rev))
    } else {
      return(list("counts_fw" = counts_fw, 
                  "counts_rev" = counts_rev))
    }
  }
  
  ACGT <- lapply(variantFiles, importSMs)
  names(ACGT) <- c("A", "C", "G", "T")
  
  # Create colData
  depth <- data.frame(importDT(depthFile))
  sdf <- merge(data.frame(sample = samplesOrder), depth, by.x = "sample", by.y = "V1")
  rownames(sdf) <- samplesOrder
  colnames(sdf) <- c("sample", "depth")
  
  # Make row Ranges for each object
  row_g_cov <- GenomicRanges::GRanges(seqnames = mitoChr,
                                      IRanges::IRanges(1:maxpos, width = 1))
  GenomicRanges::mcols(row_g_cov) <- data.frame(refAllele = toupper(ref[[2]][1:maxpos]))
  
  # Export a signac object
  refallele <- data.frame(
    pos = 1:maxpos,
    ref = as.character(toupper(ref[[2]][1:maxpos])),
    stringsAsFactors = FALSE
  )
  signac_counts <-  rbind(ACGT[["A"]][["counts_fw"]],  ACGT[["C"]][["counts_fw"]],
                          ACGT[["T"]][["counts_fw"]],  ACGT[["G"]][["counts_fw"]],
                          ACGT[["A"]][["counts_rev"]],  ACGT[["C"]][["counts_rev"]],
                          ACGT[["T"]][["counts_rev"]],  ACGT[["G"]][["counts_rev"]])
  
  signac_object <- list("counts" = signac_counts, "depth" = sdf, "refallele" = refallele)
  
  # Make summarized experiments 
  if(length(ACGT[["A"]]) > 2){
    # With base qualities
    SE <- SummarizedExperiment::SummarizedExperiment(
      assays = list(
        "A_counts_fw" = ACGT[["A"]][["counts_fw"]], "A_counts_rev" = ACGT[["A"]][["counts_rev"]],
        "A_qual_fw" = ACGT[["A"]][["qual_fw"]], "A_qual_rev" = ACGT[["A"]][["qual_rev"]],
        "C_counts_fw" = ACGT[["C"]][["counts_fw"]], "C_counts_rev" = ACGT[["C"]][["counts_rev"]],
        "C_qual_fw" = ACGT[["C"]][["qual_fw"]], "C_qual_rev" = ACGT[["C"]][["qual_rev"]],
        "G_counts_fw" = ACGT[["G"]][["counts_fw"]], "G_counts_rev" = ACGT[["G"]][["counts_rev"]],
        "G_qual_fw" = ACGT[["G"]][["qual_fw"]], "G_qual_rev" = ACGT[["G"]][["qual_rev"]],
        "T_counts_fw" = ACGT[["T"]][["counts_fw"]], "T_counts_rev" = ACGT[["T"]][["counts_rev"]],
        "T_qual_fw" = ACGT[["T"]][["qual_fw"]], "T_qual_rev" = ACGT[["T"]][["qual_rev"]],
        "coverage" =  covmat
      ),
      colData = S4Vectors::DataFrame(sdf),
      rowData = row_g_cov
    )
  } else {
    # Without base qualities
    SE <- SummarizedExperiment::SummarizedExperiment(
      assays = list(
        "A_counts_fw" = ACGT[["A"]][["counts_fw"]], "A_counts_rev" = ACGT[["A"]][["counts_rev"]], 
        "C_counts_fw" = ACGT[["C"]][["counts_fw"]], "C_counts_rev" = ACGT[["C"]][["counts_rev"]], 
        "G_counts_fw" = ACGT[["G"]][["counts_fw"]], "G_counts_rev" = ACGT[["G"]][["counts_rev"]],
        "T_counts_fw" = ACGT[["T"]][["counts_fw"]], "T_counts_rev" = ACGT[["T"]][["counts_rev"]], 
        "coverage" =  covmat
      ),
      colData = S4Vectors::DataFrame(sdf),
      rowData = row_g_cov
    )
  }
  return(list(SE, signac_object))
}

#---------------------------------------
# Function to parse the folder hierarchy
#---------------------------------------

importMito <- function(folder, ...){
  
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
  
  SElist <- importMito.explicit(Afile, Cfile, Gfile, Tfile,
                                coverageFile, depthFile, referenceAlleleFile, mitoChr, ...)
  return(SElist)
}


#-----------------
# Command line i/o
#-----------------
args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]
name <- args[2]

SElist <- importMito(folder)
saveRDS(SElist[[1]], file = paste0(folder, "/", name, ".rds"))
saveRDS(SElist[[2]], file = paste0(folder, "/", name, ".signac.rds"))
