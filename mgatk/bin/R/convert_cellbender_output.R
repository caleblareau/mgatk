library(dplyr)
library(Matrix)
library(tidyr)
library(optparse)
library(SummarizedExperiment)

option_list <- list(
  make_option(c("-i","--input"),type="character",help="input MGATK directory",metavar="file"),
  make_option(c("-c","--cellbender_h5"),type="character",help="cellbender output file (all cells!)",metavar="file"),
  make_option(c("-o","--output"),type="character",help="output file with filtered counts",metavar="file"),
  make_option("--cells",type="character",help="barcodes.tsv file with accepted cell barcodes",metavar="file")
)

opt_parser <- OptionParser(option_list = option_list)
options <- parse_args(opt_parser)

if (is.null(options$input)) {
  print_help(opt_parser)
  stop("please supply input directory",call.=FALSE)
}

if (is.null(options$output)) {
  print_help(opt_parser)
  stop("please supply output directory",call.=FALSE)
}

if (is.null(options$cellbender_h5)) {
  print_help(opt_parser)
  stop("please supply output directory",call.=FALSE)
}

ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (hdf5r::existsGroup(infile, 'matrix')) {
    # cellranger version 3
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

SparseMatrixFromBaseCounts <- function(basecounts, cells, dna.base, maxpos) {

  # Vector addition guarantee correct dimension
  fwd.mat <- sparseMatrix(
    i = c(basecounts$pos,maxpos),
    j = c(cells[basecounts$cellbarcode],1),
    x = c(basecounts$plus,0)
  )
  colnames(x = fwd.mat) <- names(x = cells)
  rownames(x = fwd.mat) <- paste(
    dna.base,
    seq_len(length.out = nrow(fwd.mat)),
    "fwd",
    sep = "-"
  )
  rev.mat <- sparseMatrix(
    i = c(basecounts$pos,maxpos),
    j = c(cells[basecounts$cellbarcode],1),
    x = c(basecounts$minus,0)
  )
  colnames(x = rev.mat) <- names(x = cells)
  rownames(x = rev.mat) <- paste(
    dna.base,
    seq_len(length.out = nrow(rev.mat)),
    "rev",
    sep = "-"
  )
  return(list(fwd.mat, rev.mat))
}

ReadMGATK <- function(dir, verbose = TRUE) {

  if (!dir.exists(paths = dir)) {
    stop("Directory not found")
  }
  a.path <- list.files(path = dir, pattern = "*.A.txt.gz", full.names = TRUE)
  c.path <- list.files(path = dir, pattern = "*.C.txt.gz", full.names = TRUE)
  t.path <- list.files(path = dir, pattern = "*.T.txt.gz", full.names = TRUE)
  g.path <- list.files(path = dir, pattern = "*.G.txt.gz", full.names = TRUE)

  refallele.path <- list.files(
    path = dir,
    pattern = "*_refAllele.txt*", # valid mtDNA contigs include chrM and MT
    full.names = TRUE
  )

  # The depth file lists all barcodes that were genotyped
  depthfile.path <- list.files(
    path = dir,
    pattern = "*.depthTable.txt",
    full.names = TRUE
  )

  if (verbose) {
    message("Reading allele counts")
  }
  column.names <- c("pos", "cellbarcode", "plus", "minus")
  a.counts <- read.table(
    file = a.path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = column.names
  )
  c.counts <- read.table(
    file = c.path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = column.names
  )
  t.counts <- read.table(
    file = t.path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = column.names
  )
  g.counts <- read.table(
    file = g.path,
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = column.names
  )

  if (verbose) {
    message("Reading metadata")
  }
  refallele <- read.table(
    file = refallele.path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("pos", "ref")
  )
  refallele$ref <- toupper(x = refallele$ref)
  depth <- read.table(
    file = depthfile.path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("cellbarcode", "mito.depth"),
    row.names = 1
  )
  cellbarcodes <- unique(x = rownames(depth))
  cb.lookup <- seq_along(along.with = cellbarcodes)
  names(cb.lookup) <- cellbarcodes

  if (verbose) {
    message("Building matrices")
  }

  maxpos <- dim(refallele)[1]
  a.mat <- SparseMatrixFromBaseCounts(
    basecounts = a.counts, cells = cb.lookup, dna.base = "A", maxpos = maxpos
  )
  c.mat <- SparseMatrixFromBaseCounts(
    basecounts = c.counts, cells = cb.lookup, dna.base = "C", maxpos = maxpos
  )
  t.mat <- SparseMatrixFromBaseCounts(
    basecounts = t.counts, cells = cb.lookup, dna.base = "T", maxpos = maxpos
  )
  g.mat <- SparseMatrixFromBaseCounts(
    basecounts = g.counts, cells = cb.lookup, dna.base = "G", maxpos = maxpos
  )

  counts <- rbind(a.mat[[1]], c.mat[[1]], t.mat[[1]], g.mat[[1]],
                  a.mat[[2]], c.mat[[2]], t.mat[[2]], g.mat[[2]])

  return(list("counts" = counts, "depth" = depth, "refallele" = refallele))
}

raw <- ReadMGATK(file.path(options$input))
rc <- raw$counts
CB <- ReadCB_h5(options$cellbender_h5)
pos <- intersect(row.names(rc),row.names(CB))
cells <- intersect(colnames(rc),colnames(CB))

# if --cells is specified, restrict to those cells
if (!is.null(options$cells)) {
   barcodes <- read.table(options$cells,header=FALSE,row.names=NULL)
   cells <- intersect(cells,barcodes[,1])
}

rc <- rc[,cells]
# replace raw values by cellbender estimates
rc[pos,cells] <- CB[pos,cells]

# at this point we probably should ensure that zero counts on one strand with non-zero counts on the other are actually zero and not NA; otherwise this triggers an error in Signac: https://github.com/timoast/signac/blob/2e4e18354297adc4b069b48a7794aa1efed35016/R/mito.R#L603
assays <- list()
for (let in c('A','C','G','T')) {
  assays[[paste0(let,'_counts_fw')]] <- rc[grepl(paste0(let,'-[0-9]*-fwd'),row.names(rc)),]
  row.names(assays[[paste0(let,'_counts_fw')]]) <- raw$refallele$pos
  assays[[paste0(let,'_counts_rev')]] <- rc[grepl(paste0(let,'-[0-9]*-rev'),row.names(rc)),]
  row.names(assays[[paste0(let,'_counts_rev')]]) <- raw$refallele$pos
}
assays[['coverage']] <- Reduce('+', assays)
tmp <- SummarizedExperiment(assays=assays,
                            colData=data.frame(sample=colnames(rc),
                              depth=Matrix::colSums(rc)/(dim(rc)[1]/8)),
                            rowData=data.frame(refAllele=raw$refallele$ref,
                              row.names=raw$refallele$pos))

saveRDS(tmp, options$output)
