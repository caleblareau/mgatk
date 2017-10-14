
context("Verify importing and blacklisting works")

path <-paste0(system.file('extdata',package='mgatk'),"/glioma/final/")
Afile <-paste0(path, "glio.A.txt")
Cfile <-paste0(path, "glio.C.txt")
Gfile <-paste0(path, "glio.G.txt")
Tfile <-paste0(path, "glio.T.txt")

coverageFile <- paste0(path, "glio.coverage.txt")
depthFile <- paste0(path, "glio.depthTable.txt")
referenceAlleleFile <- paste0(path, "chrM_refAllele.txt")
mitoMAE1 <- importMito.explicit(Afile, Cfile, Gfile, Tfile,
   coverageFile, depthFile, referenceAlleleFile)

folder <-paste0(system.file('extdata',package='mgatk'),"/glioma/final")
mitoMAE2 <- importMito(folder)

test_that("Explicit and casual imports work the same", {
 expect_equal(dim(mitoMAE1@colData)[1], dim(mitoMAE2@colData)[1])
 expect_equal(dim(mitoMAE1@colData)[2], dim(mitoMAE2@colData)[2])
})

test_that("Reference alleles and alternate alleles are the same", {
 x <- mcols(rowRanges(mitoMAE1[["alleles"]]))
 expect_equal(sum(x$refAllele == x$altAllele), 0)
})

test_that("Blacklist subsetting works", {
 mitoMAEbl <- filterKnownBlacklist(mitoMAE2, "hg19_TF1")
 expect_equal(dim(mitoMAEbl@ExperimentList[["coverage"]])[2],
              dim(mitoMAE2@ExperimentList[["coverage"]])[2])

 expect_less_than(dim(mitoMAEbl@ExperimentList[["coverage"]])[1],
              dim(mitoMAE2@ExperimentList[["coverage"]])[1])
})
