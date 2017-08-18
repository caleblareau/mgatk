
context("Verify importing and blacklisting works")

path <-paste0(system.file('extdata',package='mgatk'),"/glioma/final/")
Afile <-paste0(path, "mgatk.A.txt")
Cfile <-paste0(path, "mgatk.C.txt")
Gfile <-paste0(path, "mgatk.G.txt")
Tfile <-paste0(path, "mgatk.T.txt")

coverageFile <- paste0(path, "mgatk.coverage.txt")
depthFile <- paste0(path, "mgatk.depthTable.txt")
referenceAlleleFile <- paste0(path, "mgatk.chrM_refAllele.txt")
mitoSE1 <- importMito.explicit(Afile, Cfile, Gfile, Tfile,
   coverageFile, depthFile, referenceAlleleFile)

folder <-paste0(system.file('extdata',package='mgatk'),"/glioma/final")
mitoSE2 <- importMito(folder)


test_that("Explicit and casual imports work the same", {
 expect_equal(dim(mitoSE1)[1], dim(mitoSE2)[1])
 expect_equal(dim(mitoSE1)[2], dim(mitoSE2)[2])
})

test_that("Blacklist subsetting works", {
 mitoSEbl <- filterKnownBlacklist(mitoSE2, "hg19_TF1")
 expect_equal(dim(mitoSEbl)[2], dim(mitoSE2)[2])
})
