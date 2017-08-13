
context("distance functions work")
set.seed(14651)

# Essentially Jacob's original implementation
dist_mito <- function(df.ratio, df.coverage = NULL, dist.type = "abs") {
  # Setup output distance matrix
  dist.out <- matrix(nrow = dim(df.ratio)[2], ncol = dim(df.ratio)[2])

  # If coverage matrix is not provided, assume equal weights
  if (is.null(df.coverage)) {
    df.coverage <- matrix(1L, nrow = length(dist.rows), ncol = length(dist.rows))
  }

  # Loop through each
  for (i in seq(1:dim(df.ratio)[2])) {
    wi <- 1 / df.coverage[, i] # inverse coverage of row i
    for (j in seq(1:dim(df.ratio)[2])) {
      wj <- 1 / df.coverage[, j] #inverse coverage of row j
      w <- 1 / (wi + wj) # weight is inverse of sum of coverages
      if (dist.type == "euclidean") {
        # Calculate euclidean distance weighting each ratio
        dist.vect <- as.vector(na.omit(w * sqrt(df.ratio[, i] - df.ratio[, j])^2))
        # Calculate sqrt distance weighting each ratio
      } else if (dist.type=="sqrt") {
        dist.vect <- as.vector(na.omit(w * sqrt(abs(df.ratio[, i] - df.ratio[, j]))))
      } else {
        # Calculate absolute distance weighting each ratio
        dist.vect <- as.vector(na.omit(w * abs(df.ratio[, i] - df.ratio[, j])))
      }
      dist.vect.length <- length(dist.vect)
      # Normalize by total number of ratios compared
      dist.val <- sum(dist.vect) / dist.vect.length
      # Save as distance between i and j
      dist.out[i, j] <- dist.val
    }
  }
  return(dist.out)
}


x <- matrix(runif(25), ncol = 5) # 5 samples
y <- matrix(rpois(25,10), ncol = 5) # 5 samples

abs_d1 <- dist_mito(x, y, "abs")
abs_d2 <- mgatk::computeWeightedDistance(x, y, "abs")

tm <- microbenchmark::microbenchmark(
  dist_mito(x, y, "abs"),
  computeWeightedDistance(x, y, "abs"), times=100)
sdm <- summary(tm)

test_that("C++ and R versions give the same result", {
 expect_true(all(round(abs_d1,3)==round(abs_d2,3)))
})





