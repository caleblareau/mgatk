dist_mito <-
  function(df.ratio, df.coverage = NULL, dist.type = "abs") {
    # Setup output distance matrix
    dist.rows <- colnames(df.ratio)
    dist.out <- matrix(nrow = length(dist.rows), ncol = length(dist.rows))
    row.names(dist.out) <- dist.rows
    colnames(dist.out) <- dist.rows
    # If coverage matrix is not provided, assume equal weights
    if (is.null(df.coverage)) {
      df.coverage <- matrix(1L, nrow = length(dist.rows), ncol = length(dist.rows))
    }
    # Loop through each
    for (i in seq(1:length(dist.rows))) {
      print(i)
      wi <- 1 / df.coverage[, i] # inverse coverage of row i
      for (j in seq(1:length(dist.rows))) {
        wj <- 1 / df.coverage[, j] #inverse coverage of row j
        w <- 1 / (wi + wj) # weight is inverse of sum of coverages
        if (dist.type == "euclidean") {
          # Calculate euclidean distance weighting each ratio
          dist.vect <- as.vector(na.omit(w * (df.ratio[, i] - df.ratio[, j])^2))
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
