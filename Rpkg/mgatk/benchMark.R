library(microbenchmark)
library(mgatk)
library(BuenColors)

x <- matrix(rnorm(100000), ncol = 500) # 50 samples
y <- matrix(rpois(100000, 10), ncol = 500) # 50 samples

dim(x)
dim(y)

library(ggplot2)
tm <- microbenchmark(dist_mito(x, y), computeWeightedDistance(x, y),  times=3)
autoplot(tm) + pretty_plot()



