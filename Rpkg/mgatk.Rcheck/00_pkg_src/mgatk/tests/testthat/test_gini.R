
context("Verify Cpp implementations of Gini index works / agrees with R")
set.seed(14651)

x <- matrix(runif(200), ncol = 5) # 5 samples
v <- runif(100)

rowsR <- apply(x, 1, ineq::ineq, "rows", "Gini")
colsR <- apply(x, 2, ineq::ineq, "cols", "Gini")
rowsCpp <- giniRows(x)
colsCpp <- giniCols(x)

test_that("C++ and R versions give the same result", {
 expect_true(all(round(rowsR,3)==round(rowsCpp,3)))
 expect_true(all(round(colsR,3)==round(colsCpp,3)))
 expect_equal(ineq::ineq(v, "base", "Gini"), giniCpp(v))
})
