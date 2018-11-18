#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// Sort vector
NumericVector stl_sort(NumericVector x) {
   NumericVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}

// Remove NAs from vector
NumericVector naomit(NumericVector x){
  std::vector<double> r(x.size());
  int k=0;
    for (int i = 0; i < x.size(); ++i) {
      if (x[i]==x[i]) {
        r[k] = x[i];
        k++;
      }
    }
 r.resize(k);
 return wrap(r);
}

// Compute weighted ( by index )  sum of ordered vector
double osum (NumericVector x){
  int n = x.length();
  double total = 0;
  for (int i = 0; i < n ; ++i) {
    total += ((i+1)*x(i));
  }
  return total;
}

//' Calculate a vector of gini indices from the rows of a matrix
//'
//' This function computes the row-wise Gini indices
//' sensitive to NAs
//'
//' @param x An m x n numeric matrix
//' @return A vector of length m with the Gini Index
//' computed for each row across all of the columns
//'
//' @examples
//'
//' x <- matrix(runif(1000), nrow = 20) # 50 samples
//' gi <- giniRows(x)
//'
//' @export
// [[Rcpp::export]]
NumericVector giniRows (NumericMatrix x){

  int nrow = x.nrow();
  Rcpp::NumericVector xx(nrow);

  // https://github.com/cran/ineq/blob/master/R/ineq.R
  for (int i = 0 ; i < nrow; i++){
    NumericVector sortedrow = stl_sort(naomit(x(i,_)));
    float G = ((2*osum(sortedrow))/sum(sortedrow)) - (sortedrow.length() + 1);
    xx(i) = G/sortedrow.length();
  }

  return (xx) ;
}

//' Calculate a vector of gini indices from the rows of a matrix
//'
//' This function computes the column-wise Gini indices
//' sensitive to NAs
//'
//' @param x An m x n numeric matrix
//' @return A vector of length n with the Gini Index
//' computed for each column across all of the rows
//'
//' @examples
//'
//' x <- matrix(runif(1000), nrow = 20) # 50 samples
//' gi <- giniCols(x)
//'
//'
//' @export
// [[Rcpp::export]]
NumericVector giniCols (NumericMatrix x){

  int ncol = x.ncol();
  Rcpp::NumericVector xx(ncol);

  // https://github.com/cran/ineq/blob/master/R/ineq.R
  for (int i = 0 ; i < ncol; i++){
    NumericVector sortedrow = stl_sort(naomit(x(_,i)));
    float G = ((2*osum(sortedrow))/sum(sortedrow)) - (sortedrow.length() + 1);
    xx(i) = G/sortedrow.length();
  }

  return (xx) ;
}

//' Calculate the Gini Index using Rcpp
//'
//' Given a vector of numbers, calcuate the
//' Gini index being sensitive to NAs
//'
//' @param x A numeric vector
//' @return A float of the Gini index where values are
//' between [0,1]
//'
//' @examples
//'
//' x <- runif(1000)
//' gi <- giniCpp(x)
//'
//' @export
// [[Rcpp::export]]
float giniCpp (NumericVector x){

  int n = x.length();
  NumericVector sortedrow = stl_sort(naomit(x));
  float G = ((2*osum(sortedrow))/sum(sortedrow)) - (n + 1);
  float gini = G/n;
  return (gini);

}
