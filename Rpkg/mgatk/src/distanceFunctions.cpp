
#include <Rcpp.h>
using namespace Rcpp;

double absDist (NumericVector x, NumericVector y, NumericVector w){
  int n = y.length();
  double total = 0;
  for (int i = 0; i < n ; ++i) {
    total += w(i)*abs(x(i)-y(i));
  }
  return total;
}

double eucDist (NumericVector x, NumericVector y, NumericVector w){
  int n = y.length();
  double total = 0;
  for (int i = 0; i < n ; ++i) {
    total += w(i)*sqrt(pow(x(i)-y(i),2.0));
  }
  return total;
}

double sqrtDist (NumericVector x, NumericVector y, NumericVector w){
  int n = y.length();
  double total = 0;
  for (int i = 0; i < n ; ++i) {
    total += w(i)*sqrt(abs(x(i)-y(i)));
  }
  return total;
}

NumericMatrix calcWdist_abs (NumericMatrix x, NumericMatrix w){

  int outrows = x.ncol();
  int outcols = x.ncol();

  NumericMatrix out(outrows,outcols);

  for (int i = 0 ; i < outrows - 1; i++){
    NumericVector wi = 1 / w(_, i);

    for (int j = i + 1  ; j < outcols ; j ++) {
      NumericVector wj = 1 / w(_, j);

      NumericVector v1 = x(_, i);
      NumericVector v2 = x(_, j);

      double d = absDist(v1, v2, 1 / (wi + wj));

      out(j,i) = d;
      out(i,j) = d;
    }
  }
  return (out) ;
}

// [[Rcpp::export]]
NumericMatrix calcWdist_euclidean (NumericMatrix x, NumericMatrix w){

  int outrows = x.ncol();
  int outcols = x.ncol();

  NumericMatrix out(outrows,outcols);

  for (int i = 0 ; i < outrows - 1; i++){
    NumericVector wi = 1 / w(_, i);

    for (int j = i + 1  ; j < outcols ; j ++) {
      NumericVector wj = 1 / w(_, j);

      NumericVector v1 = x(_, i);
      NumericVector v2 = x(_, j);

      double d = eucDist(v1, v2, 1 / (wi + wj));

      out(j,i) = d;
      out(i,j) = d;
    }
  }
  return (out) ;
}


// [[Rcpp::export]]
NumericMatrix calcWdist_sqrt (NumericMatrix x, NumericMatrix w){

  int outrows = x.ncol();
  int outcols = x.ncol();

  NumericMatrix out(outrows,outcols);

  for (int i = 0 ; i < outrows - 1; i++){
    NumericVector wi = 1 / w(_, i);

    for (int j = i + 1  ; j < outcols ; j ++) {
      NumericVector wj = 1 / w(_, j);

      NumericVector v1 = x(_, i);
      NumericVector v2 = x(_, j);

      double d = sqrtDist(v1, v2, 1 / (wi + wj));

      out(j,i) = d;
      out(i,j) = d;
    }
  }
  return (out) ;
}
