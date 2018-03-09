#include <Rcpp.h>
using namespace Rcpp;


// DON'T FORGET TO compileAttributes("albatross/")!!!

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



//' Compute Euclidean distances between a specified row and all other row
//' Returns an unscaled positive numeric vector
//'
//'
//'
//' @param X A numeric matrix
//' @param x a row index
//' @export
// [[Rcpp::export]]
NumericVector dist_one(NumericMatrix X, int x) {
  x = x - 1;
  int nrow = X.nrow(), ncol = X.ncol();
  NumericVector out(nrow);

  for(int i = 0; i < nrow; ++i) {
    for(int j = 0; j < ncol; ++j) {
      out[i] += pow(X(i,j) - X(x,j), 2.0);
    }
    out[i] = sqrt(out[i]);
  }
  return out;
}


//' Compute Euclidean distances between all pairs of rows
//' Returns numeric matrix, scaled to [0,1]
//'
//'
//'
//' @param X A numeric matrix
//' @export
// [[Rcpp::export]]
NumericMatrix dist_all(NumericMatrix X) {
  int nrow = X.nrow(), ncol = X.ncol();
  NumericMatrix out(nrow, nrow);
  int i, j, k;

  for(j = 0; j < nrow; ++j) {
    for(i = 0; i < nrow; ++i) {
      for(k = 0; k < ncol; ++k) {
        out(i,j) += pow(X(i,k) - X(j,k), 2.0);
      }
      out(i,j) = sqrt(out(i,j));
    }
  }
  return out/max(out);
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
