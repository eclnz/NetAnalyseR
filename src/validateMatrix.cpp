//' @keywords internal

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat validateMatrix(const arma::mat& W) {
  // Check if the matrix is square
  if (W.n_rows != W.n_cols) {
    stop("Matrix W must be square.");
  }
  
  // Check for NA values in the matrix
  if (W.has_nan()) {
    stop("Matrix contains NA values, which are not allowed.");
  }
  
  // Check if the matrix is symmetric
  bool isSymmetric = arma::approx_equal(W, W.t(), "absdiff", 0.0001);
  if (!isSymmetric) {
    // Make the matrix symmetric by mirroring the lower triangular part onto the upper triangular part
    arma::mat symW = arma::symmatl(W);
    return symW;
  }

  // If the matrix passes all checks and is symmetric, return it as is
  return W;
}
