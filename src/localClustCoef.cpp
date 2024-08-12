//' @keywords internal

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec localClusteringCoefficientWei(const arma::mat& W_original) {
  // Clone the original matrix to avoid modifying it
  arma::mat W = W_original;
  
  // Create the adjacency matrix
  arma::mat A = arma::conv_to<arma::mat>::from(W != 0);
  
  // Compute the symmetrized weights matrix raised to the power of 1/3
  arma::mat S = arma::pow(W, 1.0 / 3.0) + arma::pow(W.t(), 1.0 / 3.0);
  
  // Calculate the total degree (sum of in-degrees and out-degrees)
  arma::vec K = arma::sum(A + A.t(), 1);
  
  // Calculate the number of 3-cycles (i.e., directed triangles)
  arma::vec cyc3 = arma::diagvec(S * S * S) / 2.0;
  
  // Set K to infinity where no 3-cycles exist to ensure C=0
  K.elem(find(cyc3 == 0)).fill(arma::datum::inf);
  
  // Calculate the number of all possible 3-cycles
  arma::vec CYC3 = K % (K - 1) - 2 * arma::diagvec(A * A);
  
  // Calculate the clustering coefficient
  arma::vec C = cyc3 / CYC3;
  
  // Replace NaN values with 0.0 for valid indices (optional, usually not needed)
  C.replace(arma::datum::nan, 0.0);

  // Return the vector of clustering coefficients
  return C;
}

