//' @keywords internal

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec localClusteringCoefficientWei(const arma::mat& W_original) {
  // Clone the original matrix to avoid modifying it
  arma::mat W = W_original;
  
  // Zero out the diagonal to ignore self-loops
  W.diag().zeros();
  
  // Compute the cube root of the weights, as per Onnela's formula
  arma::mat W_13 = arma::pow(W, 1.0 / 3.0);
  
  // Calculate the numerator of the clustering coefficient (cycle triangles)
  arma::vec cyc3 = arma::diagvec(W_13 * W_13 * W_13);
  
  // Calculate the degree of each node based on the number of non-zero connections
  arma::vec K = arma::sum(arma::conv_to<arma::mat>::from(W != 0), 1);
  
  // Initialize clustering coefficients with NaN to handle division by zero or undefined cases
  arma::vec C = arma::vec(size(K)).fill(arma::datum::nan);

  // Identify nodes with at least one triadic closure and a degree of at least 2
  arma::uvec valid_indices = find(cyc3 > 0 && K > 1);

  // Compute the local clustering coefficient only for valid nodes
  C.elem(valid_indices) = cyc3.elem(valid_indices) / (K.elem(valid_indices) % (K.elem(valid_indices) - 1));
  
  // For nodes with a degree of less than 2, set clustering coefficient to 0
  C.elem(find(K <= 1)).zeros();
    
  // Correct NaN values to 0 after computation
  C.replace(arma::datum::nan, 0.0);

  // Return the vector of clustering coefficients
  return C;
}

