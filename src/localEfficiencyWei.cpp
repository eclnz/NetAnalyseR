//' @export

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Assuming floydWarshallRcpp is defined as follows:
arma::mat floydWarshallRcpp(const arma::mat& graph);

// [[Rcpp::export]]
arma::vec localEfficiencyWei(const arma::mat& W_original) {
  arma::mat W = W_original;
  W.diag().zeros(); // Remove self-loops

  int numNodes = W.n_rows;
  arma::vec localEfficiencies(numNodes, fill::zeros); // Initialize vector for local efficiencies
  
  // Convert weights to lengths
  arma::mat lengths = 1 / W; // Assuming length_inversion() translates to this
  lengths.diag().zeros(); // Ensure diagonal elements are zero
  
  for (int node = 0; node < numNodes; ++node) {
    arma::uvec neighbors = find(W.row(node) > 0); // Identify neighbors of the current node
    
    if (neighbors.size() <= 1) continue; // Skip if node has 1 or fewer neighbors

    // Create subgraph
    arma::mat subgraphLengths = lengths(neighbors, neighbors);

    if (subgraphLengths.max() == 0) continue; // Skip if neighbors are not connected
    
    // Calculate the shortest paths in the subgraph
    arma::mat distances = 1 / floydWarshallRcpp(subgraphLengths);
    distances.diag().zeros(); // Zero out the diagonal to avoid division by zero in efficiency calculation
    
    // Calculate efficiency
    arma::mat efficiencyMatrix = distances + trans(distances);
    arma::vec numerator = sum(efficiencyMatrix, 1);
    double denominator = sum(numerator);

    if (denominator != 0) {
      localEfficiencies(node) = sum(numerator) / (denominator * 2);
    }
  }
  
  return localEfficiencies;
}
