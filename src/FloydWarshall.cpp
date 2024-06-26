//' @name floydWarshallRcpp
//' @title Floyd-Warshall Algorithm in R
//' @description This function implements the Floyd-Warshall algorithm using Rcpp for finding the shortest distances
//' between every pair of vertices in a given edge weighted directed Graph. It operates on a matrix representation of the graph.
//' @param inputMatrix A NumericMatrix representing the adjacency matrix of the graph where the element at the ith row
//' and jth column represents the weight of the edge from vertex i to vertex j. If there is no edge between vertex i and vertex j,
//' the value should be 0 for i != j.
//' @return A NumericMatrix where the element at the ith row and jth column represents the shortest distance from vertex i
//' to vertex j in the input graph.
//' @examples
//' W <- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
//' shortestPaths <- floydWarshallRcpp(W)
//' @export
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat floydWarshallRcpp(const arma::mat& inputMatrix) {
    int numVertices = inputMatrix.n_rows; // Assuming a square matrix
    arma::mat distanceMatrix = inputMatrix;

    // Initialize the distance matrix
    for (unsigned int i = 0; i < numVertices; ++i) {
        for (unsigned int j = 0; j < numVertices; ++j) {
            if (inputMatrix(i, j) == 0 && i != j) {
                // Use Armadillo's representation of infinity for unreachable vertices
                distanceMatrix(i, j) = arma::datum::inf;
            }
        }
    }

    // Apply the Floyd-Warshall algorithm
    for (unsigned int k = 0; k < numVertices; ++k) {
        for (unsigned int i = 0; i < numVertices; ++i) {
            for (unsigned int j = 0; j < numVertices; ++j) {
                if (distanceMatrix(i, k) + distanceMatrix(k, j) < distanceMatrix(i, j)) {
                    distanceMatrix(i, j) = distanceMatrix(i, k) + distanceMatrix(k, j);
                }
            }
        }
    }

    return distanceMatrix;
}
