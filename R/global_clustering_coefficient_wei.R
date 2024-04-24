#' @title Calculate the Global Clustering Coefficient for a Weighted Undirected Graph
#'
#' @description This function calculates the global clustering coefficient of a weighted undirected graph.
#' The global clustering coefficient is a measure of the degree to which nodes in a graph tend to cluster together.
#' The function first validates the input matrix to ensure it represents a proper adjacency matrix of a graph,
#' sets the diagonal elements to zero to ignore self-loops, and then computes the mean local clustering coefficient
#' across all nodes as the global clustering coefficient.
#'
#' @param W A square, symmetric matrix representing the weighted adjacency matrix of an undirected graph.
#'          The matrix should have non-negative weights, and diagonal elements are ignored in the computation.
#' @return A numeric value representing the global clustering coefficient of the graph.
#' @examples
#' # Example: Create a 3x3 weighted adjacency matrix
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' global_clustering_coefficient_wei(W)
#' @export
global_clustering_coefficient_wei <- function(W) {
  # Validate the input matrix to ensure it is a proper adjacency matrix for a graph
  validate_matrix(W) # Check if matrix is valid
  # Remove self-loops by setting diagonal elements to zero
  Cl = mean(localClusteringCoefficientWei(W))
  # Return the global clustering coefficient
  return(Cl)
}
