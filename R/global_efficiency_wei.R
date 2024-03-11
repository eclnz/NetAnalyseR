#' Calculate Global Efficiency of a Weighted Graph
#'
#' This function calculates the global efficiency of a weighted graph based on the weighted adjacency matrix.
#' Global efficiency is a measure of the average inverse shortest path length in the graph, providing an indication
#' of how efficiently information is exchanged across the entire network. The function validates the input matrix,
#' sets diagonal elements to zero to ignore self-loops, converts the weights to lengths, calculates the inverse of
#' the shortest path length for each pair of nodes, and then computes the global efficiency from these values.
#'
#' @param W A square, symmetric matrix representing the weighted adjacency matrix of an undirected graph.
#'          Weights should be non-negative, and diagonal elements (self-loops) are ignored in the calculation.
#' @return A numeric value representing the global efficiency of the graph.
#' @examples
#' # Example: Create a 3x3 weighted adjacency matrix
#' W <- matrix(c(0, 2, 1, 2, 0, 3, 1, 3, 0), nrow = 3, byrow = TRUE)
#' global_efficiency_wei(W)
#' @export
global_efficiency_wei <- function(W) {
  # Validate the input matrix to ensure it is a proper adjacency matrix for a graph
  W <- validate_matrix(W)
  
  # Remove self-loops by setting diagonal elements to zero
  diag(W) <- 0
  
  # Number of nodes in the graph
  numNodes <- nrow(W)
  
  # Conversion factor for cube root, used in distance calculations
  cubeRootFactor <- 1 / 3
  
  # Convert weights to lengths for distance calculation
  lengthsMatrix <- length_inversion(W)
  
  # Calculate the inverse distance matrix, representing the efficiency of each pair of nodes
  inverseDistanceMatrix <- 1 / shortest_distance(lengthsMatrix)
  
  # Calculate global efficiency as the sum of upper triangular part of the inverse distance matrix
  # divided by the total number of possible connections (excluding self-connections)
  globalEfficiency <- sum(inverseDistanceMatrix[upper.tri(inverseDistanceMatrix)], na.rm = TRUE) / (numNodes * (numNodes - 1) / 2)
  
  # Return the global efficiency of the graph
  return(globalEfficiency)
}
