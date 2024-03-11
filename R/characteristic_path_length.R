#' Calculate the Characteristic Path Length of a Graph
#'
#' This function computes the characteristic path length of a graph,
#' which is defined as the average length of the shortest paths
#' for all pairs of nodes in the graph. The function expects a
#' distance matrix where the element at the i-th row and j-th
#' column represents the length of the shortest path from node i
#' to node j.
#'
#' @param W A square matrix adjacency matrix, where edges represent connection
#'          strength. 
#' @return The characteristic path length of the graph as a single
#'         numerical value. This represents the average shortest path
#'         length between all pairs of nodes.
#' @examples
#' # Example usage with a simple 3-node (triangular) graph
#' D <- matrix(c(0, 1, 2,
#'               1, 0, 1,
#'               2, 1, 0), nrow = 3, byrow = TRUE)
#' characteristic_path_length(W)
#' @export
#' 
characteristic_path_length <- function(W) {
  # Ensure the input matrix is valid for the operation
  W <- validate_matrix(W) # Check if the matrix is valid
  
  # Set diagonal elements to zero to ignore self-paths in calculations
  diag(W) <- 0
  
  # Get the number of rows (nodes) in the matrix
  n <- nrow(W)
  
  # Calculate lengths
  L <- length_inversion(W)
  
  # Calculate shortest distance between nodes
  D <- shortest_distance(L)
  
  # Replace infinite distances with NA to ignore them in the mean calculation
  D[is.infinite(D)] <- NA  
  
  # Extract all valid (non-NA) distances from the matrix
  Dv <- D[!is.na(D)]  
  
  # Calculate the mean of all valid distances as the characteristic path length
  lambda <- mean(Dv)
  
  # Return the characteristic path length
  return(lambda)
}
