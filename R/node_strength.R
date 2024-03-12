#' @title Calculate the Strength of Nodes in an Undirected Network
#'
#' @description This function computes the node strength in an undirected network,
#' defined as the sum of weights of edges connected to each node.
#'
#' @param W A square, numeric matrix representing an undirected, weighted connection matrix of the network.
#' The elements of the matrix represent the weights of the edges between nodes. The matrix must be symmetric.
#' @return A numeric vector containing the strength of each node.
#' @examples
#' # Create a 3x3 undirected, weighted connection matrix
#' W <- matrix(c(0, 2, 1,
#'               2, 0, 3,
#'               1, 3, 0), nrow = 3, byrow = TRUE)
#' node_strength(W)
#' @export
node_strength <- function(W) {
  W <- validate_matrix(W) # Validate the input matrix to ensure it's appropriate for calculations
  diag(W) <- 0 # Remove self-loops by setting the diagonal elements to zero

  # Calculate node strength as the sum of weights of links connected to the node.
  # For an undirected network, column sums and row sums are equivalent.
  node_strength <- colSums(W)

  return(node_strength)
}
