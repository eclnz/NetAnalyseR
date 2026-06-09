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
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' node_strength(W)
#' @export
node_strength_ <- function(W) {
  diag(W) <- 0
  colSums(W)
}

node_strength <- function(W) {
  node_strength_(validate_matrix(W))
}
