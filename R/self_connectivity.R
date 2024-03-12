#' @title Calculate Self-Connection Weights from a Network Matrix
#'
#' @description This function returns the weights of self-connections within a network,
#' where each self-connection is a connection from a node to itself. It directly extracts
#' the diagonal of the input matrix, which represents these self-connections.
#'
#' @param W A square, numeric matrix representing a weighted connection matrix of the network.
#' Self-connections are represented by the diagonal elements of this matrix.
#' @return A numeric vector containing the weights of self-connections for each node in the network.
#' @examples
#' # Create a 3x3 matrix with self-connections
#' W <- matrix(c(1, 2, 3,
#'               4, 5, 6,
#'               7, 8, 9), byrow = TRUE, nrow = 3)
#' self_connectivity(W)
#'
#' @export
#'
self_connectivity <- function(W) {
  W <- validate_matrix(W) # Ensure the matrix is valid for analysis

  # Directly return the diagonal of the matrix, which represents self-connections
  # Each element on the diagonal is the weight of the self-connection for that node
  self_connections_weights <- diag(W)

  return(self_connections_weights)
}
