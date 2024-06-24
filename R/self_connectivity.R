#' @title Calculate Self-Connection Weights from a Network Matrix
#'
#' @description This function returns the weights of self-connections within a network,
#' where each self-connection is a connection from a node to itself. It directly extracts
#' the diagonal of the input matrix, which represents these self-connections.
#'
#' @param W A square, numeric matrix representing a weighted connection matrix of the network.
#' @param validate Whether to validate the input matrix.
#' Self-connections are represented by the diagonal elements of this matrix.
#' @return A numeric vector containing the weights of self-connections for each node in the network.
#' @examples
#' # Create a 3x3 matrix with self-connections
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' self_connectivity(W)
#' @export
self_connectivity <- function(W, validate = TRUE) {
  if(validate){validate_matrix(W)} # Ensure the matrix is valid for analysis
  # Each element on the diagonal is the weight of the self-connection for that node
  self_connections_weights <- diag(W)
  return(self_connections_weights)
}
