#' @title Calculate the Total Weight of Self-Connections in a Network
#'
#' @description This function computes the total weight of self-connections within a network,
#' where a self-connection is a connection that initiates and terminates on the same node.
#'
#' @param W A square, numeric matrix representing a weighted connection matrix of the network.
#' Self-connections are represented by the diagonal elements of the matrix.
#' @param validate Whether to validate the input matrix.
#' @return A single numeric value representing the sum of weights of all self-connections in the network.
#' @examples
#' # Create a 3x3 matrix with self-connections
#' W <- matrix(c(1, 2, 3,
#'               4, 5, 6,
#'               7, 8, 9), byrow = TRUE, nrow = 3)
#' intra_node(W)
#' @export
intra_node <- function(W, validate = TRUE) {
  if(validate){validate_matrix(W)} # Validate the input matrix

  # Extract the diagonal of the matrix to focus on self-connections
  # The diagonal contains the weights of connections that initiate and terminate on the same node
  self_connections_weights <- diag(W)

  # Calculate the total weight of self-connections by summing up the weights on the diagonal
  total_self_connections_weight <- sum(self_connections_weights)

  return(total_self_connections_weight)
}
