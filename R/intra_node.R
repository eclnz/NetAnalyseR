#' @title Calculate the Total Weight of Self-Connections in a Network
#'
#' @description This function computes the total weight of self-connections within a network,
#' where a self-connection is a connection that initiates and terminates on the same node.
#'
#' @param W A square, numeric matrix representing a weighted connection matrix of the network.
#' Self-connections are represented by the diagonal elements of the matrix.
#' @return A single numeric value representing the sum of weights of all self-connections in the network.
#' @examples
#' # Create a 3x3 matrix with self-connections
#' W <- matrix(c(1, 2, 3,
#'               4, 5, 6,
#'               7, 8, 9), byrow = TRUE, nrow = 3)
#' intra_node(W)
#' @export
intra_node_ <- function(W) {
  sum(diag(W))
}

intra_node <- function(W) {
  intra_node_(validate_matrix(W))
}
