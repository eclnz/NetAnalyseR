#' @title Calculate the Total Inter-Node Connections in a Network
#'
#' @description This function calculates the total number of connections between nodes in a network,
#' specifically focusing on undirected networks. It processes a weighted connection matrix,
#' setting self-connections to zero and considering only unique connections due to the network's undirected nature.
#'
#' @param W A square, numeric matrix representing a weighted connection matrix of the network.
#' The matrix should represent an undirected network, meaning connections are symmetric.
#' @return A single numeric value representing the sum of weights of all unique inter-node connections.
#' @examples
#' # Create a 4x4 undirected, weighted connection matrix
#' W <- matrix(c(0, 1, 3, 2,
#'               1, 0, 2, 4,
#'               3, 2, 0, 5,
#'               2, 4, 5, 0), byrow = TRUE, nrow = 4)
#' inter_node(W)
#' @export
inter_node_ <- function(W) {
  diag(W) <- 0
  sum(W[lower.tri(W)])
}

inter_node <- function(W) {
  inter_node_(validate_matrix(W))
}
