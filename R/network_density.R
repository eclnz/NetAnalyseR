#' @title Calculate Network Density for an Undirected Network
#'
#' @description This function calculates the density of an undirected network, defined as the ratio of
#' actual connections to possible connections, excluding self-connections. It transforms the input weighted
#' connection matrix into a binary adjacency matrix (where connections are marked as 1, no connections as 0)
#' and computes the density based on this binary representation.
#'
#' @param W A square, numeric matrix representing a weighted connection matrix of the network.
#' The matrix should be symmetric to represent an undirected network.
#' @return A single numeric value representing the network's density.
#' @examples
#' # Create a 4x4 undirected, weighted connection matrix
#' W <- matrix(c(0, 1, 3, 2,
#'               1, 0, 2, 4,
#'               3, 2, 0, 5,
#'               2, 4, 5, 0), byrow = TRUE, nrow = 4)
#' network_density(W)
#' @export
network_density <- function(W) {
  W <- validate_matrix(W) # Validate the input matrix
  diag(W) <- 0 # Set diagonal to zero to ignore self-connections

  # Convert the weighted connection matrix to a binary adjacency matrix
  # where 1 represents the presence of a connection and 0 represents its absence
  A <- W > 0

  # Calculate network density by averaging the values of the lower triangle of the binary adjacency matrix.
  # This avoids counting each connection twice in an undirected network.
  density <- mean(A[lower.tri(A)])

  return(density)
}
