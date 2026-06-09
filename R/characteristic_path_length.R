#' @title Calculate the Characteristic Path Length of a Graph
#'
#' @description This function computes the characteristic path length of a graph,
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
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' characteristic_path_length(W)
#' @export
#'
characteristic_path_length_ <- function(W) {
  diag(W) <- 0
  D <- shortest_distance_(length_inversion_(W))
  D[is.infinite(D)] <- NA
  mean(D[lower.tri(D)])
}

characteristic_path_length <- function(W) {
  characteristic_path_length_(validate_matrix(W))
}
