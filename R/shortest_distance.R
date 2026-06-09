#' @title Calculate the Shortest Path Between Nodes in a Graph
#'
#' @description
#' Computes the shortest paths between nodes in a graph using the Floyd-Warshall or
#' Dijkstra algorithm, selected automatically based on network density. The function
#' takes a matrix of distances (lengths) between nodes, where the distance inversely
#' represents the strength of connection between the nodes. Stronger connections have
#' shorter distances. This function is particularly useful for processing length matrices
#' derived from the function `length_inversion(W)`.
#'
#' @param L A numeric matrix representing the lengths or distances between nodes in the graph.
#' The matrix should be square, with dimensions N x N, where N is the number of nodes. Each
#' connection represents the distance from node i to node j.
#' @return A matrix of the same dimension as L, where each element represents the
#' shortest distance from node i to node j in the graph.
#' @examples
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' L <- length_inversion(W)
#' shortest_paths <- shortest_distance(L)
#'
#' @export
#'
# Floyd-Warshall O(n^3) outperforms Dijkstra O(n * E * log n) when the graph is
# dense (most node pairs directly connected); empirically the crossover is ~0.55.
shortest_distance_ <- function(L) {
  diag(L) <- 0
  if (network_density_(L) > 0.55) floydWarshallRcpp(L) else dijkstraAllPairs(L)
}

shortest_distance <- function(L) {
  shortest_distance_(validate_matrix(L))
}
