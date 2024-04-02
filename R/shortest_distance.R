#' @title Calculate the Shortest Path Between Nodes in a Graph
#'
#' @description
#' Computes the shortest paths between nodes in a graph using the Floyd Warshall Algorithm.
#' The function takes a matrix of distances (lengths) between nodes, where the distance inversely
#' represents the strength of connection between the nodes. Stronger connections have shorter
#' distances. This function is particularly useful for processing length matrices derived from the
#' function `length_inversion(W)`.
#'
#' @param L A numeric matrix representing the lengths or distances between nodes in the graph.
#' The matrix should be square, with dimensions N x N, where N is the number of nodes. Each
#' connection represents the distance from node i to node j.
#' @return A matrix of the same dimension as L, where each connection represents the
#' shortest distance from node i to node j in the graph. If a node is unreachable, the distance
#' is set to 0.
#' @examples
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' L <- length_inversion(W) # Hypothetical function to convert weights to lengths
#' shortest_paths <- shortest_distance(L)
#'
#' @export
#'
shortest_distance <- function(L) {
  # Check input matrix
  L <- validate_matrix(L)

  # Initialize the distance matrix and set diagonal distances to zero
  diag(L) <- 0 # Self-distances are always zero

  # Compute shortest distance with Floyd Warshall Algorithm in C++ call
  D <- floydWarshallRcpp(L)

  return(D)
}
