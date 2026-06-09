#' @title Calculate Global Efficiency of a Weighted Graph
#'
#' @description This function calculates the global efficiency of a weighted graph based on the weighted adjacency matrix.
#' Global efficiency is a measure of the average inverse shortest path length in the graph, providing an indication
#' of how efficiently information is exchanged across the entire network. The function validates the input matrix,
#' sets diagonal elements to zero to ignore self-loops, converts the weights to lengths, calculates the inverse of
#' the shortest path length for each pair of nodes, and then computes the global efficiency from these values.
#'
#' @param W A square, symmetric matrix representing the weighted adjacency matrix of an undirected graph.
#'          Weights should be non-negative, and diagonal elements (self-loops) are ignored in the calculation.
#' @return A numeric value representing the global efficiency of the graph.
#' @examples
#' # Example: Create a 3x3 weighted adjacency matrix
#' W <- matrix(c(0, 2, 1, 2, 0, 3, 1, 3, 0), nrow = 3, byrow = TRUE)
#' global_efficiency_wei(W)
#' @export
global_efficiency_wei <- function(W) {
  global_efficiency_wei_(validate_matrix(W))
}
