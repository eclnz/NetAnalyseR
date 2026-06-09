#' @title Calculate the Local Efficiency of a Weighted Graph
#'
#' @description This function computes the local efficiency of a weighted graph. Local efficiency is a measure of how well information
#' is exchanged by the neighbors of a given node when that node is removed. The function processes a weighted adjacency
#' matrix representing the graph, validates it, sets self-loops to zero, and applies a weight conversion to treat the weights
#' as lengths. It then calculates the local efficiency for each node by considering the efficiency of its subgraph formed by its
#' immediate neighbors.
#'
#' @param W A square, symmetric matrix representing the weighted adjacency matrix of an undirected graph.
#'          Weights should be non-negative, and diagonal elements (self-loops) are ignored.
#' @return A numeric vector of length equal to the number of nodes in the graph, where each element represents
#'         the local efficiency of the corresponding node.
#' @examples
#' # Example: Create a 4x4 weighted adjacency matrix
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' local_efficiency_wei(W)
#' @export
#'
local_efficiency_wei_ <- function(W) {
  localEfficiencyWei(W)
}

local_efficiency_wei <- function(W) {
  local_efficiency_wei_(validate_matrix(W))
}
