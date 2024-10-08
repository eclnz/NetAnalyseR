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
#' @param validate Whether to validate the input matrix.
#' @return A numeric vector of length equal to the number of nodes in the graph, where each element represents
#'         the local efficiency of the corresponding node.
#' @examples
#' # Example: Create a 4x4 weighted adjacency matrix
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' local_efficiency_wei(W)
#' @export
#'
local_efficiency_wei <- function(W, validate = TRUE) {
  # Validate the input matrix to ensure it is a proper adjacency matrix for a graph
  if(validate){validate_matrix(W)}
  # Remove self-loops by setting diagonal elements to zero
  # diag(W) <- 0
  # Determine if there is a connection (edge) between nodes
  isConnected <- W > 0
  # Number of nodes in the graph
  numNodes <- nrow(W)
  # Conversion factor for cube root
  cubeRootFactor <- 1 / 3
  # Convert weights to lengths for distance calculation
  lengths <- length_inversion(W, FALSE)
  # Initialize vector for local efficiencies
  localEfficiencies <- numeric(numNodes)
  # Calculate cube root of weights and lengths for harmonic mean calculation
  cubeRootWeights <- W ^ cubeRootFactor
  cubeRootLengths <- lengths ^ cubeRootFactor
  # Ensure diagonal elements of cube root lengths are zero
  # diag(cubeRootLengths) <- 0
  for (node in 1:numNodes) {
    # Identify neighbors of the current node
    neighbors <- which(isConnected[node, ] | isConnected[, node])
    # Nodes with one or less neighbors have an efficiency of zero.
    if (length(neighbors) <= 1) {
      localEfficiencies[node] <- 0
      next
    }
    # Subgraph weights for neighbors, using cube root weights to calculate subgraph
    subgraphWeights <- (as.matrix(cubeRootWeights[node, neighbors] + cubeRootWeights[neighbors, node]))
    # Subgraph lengths for neighbors, now using cube root lengths similar to `local_efficiency`
    subgraphLengths <- cubeRootLengths[neighbors, neighbors]
    # Calculate distance matrix for the subgraph using floydWarshallRcpp, mirroring `local_efficiency` logic
    distances <- 1 / shortest_distance(subgraphLengths, FALSE)
    # Symmetric efficiency matrix for the subgraph, derived from distances
    efficiencyMatrix <- distances + t(distances)
    # Pre-sum matrix for numerator calculation, combining weights and efficiency
    preSumMatrix <- ((subgraphWeights %*% t(subgraphWeights)) * efficiencyMatrix)
    # Numerator: sum of pre-sum matrix elements, halved and ensuring finite values
    numerator <- sum(preSumMatrix[is.finite(preSumMatrix)], na.rm = TRUE) / 2
    if (numerator != 0) {
      # Sum of adjacency values for normalization, reflecting connectivity
      adjacencySums <- isConnected[node, neighbors] + isConnected[neighbors, node]
      # Denominator for efficiency calculation, based on connectivity sums
      denominator <- sum(adjacencySums)^2 - sum(adjacencySums^2)
      # Calculate and assign local efficiency
      localEfficiencies[node] <- numerator / denominator
    }
  }
  # Return the vector of local efficiencies for each node
  return(localEfficiencies)
}


