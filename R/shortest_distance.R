#' Calculate the Shortest Path Between Nodes in a Graph
#'
#' @description
#' Computes the shortest paths between nodes in a graph using a variant of Dijkstra's algorithm.
#' The function takes a matrix of distances (lengths) between nodes, where the distance inversely
#' represents the strength of connection between the nodes. Stronger connections have shorter
#' distances. This function is particularly useful for processing matrices derived from weight
#' conversion processes, like `weight_conversion(W, 'lengths')`.
#'
#' @param L A numeric matrix representing the lengths or distances between nodes in the graph.
#' The matrix should be square, with dimensions N x N, where N is the number of nodes. Each
#' element L[i, j] represents the distance from node i to node j.
#'
#' @return A matrix of the same dimension as `L`, where each element [i, j] represents the
#' shortest distance from node i to node j in the graph. If a node is unreachable, the distance
#' is set to 0.
#'
#' @examples
#' W <- matrix(runif(16), 4, 4) # Random weights for illustration
#' L <- weight_conversion(W, 'lengths') # Hypothetical function to convert weights to lengths
#' shortest_paths <- shortest_distance(L)
#'
#' @export
#'
shortest_distance <- function(L) {
  L <- validate_matrix(L) # Check if the input matrix is valid
  if(all(L == floor(L))) {
    warning("All lengths are integers, this can be a sign you have not converted the matrix to a length matrix")
  } # Warn if the matrix doesn't seem to be converted to lengths

  # Initialize the distance matrix and set diagonal distances to zero
  diag(L) <- 0 # Self-distances are always zero
  n <- nrow(L) # Number of nodes in the graph
  D <- matrix(Inf, n, n) # Initialize the distance matrix with infinite distances
  diag(D) <- 0 # Distance to self is zero

  # Main loop to compute shortest distances using Dijkstra's algorithm
  for (u in 1:n) {
    S <- rep(TRUE, n) # Initialize all distances as temporary (TRUE)
    L1 <- L # Copy of the length matrix for manipulation
    V <- u # Start with the source node
    repeat {
      S[V] <- FALSE # Mark the current node's distance as permanent
      L1[, V] <- 0 # Remove in-edges to the current node as it's already processed

      # Update distances to neighbors of the current node
      for (v in V) {
        T <- which(L1[v, ] > 0) # Find neighbors of the current node
        if (length(T) == 0) next # Skip if no neighbors

        # Loop through neighbors to update distances
        for (t in T) {
          alt <- D[u, v] + L1[v, t] # Alternative path distance through node v
          if (alt < D[u, t]) {
            D[u, t] <- alt # Update distance if it's shorter
          }
        }
      }

      # Check for completion: if all nodes are permanent, or no reachable nodes left
      if (all(!S)) break
      minD <- min(D[u, S], na.rm = TRUE)
      if (is.infinite(minD)) break

      V <- which(D[u, ] == minD) # Find the next node(s) with the shortest tentative distance
    }
  }

  D[is.infinite(D)] <- 0 # Replace infinite distances with 0 indicating no path exists
  return(D)
}
