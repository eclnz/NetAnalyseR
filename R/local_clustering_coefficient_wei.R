#' @title Calculate Local Clustering Coefficient for Weighted Undirected Networks
#'
#' @description This function computes the local clustering coefficient for each node
#' in a weighted undirected network. The clustering coefficient is a measure
#' of the degree to which nodes in a network cluster together, taking into
#' account the weight of the edges. This implementation is based on the method
#' proposed by Onnela et al. (2005), which is suitable for networks where
#' edge weights represent the strength of connections between nodes.
#'
#' @param W A square, symmetric matrix representing the weighted undirected
#'          connection matrix of the network.
#' @return A numeric vector containing the local clustering coefficient for
#'         each node in the network.
#' @examples
#' # Example usage with a simple weighted network
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' local_clustering_coefficient_wei(W)
#' @export
local_clustering_coefficient_wei<- function(W) {
  # Ensure W is a valid matrix and convert it to matrix form if necessary
  W <- validate_matrix(W)
  W <- as.matrix(W) # Convert to matrix if not already
  diag(W) <- 0 # Zero out the diagonal to ignore self-loops

  # Compute the cube root of the weights, as per Onnela's formula
  W_13 <- W ^ (1/3) # Element-wise power for transformation

  # Calculate the numerator of the clustering coefficient (cycle triangles)
  cyc3 <- diag(W_13 %*% W_13 %*% W_13) # Triadic closure

  # Calculate the degree of each node based on the number of non-zero connections
  K <- rowSums(W != 0)

  # Adjust degrees for nodes with no triadic closure to avoid division by zero
  K[cyc3 == 0] <- Inf # Set to Inf, effectively making clustering coefficient 0

  # Compute the local clustering coefficient for each node
  C <- cyc3 / (K * (K - 1)) # Onnela's clustering coefficient formula

  # Return the vector of clustering coefficients
  return(C)
}
