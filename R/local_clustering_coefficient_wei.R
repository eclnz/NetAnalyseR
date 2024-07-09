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
#' @param validate Whether to validate the input matrix.
#' @return A numeric vector containing the local clustering coefficient for
#'         each node in the network.
#' @examples
#' # Example usage with a simple weighted network
#' W <- matrix(c(0, 2, 1, 4, 2, 0, 3, 5, 1, 3, 0, 6, 4, 5, 6, 0), nrow = 4, byrow = TRUE)
#' local_clustering_coefficient_wei(W)
#' @export
local_clustering_coefficient_wei<- function(W, validate = TRUE) {
  # Ensure W is a valid matrix and convert it to matrix form if necessary
  if(validate){validate_matrix(W)}
  # Call C++ function
  C <- localClusteringCoefficientWei(W)
  # Return the vector of clustering coefficients
  return(C)
}
