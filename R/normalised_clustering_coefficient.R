#' Calculate Normalized Clustering Coefficient for Single or List of Matrices
#'
#' @title Normalized Clustering Coefficient Calculation
#' @description Computes the normalized clustering coefficient for a given matrix.
#' The normalized clustering coefficient is calculated by dividing the global clustering coefficient of the input matrix
#' by the mean global clustering coefficient of a set of random matrices with a matched strength distribution.
#' Accepts either a numeric matrix or a list containing a matrix and array of random matrices already generated.
#' If no array of randomized matrices is provided and the input is a matrix, the function generates one.
#' If the global clustering coefficient of the input is 0, the function returns 0 to avoid division by zero errors.
#' The list elements must not be named.
#' @param mat_list A square matrix or a list containing a matrix and an optional 3D array of randomized matrices for comparison.
#' @param rand_array An optional 3D array of randomized matrices for comparison if the input is a single matrix.
#' @return The normalized clustering coefficient.
#' @examples
#' W <- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
#' norm_clust_coeff <- normalised_clustering_coefficient(W)
#' rand_matrices <- generateRewiredMatrices(W, 100)
#' rand_matrices <- abind::abind(rand_matrices, along = 3)
#' norm_clust_coeff <- normalised_clustering_coefficient(list(W, rand_matrices))
#' norm_clust_coeff <- normalised_clustering_coefficient(W, rand_matrices)
#' @importFrom abind abind
#' @export

normalised_clustering_coefficient <- function(mat_list, rand_array = NULL) {
  if (is.matrix(mat_list)) {
    W <- validate_matrix(mat_list)
    if (network_density_(W) == 1) {
      warning("Network density is equal to 1. Network cannot be rewired while maintaining degree distribution\n
              This may make results incorrect.")
    }
    if (is.null(rand_array)) {
      rand_array <- abind::abind(generateRewiredMatrices(W, 100), along = 3)
    }
    if (!is.array(rand_array)) stop("Random array specified is not in the form of an array")
    return(normalised_clustering_coefficient_(W, rand_array))
  }
  if (is.list(mat_list)) {
    W <- mat_list[[1]]
    ra <- mat_list[[2]]
    if (network_density_(W) == 1) {
      warning("Network density is equal to 1. Network cannot be rewired while maintaining degree distribution\n
              This may make results incorrect.")
    }
    if (!is.array(ra)) stop("Random array specified is not in the form of an array")
    return(normalised_clustering_coefficient_(W, ra))
  }
}
