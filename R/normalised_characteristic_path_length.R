#' Calculate Normalized Characteristic Path Length
#'
#' @title Normalized Characteristic Path Length Calculation
#' @description Computes the normalized characteristic path length for a given network.
#' Accepts either a numeric matrix or a list containing a matrix and array of random matrices already generated.
#' If no array of randomized matrices is provided, one is generated.
#' The list elements must not be named.
#' @param mat_list Either a square matrix or a list containing the network matrix and an optional 3D array of randomized matrices.
#' @param rand_array An optional 3D array of randomized matrices for comparison if the input is a single matrix.
#' @return The normalized characteristic path length.
#' @examples
#' W <- matrix(c(0, 2, 1, 0, 2, 0, 3, 5, 1, 3, 0, 6, 0, 5, 6, 0), nrow = 4, byrow = TRUE)
#' norm_char_path_length <- normalised_characteristic_path_length(W)
#' W_rand <- generateRewiredMatrices(W, 100)
#' W_rand <- abind::abind(W_rand, along = 3)
#' norm_char_path_length <- normalised_characteristic_path_length(list(W, W_rand))
#' norm_char_path_length <- normalised_characteristic_path_length(W, W_rand)
#' @export

normalised_characteristic_path_length <- function(mat_list, rand_array = NULL) {
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
    return(normalised_characteristic_path_length_(W, rand_array))
  }
  if (is.list(mat_list)) {
    W <- mat_list[[1]]
    ra <- mat_list[[2]]
    if (network_density_(W) == 1) {
      warning("Network density is equal to 1. Network cannot be rewired while maintaining degree distribution\n
              This may make results incorrect.")
    }
    if (!is.array(ra)) stop("Random array specified is not in the form of an array")
    return(normalised_characteristic_path_length_(W, ra))
  }
}
